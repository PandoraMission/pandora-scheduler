#!/usr/bin/env python3
"""Exercise all new constraint features with synthetic GMAT data.

This script builds a small synthetic GMAT ephemeris valid for a 6-hour window,
then runs visibility computation under several configurations to demonstrate:

 1. Baseline (no day/night, no ST) — legacy-equivalent behaviour
 2. Day/night Earth avoidance (110° day, 80° night)
 3. Star tracker keepout with roll sweep
 4. Full pipeline: day/night + ST + roll

Each case prints a concise summary and the script produces CSV files in a
temporary directory.

Usage:
    poetry run python scripts/run_example_cases.py
"""

from __future__ import annotations

import textwrap
import warnings
from datetime import datetime, timedelta
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.catalog import (
    _build_base_payload,
    _build_star_visibility,
)
from pandorascheduler_rework.visibility.geometry import (
    build_minute_cadence,
    interpolate_gmat_ephemeris,
)


# ── Synthetic GMAT builder ──────────────────────────────────────────────────

def _to_gmat_mod_julian(times: Time) -> np.ndarray:
    """Convert astropy Time to GMAT UTCModJulian."""
    return np.asarray(times.to_value("mjd"), dtype=float) - 29999.5


def make_synthetic_gmat(tmp_path: Path, start: datetime, end: datetime) -> Path:
    """Build a realistic-ish synthetic GMAT ephemeris.

    ~96-minute LEO, 51.6° inclination, 600 km altitude.
    Sun along +X (ecliptic), Moon along +Y.
    """
    pad = timedelta(minutes=15)
    n_pts = int((end - start + 2 * pad).total_seconds() / 30) + 1
    times = Time(
        [start - pad + timedelta(seconds=30 * i) for i in range(n_pts)],
        format="datetime",
        scale="utc",
    )

    # Orbital parameters
    R = 6971.0  # km (altitude ≈ 600 km)
    period_s = 96.0 * 60  # ~96 min
    omega = 2 * np.pi / period_s
    inc = np.deg2rad(51.6)

    t_sec = np.array([(t - times[0]).sec for t in times])
    nu = omega * t_sec  # true anomaly proxy

    # Spacecraft position in ECI (simplified circular orbit)
    sc_x = R * np.cos(nu)
    sc_y = R * np.sin(nu) * np.cos(inc)
    sc_z = R * np.sin(nu) * np.sin(inc)

    # Sun: far away along +X (simplification for fixed geometry over 6h)
    au_km = 1.496e8
    sun_x = np.full(n_pts, au_km)
    sun_y = np.zeros(n_pts)
    sun_z = np.zeros(n_pts)

    # Moon: along +Y, ~384400 km
    moon_x = np.zeros(n_pts)
    moon_y = np.full(n_pts, 3.844e5)
    moon_z = np.zeros(n_pts)

    # Sub-satellite latitude/longitude
    lat = np.rad2deg(np.arcsin(sc_z / R))
    lon = np.rad2deg(np.arctan2(sc_y, sc_x))

    df = pd.DataFrame({
        "Earth.UTCModJulian": _to_gmat_mod_julian(times),
        "Earth.EarthMJ2000Eq.X": np.zeros(n_pts),
        "Earth.EarthMJ2000Eq.Y": np.zeros(n_pts),
        "Earth.EarthMJ2000Eq.Z": np.zeros(n_pts),
        "Pandora.EarthMJ2000Eq.X": sc_x,
        "Pandora.EarthMJ2000Eq.Y": sc_y,
        "Pandora.EarthMJ2000Eq.Z": sc_z,
        "Sun.EarthMJ2000Eq.X": sun_x,
        "Sun.EarthMJ2000Eq.Y": sun_y,
        "Sun.EarthMJ2000Eq.Z": sun_z,
        "Luna.EarthMJ2000Eq.X": moon_x,
        "Luna.EarthMJ2000Eq.Y": moon_y,
        "Luna.EarthMJ2000Eq.Z": moon_z,
        "Pandora.Earth.Latitude": lat,
        "Pandora.Earth.Longitude": lon,
    })
    path = tmp_path / "synthetic_gmat.csv"
    df.to_csv(path, index=False)
    return path


# ── Helpers ──────────────────────────────────────────────────────────────────

def _run_case(
    name: str,
    gmat_path: Path,
    tmp_path: Path,
    window_start: datetime,
    window_end: datetime,
    star_coord: SkyCoord,
    star_label: str,
    **config_overrides,
) -> pd.DataFrame:
    """Run one visibility case and return the DataFrame."""
    config = PandoraSchedulerConfig(
        window_start=window_start,
        window_end=window_end,
        gmat_ephemeris=gmat_path,
        output_dir=tmp_path,
        **config_overrides,
    )
    cadence = build_minute_cadence(window_start, window_end)
    ephemeris = interpolate_gmat_ephemeris(gmat_path, cadence)
    payload = _build_base_payload(ephemeris, cadence)
    df = _build_star_visibility(payload, star_coord, config)
    return df


def _summarise(name: str, df: pd.DataFrame, star_label: str) -> str:
    """Produce a summary block for one case."""
    total = len(df)
    vis_count = int((df["Visible"] > 0).sum())
    pct = 100.0 * vis_count / total if total else 0
    earth_min = df["Earth_Sep"].min()
    earth_max = df["Earth_Sep"].max()
    sun_min = df["Sun_Sep"].min()
    moon_min = df["Moon_Sep"].min()
    roll_valid = df["Roll_Deg"].notna().sum()
    st_pass = int(df["N_ST_Pass"].sum())

    lines = [
        f"  Case: {name}",
        f"  Target: {star_label}",
        f"  Minutes: {total}   Visible: {vis_count} ({pct:.1f}%)",
        f"  Earth sep range: [{earth_min:.1f}°, {earth_max:.1f}°]",
        f"  Sun sep min: {sun_min:.1f}°   Moon sep min: {moon_min:.1f}°",
        f"  Roll valid: {roll_valid}   ST pass sum: {st_pass}",
    ]
    return "\n".join(lines)


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning)

    window_start = datetime(2026, 6, 15, 0, 0, 0)
    window_end = datetime(2026, 6, 15, 6, 0, 0)

    # Targets at different sky positions
    # Sun is at RA≈0°, so targets with RA>90° clear the 91° Sun keepout
    targets = [
        ("Anti-Sun moderate dec", SkyCoord(ra=160 * u.deg, dec=30 * u.deg, frame="icrs")),
        ("Near Sun (blocked)", SkyCoord(ra=10 * u.deg, dec=5 * u.deg, frame="icrs")),
        ("Opposite Sun", SkyCoord(ra=180 * u.deg, dec=0 * u.deg, frame="icrs")),
        ("Anti-Sun high dec", SkyCoord(ra=150 * u.deg, dec=70 * u.deg, frame="icrs")),
    ]

    cases = [
        (
            "1. Baseline (no day/night, no ST)",
            dict(
                earth_avoidance_deg=86.0,  # old default
                earth_avoidance_day_deg=None,
                earth_avoidance_night_deg=None,
                st_required=0,
            ),
        ),
        (
            "2. Day/night Earth avoidance (110°/80°)",
            dict(
                earth_avoidance_deg=110.0,
                earth_avoidance_day_deg=110.0,
                earth_avoidance_night_deg=80.0,
                st_required=0,
            ),
        ),
        (
            "3. ST keepout + roll (Sun=45°, Moon=20°, Limb=10°)",
            dict(
                earth_avoidance_deg=86.0,
                earth_avoidance_day_deg=None,
                earth_avoidance_night_deg=None,
                st_sun_min_deg=45.0,
                st_moon_min_deg=20.0,
                st_earthlimb_min_deg=10.0,
                st_required=1,
                roll_step_deg=5.0,
                min_power_frac=0.5,
            ),
        ),
        (
            "4. Full: day/night + ST + roll",
            dict(
                earth_avoidance_deg=110.0,
                earth_avoidance_day_deg=110.0,
                earth_avoidance_night_deg=80.0,
                st_sun_min_deg=45.0,
                st_moon_min_deg=20.0,
                st_earthlimb_min_deg=10.0,
                st_required=1,
                roll_step_deg=5.0,
                min_power_frac=0.5,
            ),
        ),
    ]

    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        gmat_path = make_synthetic_gmat(tmp_path, window_start, window_end)

        all_summaries: list[str] = []
        all_csvs: list[Path] = []

        for case_name, overrides in cases:
            print(f"\n{'='*72}")
            print(f"  {case_name}")
            print(f"{'='*72}")

            for star_label, star_coord in targets:
                df = _run_case(
                    case_name,
                    gmat_path,
                    tmp_path,
                    window_start,
                    window_end,
                    star_coord,
                    star_label,
                    **overrides,
                )

                summary = _summarise(case_name, df, star_label)
                all_summaries.append(summary)
                print(f"\n{summary}")

                # Save CSV
                safe_name = case_name.split(".")[0].strip() + "_" + star_label.replace(" ", "_")
                csv_path = tmp_path / f"{safe_name}.csv"
                df.to_csv(csv_path, index=False)
                all_csvs.append(csv_path)

        # ── Cross-case comparison ────────────────────────────────────
        print(f"\n\n{'='*72}")
        print("  CROSS-CASE COMPARISON (Opposite Sun target)")
        print(f"{'='*72}")

        comparison_target = targets[2]  # "Opposite Sun"
        star_label, star_coord = comparison_target
        rows = []
        for case_name, overrides in cases:
            df = _run_case(
                case_name, gmat_path, tmp_path,
                window_start, window_end, star_coord, star_label,
                **overrides,
            )
            total = len(df)
            vis = int((df["Visible"] > 0).sum())
            rows.append({
                "Case": case_name,
                "Total min": total,
                "Visible min": vis,
                "Vis %": f"{100*vis/total:.1f}",
                "Roll valid": int(df["Roll_Deg"].notna().sum()),
                "ST pass sum": int(df["N_ST_Pass"].sum()),
            })

        cmp_df = pd.DataFrame(rows)
        print(f"\nTarget: {star_label}")
        print(cmp_df.to_string(index=False))

        # ── Day/night detail for one target ──────────────────────────
        print(f"\n\n{'='*72}")
        print("  DAY vs NIGHT DETAIL (Anti-Sun moderate dec)")
        print(f"{'='*72}")

        star_label_pole, star_coord_pole = targets[0]

        # Uniform 110°
        df_uniform = _run_case(
            "Uniform 110°", gmat_path, tmp_path,
            window_start, window_end, star_coord_pole, star_label_pole,
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=None,
            earth_avoidance_night_deg=None,
            st_required=0,
        )

        # Day 110° / Night 80°
        df_daynight = _run_case(
            "Day 110° / Night 80°", gmat_path, tmp_path,
            window_start, window_end, star_coord_pole, star_label_pole,
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            st_required=0,
        )

        uniform_vis = int((df_uniform["Visible"] > 0).sum())
        daynight_vis = int((df_daynight["Visible"] > 0).sum())
        gained = daynight_vis - uniform_vis

        print(f"\nTarget: {star_label_pole}")
        print(f"  Uniform 110°:       {uniform_vis} visible minutes")
        print(f"  Day 110° / Night 80°: {daynight_vis} visible minutes")
        print(f"  Gained from night:  {gained} minutes ({100*gained/len(df_uniform):.1f}% of window)")

        # ── Roll sweep detail ────────────────────────────────────────
        print(f"\n\n{'='*72}")
        print("  ROLL SWEEP DETAIL (Opposite Sun target)")
        print(f"{'='*72}")

        star_label_opp, star_coord_opp = targets[2]
        df_roll = _run_case(
            "Roll sweep", gmat_path, tmp_path,
            window_start, window_end, star_coord_opp, star_label_opp,
            earth_avoidance_deg=86.0,
            st_sun_min_deg=45.0,
            st_moon_min_deg=20.0,
            st_earthlimb_min_deg=10.0,
            st_required=1,
            roll_step_deg=5.0,
            min_power_frac=0.5,
        )

        vis_mask = df_roll["Visible"] > 0
        roll_vals = df_roll.loc[vis_mask, "Roll_Deg"]
        if len(roll_vals) > 0:
            print(f"\nTarget: {star_label_opp}")
            print(f"  Visible minutes: {vis_mask.sum()}")
            print(f"  Roll angles used: {sorted(roll_vals.dropna().unique())}")
            print(f"  Roll angle range: [{roll_vals.min():.0f}°, {roll_vals.max():.0f}°]")
        else:
            print(f"\nTarget: {star_label_opp}")
            print("  No visible minutes with ST constraints (target may be too close to Sun/Earth)")

        # Without ST for comparison
        df_no_st = _run_case(
            "No ST", gmat_path, tmp_path,
            window_start, window_end, star_coord_opp, star_label_opp,
            earth_avoidance_deg=86.0,
            st_required=0,
        )
        no_st_vis = int((df_no_st["Visible"] > 0).sum())
        st_vis = int(vis_mask.sum())
        lost = no_st_vis - st_vis
        print(f"\n  Without ST: {no_st_vis} visible minutes")
        print(f"  With ST:    {st_vis} visible minutes")
        print(f"  Lost to ST: {lost} minutes")

        print(f"\n\n{'='*72}")
        print("  ALL EXAMPLES COMPLETED SUCCESSFULLY")
        print(f"{'='*72}\n")


if __name__ == "__main__":
    main()
