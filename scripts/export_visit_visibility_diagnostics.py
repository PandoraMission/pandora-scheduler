#!/usr/bin/env python3
"""Export minute-by-minute visit visibility diagnostics to CSV.

This script writes one row per minute for a requested visit window and includes:
- boresight Sun / Moon / Earth-center separations and pass flags
- day/night Earth threshold applied that minute
- chosen roll angle and star-tracker pass counts
- per-tracker Sun / Moon / Earth-limb separations and pass flags
- a compact blocking-reasons summary column
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyCoord as CartesianSkyCoord
from astropy.time import Time

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.constraints import (
    ST1_BODY,
    ST2_BODY,
    _R_EARTH_KM,
    _normalise,
    compute_visibility_with_constraints,
    detect_orbit_boundaries,
    effective_earth_threshold,
    evaluate_star_tracker,
    fast_limb_deg,
    fast_sep_deg,
    orbit_slices_from_boundaries,
    star_tracker_eci,
    subsatellite_is_sunlit,
)
from pandorascheduler_rework.visibility.geometry import (
    build_minute_cadence,
    compute_saa_crossings,
    interpolate_gmat_ephemeris,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export minute-by-minute visit visibility diagnostics."
    )
    parser.add_argument(
        "--start",
        required=True,
        help="Visit start in 'YYYY-MM-DD HH:MM:SS' or ISO format.",
    )
    parser.add_argument(
        "--stop",
        required=True,
        help="Visit stop in 'YYYY-MM-DD HH:MM:SS' or ISO format.",
    )
    parser.add_argument(
        "--target-name",
        required=True,
        help="Target/planet name used for lookup and output naming.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Scheduler JSON config used to infer keepouts and GMAT path.",
    )
    parser.add_argument(
        "--schedule-csv",
        type=Path,
        help="Optional schedule CSV used to look up RA/DEC from the matching row.",
    )
    parser.add_argument(
        "--ra",
        type=float,
        help="Target right ascension in degrees. Required unless lookup succeeds.",
    )
    parser.add_argument(
        "--dec",
        type=float,
        help="Target declination in degrees. Required unless lookup succeeds.",
    )
    parser.add_argument(
        "--target-definition-json",
        type=Path,
        help="Optional explicit target-definition JSON for RA/DEC lookup.",
    )
    parser.add_argument(
        "--gmat",
        type=Path,
        help="Optional explicit GMAT ephemeris path. Overrides --config.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional output CSV path. Defaults beside the config file.",
    )
    return parser.parse_args()


def _parse_datetime(value: str) -> datetime:
    return datetime.fromisoformat(value.strip().replace("T", " ").replace("Z", ""))


def _json_value(data: dict[str, Any], key: str, default: Any = None) -> Any:
    if key in data:
        return data[key]
    nested = data.get("extra_inputs", {})
    if isinstance(nested, dict) and key in nested:
        return nested[key]
    return default


def _load_config(config_path: Path) -> dict[str, Any]:
    with config_path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _lookup_radec_from_schedule(
    schedule_csv: Path,
    target_name: str,
    start: datetime,
) -> tuple[float, float] | None:
    schedule = pd.read_csv(schedule_csv)
    start_text = start.strftime("%Y-%m-%d %H:%M:%S")
    mask = (
        schedule.get("Target", pd.Series(dtype=object)).astype(str).str.strip().eq(target_name)
        & schedule.get("Observation Start", pd.Series(dtype=object)).astype(str).str.strip().eq(start_text)
    )
    match = schedule.loc[mask]
    if match.empty:
        return None
    row = match.iloc[0]
    return float(row["RA"]), float(row["DEC"])


def _lookup_radec_from_target_definition(
    target_name: str,
    target_definition_json: Path | None,
    config_data: dict[str, Any],
) -> tuple[float, float] | None:
    candidate = target_definition_json
    if candidate is None:
        base = _json_value(config_data, "target_definition_base")
        if base:
            candidate = (
                Path(str(base)).expanduser().resolve()
                / "exoplanet"
                / f"{target_name}_target_definition.json"
            )
    if candidate is None or not candidate.exists():
        return None
    with candidate.open(encoding="utf-8") as handle:
        data = json.load(handle)
    return float(data["RA"]), float(data["DEC"])


def _resolve_radec(args: argparse.Namespace, config_data: dict[str, Any]) -> tuple[float, float]:
    if args.ra is not None and args.dec is not None:
        return float(args.ra), float(args.dec)

    start = _parse_datetime(args.start)
    if args.schedule_csv is not None:
        from_schedule = _lookup_radec_from_schedule(
            args.schedule_csv.resolve(),
            args.target_name,
            start,
        )
        if from_schedule is not None:
            return from_schedule

    from_target_definition = _lookup_radec_from_target_definition(
        args.target_name,
        args.target_definition_json.resolve() if args.target_definition_json else None,
        config_data,
    )
    if from_target_definition is not None:
        return from_target_definition

    raise ValueError(
        "Unable to resolve RA/DEC. Provide --ra/--dec, or use --schedule-csv, "
        "--target-definition-json, or a config whose extra_inputs.target_definition_base "
        "contains the target definition."
    )


def _default_output(args: argparse.Namespace) -> Path:
    safe_name = args.target_name.replace("/", "_")
    stamp = args.start.replace(":", "").replace(" ", "_")
    return args.config.resolve().parent / f"{safe_name}_{stamp}_visibility_diagnostics.csv"


def _build_payload_frame(
    target_unit: np.ndarray,
    roll_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    north = np.array([0.0, 0.0, 1.0])
    north_tiled = np.broadcast_to(north, target_unit.shape).copy()
    dot_nz = np.einsum("ij,ij->i", north_tiled, target_unit)[:, None]
    north_proj = north_tiled - target_unit * dot_nz
    proj_norm = np.linalg.norm(north_proj, axis=1, keepdims=True)

    degenerate = (proj_norm < 1e-10).squeeze()
    if np.any(degenerate):
        fallback = np.array([1.0, 0.0, 0.0])
        fallback_tiled = np.broadcast_to(fallback, target_unit.shape).copy()
        dot_fz = np.einsum("ij,ij->i", fallback_tiled, target_unit)[:, None]
        fallback_proj = fallback_tiled - target_unit * dot_fz
        north_proj[degenerate] = fallback_proj[degenerate]
        proj_norm = np.linalg.norm(north_proj, axis=1, keepdims=True)

    proj_norm = np.where(proj_norm == 0, 1.0, proj_norm)
    x_ref = north_proj / proj_norm
    y_ref = np.cross(target_unit, x_ref)
    y_ref_norm = np.linalg.norm(y_ref, axis=1, keepdims=True)
    y_ref_norm = np.where(y_ref_norm == 0, 1.0, y_ref_norm)
    y_ref = y_ref / y_ref_norm

    finite_roll = np.isfinite(roll_deg)
    x_pay = np.full_like(target_unit, np.nan, dtype=float)
    y_pay = np.full_like(target_unit, np.nan, dtype=float)
    z_pay = np.full_like(target_unit, np.nan, dtype=float)
    if np.any(finite_roll):
        roll_rad = np.deg2rad(roll_deg[finite_roll])[:, None]
        cos_r = np.cos(roll_rad)
        sin_r = np.sin(roll_rad)
        x_pay[finite_roll] = cos_r * x_ref[finite_roll] + sin_r * y_ref[finite_roll]
        y_pay[finite_roll] = -sin_r * x_ref[finite_roll] + cos_r * y_ref[finite_roll]
        z_pay[finite_roll] = target_unit[finite_roll]
    return x_pay, y_pay, z_pay


def _tracker_columns(
    x_pay: np.ndarray,
    y_pay: np.ndarray,
    z_pay: np.ndarray,
    sun_unit: np.ndarray,
    moon_unit: np.ndarray,
    zenith_unit: np.ndarray,
    limb_angle_rad: np.ndarray,
    config: PandoraSchedulerConfig,
    st_body: np.ndarray,
    prefix: str,
) -> dict[str, np.ndarray]:
    valid = np.isfinite(x_pay).all(axis=1)
    sun_sep = np.full(len(x_pay), np.nan)
    moon_sep = np.full(len(x_pay), np.nan)
    earthlimb_sep = np.full(len(x_pay), np.nan)
    ok = np.zeros(len(x_pay), dtype=bool)

    if np.any(valid):
        st_eci = star_tracker_eci(x_pay[valid], y_pay[valid], z_pay[valid], st_body)
        sun_sep[valid] = fast_sep_deg(st_eci, sun_unit[valid])
        moon_sep[valid] = fast_sep_deg(st_eci, moon_unit[valid])
        earthlimb_sep[valid] = fast_limb_deg(st_eci, zenith_unit[valid], limb_angle_rad[valid])
        ok[valid] = evaluate_star_tracker(
            st_eci,
            sun_unit[valid],
            moon_unit[valid],
            zenith_unit[valid],
            limb_angle_rad[valid],
            config.st_sun_min_deg,
            config.st_moon_min_deg,
            config.st_earthlimb_min_deg,
        )

    return {
        f"{prefix}_sun_sep_deg": np.round(sun_sep, 6),
        f"{prefix}_sun_min_deg": np.full(len(x_pay), config.st_sun_min_deg, dtype=float),
        f"{prefix}_sun_ok": sun_sep >= config.st_sun_min_deg,
        f"{prefix}_moon_sep_deg": np.round(moon_sep, 6),
        f"{prefix}_moon_min_deg": np.full(len(x_pay), config.st_moon_min_deg, dtype=float),
        f"{prefix}_moon_ok": moon_sep >= config.st_moon_min_deg,
        f"{prefix}_earthlimb_sep_deg": np.round(earthlimb_sep, 6),
        f"{prefix}_earthlimb_min_deg": np.full(len(x_pay), config.st_earthlimb_min_deg, dtype=float),
        f"{prefix}_earthlimb_ok": earthlimb_sep >= config.st_earthlimb_min_deg,
        f"{prefix}_ok": ok,
    }


def _blocking_reasons(row: pd.Series) -> str:
    reasons: list[str] = []
    if not bool(row["sun_ok"]):
        reasons.append("boresight_sun")
    if not bool(row["moon_ok"]):
        reasons.append("boresight_moon")
    if not bool(row["earth_ok"]):
        reasons.append("boresight_earth")
    if bool(row["boresight_ok"]) and not bool(row["st_requirement_ok"]):
        for key in (
            "st1_sun_ok",
            "st1_moon_ok",
            "st1_earthlimb_ok",
            "st2_sun_ok",
            "st2_moon_ok",
            "st2_earthlimb_ok",
        ):
            if not bool(row[key]):
                reasons.append(key.replace("_ok", ""))
    return ",".join(reasons)


def main() -> int:
    args = _parse_args()
    config_data = _load_config(args.config.resolve())

    start = _parse_datetime(args.start)
    stop = _parse_datetime(args.stop)
    if stop <= start:
        raise ValueError("--stop must be later than --start")

    ra_deg, dec_deg = _resolve_radec(args, config_data)
    gmat_path = args.gmat.resolve() if args.gmat else Path(
        str(_json_value(config_data, "visibility_gmat"))
    ).expanduser().resolve()
    if not gmat_path.exists():
        raise FileNotFoundError(f"GMAT ephemeris not found: {gmat_path}")

    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=stop,
        sun_avoidance_deg=float(_json_value(config_data, "sun_avoidance_deg", 91.0)),
        moon_avoidance_deg=float(_json_value(config_data, "moon_avoidance_deg", 25.0)),
        earth_avoidance_deg=float(_json_value(config_data, "earth_avoidance_day_deg", 110.0)),
        earth_avoidance_day_deg=_json_value(config_data, "earth_avoidance_day_deg"),
        earth_avoidance_night_deg=_json_value(config_data, "earth_avoidance_night_deg"),
        daynight_mode=str(_json_value(config_data, "daynight_mode", "subsatellite")),
        twilight_margin_deg=float(_json_value(config_data, "twilight_margin_deg", 0.0)),
        st_sun_min_deg=float(_json_value(config_data, "st_sun_min_deg", 0.0)),
        st_moon_min_deg=float(_json_value(config_data, "st_moon_min_deg", 0.0)),
        st_earthlimb_min_deg=float(_json_value(config_data, "st_earthlimb_min_deg", 0.0)),
        st_required=int(_json_value(config_data, "st_required", 1)),
        roll_step_deg=float(_json_value(config_data, "roll_step_deg", 2.0)),
        min_power_frac=float(_json_value(config_data, "min_power_frac", 0.7)),
    )

    cadence = build_minute_cadence(start, stop)
    ephemeris = interpolate_gmat_ephemeris(gmat_path, cadence)
    mjd_array = np.asarray(cadence.mjd_utc, dtype=float)
    time_utc = pd.to_datetime(Time(mjd_array, format="mjd", scale="utc").datetime).round("s")

    earth_pc_xyz = ephemeris.earth_pc
    sun_pc_xyz = ephemeris.sun_pc
    moon_pc_xyz = ephemeris.moon_pc
    nadir_unit = _normalise(earth_pc_xyz)
    zenith_unit = -nadir_unit
    sun_unit = _normalise(sun_pc_xyz)
    moon_unit = _normalise(moon_pc_xyz)
    observer_dist_km = np.linalg.norm(earth_pc_xyz, axis=1)
    limb_angle_rad = np.arccos(np.clip(_R_EARTH_KM / observer_dist_km, -1.0, 1.0))
    orbit_boundaries = detect_orbit_boundaries(ephemeris.spacecraft_lat_deg)
    orbit_slices = orbit_slices_from_boundaries(orbit_boundaries, len(mjd_array))

    star_coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    target_cart = star_coord.icrs.cartesian
    target_unit_1 = np.array(
        [target_cart.x.value, target_cart.y.value, target_cart.z.value],
        dtype=float,
    )
    target_unit_1 = target_unit_1 / np.linalg.norm(target_unit_1)
    target_unit = np.broadcast_to(target_unit_1, (len(mjd_array), 3)).copy()

    earth_pc_sc = CartesianSkyCoord(
        ephemeris.earth_pc,
        unit=u.km,
        representation_type="cartesian",
    )
    earth_sep_deg = earth_pc_sc.separation(star_coord).deg

    results = compute_visibility_with_constraints(
        target_unit=target_unit,
        nadir_unit=nadir_unit,
        sun_unit=sun_unit,
        moon_unit=moon_unit,
        observer_dist_km=observer_dist_km,
        zenith_unit=zenith_unit,
        limb_angle_rad=limb_angle_rad,
        orbit_slices=orbit_slices,
        earth_center_sep_deg=earth_sep_deg,
        config=config,
    )

    sun_sep_deg = fast_sep_deg(target_unit, sun_unit)
    moon_sep_deg = fast_sep_deg(target_unit, moon_unit)
    sun_ok = sun_sep_deg > config.sun_avoidance_deg
    moon_ok = moon_sep_deg > config.moon_avoidance_deg
    earth_threshold = effective_earth_threshold(
        target_unit,
        nadir_unit,
        sun_unit,
        config.earth_avoidance_day_deg,
        config.earth_avoidance_night_deg,
        config.earth_avoidance_deg,
        limb_angle_rad=limb_angle_rad,
        twilight_margin_deg=config.twilight_margin_deg,
        daynight_mode=config.daynight_mode,
    )
    earth_ok = earth_sep_deg > earth_threshold
    boresight_ok = sun_ok & moon_ok & earth_ok

    roll_deg = np.asarray(results["roll_deg"], dtype=float)
    x_pay, y_pay, z_pay = _build_payload_frame(target_unit, roll_deg)
    st1 = _tracker_columns(
        x_pay,
        y_pay,
        z_pay,
        sun_unit,
        moon_unit,
        zenith_unit,
        limb_angle_rad,
        config,
        ST1_BODY,
        "st1",
    )
    st2 = _tracker_columns(
        x_pay,
        y_pay,
        z_pay,
        sun_unit,
        moon_unit,
        zenith_unit,
        limb_angle_rad,
        config,
        ST2_BODY,
        "st2",
    )

    out = pd.DataFrame(
        {
            "time_utc": time_utc,
            "visible": np.asarray(results["visible"], dtype=bool),
            "saa_crossing": compute_saa_crossings(
                ephemeris.spacecraft_lat_deg,
                ephemeris.spacecraft_lon_deg,
            ).astype(bool),
            "sunlit_subsatellite": subsatellite_is_sunlit(nadir_unit, sun_unit),
            "sun_sep_deg": np.round(sun_sep_deg, 6),
            "sun_min_deg": config.sun_avoidance_deg,
            "sun_ok": sun_ok,
            "moon_sep_deg": np.round(moon_sep_deg, 6),
            "moon_min_deg": config.moon_avoidance_deg,
            "moon_ok": moon_ok,
            "earth_sep_deg": np.round(earth_sep_deg, 6),
            "earth_threshold_deg": np.round(np.asarray(earth_threshold, dtype=float), 6),
            "earth_ok": earth_ok,
            "boresight_ok": boresight_ok,
            "roll_deg": np.round(roll_deg, 2),
            "n_st_pass": np.asarray(results["n_st_pass"], dtype=int),
            "st_required": config.st_required,
            "st_requirement_ok": np.asarray(results["n_st_pass"], dtype=int) >= config.st_required,
        }
    )
    for tracker_cols in (st1, st2):
        for key, value in tracker_cols.items():
            out[key] = value

    out["blocking_reasons"] = out.apply(_blocking_reasons, axis=1)

    output_path = args.output.resolve() if args.output else _default_output(args)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, index=False)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
