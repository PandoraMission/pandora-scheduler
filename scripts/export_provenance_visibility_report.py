#!/usr/bin/env python3
"""Export a minute-by-minute visibility report for a provenance CSV.

The report is intentionally similar to the existing external
``*-violations.txt`` artifacts:
- one line per minute inside each provenance interval
- visit / sequence / target identifiers
- day/night branch
- boresight Sun / Moon / Earth-limb angles
- ST1 / ST2 Sun / Moon / Earth-limb angles
- violation marker when the scheduler-style keepouts fail

Important note:
- This script uses the scheduler's own Earth-center day/night keepout logic.
- It also prints the Earth-limb angle for reference.
- It does *not* attempt to reproduce an external validator's separate
  boresight-limb threshold model.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
import types
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

R_EARTH_KM = 6371.0
ST1_BODY = np.array([0.6804, -0.7071, -0.1923], dtype=float)
ST1_BODY = ST1_BODY / np.linalg.norm(ST1_BODY)
ST2_BODY = np.array([0.6804, 0.7071, -0.1923], dtype=float)
ST2_BODY = ST2_BODY / np.linalg.norm(ST2_BODY)


@dataclass(frozen=True)
class Thresholds:
    sun_min: float
    moon_min: float
    earth_day_min: float
    earth_night_min: float
    st_sun_min: float
    st_moon_min: float
    st_earthlimb_min: float
    st_required: int
    roll_step_deg: float
    min_power_frac: float


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a minute-by-minute visibility report from provenance CSV."
    )
    parser.add_argument(
        "provenance_csv",
        type=Path,
        help="Path to Pandora_science_calendar_sequence_provenance.csv",
    )
    parser.add_argument(
        "--run-manifest",
        type=Path,
        help="Path to run_config_manifest.json. Defaults next to the provenance CSV.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Path to run data directory. Defaults from the provenance parent.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output text path. Defaults next to the provenance CSV.",
    )
    return parser.parse_args()


def _load_geometry_module(repo_root: Path):
    path = repo_root / "src" / "pandorascheduler_rework" / "visibility" / "geometry.py"
    spec = importlib.util.spec_from_file_location("ps_geom_report", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["ps_geom_report"] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_constraints_module(repo_root: Path):
    pkg = sys.modules.setdefault(
        "pandorascheduler_rework", types.ModuleType("pandorascheduler_rework")
    )
    config_mod = sys.modules.setdefault(
        "pandorascheduler_rework.config",
        types.ModuleType("pandorascheduler_rework.config"),
    )
    if not hasattr(config_mod, "PandoraSchedulerConfig"):
        config_mod.PandoraSchedulerConfig = object
    setattr(pkg, "config", config_mod)

    path = repo_root / "src" / "pandorascheduler_rework" / "visibility" / "constraints.py"
    spec = importlib.util.spec_from_file_location("ps_constraints_report", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["ps_constraints_report"] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _sep_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    dot = np.einsum("ij,ij->i", a, b)
    return np.rad2deg(np.arccos(np.clip(dot, -1.0, 1.0)))


def _fast_limb_deg(
    target_unit: np.ndarray,
    zenith_unit: np.ndarray,
    limb_angle_rad: np.ndarray,
) -> np.ndarray:
    dot = np.einsum("ij,ij->i", target_unit, zenith_unit)
    elev = np.arcsin(np.clip(dot, -1.0, 1.0))
    return np.rad2deg(elev + limb_angle_rad)


def _sun_constrained_attitude(
    target_unit: np.ndarray,
    sun_unit: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z_payload = target_unit
    y_payload = np.cross(sun_unit, z_payload)
    y_norm = np.linalg.norm(y_payload, axis=1, keepdims=True)
    y_norm = np.where(y_norm == 0, 1.0, y_norm)
    y_payload = y_payload / y_norm
    x_payload = np.cross(y_payload, z_payload)
    x_norm = np.linalg.norm(x_payload, axis=1, keepdims=True)
    x_norm = np.where(x_norm == 0, 1.0, x_norm)
    x_payload = x_payload / x_norm
    return x_payload, y_payload, z_payload


def _star_tracker_eci(
    x_payload: np.ndarray,
    y_payload: np.ndarray,
    z_payload: np.ndarray,
    st_body: np.ndarray,
) -> np.ndarray:
    st_eci = (
        x_payload * st_body[0]
        + y_payload * st_body[1]
        + z_payload * st_body[2]
    )
    return _unit(st_eci)


def _radec_unit(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cos_dec = np.cos(dec)
    return np.array(
        [[cos_dec * np.cos(ra), cos_dec * np.sin(ra), np.sin(dec)]],
        dtype=float,
    )


def _parse_utc(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", ""))


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _json_value(data: dict[str, Any], key: str, default: Any = None) -> Any:
    if key in data:
        return data[key]
    nested = data.get("extra_inputs", {})
    if isinstance(nested, dict) and key in nested:
        return nested[key]
    return default


def _load_thresholds(cfg: dict[str, Any]) -> Thresholds:
    return Thresholds(
        sun_min=float(_json_value(cfg, "sun_avoidance_deg", 0.0)),
        moon_min=float(_json_value(cfg, "moon_avoidance_deg", 0.0)),
        earth_day_min=float(_json_value(cfg, "earth_avoidance_day_deg", 0.0)),
        earth_night_min=float(_json_value(cfg, "earth_avoidance_night_deg", 0.0)),
        st_sun_min=float(_json_value(cfg, "st_sun_min_deg", 0.0)),
        st_moon_min=float(_json_value(cfg, "st_moon_min_deg", 0.0)),
        st_earthlimb_min=float(_json_value(cfg, "st_earthlimb_min_deg", 0.0)),
        st_required=int(_json_value(cfg, "st_required", 0)),
        roll_step_deg=float(_json_value(cfg, "roll_step_deg", 5.0)),
        min_power_frac=float(_json_value(cfg, "min_power_frac", 0.0)),
    )


def _build_target_lookup(data_dir: Path) -> dict[str, tuple[float, float]]:
    candidates = [
        data_dir / "all_targets.csv",
        data_dir / "exoplanet_targets.csv",
        data_dir / "auxiliary-standard_targets.csv",
        data_dir / "monitoring-standard_targets.csv",
        data_dir / "occultation-standard_targets.csv",
    ]
    lookup: dict[str, tuple[float, float]] = {}
    for path in candidates:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if "RA" not in df.columns or "DEC" not in df.columns:
            continue
        for _, row in df.iterrows():
            ra = row.get("RA")
            dec = row.get("DEC")
            if pd.isna(ra) or pd.isna(dec):
                continue
            for key in (
                str(row.get("Planet Name", "")).strip(),
                str(row.get("Star Name", "")).strip(),
            ):
                if key:
                    lookup.setdefault(key, (float(ra), float(dec)))
    return lookup


def _resolve_target_radec(
    target_name: str,
    target_lookup: dict[str, tuple[float, float]],
) -> tuple[float, float]:
    direct = target_lookup.get(target_name)
    if direct is not None:
        return direct
    stripped = target_name.removesuffix(" STD").strip()
    direct = target_lookup.get(stripped)
    if direct is not None:
        return direct
    raise KeyError(f"Could not resolve RA/DEC for target {target_name!r}")


def _default_output(provenance_csv: Path) -> Path:
    return provenance_csv.with_name(f"{provenance_csv.stem}_visibility_report.txt")


def _choose_roll_per_orbit(
    constraints,
    target_unit: np.ndarray,
    zenith_unit: np.ndarray,
    sun_unit: np.ndarray,
    moon_unit: np.ndarray,
    limb_angle_rad: np.ndarray,
    orbit_slices: list[slice],
    boresight_visible: np.ndarray,
    thresholds: Thresholds,
) -> np.ndarray:
    chosen_roll = np.full(target_unit.shape[0], np.nan)
    roll_angles = np.arange(0, 360, thresholds.roll_step_deg)

    for orb_slice in orbit_slices:
        bvis = boresight_visible[orb_slice]
        if not np.any(bvis):
            continue

        tgt = target_unit[orb_slice]
        zen = zenith_unit[orb_slice]
        sun = sun_unit[orb_slice]
        moon = moon_unit[orb_slice]
        limb = limb_angle_rad[orb_slice]

        best_roll_deg = np.nan
        best_vis_count = -1
        best_power_mean = -1.0

        for roll_deg in roll_angles:
            x_pay, y_pay, z_pay = constraints.fixed_roll_attitude(
                tgt, np.deg2rad(float(roll_deg))
            )
            power = constraints.solar_power_fraction(y_pay, sun)
            mean_power = float(np.mean(power[bvis]))
            if mean_power < thresholds.min_power_frac:
                continue

            st1_eci = constraints.star_tracker_eci(
                x_pay, y_pay, z_pay, constraints.ST1_BODY
            )
            st2_eci = constraints.star_tracker_eci(
                x_pay, y_pay, z_pay, constraints.ST2_BODY
            )
            st1_ok = constraints.evaluate_star_tracker(
                st1_eci,
                sun,
                moon,
                zen,
                limb,
                thresholds.st_sun_min,
                thresholds.st_moon_min,
                thresholds.st_earthlimb_min,
            )
            st2_ok = constraints.evaluate_star_tracker(
                st2_eci,
                sun,
                moon,
                zen,
                limb,
                thresholds.st_sun_min,
                thresholds.st_moon_min,
                thresholds.st_earthlimb_min,
            )
            combined = st1_ok | st2_ok if thresholds.st_required == 1 else st1_ok & st2_ok
            vis_count = int(np.sum(bvis & combined))

            if vis_count > best_vis_count or (
                vis_count == best_vis_count and mean_power > best_power_mean
            ):
                best_roll_deg = float(roll_deg)
                best_vis_count = vis_count
                best_power_mean = mean_power

        if np.isfinite(best_roll_deg):
            best_roll_deg = ((best_roll_deg + 180.0) % 360.0) - 180.0
            chosen_roll[orb_slice] = best_roll_deg

    return chosen_roll


def main() -> int:
    args = _parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    geometry = _load_geometry_module(repo_root)
    constraints = _load_constraints_module(repo_root)

    provenance_csv = args.provenance_csv.resolve()
    run_dir = provenance_csv.parent
    manifest_path = args.run_manifest.resolve() if args.run_manifest else run_dir / "run_config_manifest.json"
    data_dir = args.data_dir.resolve() if args.data_dir else next(
        (p for p in run_dir.iterdir() if p.is_dir() and p.name.startswith("data_")),
        None,
    )
    if data_dir is None:
        raise ValueError("Could not infer data directory; provide --data-dir")

    manifest = _load_json(manifest_path)
    json_cfg = manifest.get("json_config", {})
    source_config_path = manifest.get("source_config_path")
    source_cfg = _load_json(Path(source_config_path)) if source_config_path else {}
    visibility_gmat = _json_value(json_cfg, "visibility_gmat", _json_value(source_cfg, "visibility_gmat"))
    if not visibility_gmat:
        raise ValueError("Could not determine visibility_gmat from run manifest/config")
    gmat_path = Path(str(visibility_gmat)).expanduser().resolve()

    thresholds = _load_thresholds(json_cfg if json_cfg else source_cfg)
    target_lookup = _build_target_lookup(data_dir)
    prov_df = pd.read_csv(provenance_csv, dtype={"visit_id": str, "sequence_id": str})
    if prov_df.empty:
        raise ValueError(f"No rows in provenance CSV: {provenance_csv}")

    start_all = min(_parse_utc(v) for v in prov_df["start_utc"])
    stop_all = max(_parse_utc(v) for v in prov_df["stop_utc"])
    cadence = geometry.build_minute_cadence(start_all, stop_all)
    eph = geometry.interpolate_gmat_ephemeris(gmat_path, cadence)

    ephemeris_index = pd.DatetimeIndex(
        pd.date_range(start_all, stop_all, freq="min").tz_localize(None)
    )
    eph_lookup = {ts.to_pydatetime(): idx for idx, ts in enumerate(ephemeris_index)}
    orbit_boundaries = constraints.detect_orbit_boundaries(eph.spacecraft_lat_deg)
    orbit_slices = constraints.orbit_slices_from_boundaries(
        orbit_boundaries, len(ephemeris_index)
    )

    out_path = args.output.resolve() if args.output else _default_output(provenance_csv)
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(
            f"Calendar: {start_all.isoformat()}.000 → {stop_all.isoformat()}.000\n"
        )
        handle.write(f"  Source: {provenance_csv.name}\n")
        handle.write(f"  Rows: {len(prov_df)}\n\n")
        handle.write(
            "Keep-out thresholds: "
            + str(
                {
                    "bs_sun_min": thresholds.sun_min,
                    "bs_moon_min": thresholds.moon_min,
                    "bs_earth_day_center_min": thresholds.earth_day_min,
                    "bs_earth_night_center_min": thresholds.earth_night_min,
                    "st_sun_min": thresholds.st_sun_min,
                    "st_moon_min": thresholds.st_moon_min,
                    "st_earthlimb_min": thresholds.st_earthlimb_min,
                }
            )
            + "\n"
        )
        handle.write(f"Star trackers required: {thresholds.st_required}\n\n")
        handle.write(
            f"{'Visit':>8} {'Seq':>4} {'Target':<28} {'D/N':>3} "
            f"{'BS Sun':>7} {'BS Moon':>7} {'BS Ctr':>7} {'BS Limb':>8} "
            f"{'ST1 Sun':>8} {'ST1 Moon':>8} {'ST1 Limb':>8} "
            f"{'ST2 Sun':>8} {'ST2 Moon':>8} {'ST2 Limb':>8}\n"
        )
        handle.write("-" * 140 + "\n")

        for _, row in prov_df.iterrows():
            visit_id = str(row["visit_id"])
            seq_id = str(row["sequence_id"])
            target_name = str(row["target"]).strip()
            ra_deg, dec_deg = _resolve_target_radec(target_name, target_lookup)
            start = _parse_utc(str(row["start_utc"]))
            stop = _parse_utc(str(row["stop_utc"]))
            n_minutes = int(round((stop - start).total_seconds() / 60.0))
            minute_index = pd.date_range(start, periods=n_minutes, freq="min").tz_localize(None)
            target_unit = np.repeat(_radec_unit(ra_deg, dec_deg), len(minute_index), axis=0)

            indices = np.array([eph_lookup[ts.to_pydatetime()] for ts in minute_index], dtype=int)
            earth_pc = eph.earth_pc[indices]
            sun_pc = eph.sun_pc[indices]
            moon_pc = eph.moon_pc[indices]
            nadir_unit = _unit(earth_pc)
            zenith_unit = -nadir_unit
            sun_unit = _unit(sun_pc)
            moon_unit = _unit(moon_pc)

            obs_dist = np.linalg.norm(earth_pc, axis=1)
            limb_angle_rad = np.arccos(np.clip(R_EARTH_KM / obs_dist, -1.0, 1.0))

            bs_sun = _sep_deg(target_unit, sun_unit)
            bs_moon = _sep_deg(target_unit, moon_unit)
            bs_ctr = _sep_deg(target_unit, nadir_unit)
            bs_limb = _fast_limb_deg(target_unit, zenith_unit, limb_angle_rad)

            day = np.einsum("ij,ij->i", zenith_unit, sun_unit) > 0.0
            earth_threshold = np.where(day, thresholds.earth_day_min, thresholds.earth_night_min)
            boresight_ok = (
                (bs_sun >= thresholds.sun_min)
                & (bs_moon >= thresholds.moon_min)
                & (bs_ctr > earth_threshold)
            )

            chosen_roll = _choose_roll_per_orbit(
                constraints,
                target_unit,
                zenith_unit,
                sun_unit,
                moon_unit,
                limb_angle_rad,
                orbit_slices,
                boresight_ok,
                thresholds,
            )

            st1_eci = np.full_like(target_unit, np.nan)
            st2_eci = np.full_like(target_unit, np.nan)
            finite_roll = np.isfinite(chosen_roll)
            if np.any(finite_roll):
                for roll_deg in np.unique(chosen_roll[finite_roll]):
                    mask = finite_roll & (chosen_roll == roll_deg)
                    x_pay, y_pay, z_pay = constraints.fixed_roll_attitude(
                        target_unit[mask], np.deg2rad(float(roll_deg))
                    )
                    st1_eci[mask] = constraints.star_tracker_eci(
                        x_pay, y_pay, z_pay, constraints.ST1_BODY
                    )
                    st2_eci[mask] = constraints.star_tracker_eci(
                        x_pay, y_pay, z_pay, constraints.ST2_BODY
                    )

            if np.any(~finite_roll):
                x_pay, y_pay, z_pay = _sun_constrained_attitude(
                    target_unit[~finite_roll], sun_unit[~finite_roll]
                )
                st1_eci[~finite_roll] = _star_tracker_eci(
                    x_pay, y_pay, z_pay, ST1_BODY
                )
                st2_eci[~finite_roll] = _star_tracker_eci(
                    x_pay, y_pay, z_pay, ST2_BODY
                )

            st1_sun = _sep_deg(st1_eci, sun_unit)
            st1_moon = _sep_deg(st1_eci, moon_unit)
            st1_limb = _fast_limb_deg(st1_eci, zenith_unit, limb_angle_rad)
            st2_sun = _sep_deg(st2_eci, sun_unit)
            st2_moon = _sep_deg(st2_eci, moon_unit)
            st2_limb = _fast_limb_deg(st2_eci, zenith_unit, limb_angle_rad)
            st1_ok = (
                (st1_sun >= thresholds.st_sun_min)
                & (st1_moon >= thresholds.st_moon_min)
                & (st1_limb >= thresholds.st_earthlimb_min)
            )
            st2_ok = (
                (st2_sun >= thresholds.st_sun_min)
                & (st2_moon >= thresholds.st_moon_min)
                & (st2_limb >= thresholds.st_earthlimb_min)
            )
            if thresholds.st_required <= 0:
                st_ok = np.ones(len(minute_index), dtype=bool)
            elif thresholds.st_required == 1:
                st_ok = st1_ok | st2_ok
            else:
                st_ok = st1_ok & st2_ok

            overall_ok = boresight_ok & st_ok

            for idx, timestamp in enumerate(minute_index):
                line = (
                    f"{visit_id:>8} {seq_id:>4} {target_name:<28} "
                    f"{('D' if bool(day[idx]) else 'N'):>3} "
                    f"{bs_sun[idx]:7.1f} {bs_moon[idx]:7.1f} {bs_ctr[idx]:7.1f} {bs_limb[idx]:8.1f} "
                    f"{st1_sun[idx]:8.1f} {st1_moon[idx]:8.1f} {st1_limb[idx]:8.1f} "
                    f"{st2_sun[idx]:8.1f} {st2_moon[idx]:8.1f} {st2_limb[idx]:8.1f}"
                )
                if not bool(overall_ok[idx]):
                    line += " *** VIOLATION"
                handle.write(line + "\n")

    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
