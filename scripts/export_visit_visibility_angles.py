#!/usr/bin/env python3
"""Export per-minute visibility geometry for a scheduled visit to CSV.

This script computes the same boresight geometry we inspected manually:
- target-Sun separation
- target-Moon separation
- target-Earth-center separation
- target-Earth-limb angle
- subsatellite day/night branch
- scheduler Earth-center threshold used that minute
- whether the scheduler Earth-center keepout passes that minute

It intentionally avoids importing the full ``pandorascheduler_rework`` package,
so it can run even in environments where parquet-related extras are absent.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

R_EARTH_KM = 6371.0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export per-minute Sun/Moon/Earth geometry for a visit."
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
        help="Target label for the output and optional schedule lookup.",
    )
    parser.add_argument(
        "--ra",
        type=float,
        help="Target right ascension in degrees. Required unless --schedule-csv is provided.",
    )
    parser.add_argument(
        "--dec",
        type=float,
        help="Target declination in degrees. Required unless --schedule-csv is provided.",
    )
    parser.add_argument(
        "--schedule-csv",
        type=Path,
        help="Optional schedule CSV used to look up RA/DEC from the matching row.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Optional scheduler JSON config used to infer GMAT path and Earth thresholds.",
    )
    parser.add_argument(
        "--gmat",
        type=Path,
        help="Path to GMAT ephemeris file. Overrides values inferred from --config.",
    )
    parser.add_argument(
        "--earth-day-deg",
        type=float,
        help="Scheduler Earth day threshold in degrees. Overrides --config.",
    )
    parser.add_argument(
        "--earth-night-deg",
        type=float,
        help="Scheduler Earth night threshold in degrees. Overrides --config.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output CSV path. Defaults next to the schedule/config context.",
    )
    return parser.parse_args()


def _parse_datetime(value: str) -> datetime:
    value = value.strip().replace("T", " ").replace("Z", "")
    return datetime.fromisoformat(value)


def _load_geometry_module(repo_root: Path):
    path = repo_root / "src" / "pandorascheduler_rework" / "visibility" / "geometry.py"
    spec = importlib.util.spec_from_file_location("ps_geom_export", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["ps_geom_export"] = module
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


def _radec_unit(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cos_dec = np.cos(dec)
    return np.array(
        [[cos_dec * np.cos(ra), cos_dec * np.sin(ra), np.sin(dec)]],
        dtype=float,
    )


def _json_value(data: dict[str, Any], key: str, default: Any = None) -> Any:
    if key in data:
        return data[key]
    nested = data.get("extra_inputs", {})
    if isinstance(nested, dict) and key in nested:
        return nested[key]
    return default


def _load_config_defaults(config_path: Path | None) -> dict[str, Any]:
    if config_path is None:
        return {}
    with config_path.open(encoding="utf-8") as handle:
        raw = json.load(handle)
    return {
        "gmat": _json_value(raw, "visibility_gmat"),
        "earth_day_deg": _json_value(raw, "earth_avoidance_day_deg"),
        "earth_night_deg": _json_value(raw, "earth_avoidance_night_deg"),
    }


def _lookup_radec(
    schedule_csv: Path,
    target_name: str,
    start: datetime,
) -> tuple[float, float]:
    with schedule_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row_target = str(row.get("Target", "")).strip()
            row_start = str(row.get("Observation Start", "")).strip()
            if row_target != target_name:
                continue
            if row_start != start.strftime("%Y-%m-%d %H:%M:%S"):
                continue
            return float(row["RA"]), float(row["DEC"])
    raise ValueError(
        f"Could not find schedule row for target={target_name!r} "
        f"start={start.strftime('%Y-%m-%d %H:%M:%S')} in {schedule_csv}"
    )


def _default_output(args: argparse.Namespace, repo_root: Path) -> Path:
    safe_name = args.target_name.replace("/", "_")
    stem = f"{safe_name}_{args.start.replace(':', '').replace(' ', '_')}_angles.csv"
    if args.schedule_csv is not None:
        return args.schedule_csv.resolve().parent / stem
    if args.config is not None:
        return args.config.resolve().parent / stem
    return repo_root / stem


def main() -> int:
    args = _parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    geometry = _load_geometry_module(repo_root)

    start = _parse_datetime(args.start)
    stop = _parse_datetime(args.stop)
    if stop <= start:
        raise ValueError("stop must be later than start")

    cfg_defaults = _load_config_defaults(args.config)

    gmat_path = args.gmat or (
        Path(str(cfg_defaults["gmat"])).expanduser().resolve()
        if cfg_defaults.get("gmat")
        else None
    )
    if gmat_path is None:
        raise ValueError("Provide --gmat or --config with visibility_gmat")

    earth_day_deg = (
        args.earth_day_deg
        if args.earth_day_deg is not None
        else cfg_defaults.get("earth_day_deg")
    )
    earth_night_deg = (
        args.earth_night_deg
        if args.earth_night_deg is not None
        else cfg_defaults.get("earth_night_deg")
    )
    if earth_day_deg is None or earth_night_deg is None:
        raise ValueError(
            "Provide Earth thresholds via --earth-day-deg/--earth-night-deg or --config"
        )

    if args.ra is None or args.dec is None:
        if args.schedule_csv is None:
            raise ValueError("Provide --ra/--dec or use --schedule-csv for lookup")
        ra_deg, dec_deg = _lookup_radec(args.schedule_csv, args.target_name, start)
    else:
        ra_deg = float(args.ra)
        dec_deg = float(args.dec)

    n_minutes = int(round((stop - start).total_seconds() / 60.0))
    minute_index = pd.date_range(start, periods=n_minutes, freq="min")
    window_end = minute_index[-1].to_pydatetime()

    cadence = geometry.build_minute_cadence(start, window_end)
    eph = geometry.interpolate_gmat_ephemeris(gmat_path, cadence)

    target = np.repeat(_radec_unit(ra_deg, dec_deg), len(minute_index), axis=0)
    sun_unit = _unit(eph.sun_pc)
    moon_unit = _unit(eph.moon_pc)
    nadir_unit = _unit(eph.earth_pc)
    zenith_unit = -nadir_unit

    obs_dist = np.linalg.norm(eph.earth_pc, axis=1)
    limb_angle_rad = np.arccos(np.clip(R_EARTH_KM / obs_dist, -1.0, 1.0))

    sun_sep = _sep_deg(target, sun_unit)
    moon_sep = _sep_deg(target, moon_unit)
    earth_center_sep = _sep_deg(target, nadir_unit)
    earth_limb_deg = _fast_limb_deg(target, zenith_unit, limb_angle_rad)

    dot_zenith_sun = np.einsum("ij,ij->i", zenith_unit, sun_unit)
    subsatellite_day = dot_zenith_sun > 0.0
    earth_threshold = np.where(subsatellite_day, earth_day_deg, earth_night_deg)
    scheduler_earth_visible = earth_center_sep > earth_threshold

    output_path = args.output.resolve() if args.output else _default_output(args, repo_root)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "target_name",
                "utc",
                "sun_deg",
                "moon_deg",
                "earth_center_deg",
                "earth_limb_deg",
                "subsatellite_branch",
                "scheduler_earth_threshold_deg",
                "scheduler_earth_visible",
            ]
        )
        for idx, timestamp in enumerate(minute_index):
            writer.writerow(
                [
                    args.target_name,
                    timestamp.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    round(float(sun_sep[idx]), 2),
                    round(float(moon_sep[idx]), 2),
                    round(float(earth_center_sep[idx]), 2),
                    round(float(earth_limb_deg[idx]), 2),
                    "day" if bool(subsatellite_day[idx]) else "night",
                    round(float(earth_threshold[idx]), 1),
                    "yes" if bool(scheduler_earth_visible[idx]) else "no",
                ]
            )

    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
