#!/usr/bin/env python3
"""Verify that scheduled transit observations overlap actual transit windows.

This script is intended as an independent consistency check on produced outputs.
It compares each scheduled transit observation window in the schedule CSV
against the transit start/stop times in the per-planet visibility parquet.

It reports observations that:
- do not overlap any transit window for the scheduled planet, or
- have a transit-coverage mismatch vs. the matched transit's Transit_Coverage.

Note: the schedule CSV's "Transit Coverage" is the *visibility coverage of the transit*
(computed during visibility generation), not the fraction of the scheduled visit spent
inside the transit window.

Usage examples:
  poetry run python scripts/verify_scheduled_transits.py --output-dir output_standalone
  poetry run python scripts/verify_scheduled_transits.py --schedule-csv output_standalone/Pandora_Schedule_*.csv --data-dir output_standalone/data

Exit status:
  0 if all checks pass (within tolerances)
  2 if any mismatches are found
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy import units as u


@dataclass(frozen=True)
class TransitWindow:
    start: datetime
    stop: datetime
    transit_coverage: Optional[float]


def _as_float(value: object) -> Optional[float]:
    try:
        val = float(value)
    except Exception:
        return None
    if np.isnan(val):
        return None
    return val


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify schedule primary windows overlap transit windows"
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory for a run (e.g. output_standalone). If provided, the script "
            "will infer --data-dir and auto-pick the newest Pandora_Schedule_*.csv unless --schedule-csv is set."
        ),
    )
    parser.add_argument(
        "--schedule-csv",
        type=Path,
        default=None,
        help="Path to schedule CSV. If omitted, inferred from --output-dir.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help=(
            "Path to the run data directory (contains exoplanet_targets.csv and targets/). "
            "If omitted, inferred as <output-dir>/data."
        ),
    )

    parser.add_argument(
        "--transit-source",
        choices=["parquet", "ephemeris"],
        default="parquet",
        help=(
            "Source of transit windows used for verification. "
            "'parquet' reads Transit_Start/Transit_Stop from the per-planet visibility parquet. "
            "'ephemeris' regenerates transit windows from exoplanet_targets.csv period/epoch/duration."
        ),
    )

    parser.add_argument(
        "--coverage-tol",
        type=float,
        default=0.05,
        help="Allowed absolute difference between computed and recorded Transit Coverage (default: 0.05).",
    )
    parser.add_argument(
        "--min-overlap-minutes",
        "--min-overlap",
        type=float,
        default=1.0,
        help=(
            "Minimum required overlap, in minutes, to count as matching a transit (default: 1.0). "
            "This avoids edge-only overlaps due to rounding."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of primary observations checked (for quick runs).",
    )
    parser.add_argument(
        "--write-report",
        type=Path,
        default=None,
        help="Optional path to write a CSV report of mismatches.",
    )

    return parser.parse_args()


def _pick_latest_schedule_csv(output_dir: Path) -> Path:
    candidates = sorted(output_dir.glob("Pandora_Schedule_*.csv"), key=lambda p: p.stat().st_mtime)
    if not candidates:
        raise FileNotFoundError(f"No Pandora_Schedule_*.csv found in {output_dir}")
    return candidates[-1]


def _coerce_datetime_col(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors="coerce", utc=False)


def _load_transit_windows(
    targets_dir: Path, planet_name: str, star_name: str
) -> list[TransitWindow]:
    # Parquet is produced under <data>/targets/<Star Name>/<Planet Name>/Visibility for <Planet>.parquet
    parquet_path = (
        targets_dir
        / str(star_name)
        / str(planet_name)
        / f"Visibility for {planet_name}.parquet"
    )
    if not parquet_path.exists():
        # Some generators might use a different naming; fail loudly.
        raise FileNotFoundError(
            f"Visibility parquet missing for planet={planet_name!r}, star={star_name!r}: {parquet_path}"
        )

    df = pd.read_parquet(
        parquet_path, columns=["Transit_Start", "Transit_Stop", "Transit_Coverage"]
    )
    if df.empty:
        return []

    starts_dt = Time(df["Transit_Start"].to_numpy(), format="mjd", scale="utc").to_value(
        "datetime"
    )
    stops_dt = Time(df["Transit_Stop"].to_numpy(), format="mjd", scale="utc").to_value(
        "datetime"
    )

    windows: list[TransitWindow] = []
    cov_values = pd.to_numeric(df.get("Transit_Coverage"), errors="coerce")
    for start, stop, cov in zip(starts_dt, stops_dt, cov_values.to_numpy()):
        if start is None or stop is None:
            continue
        transit_cov: Optional[float]
        if cov is None or (isinstance(cov, float) and np.isnan(cov)):
            transit_cov = None
        else:
            transit_cov = float(cov)

        windows.append(
            TransitWindow(
                start=pd.Timestamp(start).to_pydatetime(),
                stop=pd.Timestamp(stop).to_pydatetime(),
                transit_coverage=transit_cov,
            )
        )

    return windows


def _load_transit_windows_from_ephemeris(
    planet_row: pd.Series,
    schedule_start_utc: datetime,
    schedule_stop_utc: datetime,
) -> list[TransitWindow]:
    period_days = _as_float(planet_row.get("Period (days)"))
    duration_hours = _as_float(planet_row.get("Transit Duration (hrs)"))

    # The manifest may include either full JD or JD-2400000.5.
    epoch_jd = _as_float(planet_row.get("Transit Epoch (BJD_TDB)"))
    if epoch_jd is None:
        epoch_reduced = _as_float(planet_row.get("Transit Epoch (BJD_TDB-2400000.5)"))
        if epoch_reduced is not None:
            epoch_jd = epoch_reduced + 2400000.5

    if period_days is None or duration_hours is None or epoch_jd is None:
        return []

    # Treat the epoch as a JD-like value in TDB for generating a consistent sequence.
    # Note: If the epoch is a true BJD, this ignores barycentric correction details.
    epoch_tdb = Time(epoch_jd, format="jd", scale="tdb")

    start_tdb_jd = Time(schedule_start_utc, scale="utc").tdb.jd
    stop_tdb_jd = Time(schedule_stop_utc, scale="utc").tdb.jd

    # Expand bounds slightly so edge cases still find a match.
    start_tdb_jd -= period_days
    stop_tdb_jd += period_days

    k_start = int(np.floor((start_tdb_jd - epoch_tdb.jd) / period_days))
    k_stop = int(np.ceil((stop_tdb_jd - epoch_tdb.jd) / period_days))

    half_duration = timedelta(seconds=float(duration_hours) * 3600.0 / 2.0)

    windows: list[TransitWindow] = []
    for k in range(k_start, k_stop + 1):
        center_tdb = epoch_tdb + float(k) * float(period_days) * u.day  # type: ignore[name-defined]
        center_utc_dt = pd.Timestamp(center_tdb.utc.to_datetime()).to_pydatetime()
        windows.append(
            TransitWindow(
                start=center_utc_dt - half_duration,
                stop=center_utc_dt + half_duration,
                transit_coverage=None,
            )
        )
    return windows


def _overlap_seconds(obs_start: datetime, obs_stop: datetime, tw: TransitWindow) -> float:
    latest_start = max(obs_start, tw.start)
    earliest_stop = min(obs_stop, tw.stop)
    return float(max((earliest_stop - latest_start).total_seconds(), 0.0))


def main() -> int:
    args = _parse_args()

    if args.output_dir is None and (args.schedule_csv is None or args.data_dir is None):
        raise SystemExit(
            "Provide either --output-dir, or both --schedule-csv and --data-dir."
        )

    output_dir = args.output_dir
    if output_dir is not None:
        output_dir = output_dir.expanduser().resolve()

    schedule_csv = args.schedule_csv
    if schedule_csv is None:
        if output_dir is None:
            raise SystemExit("--schedule-csv is required if --output-dir is not provided")
        schedule_csv = _pick_latest_schedule_csv(output_dir)
    schedule_csv = schedule_csv.expanduser().resolve()

    data_dir = args.data_dir
    if data_dir is None:
        if output_dir is None:
            raise SystemExit("--data-dir is required if --output-dir is not provided")
        data_dir = output_dir / "data"
    data_dir = data_dir.expanduser().resolve()

    targets_dir = data_dir / "targets"
    manifest_path = data_dir / "exoplanet_targets.csv"

    if not schedule_csv.exists():
        raise FileNotFoundError(f"Schedule CSV not found: {schedule_csv}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Primary manifest not found: {manifest_path}")
    if not targets_dir.exists():
        raise FileNotFoundError(f"Targets directory not found: {targets_dir}")

    schedule = pd.read_csv(schedule_csv)

    required_cols = {"Target", "Observation Start", "Observation Stop"}
    missing = required_cols.difference(schedule.columns)
    if missing:
        raise ValueError(f"Schedule CSV missing required columns: {sorted(missing)}")

    schedule["Target"] = schedule["Target"].astype(str).str.strip()
    schedule["Observation Start"] = _coerce_datetime_col(schedule["Observation Start"])
    schedule["Observation Stop"] = _coerce_datetime_col(schedule["Observation Stop"])

    # Transit observations: rows with numeric Transit Coverage (aux/STD/free time have blanks)
    if "Transit Coverage" not in schedule.columns:
        raise ValueError("Schedule CSV missing required column: 'Transit Coverage'")
    transit_cov_series = pd.to_numeric(schedule["Transit Coverage"], errors="coerce")
    is_transit_obs = transit_cov_series.notna()

    transit_obs = schedule.loc[is_transit_obs].copy()
    if args.limit is not None:
        transit_obs = transit_obs.head(int(args.limit))

    manifest = pd.read_csv(manifest_path)
    if "Planet Name" not in manifest.columns or "Star Name" not in manifest.columns:
        raise ValueError("exoplanet_targets.csv must contain 'Planet Name' and 'Star Name'")

    manifest = manifest.dropna(subset=["Planet Name"]).copy()
    manifest["Planet Name"] = manifest["Planet Name"].astype(str)
    planet_to_star = (
        manifest[["Planet Name", "Star Name"]]
        .dropna()
        .astype({"Planet Name": str, "Star Name": str})
        .set_index("Planet Name")["Star Name"]
        .to_dict()
    )

    planet_rows = manifest.set_index("Planet Name").to_dict(orient="index")

    transit_cache: dict[str, list[TransitWindow]] = {}

    mismatches: list[dict] = []
    checked = 0

    schedule_min_start = pd.Timestamp(transit_obs["Observation Start"].min()).to_pydatetime()
    schedule_max_stop = pd.Timestamp(transit_obs["Observation Stop"].max()).to_pydatetime()

    for _, row in transit_obs.iterrows():
        planet = str(row["Target"])
        obs_start = row["Observation Start"]
        obs_stop = row["Observation Stop"]

        if pd.isna(obs_start) or pd.isna(obs_stop):
            mismatches.append(
                {
                    "Target": planet,
                    "Observation Start": row.get("Observation Start"),
                    "Observation Stop": row.get("Observation Stop"),
                    "Issue": "invalid observation timestamps",
                }
            )
            continue

        obs_start_dt = pd.Timestamp(obs_start).to_pydatetime()
        obs_stop_dt = pd.Timestamp(obs_stop).to_pydatetime()

        if planet not in planet_to_star:
            mismatches.append(
                {
                    "Target": planet,
                    "Observation Start": obs_start_dt,
                    "Observation Stop": obs_stop_dt,
                    "Issue": "planet not found in exoplanet_targets.csv",
                }
            )
            continue

        if planet not in transit_cache:
            if args.transit_source == "parquet":
                star = planet_to_star.get(planet)
                if star is None:
                    transit_cache[planet] = []
                else:
                    transit_cache[planet] = _load_transit_windows(targets_dir, planet, star)
            else:
                row_dict = planet_rows.get(planet)
                if row_dict is None:
                    transit_cache[planet] = []
                else:
                    transit_cache[planet] = _load_transit_windows_from_ephemeris(
                        pd.Series(row_dict), schedule_min_start, schedule_max_stop
                    )

        windows = transit_cache[planet]
        if not windows:
            mismatches.append(
                {
                    "Target": planet,
                    "Observation Start": obs_start_dt,
                    "Observation Stop": obs_stop_dt,
                    "Issue": "no transits found in visibility parquet",
                }
            )
            continue

        best_overlap_sec = 0.0
        best_tw: Optional[TransitWindow] = None
        for tw in windows:
            overlap_sec = _overlap_seconds(obs_start_dt, obs_stop_dt, tw)
            if overlap_sec > best_overlap_sec:
                best_overlap_sec = overlap_sec
                best_tw = tw

        checked += 1

        best_overlap_min = best_overlap_sec / 60.0
        if best_tw is None or best_overlap_min < float(args.min_overlap_minutes):
            mismatches.append(
                {
                    "Target": planet,
                    "Observation Start": obs_start_dt,
                    "Observation Stop": obs_stop_dt,
                    "Issue": "no matching transit overlap",
                    "Best Overlap Minutes": best_overlap_min,
                }
            )
            continue

        # Compare schedule Transit Coverage against the matched transit's Transit_Coverage
        # (only available/meaningful when transit windows come from the visibility parquet).
        if args.transit_source == "parquet":
            recorded_cov = _as_float(row.get("Transit Coverage"))
            if recorded_cov is not None:
                parquet_cov = best_tw.transit_coverage
                if parquet_cov is not None and abs(recorded_cov - parquet_cov) > float(
                    args.coverage_tol
                ):
                    mismatches.append(
                        {
                            "Target": planet,
                            "Observation Start": obs_start_dt,
                            "Observation Stop": obs_stop_dt,
                            "Issue": "coverage mismatch",
                            "Recorded Transit Coverage": recorded_cov,
                            "Parquet Transit_Coverage": parquet_cov,
                            "Best Overlap Minutes": best_overlap_min,
                            "Matched Transit Start": best_tw.start,
                            "Matched Transit Stop": best_tw.stop,
                        }
                    )

    total_primary = int(len(transit_obs))
    mismatch_count = int(len(mismatches))

    print("=")
    print("Scheduled transit verification")
    print("=")
    print(f"Schedule CSV: {schedule_csv}")
    print(f"Data dir:     {data_dir}")
    print(f"Transit src:  {args.transit_source}")
    print(f"Primary rows: {total_primary}")
    print(f"Checked:      {checked}")
    print(f"Mismatches:   {mismatch_count}")

    if mismatch_count:
        sample = pd.DataFrame(mismatches)
        with pd.option_context("display.max_rows", 20, "display.max_columns", 50):
            print("\nSample mismatches (up to 20):")
            print(sample.head(20).to_string(index=False))

        if args.write_report is not None:
            report_path = args.write_report.expanduser().resolve()
        elif output_dir is not None:
            report_path = output_dir / "transit_verification_mismatches.csv"
        else:
            report_path = None

        if report_path is not None:
            sample.to_csv(report_path, index=False)
            print(f"\nWrote mismatch report: {report_path}")

        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
