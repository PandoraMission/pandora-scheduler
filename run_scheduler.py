#!/usr/bin/env python3
"""
Pandora Scheduler Driver Script
===============================

Pandora Scheduler - Main Entry Point

This script runs the complete Pandora observation scheduling pipeline from start
to finish, generating all necessary output files including:
  - Target manifests (from target definition files)
  - Visibility catalogs (if requested)
  - Observation schedule CSV
  - Science calendar XML
  - Observation time reports
  - Tracker files (CSV and pickle)

Usage:
 # Basic run with default configuration
    poetry run python run_scheduler.py \\
        --start "2026-02-05" \\
        --end "2026-02-12" \\
        --output ./output

    # Run with custom configuration
    poetry run python run_scheduler.py \\
        --start "2026-02-05" \\
        --end "2026-02-12" \\
        --output ./output \\
        --config config.json

    # Generate visibility data as part of the run
    poetry run python run_scheduler.py \\
        --start "2026-02-05" \\
        --end "2026-02-12" \\
        --output ./output \\
        --generate-visibility

    # Use custom target definition files
    poetry run python run_scheduler.py \\
        --start "2026-02-05" \\
        --end "2026-02-12" \\
        --output ./output \\
        --target-definitions ./custom_targets
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, Optional

from pandorascheduler_rework.config import PandoraSchedulerConfig, resolve_data_subdir
from pandorascheduler_rework.pipeline import SchedulerResult, build_schedule
from pandorascheduler_rework.science_calendar import (
    ScienceCalendarInputs,
    generate_science_calendar,
)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the application."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run the Pandora observation scheduler pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Required arguments
    parser.add_argument(
        "--start",
        type=str,
        required=True,
        help="Schedule window start date (YYYY-MM-DD or YYYY-MM-DD HH:MM:SS)",
    )
    parser.add_argument(
        "--end",
        type=str,
        required=True,
        help="Schedule window end date (YYYY-MM-DD or YYYY-MM-DD HH:MM:SS)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output directory for generated files",
    )

    # Optional configuration
    parser.add_argument(
        "--config",
        type=Path,
        help="Path to JSON configuration file (overrides defaults)",
    )
    parser.add_argument(
        "--target-definitions",
        type=Path,
        help=(
            "Base directory containing target definition files (e.g., PandoraTargetList/target_definition_files/). "
            "When provided, generates *_targets.csv manifests from JSON definitions."
        ),
    )
    parser.add_argument(
        "--generate-visibility",
        action="store_true",
        help="Generate visibility catalogs before scheduling (requires target manifests)",
    )

    # Visibility configuration (if generating)
    parser.add_argument(
        "--gmat-ephemeris",
        type=Path,
        help="Path to GMAT ephemeris file (for visibility generation)",
    )
    parser.add_argument(
        "--sun-avoidance",
        type=float,
        default=91.0,
        help="Sun avoidance angle in degrees (default: 91.0)",
    )
    parser.add_argument(
        "--moon-avoidance",
        type=float,
        default=25.0,
        help="Moon avoidance angle in degrees (default: 25.0)",
    )
    parser.add_argument(
        "--earth-avoidance",
        type=float,
        default=110.0,
        help="Earth avoidance angle in degrees (default: 110.0)",
    )
    parser.add_argument(
        "--earth-avoidance-day",
        type=float,
        default=None,
        help="Earth avoidance when nearest limb is sunlit (degrees). None = use --earth-avoidance uniformly.",
    )
    parser.add_argument(
        "--earth-avoidance-night",
        type=float,
        default=None,
        help="Earth avoidance when nearest limb is in shadow (degrees). None = use --earth-avoidance uniformly.",
    )
    parser.add_argument(
        "--twilight-margin",
        type=float,
        default=0.0,
        help="Degrees past geometric terminator to classify as sunlit for day/night keepout (default: 0 = sharp terminator).",
    )

    parser.add_argument(
        "--daynight-mode",
        type=str,
        default="subsatellite",
        dest="daynight_mode",
        help=(
            "Day/night mode for visibility calculations: "
            "'subsatellite' = classify by the subsatellite point, "
            "'limb' = classify by the nearest Earth limb point in the target direction "
            "(default: 'subsatellite')"
        ),
    )

    # Star tracker keepout configuration
    parser.add_argument(
        "--st-sun-min",
        type=float,
        default=0.0,
        help="Star tracker Sun keepout angle in degrees (default: 0 = disabled)",
    )
    parser.add_argument(
        "--st-moon-min",
        type=float,
        default=0.0,
        help="Star tracker Moon keepout angle in degrees (default: 0 = disabled)",
    )
    parser.add_argument(
        "--st-earthlimb-min",
        type=float,
        default=0.0,
        help="Star tracker Earth-limb keepout angle in degrees (default: 0 = disabled)",
    )
    parser.add_argument(
        "--st1-earthlimb-min",
        type=float,
        default=None,
        help="ST1 Earth-limb keepout override (degrees). None = use --st-earthlimb-min.",
    )
    parser.add_argument(
        "--st2-earthlimb-min",
        type=float,
        default=None,
        help="ST2 Earth-limb keepout override (degrees). None = use --st-earthlimb-min.",
    )
    parser.add_argument(
        "--st-required",
        type=int,
        default=1,
        help="Number of star trackers required: 0 (skip), 1 (OR), 2 (AND) (default: 1)",
    )

    # Roll sweep configuration
    parser.add_argument(
        "--roll-step",
        type=float,
        default=2.0,
        help="Roll sweep step size in degrees (default: 2.0)",
    )
    parser.add_argument(
        "--min-power-frac",
        type=float,
        default=0.7,
        help="Minimum solar power fraction to accept a roll angle (default: 0.7)",
    )
    parser.add_argument(
        "--min-sequence-minutes",
        type=int,
        default=None,
        help="Minimum contiguous sequence duration in minutes (default: 8)",
    )

    # Scheduling configuration
    parser.add_argument(
        "--schedule-step-hours",
        type=float,
        default=24.0,
        help=(
            "Scheduler rolling window step size in hours (default: 24.0). "
            "Per-target visit duration comes from target manifests (Obs Window (hrs))."
        ),
    )
    parser.add_argument(
        "--transit-coverage",
        type=float,
        default=0.4,
        help="Minimum transit coverage fraction (default: 0.4)",
    )
    parser.add_argument(
        "--weights",
        type=str,
        default="0.8,0.0,0.2",
        help="Schedule weights as comma-separated values: coverage,saa,schedule (default: 0.8,0.0,0.2)",
    )
    parser.add_argument(
        "--min-visibility",
        type=float,
        default=0.5,
        help="Minimum visibility fraction for non-transit observations (default: 0.5)",
    )

    # Flags
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--skip-xml",
        action="store_true",
        help="Skip generation of the science calendar XML",
    )
    parser.add_argument(
        "--primary-only",
        action="store_true",
        help="Only schedule primary science targets; skip non-primary gap filling",
    )
    parser.add_argument(
        "--use-target-list-for-occultations",
        action="store_true",
        help="Use the target list for occultation scheduling instead of a separate list",
    )
    parser.add_argument(
        "--prioritise-occultations-by-slew",
        action="store_true",
        help="Prioritise occultation targets based on slew cost",
    )
    parser.add_argument(
        "--no-break-occultation-sequences",
        action="store_true",
        help="Disable splitting long occultation sequences into chunks",
    )
    parser.add_argument(
        "--no-occultation-xml",
        action="store_true",
        help="Skip occultation-target calculations during XML generation",
    )
    parser.add_argument(
        "--skip-occultation-pass1",
        action="store_true",
        help="Skip Pass 1 in occultation assignment (single target must cover all intervals)",
    )
    parser.add_argument(
        "--requested-occ-time-override",
        action="store_true",
        help="Allow occultation scheduling to continue when requested-hours has been met",
    )
    parser.add_argument(
        "--allow-occ-st-violation",
        action="store_true",
        help="Allow occultation targets that violate only star-tracker keepout (prefer shortest loss)",
    )

    # Profiling configuration
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Enable profiling",
    )
    parser.add_argument(
        "--profile-output",
        default="profile_output.prof",
        help="Profile output file",
    )

    parser.add_argument(
        "--show-progress",
        action="store_true",
        help="Show progress bars during execution",
    )
    parser.add_argument(
        "--skip-manifests",
        action="store_true",
        help="Skip regenerating target manifests from target definition files",
    )
    parser.add_argument(
        "--legacy-mode",
        action="store_true",
        help=(
            "Use legacy scheduling algorithms for validation against historical outputs. "
            "When enabled, uses MJD-based visibility filtering which matches the original "
            "scheduler exactly. Default (disabled) uses improved datetime-based filtering."
        ),
    )
    parser.add_argument(
        "--parallel-workers",
        type=int,
        default=0,
        help=(
            "Number of parallel workers for visibility generation. "
            "0 = auto (all CPUs), 1 = serial. (default: 0)"
        ),
    )

    return parser.parse_args()


def parse_datetime(date_str: str) -> datetime:
    """Parse a date string into a datetime object."""
    try:
        # Try full datetime format first
        return datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        try:
            # Fallback to date only
            return datetime.strptime(date_str, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"Invalid date format: {date_str}. Expected YYYY-MM-DD or YYYY-MM-DD HH:MM:SS"
            )


def print_summary(result: SchedulerResult, xml_path: Optional[Path]) -> None:
    """Print a summary of the scheduling run."""
    try:

        def _hours(td: pd.Timedelta) -> float:
            return float(td.total_seconds() / 3600.0)

        def _fmt_hours(value: float) -> str:
            return f"{value:,.2f} h"

        def _safe_to_datetime(series: pd.Series) -> pd.Series:
            # Schedule CSV can contain strings or timestamps; normalize defensively.
            return pd.to_datetime(series, errors="coerce")

        import pandas as pd

        print("\n" + "=" * 80)
        print("SCHEDULING PIPELINE COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print("\nGenerated Files:")
        print("-" * 80)

        if getattr(result, "schedule_csv", None):
            print(f"  📄 Schedule CSV:      {result.schedule_csv}")
        reports = getattr(result, "reports", {}) or {}
        if reports.get("observation_time"):
            print(f"  📊 Observation Time   {reports.get('observation_time')}")
        if reports.get("tracker_csv"):
            print(f"  📊 Tracker Csv        {reports.get('tracker_csv')}")
        if reports.get("tracker_pickle"):
            print(f"  📊 Tracker Pickle     {reports.get('tracker_pickle')}")
        if xml_path:
            print(f"  📑 Science Calendar:  {xml_path}")

        print("-" * 80)
        print("\nSchedule Statistics:")
        schedule_df = None
        diagnostics = getattr(result, "diagnostics", {}) or {}
        schedule_df = diagnostics.get("schedule_dataframe")

        if schedule_df is not None and len(schedule_df) > 0:
            try:
                df = schedule_df.copy()
                if "Target" in df.columns:
                    df["Target"] = df["Target"].astype(str).str.strip()

                if "Observation Start" in df.columns:
                    df["Observation Start"] = _safe_to_datetime(df["Observation Start"])
                if "Observation Stop" in df.columns:
                    df["Observation Stop"] = _safe_to_datetime(df["Observation Stop"])

                total_obs = int(len(df))
                unique_targets = (
                    int(df["Target"].nunique()) if "Target" in df.columns else 0
                )

                # Categorize
                if "Target" in df.columns:
                    is_free = df["Target"].eq("Free Time")
                    is_std = df["Target"].str.contains(r"\bSTD\b", na=False)
                    is_primary = (~is_free) & (~is_std)
                else:
                    is_free = pd.Series([False] * len(df))
                    is_std = pd.Series([False] * len(df))
                    is_primary = pd.Series([True] * len(df))

                primary_count = int(is_primary.sum())
                std_count = int(is_std.sum())
                free_count = int(is_free.sum())

                # Durations
                durations = None
                if (
                    "Observation Start" in df.columns
                    and "Observation Stop" in df.columns
                ):
                    durations = df["Observation Stop"] - df["Observation Start"]
                    durations = durations.where(durations.notna(), pd.Timedelta(0))
                    # Guard against negative durations due to bad parsing
                    durations = durations.clip(lower=pd.Timedelta(0))

                # Schedule span
                total_span_hours = None
                if (
                    "Observation Start" in df.columns
                    and "Observation Stop" in df.columns
                    and df["Observation Start"].notna().any()
                    and df["Observation Stop"].notna().any()
                ):
                    sched_start = df["Observation Start"].min()
                    sched_stop = df["Observation Stop"].max()
                    if pd.notna(sched_start) and pd.notna(sched_stop):
                        total_span_hours = _hours(sched_stop - sched_start)

                # Time totals by category
                primary_hours = std_hours = free_hours = None
                if durations is not None:
                    primary_hours = _hours(durations[is_primary].sum())
                    std_hours = _hours(durations[is_std].sum())
                    free_hours = _hours(durations[is_free].sum())

                print(f"  Total observations:        {total_obs}")
                print(f"  Unique targets:            {unique_targets}")
                print(f"  Primary observations:      {primary_count}")
                print(f"  Standard (STD) blocks:     {std_count}")
                print(f"  Free Time blocks:          {free_count}")

                if total_span_hours is not None:
                    print(
                        f"  Schedule span:             {_fmt_hours(total_span_hours)}"
                    )

                if (
                    primary_hours is not None
                    and std_hours is not None
                    and free_hours is not None
                ):
                    used_hours = primary_hours + std_hours
                    print(f"  Primary time:              {_fmt_hours(primary_hours)}")
                    print(f"  STD time:                  {_fmt_hours(std_hours)}")
                    print(f"  Free time:                 {_fmt_hours(free_hours)}")
                    if total_span_hours is not None and total_span_hours > 0:
                        util = 100.0 * (used_hours / total_span_hours)
                        print(f"  Utilization (non-free):    {util:,.1f}%")

                if (
                    durations is not None
                    and durations.astype("timedelta64[ns]").notna().any()
                ):

                    def _dur_stats(
                        mask: pd.Series,
                    ) -> tuple[float | None, float | None]:
                        subset = durations[mask]
                        if subset.empty:
                            return None, None
                        return _hours(subset.mean()), _hours(subset.median())

                    p_mean, p_med = _dur_stats(is_primary)
                    if p_mean is not None:
                        print(f"  Primary duration (mean):   {_fmt_hours(p_mean)}")
                        print(
                            f"  Primary duration (median): {_fmt_hours(p_med or 0.0)}"
                        )

                # Transit-quality stats (only for rows that have Transit Coverage)
                if "Transit Coverage" in df.columns and is_primary.any():
                    cov = pd.to_numeric(
                        df.loc[is_primary, "Transit Coverage"], errors="coerce"
                    )
                    cov = cov.dropna()
                    if len(cov) > 0:
                        cov_mean = float(cov.mean())
                        cov_median = float(cov.median())
                        cov_min = float(cov.min())
                        cov_p10 = float(cov.quantile(0.10))
                        cov_p90 = float(cov.quantile(0.90))
                        full = int((cov >= 0.999999).sum())
                        good = int((cov >= 0.90).sum())
                        print("\n  Primary transit coverage:")
                        print(
                            f"    mean/median:             {cov_mean:.3f} / {cov_median:.3f}"
                        )
                        print(
                            f"    min / p10 / p90:         {cov_min:.3f} / {cov_p10:.3f} / {cov_p90:.3f}"
                        )
                        print(
                            f"    >= 0.90:                 {good} ({100.0 * good / len(cov):.1f}%)"
                        )
                        print(
                            f"    ~100% (>= 0.999999):     {full} ({100.0 * full / len(cov):.1f}%)"
                        )

                # SAA overlap summary
                if "SAA Overlap" in df.columns and is_primary.any():
                    saa = pd.to_numeric(
                        df.loc[is_primary, "SAA Overlap"], errors="coerce"
                    )
                    saa = saa.dropna()
                    if len(saa) > 0:
                        print("\n  Primary SAA overlap:")
                        print(
                            f"    mean/median:             {float(saa.mean()):.3f} / {float(saa.median()):.3f}"
                        )
            except Exception:
                print("  (Unable to compute detailed schedule statistics)")
        else:
            print("  (No schedule dataframe available for statistics)")

        print("\n" + "=" * 80)
    except AttributeError as exc:
        # Be defensive: if the SchedulerResult shape is unexpected, log and provide minimal info.
        logger = logging.getLogger(__name__)
        logger.error("Failed to generate summary: %s", exc)
        try:
            # Best-effort fallback: print whatever schedule CSV path exists.
            if hasattr(result, "schedule_csv") and result.schedule_csv:
                print(f"Schedule generated: {result.schedule_csv}")
        except Exception:
            # Give up cleanly; avoid raising further exceptions from summary printing.
            pass


def main() -> int:
    """Main execution function."""
    args = parse_args()
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Initialize profiler if requested
    profiler = None
    if args.profile:
        import cProfile
        import pstats

        profiler = cProfile.Profile()
        profiler.enable()

    try:
        logger.info(f"Scheduling window: {args.start} to {args.end}")
        logger.info(f"Output directory: {args.output}")

        # 1. Load Configuration
        json_config = {}
        if args.config:
            with open(args.config, "r") as f:
                json_config = json.load(f)

        def _get_val(key: str, cli_value: Any, default: Any) -> Any:
            """Prioritize CLI-provided value > JSON config > default.

            Because many CLI args have non-None defaults, we treat "CLI provided" as
            "CLI value differs from the default".
            """

            if cli_value is not None and cli_value != default:
                return cli_value
            if key in json_config and json_config[key] is not None:
                return json_config[key]
            return default

        def _get_any(keys: list[str], cli_value: Any, default: Any) -> Any:
            if cli_value is not None and cli_value != default:
                return cli_value
            for key in keys:
                if key in json_config and json_config[key] is not None:
                    return json_config[key]
            return default

        def _as_bool(value: Any, default: bool = False) -> bool:
            if value is None:
                return default
            if isinstance(value, bool):
                return value
            if isinstance(value, str):
                normalized = value.strip().lower()
                if normalized == "":
                    return default
                return normalized in {"1", "true", "yes", "y", "on"}
            if isinstance(value, (int, float)):
                return bool(value)
            return default

        # Default weights if not provided
        transit_scheduling_weights = (0.8, 0.0, 0.2)

        # Start with any pipeline extra inputs provided via JSON.
        extra_inputs: Dict[str, Any] = {}
        if isinstance(json_config.get("extra_inputs"), dict):
            extra_inputs.update(json_config.get("extra_inputs") or {})

        # Resolve target definition base (CLI overrides JSON)
        target_def_base = args.target_definitions or extra_inputs.get(
            "target_definition_base"
        )

        # Resolve visibility GMAT file (CLI overrides JSON)
        visibility_gmat = args.gmat_ephemeris or extra_inputs.get("visibility_gmat")

        # 2. Validate Inputs
        if args.generate_visibility and not target_def_base:
            logger.error(
                "Visibility generation requires target definitions. "
                "Please provide target definitions via --target-definitions"
            )
            return 1

        # Determine whether visibility generation was requested (CLI or JSON)
        generate_visibility = bool(args.generate_visibility) or (
            str(
                extra_inputs.get(
                    "generate_visibility", json_config.get("generate_visibility", "")
                )
            ).lower()
            in {"1", "true", "yes", "y"}
        )
        # `config` is not yet constructed here, so validate the CLI/JSON inputs now
        if generate_visibility and visibility_gmat is None:
            logger.error(
                "Visibility generation requested but no GMAT ephemeris provided. "
                "Supply --gmat-ephemeris or set extra_inputs.visibility_gmat in the JSON config."
            )
            return 1

        # Build PandoraSchedulerConfig with the dataclass field names and types
        schedule_step_hours = float(
            _get_val("schedule_step_hours", args.schedule_step_hours, 24.0)
        )
        transit_cov = float(
            _get_val("transit_coverage_min", args.transit_coverage, 0.4)
        )
        min_vis = float(_get_val("min_visibility", args.min_visibility, 0.5))

        # Coerce unified transit_scheduling_weights from JSON or CLI into a 3-tuple
        raw_transit_weights = _get_val(
            "transit_scheduling_weights", args.weights, transit_scheduling_weights
        )
        if isinstance(raw_transit_weights, str):
            raw_transit_weights = tuple(
                float(x.strip()) for x in raw_transit_weights.split(",")
            )
        transit_weights_tuple = tuple(float(x) for x in raw_transit_weights)
        if len(transit_weights_tuple) != 3:
            raise ValueError(
                "transit_scheduling_weights must contain exactly 3 values"
            )

        if target_def_base:
            extra_inputs["target_definition_base"] = Path(str(target_def_base))
            # If the JSON doesn't provide an explicit list, use the standard set.
            extra_inputs.setdefault(
                "target_definition_files",
                [
                    "exoplanet",
                    "auxiliary-standard",
                    "monitoring-standard",
                    "occultation-standard",
                ],
            )

        if args.skip_manifests:
            extra_inputs["skip_manifests"] = True

        # Allow JSON-only control of skip_manifests
        if str(extra_inputs.get("skip_manifests", "")).lower() in {"1", "true", "yes", "y"}:
            extra_inputs["skip_manifests"] = True

        # Visibility GMAT goes into the typed field `gmat_ephemeris` on the config
        if visibility_gmat is None:
            gmat_path = None
        else:
            gmat_path = Path(str(visibility_gmat)).expanduser().resolve()

        sun_avoid = float(
            _get_any(
                ["sun_avoidance_deg", "visibility_sun_deg"], args.sun_avoidance, 91.0
            )
        )
        moon_avoid = float(
            _get_any(
                ["moon_avoidance_deg", "visibility_moon_deg"], args.moon_avoidance, 25.0
            )
        )
        earth_avoid = float(
            _get_any(
                ["earth_avoidance_deg", "visibility_earth_deg"], args.earth_avoidance, 110.0
            )
        )
        # Day/night Earth avoidance (None = use uniform earth_avoid)
        _raw_day = _get_val("earth_avoidance_day_deg", args.earth_avoidance_day, None)
        earth_avoid_day = float(_raw_day) if _raw_day is not None else None
        _raw_night = _get_val("earth_avoidance_night_deg", args.earth_avoidance_night, None)
        earth_avoid_night = float(_raw_night) if _raw_night is not None else None
        data_subdir = resolve_data_subdir(
            extra_inputs,
            sun_avoidance_deg=sun_avoid,
            moon_avoidance_deg=moon_avoid,
            earth_avoidance_deg=earth_avoid,
            earth_avoidance_day_deg=earth_avoid_day,
        )
        extra_inputs["data_subdir"] = data_subdir

        twilight_margin = float(
            _get_val("twilight_margin_deg", args.twilight_margin, 0.0)
        )

        daynight_mode = str(
            _get_val("daynight_mode", args.daynight_mode, "subsatellite")
        ).lower()

        # Star tracker keepouts
        st_sun_min = float(
            _get_val("st_sun_min_deg", args.st_sun_min, 0.0)
        )
        st_moon_min = float(
            _get_val("st_moon_min_deg", args.st_moon_min, 0.0)
        )
        st_earthlimb_min = float(
            _get_val("st_earthlimb_min_deg", args.st_earthlimb_min, 0.0)
        )
        _raw_st1_el = _get_any(
            ["st1_earthlimb_min_deg"],
            getattr(args, "st1_earthlimb_min", None),
            None,
        )
        st1_earthlimb_min = float(_raw_st1_el) if _raw_st1_el is not None else None
        _raw_st2_el = _get_any(
            ["st2_earthlimb_min_deg"],
            getattr(args, "st2_earthlimb_min", None),
            None,
        )
        st2_earthlimb_min = float(_raw_st2_el) if _raw_st2_el is not None else None
        st_required = int(
            _get_val("st_required", args.st_required, 1)
        )

        # Roll sweep
        roll_step = float(
            _get_val("roll_step_deg", args.roll_step, 2.0)
        )
        min_power_frac = float(
            _get_val("min_power_frac", args.min_power_frac, 0.7)
        )

        short_visit_threshold_hours = float(
            _get_val("short_visit_threshold_hours", None, 12.0)
        )
        short_visit_edge_buffer_hours = float(
            _get_val("short_visit_edge_buffer_hours", None, 1.5)
        )
        long_visit_edge_buffer_hours = float(
            _get_val("long_visit_edge_buffer_hours", None, 4.0)
        )

        obs_sequence_duration_min = int(_get_val("obs_sequence_duration_min", None, 90))
        occ_sequence_limit_min = int(_get_val("occ_sequence_limit_min", None, 40))
        min_sequence_minutes = int(_get_val("min_sequence_minutes", args.min_sequence_minutes, 8))

        std_obs_duration_hours = float(_get_val("std_obs_duration_hours", None, 0.5))
        std_obs_frequency_days = float(_get_val("std_obs_frequency_days", None, 3.0))

        force_regenerate = bool(_get_val("force_regenerate", None, False))
        primary_only_mode = bool(
            _get_val("primary_only_mode", args.primary_only, False)
        )
        use_target_list_for_occultations = _as_bool(
            _get_val("use_target_list_for_occultations", args.use_target_list_for_occultations, False),
            False,
        )
        prioritise_occultations_by_slew = _as_bool(
            _get_val("prioritise_occultations_by_slew", args.prioritise_occultations_by_slew, False),
            False,
        )
        # CLI --no-break-occultation-sequences inverts the default-True field
        break_occ_cli = not args.no_break_occultation_sequences if args.no_break_occultation_sequences else None
        break_occultation_sequences = _as_bool(
            _get_val("break_occultation_sequences", break_occ_cli, True), True
        )
        enable_occultation_xml = _as_bool(
            _get_any(
                ["generate_occultation_xml", "enable_occultation_xml"],
                not args.no_occultation_xml if args.no_occultation_xml else None,
                True,
            ),
            True,
        )
        enable_occultation_pass1 = _as_bool(
            _get_any(
                ["one_occultation_target", "enable_occultation_pass1"],
                not args.skip_occultation_pass1 if args.skip_occultation_pass1 else None,
                True,
            ),
            True,
        )
        requested_occ_time_override = _as_bool(
            _get_val(
                "requested_occ_time_override",
                True if args.requested_occ_time_override else None,
                None,
            ),
            False,
        )
        allow_occ_startracker_violation = _as_bool(
            _get_val(
                "allow_occ_startracker_violation",
                True if args.allow_occ_st_violation else None,
                None,
            ),
            False,
        )

        commissioning_days = int(_get_val("commissioning_days", None, 0))

        show_progress = bool(_get_val("show_progress", args.show_progress, False))
        use_legacy_mode = bool(_get_any(["use_legacy_mode", "legacy_mode"], args.legacy_mode, False))
        parallel_workers = int(
            _get_val("parallel_workers", args.parallel_workers, 0)
        )

        aux_sort_key = str(_get_val("aux_sort_key", None, "sort_by_tdf_priority"))
        author = _get_val("author", None, None)
        created_timestamp = _get_val("created_timestamp", None, None)
        visit_limit = _get_val("visit_limit", None, None)
        target_filters = _get_val("target_filters", None, ())
        if target_filters is None:
            target_filters = ()

        output_dir = args.output

        config = PandoraSchedulerConfig(
            window_start=parse_datetime(args.start),
            window_end=parse_datetime(args.end),
            schedule_step=timedelta(hours=schedule_step_hours),
            targets_manifest=output_dir / data_subdir,
            gmat_ephemeris=gmat_path,
            output_dir=output_dir,
            # Scheduling Thresholds
            transit_coverage_min=transit_cov,
            min_visibility=min_vis,
            commissioning_days=commissioning_days,
            # Transit edge buffers
            short_visit_threshold_hours=short_visit_threshold_hours,
            short_visit_edge_buffer_hours=short_visit_edge_buffer_hours,
            long_visit_edge_buffer_hours=long_visit_edge_buffer_hours,
            # Weights
            transit_scheduling_weights=transit_weights_tuple,
            # Keepout angles
            sun_avoidance_deg=sun_avoid,
            moon_avoidance_deg=moon_avoid,
            earth_avoidance_deg=earth_avoid,
            earth_avoidance_day_deg=earth_avoid_day,
            earth_avoidance_night_deg=earth_avoid_night,
            twilight_margin_deg=twilight_margin,
            daynight_mode=daynight_mode,
            # Star tracker keepouts
            st_sun_min_deg=st_sun_min,
            st_moon_min_deg=st_moon_min,
            st_earthlimb_min_deg=st_earthlimb_min,
            st1_earthlimb_min_deg=st1_earthlimb_min,
            st2_earthlimb_min_deg=st2_earthlimb_min,
            st_required=st_required,
            # Roll sweep
            roll_step_deg=roll_step,
            min_power_frac=min_power_frac,
            # XML / sequence generation
            obs_sequence_duration_min=obs_sequence_duration_min,
            occ_sequence_limit_min=occ_sequence_limit_min,
            min_sequence_minutes=min_sequence_minutes,
            break_occultation_sequences=break_occultation_sequences,
            # Standard observations
            std_obs_duration_hours=std_obs_duration_hours,
            std_obs_frequency_days=std_obs_frequency_days,
            # Extra inputs for pipeline
            extra_inputs=extra_inputs,
            # Flags
            show_progress=show_progress,
            force_regenerate=force_regenerate,
            primary_only_mode=primary_only_mode,
            use_target_list_for_occultations=use_target_list_for_occultations,
            prioritise_occultations_by_slew=prioritise_occultations_by_slew,
            enable_occultation_xml=enable_occultation_xml,
            enable_occultation_pass1=enable_occultation_pass1,
            requested_occ_time_override=requested_occ_time_override,
            allow_occ_startracker_violation=allow_occ_startracker_violation,
            use_legacy_mode=use_legacy_mode,
            parallel_workers=parallel_workers,
            # Sorting / metadata
            aux_sort_key=aux_sort_key,
            author=author,
            created_timestamp=created_timestamp,
            visit_limit=visit_limit,
            target_filters=target_filters,
        )

        # 3. Ensure targets manifest location exists (may be an output/data dir)
        if config.targets_manifest and not config.targets_manifest.exists():
            try:
                config.targets_manifest.mkdir(parents=True, exist_ok=True)
                logger.info(
                    "Created targets manifest directory: %s", config.targets_manifest
                )
            except Exception:
                logger.debug(
                    "Unable to create targets manifest directory: %s",
                    config.targets_manifest,
                )

        # 4. Run Scheduler (using new API)
        logger.info("Run data directory: %s", output_dir / data_subdir)
        logger.info("PRIMARY_ONLY_MODE=%s", str(primary_only_mode).upper())
        logger.info("GENERATE_OCCULTATION_XML=%s", str(enable_occultation_xml).upper())
        logger.info("OCCULTATION_PASS1=%s", str(enable_occultation_pass1).upper())
        logger.info("REQUESTED_OCC_TIME_OVERRIDE=%s", str(requested_occ_time_override).upper())
        logger.info("Starting scheduler pipeline...")
        if args.legacy_mode:
            logger.info("Legacy mode enabled - using MJD-based visibility filtering")
        result = build_schedule(config)

        # 4. Generate Science Calendar XML
        xml_path = None
        if not args.skip_xml and result.schedule_csv:
            # Use the same config object to create calendar settings
            # Ensure we point to the correct data directory (where manifests are)
            data_dir = output_dir / data_subdir

            inputs = ScienceCalendarInputs(
                schedule_csv=result.schedule_csv,
                data_dir=data_dir,
            )

            xml_path = generate_science_calendar(
                inputs=inputs,
                config=config,
            )
            logger.info(f"Science calendar written to: {xml_path}")

        # 5. Print Summary
        print_summary(result, xml_path)

        return 0

    except KeyboardInterrupt:
        logger.error("\nInterrupted by user")
        return 130
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=args.verbose)
        return 1
    finally:
        if profiler:
            profiler.disable()
            stats = pstats.Stats(profiler).sort_stats("cumulative")
            stats.dump_stats(args.profile_output)
            print(f"\nProfiling results written to {args.profile_output}")
            stats.print_stats(30)


if __name__ == "__main__":
    sys.exit(main())
