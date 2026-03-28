"""Reworked scheduler implementation mirroring the legacy behaviour."""

from __future__ import annotations

import logging
import pickle
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Optional, cast

import numpy as np
import pandas as pd
from astropy.time import Time
from tqdm import tqdm

from pandorascheduler_rework import observation_utils
from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.utils.io import (
    build_star_visibility_path,
    build_visibility_path,
    read_csv_cached,
    read_parquet_cached,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# ---------------------------------------------------------------------------
# Transit Factor thresholds — intentionally distinct, do not unify.
# ---------------------------------------------------------------------------
# TOO_OVERRIDE_THRESHOLD: transit_ratio = transits_left_lifetime / transits_needed.
#   When this ratio ≤ 1.0 the spacecraft has no margin — every remaining
#   observable transit must be captured.  Triggers a forced Target-of-Opportunity
#   override that bypasses normal scheduling priority.
# URGENCY_THRESHOLD: when ANY candidate's Transit Factor ≤ 1.2 the window sort
#   switches to urgency-first (ascending TF) so the most time-pressured target
#   is preferred over pure quality-factor ranking.
_TOO_OVERRIDE_THRESHOLD = 1.0
_URGENCY_THRESHOLD = 1.2


def _filter_visibility_by_time(
    vis: pd.DataFrame,
    start: datetime,
    stop: datetime,
    use_legacy_mode: bool,
) -> pd.DataFrame:
    """Filter visibility data to a time window.

    Args:
        vis: Visibility DataFrame with Time(MJD_UTC) and optionally Time_UTC columns
        start: Start of time window
        stop: End of time window
        use_legacy_mode: If True, use MJD-based filtering (matches legacy exactly).
                        If False, use datetime-based filtering (more precise at boundaries).

    Returns:
        Filtered visibility DataFrame

    Note:
        Legacy mode uses MJD filtering which can exclude boundary points due to
        floating-point precision. For example, if start converts to MJD 61079.6048611111
        but the data has 61079.6048610000, the row is excluded. Datetime filtering
        includes such boundary rows correctly.

        The practical difference is typically 0-2 rows at window boundaries, which
        can affect visibility percentage calculations and thus target selection
        when multiple targets have similar visibility.
    """
    if use_legacy_mode:
        # Legacy: MJD-based filtering (matches original scheduler exactly)
        start_mjd = Time(start).mjd
        stop_mjd = Time(stop).mjd
        mask = (vis["Time(MJD_UTC)"] >= start_mjd) & (vis["Time(MJD_UTC)"] <= stop_mjd)
    else:
        # Modern: Datetime-based filtering (more precise at boundaries)
        if "Time_UTC" in vis.columns:
            parsed_col = "_Time_UTC_dt64"
            time_utc = vis["Time_UTC"]
            if pd.api.types.is_datetime64_any_dtype(time_utc):
                vis_times = time_utc
            else:
                if parsed_col in vis.columns and pd.api.types.is_datetime64_any_dtype(
                    vis[parsed_col]
                ):
                    vis_times = vis[parsed_col]
                else:
                    # Parse once; store on the cached DataFrame to reuse in subsequent calls.
                    vis[parsed_col] = pd.to_datetime(time_utc, errors="coerce")
                    vis_times = vis[parsed_col]
        else:
            # Fallback to MJD conversion
            vis_times = Time(
                vis["Time(MJD_UTC)"].to_numpy(), format="mjd", scale="utc"
            ).to_datetime()
            vis_times = pd.to_datetime(vis_times)
        mask = (vis_times >= start) & (vis_times <= stop)

    return vis.loc[mask]


@dataclass(frozen=True)
class SchedulerPaths:
    """Filesystem layout for a run."""

    package_dir: Path
    data_dir: Path
    targets_dir: Path
    aux_targets_dir: Path
    baseline_dir: Path

    @classmethod
    def from_package_root(
        cls, package_dir: Path, data_dir_name: str = "data"
    ) -> "SchedulerPaths":
        data_dir = package_dir / data_dir_name
        return cls(
            package_dir=package_dir,
            data_dir=data_dir,
            targets_dir=data_dir / "targets",
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )


@dataclass
class AuxiliaryObservationStats:
    total_time: timedelta = timedelta()
    last_priority: float = 0.0


@dataclass
class SchedulerState:
    tracker: pd.DataFrame
    all_target_obs_time: Dict[str, timedelta] = field(default_factory=dict)
    non_primary_obs_time: Dict[str, "AuxiliaryObservationStats"] = field(
        default_factory=dict
    )
    last_std_obs: datetime = datetime(2025, 12, 1)


@dataclass(frozen=True)
class SchedulerInputs:
    pandora_start: datetime
    pandora_stop: datetime
    sched_start: datetime
    sched_stop: datetime
    target_list: pd.DataFrame
    paths: SchedulerPaths
    target_definition_files: list[str]
    primary_target_csv: Path
    auxiliary_target_csv: Path
    occultation_target_csv: Path
    output_dir: Path
    tracker_pickle_path: Optional[Path] = None


@dataclass
class SchedulerOutputs:
    schedule: pd.DataFrame
    tracker: pd.DataFrame
    observation_report_path: Optional[Path]
    schedule_path: Optional[Path]
    tracker_csv_path: Optional[Path]
    tracker_pickle_path: Optional[Path]


def run_scheduler(
    inputs: SchedulerInputs, config: PandoraSchedulerConfig
) -> SchedulerOutputs:
    """Execute the scheduling loop to mirror the legacy Schedule function."""

    commissioning_offset = timedelta(days=config.commissioning_days)
    pandora_start = inputs.pandora_start + commissioning_offset
    pandora_stop = inputs.pandora_stop
    sched_start = inputs.sched_start + commissioning_offset
    sched_stop = inputs.sched_stop

    state = SchedulerState(
        tracker=_initialize_tracker(
            inputs,
            config,
            pandora_start,
            pandora_stop,
            sched_start,
            sched_stop,
        ),
        all_target_obs_time={},
        non_primary_obs_time={},
        last_std_obs=datetime(2025, 12, 1),
    )

    inputs.output_dir.mkdir(parents=True, exist_ok=True)

    schedule_rows: list[pd.DataFrame] = []

    start = sched_start
    stop = start + config.schedule_step

    progress_bar = None
    previous_start = start
    if config.show_progress:
        total_minutes = max((sched_stop - sched_start).total_seconds() / 60.0, 1.0)
        progress_bar = tqdm(
            total=total_minutes,
            desc="Scheduling",
            unit="min",
            dynamic_ncols=True,
        )

    def advance_progress(new_start: datetime) -> None:
        nonlocal previous_start
        if not progress_bar:
            previous_start = new_start
            return
        delta_minutes = (new_start - previous_start).total_seconds() / 60.0
        if delta_minutes > 0:
            remaining = progress_bar.total - progress_bar.n
            progress_bar.update(min(delta_minutes, remaining))
        previous_start = new_start

    too_targets, too_starts, too_stops = _load_too_table(inputs.paths.data_dir)

    # Pre-load all transit windows once
    all_planet_names = state.tracker["Planet Name"].dropna().unique()
    transit_windows = _load_planet_transit_windows(pd.Index(all_planet_names), inputs)

    while stop <= sched_stop:
        state.tracker = state.tracker.sort_values(
            by=["Primary Target", "Transit Priority"],
            ascending=[False, True],
        ).reset_index(drop=True)
        obs_range = pd.date_range(start, stop, freq="min")
        temp_df = observation_utils.check_if_transits_in_obs_window(
            state.tracker,
            pd.DataFrame(
                [],
                columns=[
                    "Planet Name",
                    "Obs Start",
                    "Obs Gap Time",
                    "Visit Duration",
                    "Transit Coverage",
                    "SAA Overlap",
                    "Schedule Factor",
                    "Transit Factor",
                    "Quality Factor",
                    "Comments",
                ],
            ),
            inputs.target_list,
            start,
            pandora_start,
            pandora_stop,
            sched_start,
            sched_stop,
            obs_range,
            list(config.transit_scheduling_weights),
            config.transit_coverage_min,
            inputs.paths.targets_dir,
            short_visit_threshold_hours=config.short_visit_threshold_hours,
            short_visit_edge_buffer_hours=config.short_visit_edge_buffer_hours,
            long_visit_edge_buffer_hours=config.long_visit_edge_buffer_hours,
        )

        too_result = _handle_targets_of_opportunity(
            start,
            stop,
            obs_range,
            too_targets,
            too_starts,
            too_stops,
            state,
            inputs,
            config,
            transit_windows,
        )
        if too_result is not None:
            too_df, new_start = too_result
            if not too_df.empty:
                schedule_rows.append(too_df)
            start = new_start
            stop = start + config.schedule_step
            advance_progress(start)
            continue

        if temp_df.empty:
            aux_df, log_info = _schedule_auxiliary_target(
                start,
                stop,
                config,
                state,
                inputs,
            )
            if not aux_df.empty:
                schedule_rows.append(aux_df)
            logger.info(f"{log_info}; window {start} to {stop}")
            start = stop
            stop = start + config.schedule_step
            advance_progress(start)
            continue

        scheduled_visit = _schedule_primary_target(
            temp_df,
            state,
            inputs,
            config,
            start,
            obs_range,
        )
        schedule_rows.append(scheduled_visit)
        start = scheduled_visit["Observation Stop"].iloc[-1]
        stop = start + config.schedule_step
        advance_progress(start)

    schedule = (
        pd.concat(schedule_rows, ignore_index=True)
        if schedule_rows
        else pd.DataFrame(
            columns=[
                "Target",
                "Observation Start",
                "Observation Stop",
                "RA",
                "DEC",
                "Transit Coverage",
                "SAA Overlap",
                "Schedule Factor",
                "Quality Factor",
                "Comments",
            ]
        )
    )

    schedule = schedule.sort_values(by=["Observation Start"]).reset_index(drop=True)
    column_order = [
        "Target",
        "Observation Start",
        "Observation Stop",
        "RA",
        "DEC",
        "Transit Coverage",
        "SAA Overlap",
        "Schedule Factor",
        "Quality Factor",
        "Comments",
    ]
    schedule = schedule.reindex(columns=column_order)

    observation_report_path = _write_observation_report(
        inputs,
        state,
        pandora_start,
    )
    schedule_path, tracker_csv_path, tracker_pickle_path = _persist_outputs(
        schedule,
        state.tracker,
        inputs,
        config,
        pandora_start,
        pandora_stop,
    )

    if progress_bar:
        advance_progress(sched_stop)
        progress_bar.close()

    return SchedulerOutputs(
        schedule=schedule,
        tracker=state.tracker,
        observation_report_path=observation_report_path,
        schedule_path=schedule_path,
        tracker_csv_path=tracker_csv_path,
        tracker_pickle_path=tracker_pickle_path,
    )


def _initialize_tracker(
    inputs: SchedulerInputs,
    config: PandoraSchedulerConfig,
    pandora_start: datetime,
    pandora_stop: datetime,
    sched_start: datetime,
    sched_stop: datetime,
) -> pd.DataFrame:
    target_list = inputs.target_list.reset_index(drop=True)

    tracker = pd.DataFrame(
        {
            "Planet Name": target_list["Planet Name"].to_numpy(),
            "Primary Target": target_list["Primary Target"].to_numpy(),
            "RA": target_list["RA"].to_numpy(),
            "DEC": target_list["DEC"].to_numpy(),
            "Transits Needed": target_list["Number of Transits to Capture"].to_numpy(),
        }
    )
    tracker["Transits Acquired"] = np.zeros(len(target_list), dtype=float)

    archive_path = inputs.paths.data_dir / "Pandora_archive.csv"
    if archive_path.exists():
        try:
            archive = pd.read_csv(archive_path)
        except Exception:
            archive = None
        if archive is not None:
            for _, row in archive.iterrows():
                mask = tracker["Planet Name"] == row["Target"]
                tracker.loc[mask, "Transits Needed"] -= 1
                tracker.loc[mask, "Transits Acquired"] += 1

    transits_left_lifetime: list[int] = []
    transits_left_schedule: list[int] = []

    for _, row in target_list.iterrows():
        planet_name = str(row["Planet Name"])
        star_name = str(row["Star Name"])

        visibility_path = build_visibility_path(
            inputs.paths.targets_dir, star_name, planet_name
        )
        planet_data = None
        try:
            planet_data = read_parquet_cached(
                str(visibility_path),
                columns=["Transit_Start", "Transit_Stop", "Transit_Coverage"],
            )
        except Exception:
            planet_data = None

        if planet_data is None:
            raise FileNotFoundError(
                "Visibility file missing or unreadable for planet %s (expected at %s)"
                % (planet_name, visibility_path)
            )

        planet_data = planet_data.drop(
            planet_data.index[
                planet_data["Transit_Coverage"] < config.transit_coverage_min
            ]
        ).reset_index(drop=True)

        # If no transits remain after filtering, skip this planet
        if planet_data.empty or len(planet_data) == 0:
            transits_left_lifetime.append(0)
            transits_left_schedule.append(0)
            continue

        # Use pre-converted datetime if available (performance optimization)
        if "Transit_Start_UTC" in planet_data.columns:
            try:
                start_transits = pd.to_datetime(
                    planet_data["Transit_Start_UTC"]
                ).to_numpy()
            except Exception as e:
                raise ValueError(
                    f"Failed to parse Transit_Start_UTC for {planet_name}: {e}"
                )
        else:
            # Fallback to MJD conversion for backward compatibility
            try:
                start_transits = pd.to_datetime(
                    Time(
                        planet_data["Transit_Start"], format="mjd", scale="utc"
                    ).to_datetime()
                ).to_numpy()
            except Exception as e:
                raise ValueError(
                    f"Failed to convert Transit_Start MJD for {planet_name}: {e}"
                )

        if "Transit_Stop_UTC" in planet_data.columns:
            try:
                end_transits = pd.to_datetime(
                    planet_data["Transit_Stop_UTC"]
                ).to_numpy()
            except Exception as e:
                raise ValueError(
                    f"Failed to parse Transit_Stop_UTC for {planet_name}: {e}"
                )
        else:
            # Fallback to MJD conversion for backward compatibility
            try:
                end_transits = pd.to_datetime(
                    Time(
                        planet_data["Transit_Stop"], format="mjd", scale="utc"
                    ).to_datetime()
                ).to_numpy()
            except Exception as e:
                raise ValueError(
                    f"Failed to convert Transit_Stop MJD for {planet_name}: {e}"
                )

        # Ensure pandora_start/stop are compatible with numpy datetime64 arrays
        p_start = pd.to_datetime(pandora_start).to_numpy()
        p_stop = pd.to_datetime(pandora_stop).to_numpy()
        s_start = pd.to_datetime(sched_start).to_numpy()
        s_stop = pd.to_datetime(sched_stop).to_numpy()

        lifetime_mask = (p_start <= start_transits) & (end_transits <= p_stop)
        schedule_mask = (s_start <= start_transits) & (end_transits <= s_stop)

        transits_left_lifetime.append(int(np.count_nonzero(lifetime_mask)))
        transits_left_schedule.append(int(np.count_nonzero(schedule_mask)))

    tracker["Transits Left in Lifetime"] = transits_left_lifetime
    tracker["Transits Left in Schedule"] = transits_left_schedule
    tracker["Transit Priority"] = (
        tracker["Transits Left in Lifetime"] - tracker["Transits Needed"]
    )

    return tracker


def _write_observation_report(
    inputs: SchedulerInputs, state: SchedulerState, pandora_start: datetime
) -> Path:
    report_name = f"Observation_Time_Report_{pandora_start}.csv"
    report_path = inputs.output_dir / report_name
    requested_hours_catalogs: list[Path] = []
    # Prefer the canonical per-category manifests written to the output data dir.
    for target_def in inputs.target_definition_files[1:]:
        manifest_path = inputs.paths.data_dir / f"{target_def}_targets.csv"
        if manifest_path.exists():
            requested_hours_catalogs.append(manifest_path)

    observation_utils.save_observation_time_report(
        state.all_target_obs_time,
        inputs.target_list,
        str(report_path),
        requested_hours_catalogs=requested_hours_catalogs,
    )
    return report_path


def _persist_outputs(
    schedule: pd.DataFrame,
    tracker: pd.DataFrame,
    inputs: SchedulerInputs,
    config: PandoraSchedulerConfig,
    pandora_start: datetime,
    pandora_stop: datetime,
) -> tuple[Path, Path, Path]:
    schedule_name = (
        f"Pandora_Schedule_{config.transit_scheduling_weights[0]}_"
        f"{config.transit_scheduling_weights[1]}_{config.transit_scheduling_weights[2]}_"
        f"{pandora_start.strftime('%Y-%m-%d')}_to_{pandora_stop.strftime('%Y-%m-%d')}.csv"
    )
    schedule_path = inputs.output_dir / schedule_name
    schedule.to_csv(schedule_path, index=False)

    tracker_csv_path = inputs.output_dir / "tracker.csv"
    tracker.to_csv(tracker_csv_path, index=False)

    if inputs.tracker_pickle_path:
        tracker_pickle_path = inputs.tracker_pickle_path
    else:
        start_name = pandora_start.strftime("%Y-%m-%d")
        stop_name = pandora_stop.strftime("%Y-%m-%d")
        filename = f"Tracker_{start_name}_to_{stop_name}.pkl"
        tracker_pickle_path = inputs.output_dir / filename
    with tracker_pickle_path.open("wb") as handle:
        pickle.dump(tracker, handle)

    return schedule_path, tracker_csv_path, tracker_pickle_path


def _load_too_table(data_dir: Path) -> tuple[list[str], list[datetime], list[datetime]]:
    too_path = data_dir / "ToO_list.csv"
    if not too_path.exists():
        return [], [], []

    table = pd.read_csv(too_path)
    return (
        table["Target"].tolist(),
        [
            datetime.strptime(value, "%Y-%m-%d %H:%M:%S")
            for value in table["Obs Window Start"]
        ],
        [
            datetime.strptime(value, "%Y-%m-%d %H:%M:%S")
            for value in table["Obs Window Stop"]
        ],
    )


def _load_planet_transit_windows(
    planet_names: pd.Index,
    inputs: SchedulerInputs,
) -> dict[str, tuple[datetime, datetime]]:
    """Batch-load transit windows for multiple planets to avoid repeated CSV reads."""
    transit_windows = {}

    # Build planet -> star name lookup from target list
    planet_to_star = dict(
        zip(inputs.target_list["Planet Name"], inputs.target_list["Star Name"])
    )

    for planet_name in planet_names:
        planet_str = str(planet_name)
        star_name = planet_to_star.get(planet_str)
        if star_name is None:
            raise ValueError(
                f"Planet {planet_str} not found in target list. "
                f"Cannot determine star name."
            )
        star_name = str(star_name)

        vis_path = build_visibility_path(
            inputs.paths.targets_dir, star_name, planet_str
        )
        if not vis_path.is_file():
            raise FileNotFoundError(
                f"Planet visibility file not found: {vis_path}\n"
                f"  Planet: {planet_str}, Star: {star_name}\n"
                f"  Expected at: {vis_path}"
            )

        planet_visibility = read_parquet_cached(
            str(vis_path),
            columns=["Transit_Start", "Transit_Stop"],
        )

        if planet_visibility is None:
            raise ValueError(
                f"Planet visibility file exists but is unreadable: {vis_path}\n"
                f"  Planet: {planet_str}, Star: {star_name}"
            )

        # Skip planets with no transits in the scheduling window (valid case)
        if planet_visibility.empty or len(planet_visibility) == 0:
            continue

        start_transit = Time(
            planet_visibility["Transit_Start"].iloc[-1], format="mjd", scale="utc"
        ).to_datetime()
        end_transit = Time(
            planet_visibility["Transit_Stop"].iloc[-1], format="mjd", scale="utc"
        ).to_datetime()

        start_transit = start_transit.replace(second=0, microsecond=0)
        end_transit = end_transit.replace(second=0, microsecond=0)

        transit_windows[str(planet_name)] = (start_transit, end_transit)

    return transit_windows


def _handle_targets_of_opportunity(
    start: datetime,
    stop: datetime,
    obs_range: pd.DatetimeIndex,
    targets: list[str],
    starts: list[datetime],
    stops: list[datetime],
    state: SchedulerState,
    inputs: SchedulerInputs,
    config: PandoraSchedulerConfig,
    transit_windows: dict[str, tuple[datetime, datetime]],
) -> Optional[tuple[pd.DataFrame, datetime]]:
    if not targets:
        return None

    overlap = obs_range.intersection(pd.DatetimeIndex(starts))
    if len(overlap) == 0:
        return None

    first_overlap = overlap[0].to_pydatetime()
    idx = starts.index(first_overlap)
    obs_start = starts[idx]
    obs_stop = stops[idx]
    target_name = targets[idx]

    logger.info("Attempting to schedule Target of Opportunity")

    tracker = state.tracker
    positive_needed = tracker["Transits Needed"] > 0
    transit_ratio = pd.Series(np.inf, index=tracker.index)
    transit_ratio.loc[positive_needed] = (
        tracker.loc[positive_needed, "Transits Left in Lifetime"]
        / tracker.loc[positive_needed, "Transits Needed"]
    )
    critical_planets = tracker.loc[positive_needed & (transit_ratio <= _TOO_OVERRIDE_THRESHOLD)].copy()

    schedule_parts: list[pd.DataFrame] = []
    forced_observation = False

    # Batch-load all transit windows once
    # all_active_names = tracker.loc[positive_needed].index
    # transit_windows = _load_planet_transit_windows(all_active_names, inputs)

    for planet_name in critical_planets.index:
        if planet_name not in transit_windows:
            continue

        start_transit, end_transit = transit_windows[planet_name]

        # Use per-target visit duration and edge buffer for ToO overlap check
        planet_name_str = str(planet_name)
        visit_duration = observation_utils.get_target_visit_duration(
            planet_name_str, inputs.target_list
        )
        edge_buffer = observation_utils.compute_edge_buffer(
            visit_duration,
            config.short_visit_threshold_hours,
            config.short_visit_edge_buffer_hours,
            config.long_visit_edge_buffer_hours,
        )
        early_start, late_start = observation_utils.compute_transit_start_bounds(
            start_transit, end_transit, visit_duration, edge_buffer
        )

        # Skip if no valid start window
        if early_start > late_start:
            continue

        start_range = pd.date_range(early_start, late_start, freq="min")
        overlap_times = obs_range.intersection(start_range)

        if len(overlap_times) == 0:
            continue

        forced_observation = True
        forced_start = overlap_times[0].to_pydatetime()

        if obs_range[0].to_pydatetime() < forced_start:
            free_df = pd.DataFrame(
                [
                    [
                        "FREE PRE-TOO, REPLACE WITH AUX",
                        obs_range[0].to_pydatetime(),
                        forced_start,
                    ]
                ],
                columns=["Target", "Observation Start", "Observation Stop"],
            )
            schedule_parts.append(free_df)

        forced_df = pd.DataFrame(
            [
                [
                    planet_name,
                    forced_start,
                    obs_stop,
                    f"{planet_name} forced over ToO due to transit factor <=1",
                ]
            ],
            columns=["Target", "Observation Start", "Observation Stop", "Comments"],
        )
        schedule_parts.append(forced_df)
        logger.info(
            "Forced observation of %s over ToO due to critical transit factor",
            planet_name,
        )
        break

    if forced_observation:
        combined = (
            pd.concat(schedule_parts, ignore_index=True)
            if schedule_parts
            else pd.DataFrame()
        )
        return combined, obs_stop

    tf_warning_messages: list[str] = []
    active_planets = tracker.loc[positive_needed]

    # Reuse the already-loaded transit windows
    for planet_name in active_planets.index:
        if planet_name not in transit_windows:
            continue

        tf_ratio = transit_ratio.loc[planet_name]
        start_transit, end_transit = transit_windows[planet_name]

        # Use per-target visit duration and edge buffer
        planet_name_str = str(planet_name)
        visit_duration = observation_utils.get_target_visit_duration(
            planet_name_str, inputs.target_list
        )
        edge_buffer = observation_utils.compute_edge_buffer(
            visit_duration,
            config.short_visit_threshold_hours,
            config.short_visit_edge_buffer_hours,
            config.long_visit_edge_buffer_hours,
        )
        early_start, late_start = observation_utils.compute_transit_start_bounds(
            start_transit, end_transit, visit_duration, edge_buffer
        )

        # Skip if no valid start window
        if early_start > late_start:
            continue

        start_range = pd.date_range(early_start, late_start, freq="min")
        overlap_times = obs_range.intersection(start_range)

        if len(overlap_times) > 0 and tf_ratio > 1:
            tf_warning_messages.append(
                f"Warning: {planet_name} has MTRM > 1 and is transiting during ToO."
            )

    if obs_range[0].to_pydatetime() < obs_start:
        free_df = pd.DataFrame(
            [
                [
                    "FREE PRE-TOO, REPLACE WITH AUX",
                    obs_range[0].to_pydatetime(),
                    obs_start,
                ]
            ],
            columns=["Target", "Observation Start", "Observation Stop"],
        )
        schedule_parts.append(free_df)

    too_df = pd.DataFrame(
        [
            [
                target_name,
                obs_start,
                obs_stop,
                " ".join(tf_warning_messages).strip(),
            ]
        ],
        columns=["Target", "Observation Start", "Observation Stop", "Comments"],
    )
    schedule_parts.append(too_df)

    logger.info(f"Scheduled Target of Opportunity: {target_name}")
    for message in tf_warning_messages:
        logger.warning(message)

    combined = pd.concat(schedule_parts, ignore_index=True)
    return combined, obs_stop


def _schedule_auxiliary_target(
    start: datetime,
    stop: datetime,
    config: PandoraSchedulerConfig,
    state: SchedulerState,
    inputs: SchedulerInputs,
) -> tuple[pd.DataFrame, str]:
    # `obs_range` was previously created here but not used; remove to satisfy lint
    active_start = start
    scheduled_rows: list[list] = []
    row_columns = ["Target", "Observation Start", "Observation Stop", "RA", "DEC", "Comments"]
    std_visible_duration_by_target: dict[str, timedelta] = {}

    def _accumulate_observation_time(result_df: pd.DataFrame) -> None:
        for record in result_df.to_dict(orient="records"):
            target_label = str(record["Target"])
            if target_label == "Free Time":
                continue
            duration = record["Observation Stop"] - record["Observation Start"]
            if isinstance(duration, pd.Timedelta):
                duration = duration.to_pytimedelta()
            if target_label.endswith(" STD") and target_label in std_visible_duration_by_target:
                duration = std_visible_duration_by_target[target_label]
            state.all_target_obs_time[target_label] = (
                state.all_target_obs_time.get(target_label, timedelta()) + duration
            )

    if config.primary_only_mode:
        result = pd.DataFrame(
            [["Free Time", start, stop, float("nan"), float("nan"), ""]],
            columns=row_columns,
        )
        return result, "Primary-only mode enabled, skipping non-primary scheduling."

    obs_std_duration = timedelta(hours=config.std_obs_duration_hours)
    if (
        active_start - state.last_std_obs
        > timedelta(days=config.std_obs_frequency_days)
        and (stop - active_start).total_seconds() / 60.0 > config.effective_min_science_sequence_minutes
    ):
        std_path = inputs.paths.data_dir / "monitoring-standard_targets.csv"
        std_df = read_csv_cached(str(std_path))
        if std_df is not None:
            std_df = std_df.sort_values("Priority", ascending=False, ignore_index=True)
            std_records = std_df.to_dict(orient="records")
        else:
            std_records = []

        # Pre-load visibility DataFrames for all standard stars so they can
        # be re-used when we slide the search window.
        std_vis_cache: list[tuple[dict, pd.DataFrame]] = []
        for std_row in std_records:
            std_name = str(std_row["Star Name"])
            vis_file = build_star_visibility_path(
                inputs.paths.aux_targets_dir, std_name
            )
            if not vis_file.is_file():
                raise FileNotFoundError(
                    f"Standard target visibility file not found: {vis_file}\n"
                    f"  Target: {std_name}\n"
                    f"  Expected at: {vis_file}"
                )

            vis = read_parquet_cached(
                str(vis_file),
                columns=["Time(MJD_UTC)", "Time_UTC", "Visible"],
            )
            if vis is None:
                vis = read_parquet_cached(
                    str(vis_file),
                    columns=["Time(MJD_UTC)", "Visible"],
                )

            if vis is None or vis.empty or len(vis) == 0:
                raise ValueError(
                    f"Standard target visibility file exists but is empty: {vis_file}\n"
                    f"  Target: {std_name}"
                )
            std_vis_cache.append((std_row, vis))

        # ── Adaptive-duration + sliding-window STD search ──────────
        # Try the configured duration first, then progressively shorter
        # durations (down to the effective science minimum). For each duration,
        # slide through the gap in 1-minute increments looking for a
        # window where at least one standard star has full visibility.
        #
        # When the best placement starts after the gap beginning, the
        # STD row is extended backwards to fill the pre-STD time: the
        # science-calendar builder will emit occultation-target
        # observations during the non-visible prefix and STD science
        # observations during the visible portion.  This avoids
        # inserting Free Time into the schedule.
        std_candidate: Optional[tuple[str, float, float, float]] = None
        best_offset: timedelta = timedelta(0)
        best_duration: timedelta = obs_std_duration
        min_std_td = timedelta(minutes=config.effective_min_science_sequence_minutes)
        slide_step = timedelta(minutes=1)

        trial_duration = obs_std_duration
        while trial_duration >= min_std_td and std_candidate is None:
            latest_start = stop - trial_duration
            trial_offset = timedelta(0)
            while active_start + trial_offset <= latest_start:
                # Skip offsets that would create a pre-STD segment
                # shorter than the science minimum: the calendar
                # builder would emit a too-short occultation sequence.
                # After trying offset=0, jump straight to min_std_td.
                if timedelta(0) < trial_offset < min_std_td:
                    trial_offset = min_std_td
                    continue
                trial_start = active_start + trial_offset
                trial_stop = trial_start + trial_duration
                for std_row, vis in std_vis_cache:
                    vis_filtered = _filter_visibility_by_time(
                        vis,
                        trial_start,
                        trial_stop,
                        use_legacy_mode=config.use_legacy_mode,
                    )
                    if not vis_filtered.empty and vis_filtered["Visible"].all():
                        std_candidate = (
                            str(std_row["Star Name"]),
                            float(std_row["RA"]),
                            float(std_row["DEC"]),
                            float(std_row["Priority"]),
                        )
                        best_offset = trial_offset
                        best_duration = trial_duration
                        if trial_duration < obs_std_duration:
                            logger.warning(
                                "%s scheduled for STD observations "
                                "with SHORTENED duration "
                                "(%d min instead of %d min, "
                                "offset=+%d min into gap)",
                                std_row["Star Name"],
                                int(trial_duration.total_seconds() / 60),
                                int(obs_std_duration.total_seconds() / 60),
                                int(trial_offset.total_seconds() / 60),
                            )
                        elif trial_offset > timedelta(0):
                            logger.info(
                                "%s scheduled for STD observations "
                                "(offset=+%d min into gap)",
                                std_row["Star Name"],
                                int(trial_offset.total_seconds() / 60),
                            )
                        else:
                            logger.info(
                                "%s scheduled for STD observations "
                                "with full visibility",
                                std_row["Star Name"],
                            )
                        break

                if std_candidate is not None:
                    break
                trial_offset += slide_step

            if std_candidate is None:
                trial_duration -= slide_step

        if std_candidate is not None:
            std_name, std_ra, std_dec, priority_std = std_candidate

            # The STD row always starts at active_start.  When the
            # visible window begins later (best_offset > 0) the
            # science-calendar builder handles the non-visible prefix
            # via occultation-target observations.
            std_visible_end = active_start + best_offset + best_duration

            # If the post-STD remainder is shorter than the science minimum,
            # extend the STD row to fill the entire gap so no Free Time
            # or too-short auxiliary slot is created.
            post_remainder = (stop - std_visible_end).total_seconds() / 60.0
            if 0 < post_remainder <= config.effective_min_science_sequence_minutes:
                logger.warning(
                    "Extending STD row to absorb short post-STD "
                    "remainder (%.1f min < %d min minimum)",
                    post_remainder,
                    config.effective_min_science_sequence_minutes,
                )
                std_row_end = stop
            else:
                std_row_end = std_visible_end

            std_label = f"{std_name} STD"
            scheduled_rows.append(
                [std_label, active_start, std_row_end, std_ra, std_dec, ""]
            )
            std_visible_duration_by_target[std_label] = (
                std_visible_duration_by_target.get(std_label, timedelta())
                + best_duration
            )
            active_start = std_row_end
            state.last_std_obs = active_start

            stats = state.non_primary_obs_time.setdefault(
                "STD", AuxiliaryObservationStats()
            )
            stats.total_time += best_duration  # visible portion only
            stats.last_priority = priority_std
        else:
            # No standard star visible at any duration or position in
            # this gap.  Skip the STD entirely so the auxiliary target
            # gets the full window — do NOT emit a placeholder row that
            # would consume schedule time.
            logger.warning(
                "No visible standard star between %s and %s; "
                "skipping STD, full gap available for auxiliary target",
                start,
                stop,
            )

    # Check if any window remains after STD for auxiliary observation.
    # Even short windows are worth attempting; only fall back to Free Time if
    # no non-primary target can actually be scheduled.
    remaining_minutes = (stop - active_start).total_seconds() / 60.0
    if remaining_minutes <= 0:
        # Gap fully consumed (e.g. extended STD) — nothing more to schedule
        result = pd.DataFrame(scheduled_rows, columns=row_columns)
        _accumulate_observation_time(result)
        return result, "Gap fully scheduled by STD observation."

    if config.aux_sort_key is None:
        scheduled_rows.append(
            ["Free Time", active_start, stop, float("nan"), float("nan"), ""]
        )
        result = pd.DataFrame(scheduled_rows, columns=row_columns)
        _accumulate_observation_time(result)
        return result, "Free time, no observation scheduled."

    selected_row: Optional[list] = None
    priority_val = 0.0
    log_info = "No fuly or partially visible non-primary targets, Free Time..."
    selected_requested_hours: float | None = None
    selected_observed_hours: float | None = None
    fallback_over_requested_used = False
    non_primary_priorities = {
        name: stats.last_priority
        for name, stats in state.non_primary_obs_time.items()
        if not name.endswith("STD")
    }

    for allow_over_requested in (False, True):
        if selected_row is not None:
            break
        if allow_over_requested:
            logger.warning(
                "No eligible non-primary targets visible between %s and %s; "
                "allowing targets that have met requested hours",
                active_start,
                stop,
            )

        for target_def in inputs.target_definition_files[1:]:
            aux_csv = inputs.paths.data_dir / f"{target_def}_targets.csv"
            if not aux_csv.exists():
                continue

            aux_targets = read_csv_cached(str(aux_csv))
            if aux_targets is None:
                continue
            aux_targets = aux_targets.reset_index(drop=True)

            if "Number of Hours Requested" not in aux_targets.columns:
                raise ValueError(
                    "Auxiliary targets catalog is missing required 'Number of Hours Requested' "
                    f"column: {aux_csv}"
                )

            requested_hours = pd.to_numeric(
                aux_targets["Number of Hours Requested"], errors="coerce"
            )
            if requested_hours.isna().any():
                bad = (
                    aux_targets.loc[requested_hours.isna(), "Star Name"]
                    .astype(str)
                    .head(10)
                    .tolist()
                )
                raise ValueError(
                    "Invalid 'Number of Hours Requested' values in auxiliary targets catalog "
                    f"{aux_csv}; examples: {bad}"
                )
            if (requested_hours <= 0).any():
                bad = (
                    aux_targets.loc[
                        requested_hours <= 0, ["Star Name", "Number of Hours Requested"]
                    ]
                    .head(10)
                    .to_dict(orient="records")
                )
                raise ValueError(
                    "Non-positive 'Number of Hours Requested' values in auxiliary targets catalog "
                    f"{aux_csv}; examples: {bad}"
                )

            aux_targets = aux_targets.assign(_requested_hours=requested_hours)
            aux_targets = aux_targets.assign(
                _observed_hours=aux_targets["Star Name"]
                .astype(str)
                .map(
                    lambda label: (
                        state.non_primary_obs_time.get(
                            str(label)
                        ).total_time.total_seconds()
                        / 3600.0
                        if state.non_primary_obs_time.get(str(label)) is not None
                        else 0.0
                    )
                )
            )
            if not allow_over_requested:
                eligible = (
                    aux_targets["_observed_hours"] < aux_targets["_requested_hours"]
                )
                aux_targets = aux_targets.loc[eligible].reset_index(drop=True)

            if aux_targets.empty:
                continue

            mask = aux_targets["Star Name"].isin(non_primary_priorities.keys())
            if mask.any():
                mapped_priorities = aux_targets.loc[mask, "Star Name"].apply(
                    lambda value: non_primary_priorities.get(str(value), np.nan)
                )
                aux_targets.loc[mask, "Priority"] = mapped_priorities

            if config.aux_sort_key in {"sort_by_tdf_priority", "closest"}:
                aux_targets = aux_targets.sort_values(
                    "Priority", ascending=False, ignore_index=True
                )
            names = aux_targets["Star Name"].tolist()
            ras = aux_targets["RA"].tolist()
            decs = aux_targets["DEC"].tolist()
            priorities = aux_targets["Priority"].tolist()

            vis_all: list[int] = []
            vis_any: list[int] = []
            vis_percentages: list[float] = []

            for idx, name in enumerate(names):
                vis_file = build_star_visibility_path(
                    inputs.paths.aux_targets_dir, name
                )
                if not vis_file.is_file():
                    raise FileNotFoundError(
                        f"Auxiliary target visibility file not found: {vis_file}\n"
                        f"  Target: {name}\n"
                        f"  Expected at: {vis_file}"
                    )

                vis = read_parquet_cached(
                    str(vis_file),
                    columns=["Time(MJD_UTC)", "Time_UTC", "Visible"],
                )
                if vis is None:
                    vis = read_parquet_cached(
                        str(vis_file),
                        columns=["Time(MJD_UTC)", "Visible"],
                    )

                if vis is None or vis.empty or len(vis) == 0:
                    raise ValueError(
                        f"Auxiliary target visibility file exists but is empty: {vis_file}\n"
                        f"  Target: {name}"
                    )

                vis_filtered = _filter_visibility_by_time(
                    vis,
                    active_start,
                    stop,
                    use_legacy_mode=config.use_legacy_mode,
                )

                if not vis_filtered.empty and vis_filtered["Visible"].all():
                    vis_all.append(idx)
                    break
                if not vis_filtered.empty and vis_filtered["Visible"].any():
                    vis_any.append(idx)
                    visibility = 100 * (
                        vis_filtered["Visible"].sum() / len(vis_filtered)
                    )
                    vis_percentages.append(visibility)

            chosen_idx: int | None = None
            best_visibility: float | None = None

            if vis_all:
                if config.aux_sort_key in {"sort_by_tdf_priority", "closest"}:
                    chosen_idx = vis_all[0]
                else:
                    chosen_idx = int(np.random.randint(0, len(vis_all)))
                best_visibility = 100.0
            elif vis_any:
                any_idx = int(np.asarray(vis_percentages).argmax())
                best_visibility = vis_percentages[any_idx]
                if best_visibility >= 100 * config.min_visibility:
                    chosen_idx = vis_any[any_idx]
                else:
                    # Schedule the best available target anyway but flag it
                    chosen_idx = vis_any[any_idx]
                    logger.warning(
                        "No non-primary target with visibility >= %.2f%% from %s; "
                        "scheduling %s with %.2f%% visibility",
                        100 * config.min_visibility,
                        target_def,
                        names[vis_any[any_idx]],
                        best_visibility,
                    )

            if chosen_idx is None:
                continue

            name = names[chosen_idx]
            ra_val = ras[chosen_idx]
            dec_val = decs[chosen_idx]
            priority_val = priorities[chosen_idx]
            selected_requested_hours = float(
                aux_targets.loc[chosen_idx, "_requested_hours"]
            )
            selected_observed_hours = float(
                aux_targets.loc[chosen_idx, "_observed_hours"]
            )
            fallback_over_requested_used = (
                selected_observed_hours >= selected_requested_hours
            )

            if best_visibility is not None and best_visibility >= 100.0:
                log_info = f"{name} scheduled for non-primary observation with full visibility."
                comment = ""
            elif best_visibility is not None and best_visibility < 100 * config.min_visibility:
                log_info = (
                    f"{name} scheduled with {best_visibility:.2f}% visibility "
                    f"(below threshold {100 * config.min_visibility:.2f}%)"
                )
                comment = f"Visibility = {best_visibility:.2f}%"
            else:
                log_info = (
                    "No non-primary target with full visibility; "
                    f"{name} scheduled with {float(best_visibility or 0.0):.2f}% visibility"
                )
                comment = ""

            scheduled_rows.append([name, active_start, stop, ra_val, dec_val, comment])
            selected_row = scheduled_rows[-1]
            break

    if selected_row is None:
        scheduled_rows.append(
            ["Free Time", active_start, stop, float("nan"), float("nan"), ""]
        )
    else:
        name = selected_row[0]
        stats = state.non_primary_obs_time.setdefault(name, AuxiliaryObservationStats())
        stats.total_time += stop - active_start
        stats.last_priority = priority_val

        if (
            fallback_over_requested_used
            and selected_requested_hours is not None
            and selected_observed_hours is not None
        ):
            logger.warning(
                "Non-primary target %s has met requested observation time "
                "(observed %.2f h >= requested %.2f h) but is the only visible option; "
                "scheduling anyway",
                name,
                selected_observed_hours,
                selected_requested_hours,
            )

    result = pd.DataFrame(scheduled_rows, columns=row_columns)
    _accumulate_observation_time(result)

    return result, log_info


def _primary_transit_comment(target_list: pd.DataFrame, planet_name: str) -> str:
    """Return a comment labeling a transit as primary or secondary.

    Uses the ``Primary Target`` column if present; otherwise falls back to
    ``Number of Transits to Capture == 10`` as a heuristic for primary targets
    (matching the convention used in PR #6).
    """
    if target_list.empty or "Planet Name" not in target_list.columns:
        return ""
    match = target_list.loc[target_list["Planet Name"] == planet_name]
    if match.empty:
        return ""
    row = match.iloc[0]
    if "Primary Target" in target_list.columns:
        val = row.get("Primary Target")
        if pd.isna(val):
            return "secondary exoplanet transit"
        try:
            if bool(val):
                return "primary exoplanet transit"
        except (TypeError, ValueError):
            pass
        return "secondary exoplanet transit"
    if "Number of Transits to Capture" in target_list.columns:
        try:
            n_transits = int(row["Number of Transits to Capture"])
        except (TypeError, ValueError):
            return ""
        if n_transits >= 10:
            return "primary exoplanet transit"
        return "secondary exoplanet transit"
    return ""


def _schedule_primary_target(
    temp_df: pd.DataFrame,
    state: SchedulerState,
    inputs: SchedulerInputs,
    config: PandoraSchedulerConfig,
    start: datetime,
    obs_range: pd.DatetimeIndex,
) -> pd.DataFrame:
    if "Transit Coverage" not in temp_df.columns:
        raise ValueError("Primary candidate table is missing 'Transit Coverage'")

    # Treat (near) 100% transit coverage as a first-class priority signal.
    # We intentionally make this dominate over quality in the non-critical case.
    full_coverage = temp_df["Transit Coverage"].astype(float) >= 0.999999

    if (temp_df["Transit Factor"] <= _URGENCY_THRESHOLD).any():
        # Critical case: urgency still wins, but prefer full coverage when urgency ties.
        ranked = (
            temp_df.assign(_full_coverage=full_coverage)
            .sort_values(
                by=[
                    "Transit Factor",
                    "_full_coverage",
                    "Transit Coverage",
                    "Quality Factor",
                ],
                ascending=[True, False, False, False],
            )
            .reset_index(drop=True)
        )
    else:
        # Non-critical case: always take full-coverage transits over partial ones.
        ranked = (
            temp_df.assign(_full_coverage=full_coverage)
            .sort_values(
                by=[
                    "_full_coverage",
                    "Quality Factor",
                    "Transit Coverage",
                    "Transit Factor",
                ],
                ascending=[False, False, False, True],
            )
            .reset_index(drop=True)
        )

    if "_full_coverage" in ranked.columns:
        ranked = ranked.drop(columns=["_full_coverage"])

    first_row = ranked.iloc[0]
    planet_name = str(first_row["Planet Name"])
    ra_value = first_row["RA"]
    dec_value = first_row["DEC"]
    obs_start = pd.Timestamp(first_row["Obs Start"]).to_pydatetime()

    # Use per-target visit duration from candidate (already computed in check_if_transits_in_obs_window)
    visit_duration = first_row.get("Visit Duration")
    if visit_duration is None or pd.isna(visit_duration):
        raise ValueError(
            f"Missing Visit Duration for scheduled target '{planet_name}'. "
            "Expected per-target 'Obs Window (hrs)' to be present and valid."
        )
    elif isinstance(visit_duration, pd.Timedelta):
        visit_duration = visit_duration.to_pytimedelta()

    obs_stop = obs_start + visit_duration
    trans_cover = first_row["Transit Coverage"]
    saa_cover = first_row["SAA Overlap"]
    s_factor = first_row["Schedule Factor"]
    q_factor = first_row["Quality Factor"]

    dfs: list[pd.DataFrame] = []

    if obs_range[0].to_pydatetime() < obs_start:
        gap = obs_start - start
        gap_minutes = gap.total_seconds() / 60.0
        aux_df, aux_log = _schedule_auxiliary_target(
            start,
            obs_start,
            config,
            state,
            inputs,
        )
        if not aux_df.empty:
            dfs.append(aux_df)
        logger.info(f"{aux_log}; window {start} to {obs_start}")

    transit_comment = _primary_transit_comment(inputs.target_list, planet_name)
    main_schedule = pd.DataFrame(
        [
            [
                planet_name,
                obs_start,
                obs_stop,
                ra_value,
                dec_value,
                trans_cover,
                saa_cover,
                s_factor,
                q_factor,
                transit_comment,
            ]
        ],
        columns=[
            "Target",
            "Observation Start",
            "Observation Stop",
            "RA",
            "DEC",
            "Transit Coverage",
            "SAA Overlap",
            "Schedule Factor",
            "Quality Factor",
            "Comments",
        ],
    )
    dfs.append(main_schedule)

    state.all_target_obs_time[planet_name] = state.all_target_obs_time.get(
        planet_name, timedelta()
    ) + (obs_stop - obs_start)

    tracker_mask = state.tracker["Planet Name"] == planet_name
    state.tracker.loc[tracker_mask, "Transits Needed"] -= 1
    state.tracker.loc[tracker_mask, "Transits Acquired"] += 1
    state.tracker.loc[tracker_mask, "Transits Needed"] = state.tracker.loc[
        tracker_mask, "Transits Needed"
    ].clip(lower=0)
    state.tracker.loc[tracker_mask, "Transit Priority"] = (
        state.tracker.loc[tracker_mask, "Transits Left in Lifetime"]
        - state.tracker.loc[tracker_mask, "Transits Needed"]
    )

    acquired_value = state.tracker.loc[tracker_mask, "Transits Acquired"].iloc[0]
    acquired = int(cast(float, acquired_value))
    logger.info(
        "Scheduled transit %d of %s. Transit coverage: %.2f",
        acquired,
        planet_name,
        trans_cover,
    )

    return pd.concat(dfs, ignore_index=True)
