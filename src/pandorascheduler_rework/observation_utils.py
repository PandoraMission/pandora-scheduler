"""Utilities for scheduling observations and building observation sequences.

This module provides functions for:
- Creating observation sequence XML elements
- Processing visibility windows
- Managing target parameters and observational constraints
- Handling NIRDA and VDA payload parameters
- Building and manipulating observation schedules

This is a refactored version of the legacy helper_codes module with improved
naming and documentation while maintaining compatibility with the existing
scheduling pipeline.
"""

from __future__ import annotations

import logging
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Sequence, cast

import numpy as np
import pandas as pd
from astropy.time import Time
from pandas.errors import EmptyDataError
from tqdm import tqdm

from pandorascheduler_rework.targets.manifest import build_target_manifest
from pandorascheduler_rework.utils.io import (
    load_visibility_arrays_cached,
    read_parquet_cached,
    require_planet_visibility_file,
    resolve_star_visibility_file,
)

LOGGER = logging.getLogger(__name__)

_PLACEHOLDER_MARKERS = {"SET_BY_TARGET_DEFINITION_FILE", "SET_BY_SCHEDULER"}
_OCCULTATION_FULL_VISIBILITY_TOLERANCE_SAMPLES = 1

# Default column name for per-target visit duration
_OBS_WINDOW_COLUMN = "Obs Window (hrs)"


class TransitUnschedulableError(ValueError):
    """Raised when a transit cannot be scheduled with the required edge buffers."""

    pass


class MissingObsWindowError(ValueError):
    """Raised when a required per-target observation window is missing."""

    pass


def get_target_visit_duration(
    planet_name: str,
    target_list: pd.DataFrame,
) -> timedelta:
    """Return visit duration from the required 'Obs Window (hrs)' column.

    Parameters
    ----------
    planet_name
        Name of the planet to look up.
    target_list
        DataFrame containing target definitions with 'Planet Name' column.
    Returns
    -------
    timedelta
        The visit duration for this target.
    """
    if _OBS_WINDOW_COLUMN not in target_list.columns:
        raise MissingObsWindowError(
            f"Target list is missing required column '{_OBS_WINDOW_COLUMN}'"
        )

    mask = target_list["Planet Name"] == planet_name
    if not mask.any():
        raise MissingObsWindowError(
            f"Target '{planet_name}' not present in target list"
        )

    value = target_list.loc[mask, _OBS_WINDOW_COLUMN].iloc[0]
    if pd.isna(value):
        raise MissingObsWindowError(
            f"Target '{planet_name}' has missing '{_OBS_WINDOW_COLUMN}' value"
        )

    try:
        hours = float(value)
    except (TypeError, ValueError) as exc:
        raise MissingObsWindowError(
            f"Target '{planet_name}' has unparseable '{_OBS_WINDOW_COLUMN}' value: {value}"
        ) from exc

    if hours <= 0:
        raise MissingObsWindowError(
            f"Target '{planet_name}' has invalid '{_OBS_WINDOW_COLUMN}' value: {value}"
        )

    return timedelta(hours=hours)


def compute_edge_buffer(
    visit_duration: timedelta,
    short_visit_threshold_hours: float = 12.0,
    short_visit_edge_buffer_hours: float = 1.5,
    long_visit_edge_buffer_hours: float = 4.0,
) -> timedelta:
    """Compute the edge buffer based on visit duration.

    Parameters
    ----------
    visit_duration
        The total visit duration.
    short_visit_threshold_hours
        Visits shorter than this use the short buffer.
    short_visit_edge_buffer_hours
        Edge buffer for short visits.
    long_visit_edge_buffer_hours
        Edge buffer for long visits.

    Returns
    -------
    timedelta
        The edge buffer to use (same for pre and post transit).
    """
    visit_hours = visit_duration.total_seconds() / 3600.0
    if visit_hours < short_visit_threshold_hours:
        return timedelta(hours=short_visit_edge_buffer_hours)
    return timedelta(hours=long_visit_edge_buffer_hours)


def validate_transit_schedulable(
    planet_name: str,
    transit_duration: timedelta,
    visit_duration: timedelta,
    edge_buffer: timedelta,
) -> None:
    """Validate that a transit can fit within the visit with required buffers.

    Parameters
    ----------
    planet_name
        Name of the planet (for error messages).
    transit_duration
        Duration of the transit.
    visit_duration
        Duration of the observation visit.
    edge_buffer
        Required buffer before and after transit.

    Raises
    ------
    TransitUnschedulableError
        If the transit cannot fit with required edge buffers.
    """
    required_duration = transit_duration + 2 * edge_buffer
    if required_duration > visit_duration:
        raise TransitUnschedulableError(
            f"Transit for {planet_name} cannot be scheduled: "
            f"transit duration ({transit_duration}) + 2 × edge buffer ({edge_buffer}) "
            f"= {required_duration} exceeds visit duration ({visit_duration}). "
            f"Either increase visit duration or reduce edge buffer requirements."
        )


def compute_transit_start_bounds(
    transit_start: datetime,
    transit_stop: datetime,
    visit_duration: timedelta,
    edge_buffer: timedelta,
) -> tuple[datetime, datetime]:
    """Compute earliest and latest valid visit start times for a transit.

    The visit must:
    - Start at least `edge_buffer` before transit_start
    - End at least `edge_buffer` after transit_stop

    Parameters
    ----------
    transit_start
        Start time of the transit.
    transit_stop
        End time of the transit.
    visit_duration
        Total duration of the observation visit.
    edge_buffer
        Required buffer before and after transit.

    Returns
    -------
    tuple[datetime, datetime]
        (earliest_start, latest_start) for the visit.
        If earliest_start > latest_start, no valid start window exists.
    """
    # Latest start: transit must start at least edge_buffer after visit start
    latest_start = transit_start - edge_buffer

    # Earliest start: transit must end at least edge_buffer before visit end
    # So: visit_start + visit_duration >= transit_stop + edge_buffer
    # => visit_start >= transit_stop + edge_buffer - visit_duration
    earliest_start = transit_stop + edge_buffer - visit_duration

    return earliest_start, latest_start


def general_parameters(
    obs_sequence_duration: int = 90,
    occ_sequence_limit: int = 50,
) -> tuple[timedelta, int]:
    return obs_sequence_duration, occ_sequence_limit


# observation_sequence imported from xml module


# Array utilities moved to utils.array_ops
# remove_short_sequences - imported from utils.array_ops
# break_long_sequences - imported from utils.array_ops


# Helper functions moved to xml module


# _populate_nirda_parameters and _populate_vda_parameters imported from xml module


# _target_identifier moved to utils.string_ops (imported above)


def schedule_occultation_targets(
    v_names,
    starts,
    stops,
    visit_start,
    visit_stop,
    path,
    o_df,
    o_list,
    try_occ_targets: str,
    show_progress: bool = False,
    use_pass1: bool = True,
):
    starts_array = np.asarray(starts, dtype=float)
    stops_array = np.asarray(stops, dtype=float)

    schedule = pd.DataFrame(
        {
            "Stop": stops_array,
            "Target": pd.Series([None] * len(starts_array), dtype=object),
            "Visibility": np.nan,
        },
        index=pd.Index(starts_array, name="Start"),
    )

    if "Target" in o_df.columns:
        o_df["Target"] = o_df["Target"].astype(object)

    if "Visibility" not in o_df.columns:
        o_df["Visibility"] = np.nan
    if "Occultation Pass" not in o_df.columns:
        o_df["Occultation Pass"] = pd.Series([None] * len(o_df), dtype=object)

    if visit_start is not None and visit_stop is not None:
        description = "%s to %s: Searching for occultation target from %s" % (
            visit_start,
            visit_stop,
            try_occ_targets,
        )
    else:
        description = "Searching for occultation target from %s" % (try_occ_targets,)

    base_path = Path(path) if path is not None else None

    # Cache visibility data to avoid re-reading files in the second pass
    # Key: v_name, Value: (vis_times, visibility)
    visibility_cache: Dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def _get_visibility(name: str) -> Optional[tuple[np.ndarray, np.ndarray]]:
        if name in visibility_cache:
            return visibility_cache[name]

        f_path = resolve_star_visibility_file(base_path, name)
        if f_path is None:
            LOGGER.debug("Skipping %s; visibility file not found", name)
            return None

        data = load_visibility_arrays_cached(f_path)
        if data is None:
            LOGGER.debug("Skipping %s; visibility file unreadable", name)
            return None
        visibility_cache[name] = data
        return data

    def _interval_visible_with_tolerance(
        visibility: np.ndarray,
        interval_mask: np.ndarray,
        tolerance_samples: int = _OCCULTATION_FULL_VISIBILITY_TOLERANCE_SAMPLES,
    ) -> bool:
        sample_count = int(interval_mask.sum())
        if sample_count == 0:
            return False
        nonvisible = int((visibility[interval_mask] != 1).sum())
        return nonvisible <= tolerance_samples

    # PASS 1: Search for a single target that covers ALL intervals
    if not use_pass1:
        LOGGER.info("Pass 1 skipped (enable_occultation_pass1=False)")
    for v_name in tqdm(v_names, desc=f"{description} (Pass 1)", leave=False, disable=not show_progress or not use_pass1):
        vis_data = _get_visibility(v_name)
        if vis_data is None:
            continue

        vis_times, visibility = vis_data

        if not use_pass1:
            break

        # Check if visible for ALL intervals individually (less strict)
        all_visible = True
        has_nonzero_interval_samples = False
        for start, stop in zip(starts_array, stops_array):
            interval_mask = (vis_times >= start) & (vis_times <= stop)
            if stop > start and interval_mask.sum() > 0:
                has_nonzero_interval_samples = True
            if not _interval_visible_with_tolerance(visibility, interval_mask):
                all_visible = False
                break

        # Guard against false positives when a candidate has no data at all
        # across the substantive (non-zero duration) intervals.
        has_nonzero_interval = bool(np.any(stops_array > starts_array))
        if has_nonzero_interval and not has_nonzero_interval_samples:
            all_visible = False

        if not all_visible:
            continue
        # Apply this target to all intervals
        for idx, start in enumerate(starts_array):
            schedule.loc[start, "Target"] = v_name
            schedule.loc[start, "Visibility"] = 1

            match = o_list.loc[o_list["Star Name"] == v_name]
            if not match.empty:
                match_row = match.iloc[0]
                o_df.loc[idx, "Target"] = v_name
                o_df.loc[idx, "RA"] = match_row["RA"]
                o_df.loc[idx, "DEC"] = match_row["DEC"]
                o_df.loc[idx, "Visibility"] = 1
                o_df.loc[idx, "Occultation Pass"] = "Pass 1"

        return o_df, True

    # Pass 1 summary
    assigned_after_p1 = int(schedule["Target"].notna().sum())
    total_intervals = len(starts_array)
    LOGGER.debug(
        "Occultation Pass 1: %d/%d intervals assigned",
        assigned_after_p1,
        total_intervals,
    )

    # PASS 2: Fill gaps with multiple targets (Greedy approach)
    for v_name in tqdm(v_names, desc=f"{description} (Pass 2)", leave=False, disable=not show_progress):
        # If schedule is full, we are done
        if not schedule["Target"].isna().any():
            return o_df, True

        vis_data = _get_visibility(v_name)
        if vis_data is None:
            continue

        vis_times, visibility = vis_data

        # If the candidate has no samples across all non-zero intervals in this
        # visit, skip it (prevents empty-mask full-visibility false positives).
        has_nonzero_interval_samples = False
        for start, stop in zip(starts_array, stops_array):
            if stop <= start:
                continue
            interval_mask = (vis_times >= start) & (vis_times <= stop)
            if interval_mask.sum() > 0:
                has_nonzero_interval_samples = True
                break
        if not has_nonzero_interval_samples and bool(np.any(stops_array > starts_array)):
            continue

        for idx, (start, stop) in enumerate(zip(starts_array, stops_array)):
            if pd.isna(schedule.loc[start, "Target"]):
                interval_mask = (vis_times >= start) & (vis_times <= stop)

                # Guard: empty mask means no visibility data for this interval.
                if interval_mask.sum() == 0:
                    if pd.isna(schedule.loc[start, "Visibility"]):
                        schedule.loc[start, "Visibility"] = 0
                        o_df.loc[idx, "Visibility"] = 0
                    continue

                if _interval_visible_with_tolerance(visibility, interval_mask):
                    schedule.loc[start, "Target"] = v_name
                    schedule.loc[start, "Visibility"] = 1

                    match = o_list.loc[o_list["Star Name"] == v_name]
                    if match.empty:
                        continue
                    match_row = match.iloc[0]
                    o_df.loc[idx, "Target"] = v_name
                    o_df.loc[idx, "RA"] = match_row["RA"]
                    o_df.loc[idx, "DEC"] = match_row["DEC"]
                    o_df.loc[idx, "Visibility"] = 1
                    o_df.loc[idx, "Occultation Pass"] = "Pass 2"
                else:
                    if pd.isna(schedule.loc[start, "Visibility"]):
                        schedule.loc[start, "Visibility"] = 0
                        o_df.loc[idx, "Visibility"] = 0

    if not schedule["Target"].isna().any():
        return o_df, True

    # Pass 2 summary
    assigned_after_p2 = int(schedule["Target"].notna().sum())
    LOGGER.debug(
        "Occultation Pass 2: %d/%d intervals assigned (%d remaining)",
        assigned_after_p2,
        total_intervals,
        total_intervals - assigned_after_p2,
    )

    # PASS 3: Best-effort assignment for intervals not fully covered by any
    # single occultation target. For each remaining interval, choose the
    # candidate with the highest minute-coverage fraction (if > 0) and
    # assign it. This enables splitting long occultations across multiple
    # targets when no single candidate covers the whole interval.
    remaining_before_p3 = int(schedule["Target"].isna().sum())
    p3_progress_step = max(1, total_intervals // 10)
    LOGGER.info(
        "Occultation Pass 3: evaluating %d remaining interval(s) with best-effort coverage",
        remaining_before_p3,
    )
    for idx, (start, stop) in enumerate(zip(starts_array, stops_array)):
        if (idx + 1) % p3_progress_step == 0 or idx + 1 == total_intervals:
            remaining_now = int(schedule["Target"].isna().sum())
            LOGGER.info(
                "Occultation Pass 3 progress: %d/%d intervals checked (%d remaining)",
                idx + 1,
                total_intervals,
                remaining_now,
            )
        if not pd.isna(schedule.loc[start, "Target"]):
            continue

        best_name = None
        best_coverage = 0.0
        for v_name in v_names:
            vis = _get_visibility(v_name)
            if vis is None:
                continue
            vis_times, visibility = vis
            interval_mask = (vis_times >= start) & (vis_times <= stop)
            if interval_mask.sum() == 0:
                continue
            coverage_fraction = float((visibility[interval_mask] == 1).sum()) / float(
                interval_mask.sum()
            )
            if coverage_fraction > best_coverage:
                best_coverage = coverage_fraction
                best_name = v_name

        if best_name is not None and best_coverage > 0.0:
            schedule.loc[start, "Target"] = best_name
            schedule.loc[start, "Visibility"] = 1
            match = o_list.loc[o_list["Star Name"] == best_name]
            if match.empty:
                continue
            match_row = match.iloc[0]
            o_df.loc[idx, "Target"] = best_name
            o_df.loc[idx, "RA"] = match_row["RA"]
            o_df.loc[idx, "DEC"] = match_row["DEC"]
            o_df.loc[idx, "Visibility"] = 1
            o_df.loc[idx, "Occultation Pass"] = "Pass 3"

    # If PASS 3 produced any assignments into o_df, treat this as a valid
    # partial schedule and return it.
    assigned_after_p3 = int(schedule["Target"].notna().sum())
    LOGGER.debug(
        "Occultation Pass 3: %d/%d intervals assigned (%d remaining)",
        assigned_after_p3,
        total_intervals,
        total_intervals - assigned_after_p3,
    )
    if "Target" in o_df.columns:
        assigned_values = o_df["Target"]
        assigned_mask = assigned_values.notna() & (
            assigned_values.astype(str).str.strip() != ""
        )
        if assigned_mask.any():
            return o_df, True

    # PASS 4: For intervals still unassigned, split the interval into minute
    # resolution segments and greedily assign contiguous covered runs to the
    # candidate that covers the longest run. Build a new occ_df with one row
    # per assigned segment and return it if any assignments were made.
    result_rows: list[dict] = []
    uncovered_minutes = 0
    minute_scale = 1440.0
    remaining_before_p4 = int(schedule["Target"].isna().sum())
    p4_progress_step = max(1, total_intervals // 10)
    LOGGER.info(
        "Occultation Pass 4: minute-resolution fallback for %d remaining interval(s)",
        remaining_before_p4,
    )
    for idx, (start, stop) in enumerate(zip(starts_array, stops_array)):
        if (idx + 1) % p4_progress_step == 0 or idx + 1 == total_intervals:
            LOGGER.info(
                "Occultation Pass 4 progress: %d/%d intervals scanned (%d segment(s) built so far)",
                idx + 1,
                total_intervals,
                len(result_rows),
            )
        if not pd.isna(schedule.loc[start, "Target"]):
            # Already assigned by earlier passes
            continue

        # Build integer minute indices for the interval
        start_idx = int(np.floor(start * minute_scale))
        stop_idx = int(np.ceil(stop * minute_scale))
        if stop_idx <= start_idx:
            continue
        minutes_idx = np.arange(start_idx, stop_idx)
        if minutes_idx.size == 0:
            continue

        # Candidate coverage arrays
        candidate_coverages: dict[str, np.ndarray] = {}
        for v_name in v_names:
            vis = _get_visibility(v_name)
            if vis is None:
                continue
            vis_times, vis_flags = vis
            if vis_times.size == 0:
                continue
            vis_min_idx = np.round(vis_times * minute_scale).astype(int)
            visible_min_idx = set(vis_min_idx[vis_flags == 1])
            candidate_coverages[v_name] = np.isin(
                minutes_idx, np.fromiter(visible_min_idx, dtype=int)
            )

        i = 0
        interval_uncovered = 0
        while i < minutes_idx.size:
            # Find candidates that cover this minute
            available = [name for name, arr in candidate_coverages.items() if arr[i]]
            if not available:
                interval_uncovered += 1
                i += 1
                continue

            # For each available candidate, compute run length of consecutive True
            best_name = None
            best_len = 0
            for name in available:
                arr = candidate_coverages[name]
                # find first False after i
                tail = arr[i:]
                false_pos = np.argmax(~tail) if np.any(~tail) else tail.size
                run_len = false_pos if false_pos > 0 else tail.size
                if run_len > best_len:
                    best_len = run_len
                    best_name = name

            if best_name is None or best_len == 0:
                i += 1
                continue

            # Compute start/end mjd for this run
            run_start_idx = minutes_idx[i]
            run_end_idx = minutes_idx[i + best_len - 1] + 1
            run_start_mjd = run_start_idx / minute_scale
            run_end_mjd = run_end_idx / minute_scale

            # Format start/stop as ISO UTC strings
            try:
                start_dt = Time(run_start_mjd, format="mjd", scale="utc").to_datetime()
                stop_dt = Time(run_end_mjd, format="mjd", scale="utc").to_datetime()
                start_str = start_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
                stop_str = stop_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
            except Exception:
                start_str = str(run_start_mjd)
                stop_str = str(run_end_mjd)

            match = o_list.loc[o_list["Star Name"] == best_name]
            ra_val = (
                match.iloc[0]["RA"]
                if not match.empty and "RA" in match.columns
                else float("nan")
            )
            dec_val = (
                match.iloc[0]["DEC"]
                if not match.empty and "DEC" in match.columns
                else float("nan")
            )

            result_rows.append(
                {
                    "Target": best_name,
                    "start": start_str,
                    "stop": stop_str,
                    "RA": ra_val,
                    "DEC": dec_val,
                    "Visibility": 1.0,
                    "Occultation Pass": "Pass 4",
                }
            )

            i += best_len

        uncovered_minutes += interval_uncovered

    if result_rows:
        LOGGER.debug(
            "Occultation Pass 4: produced %d minute-resolution segments",
            len(result_rows),
        )
        if uncovered_minutes > 0:
            LOGGER.warning(
                "Occultation Pass 4: %d uncovered minute(s) detected across intervals",
                uncovered_minutes,
            )
        result_df = pd.DataFrame(result_rows)
        return result_df, True

    LOGGER.debug("Occultation Pass 4: no segments assigned — all intervals uncovered")

    mask = schedule["Target"].isna()
    schedule.loc[mask, "Target"] = "No target"
    schedule.loc[mask, "Visibility"] = 0

    o_df.loc[o_df["Target"].isna(), "Target"] = "No target"
    o_df.loc[o_df["Visibility"].isna(), "Visibility"] = 0

    return o_df, False


def save_observation_time_report(
    all_target_obs_time: Dict[str, timedelta],
    target_list: pd.DataFrame,
    output_path,
    requested_hours_catalogs: Sequence[Path] | None = None,
):
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    primary_targets = {
        str(name)
        for name in target_list.get("Planet Name", pd.Series(dtype=object)).dropna()
    }

    requested_hours_by_target: dict[str, float] = {}
    requested_hours_sources_by_target: dict[str, set[str]] = {}
    conflicting_requested_hours_targets: set[str] = set()
    if requested_hours_catalogs:
        for catalog_path in requested_hours_catalogs:
            path = Path(catalog_path)
            if not path.exists():
                raise FileNotFoundError(f"Requested-hours catalog missing: {path}")
            try:
                df = pd.read_csv(path)
            except EmptyDataError:
                # Treat a truly empty file as an empty catalog.
                continue
            if df.empty:
                continue
            if "Star Name" not in df.columns:
                raise ValueError(
                    f"Requested-hours catalog is missing 'Star Name' column: {path}"
                )
            if "Number of Hours Requested" not in df.columns:
                raise ValueError(
                    "Requested-hours catalog is missing 'Number of Hours Requested' "
                    f"column: {path}"
                )

            for _, row in df.iterrows():
                name = row.get("Star Name")
                if pd.isna(name):
                    continue
                key = str(name)

                requested_hours_sources_by_target.setdefault(key, set()).add(path.name)

                hours_req = row.get("Number of Hours Requested")
                if pd.isna(hours_req):
                    continue
                try:
                    value = float(hours_req)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid 'Number of Hours Requested' for {key!r} in {path}: {hours_req!r}"
                    ) from exc

                if (
                    key in requested_hours_by_target
                    and requested_hours_by_target[key] != value
                ):
                    chosen = max(requested_hours_by_target[key], value)
                    conflicting_requested_hours_targets.add(key)
                    requested_hours_by_target[key] = chosen
                else:
                    requested_hours_by_target[key] = value

        if conflicting_requested_hours_targets:
            example_targets = ", ".join(
                sorted(conflicting_requested_hours_targets)[:5]
            )
            LOGGER.info(
                "%d targets had conflicting 'Number of Hours Requested' values; "
                "using max values%s",
                len(conflicting_requested_hours_targets),
                f" (e.g. {example_targets})" if example_targets else "",
            )

        # Overlap is expected (e.g., a target appears in both auxiliary and
        # occultation lists).  Log at debug level for traceability.
        for target_name, source_names in requested_hours_sources_by_target.items():
            if len(source_names) > 1:
                sources_str = ", ".join(sorted(source_names))
                LOGGER.debug(
                    "Target %r appears in multiple requested-hours catalogs: %s",
                    target_name,
                    sources_str,
                )

    with output_file.open("w", encoding="utf-8") as handle:
        handle.write("Target,Is Primary,Hours Requested,Hours Scheduled,Hours Delta\n")
        for target, duration in all_target_obs_time.items():
            label = str(target)
            is_primary = "Yes" if label in primary_targets else "No"
            hours = duration.total_seconds() / 3600

            if is_primary == "Yes":
                handle.write(f"{label},{is_primary},,{hours:.2f},\n")
                continue

            if label == "STD" or label.endswith(" STD"):
                handle.write(f"{label},{is_primary},,{hours:.2f},\n")
                continue

            if not requested_hours_catalogs:
                raise ValueError(
                    "Cannot report requested-vs-scheduled hours: "
                    "requested_hours_catalogs not provided"
                )

            if label not in requested_hours_by_target:
                raise ValueError(
                    f"Missing 'Number of Hours Requested' for non-primary target {label!r}"
                )
            requested = requested_hours_by_target[label]
            delta = hours - requested
            handle.write(
                f"{label},{is_primary},{requested:.2f},{hours:.2f},{delta:.2f}\n"
            )

    return output_file


def check_if_transits_in_obs_window(
    tracker: pd.DataFrame,
    temp_df: pd.DataFrame,
    target_list: pd.DataFrame,
    start: datetime,
    pandora_start: datetime,
    pandora_stop: datetime,
    sched_start: datetime,
    sched_stop: datetime,
    obs_rng: pd.DatetimeIndex,
    transit_scheduling_weights: Sequence[float],
    transit_coverage_min: float,
    targets_dir: Path,
    *,
    short_visit_threshold_hours: float = 12.0,
    short_visit_edge_buffer_hours: float = 1.5,
    long_visit_edge_buffer_hours: float = 4.0,
):

    result_df = temp_df.copy()
    columns = [
        "Planet Name",
        "RA",
        "DEC",
        "Obs Start",
        "Obs Gap Time",
        "Visit Duration",
        "Transit Coverage",
        "SAA Overlap",
        "Schedule Factor",
        "Transit Factor",
        "Quality Factor",
        "Comments",
    ]

    planet_lookup = target_list.set_index("Planet Name")

    for i, row in tracker.iterrows():
        planet_name = row.get("Planet Name")
        if pd.isna(planet_name):
            continue

        mask = tracker["Planet Name"] == planet_name
        needed_series = cast(pd.Series, tracker.loc[mask, "Transits Needed"])
        if needed_series.empty:
            continue
        transits_needed = float(needed_series.iloc[0])

        if transits_needed == 0:
            continue

        if planet_name not in planet_lookup.index:
            LOGGER.debug("Planet %s missing from target list; skipping", planet_name)
            continue

        star_name = str(planet_lookup.loc[planet_name, "Star Name"])

        try:
            visibility_file = require_planet_visibility_file(
                targets_dir, star_name, str(planet_name)
            )
        except FileNotFoundError as e:
            LOGGER.error(str(e))
            raise

        planet_data = read_parquet_cached(
            str(visibility_file),
            columns=[
                "Transit_Start",
                "Transit_Stop",
                "Transit_Coverage",
                "SAA_Overlap",
            ],
        )
        if planet_data is None:
            raise ValueError(
                f"Planet visibility file exists but is unreadable: {visibility_file}"
            )

        # Skip planets with no transits in the scheduling window (valid case)
        if planet_data.empty:
            continue

        planet_data = planet_data.drop(
            planet_data.index[planet_data["Transit_Coverage"] < transit_coverage_min]
        ).reset_index(drop=True)
        if planet_data.empty:
            continue

        start_times = Time(
            planet_data["Transit_Start"].to_numpy(), format="mjd", scale="utc"
        ).to_value("datetime")
        stop_times = Time(
            planet_data["Transit_Stop"].to_numpy(), format="mjd", scale="utc"
        ).to_value("datetime")

        start_series = pd.Series(list(cast(Sequence[datetime], start_times)))
        stop_series = pd.Series(list(cast(Sequence[datetime], stop_times)))

        valid_mask = start_series >= start
        start_series = start_series.loc[valid_mask].reset_index(drop=True)
        stop_series = stop_series.loc[valid_mask].reset_index(drop=True)
        planet_data = planet_data.loc[valid_mask].reset_index(drop=True)

        if start_series.empty:
            continue

        lifetime_mask = (pandora_start <= start_series) & (stop_series <= pandora_stop)
        schedule_mask = (sched_start <= start_series) & (stop_series <= sched_stop)
        lifetime_count = int(lifetime_mask.sum())
        schedule_count = int(schedule_mask.sum())

        tracker.loc[mask, "Transits Left in Lifetime"] = lifetime_count
        tracker.loc[mask, "Transits Left in Schedule"] = schedule_count
        tracker.loc[mask, "Transit Priority"] = lifetime_count - transits_needed

        start_series = start_series.dt.floor("min")
        stop_series = stop_series.dt.floor("min")

        # Get per-target visit duration and compute edge buffer
        visit_duration = get_target_visit_duration(str(planet_name), target_list)
        edge_buffer = compute_edge_buffer(
            visit_duration,
            short_visit_threshold_hours,
            short_visit_edge_buffer_hours,
            long_visit_edge_buffer_hours,
        )

        # Validate that transits can be scheduled with the required edge buffers
        # Use the first transit as representative (they should all have similar duration)
        if len(start_series) > 0 and len(stop_series) > 0:
            sample_transit_duration = stop_series.iloc[0] - start_series.iloc[0]
            try:
                validate_transit_schedulable(
                    str(planet_name),
                    sample_transit_duration,
                    visit_duration,
                    edge_buffer,
                )
            except TransitUnschedulableError as e:
                LOGGER.error(str(e))
                raise

        # Compute earliest and latest valid start times for each transit
        early_start_list = []
        late_start_list = []
        for ts, te in zip(start_series, stop_series):
            earliest, latest = compute_transit_start_bounds(
                ts, te, visit_duration, edge_buffer
            )
            early_start_list.append(earliest)
            late_start_list.append(latest)
        early_start = pd.Series(early_start_list)
        late_start = pd.Series(late_start_list)

        try:
            ra_tar = float(row["RA"])
            dec_tar = float(row["DEC"])
        except (TypeError, ValueError):
            ra_tar = row.get("RA")
            dec_tar = row.get("DEC")

        transits_left = float(lifetime_count)
        if transits_needed == 0:
            continue

        coverage_values = planet_data["Transit_Coverage"].to_numpy(
            dtype=float, copy=False
        )
        saa_values = planet_data["SAA_Overlap"].to_numpy(dtype=float, copy=False)

        for j, (window_start, window_stop) in enumerate(zip(early_start, late_start)):
            if window_start > window_stop:
                continue

            start_range = pd.date_range(window_start, window_stop, freq="min")
            overlap_times = obs_rng.intersection(start_range)
            if overlap_times.empty:
                continue

            obs_start = overlap_times[0]
            gap_time = obs_start - obs_rng[0]
            # Use per-target visit duration for schedule factor normalization
            schedule_factor = 1 - (gap_time / visit_duration)
            transit_coverage = float(coverage_values[j])
            saa_overlap = float(saa_values[j])
            quality_factor = (
                (transit_scheduling_weights[0] * transit_coverage)
                + (transit_scheduling_weights[1] * (1 - saa_overlap))
                + (transit_scheduling_weights[2] * schedule_factor)
            )

            candidate = pd.DataFrame(
                [
                    [
                        planet_name,
                        ra_tar,
                        dec_tar,
                        obs_start,
                        gap_time,
                        visit_duration,
                        transit_coverage,
                        saa_overlap,
                        schedule_factor,
                        transits_left / transits_needed,
                        quality_factor,
                        np.nan,
                    ]
                ],
                columns=columns,
            )

            if result_df.empty:
                result_df = candidate
            else:
                result_df = pd.concat([result_df, candidate], ignore_index=True)

    return result_df


def process_target_files(keyword: str, *, base_path: Path) -> pd.DataFrame:
    """Return the scheduler manifest for *keyword* targets.

    Parameters
    ----------
    keyword
        Category name matching a directory in the target definition
        repository, e.g. ``"exoplanet"``.
    base_path
        Path to the target definition repository root (e.g. PandoraTargetList).
    """
    return build_target_manifest(keyword, base_path)


def create_aux_list(target_definition_files: Sequence[str], package_dir):
    """Build ``all_targets.csv`` from the provided target manifest CSVs.

    The legacy helper concatenated the partner target lists that share a common
    column set.  Reimplement the same behaviour explicitly to avoid importing
    the legacy module.
    """

    if not target_definition_files:
        raise ValueError("target_definition_files must contain at least one entry")

    base_dir = Path(package_dir)
    candidate_paths = [base_dir / f"{name}_targets.csv" for name in target_definition_files]
    if any(candidate.exists() for candidate in candidate_paths):
        data_dir = base_dir
    else:
        data_dir = base_dir / "data"
    csv_paths: List[Path] = []
    for name in target_definition_files:
        candidate = data_dir / f"{name}_targets.csv"
        if candidate.exists():
            csv_paths.append(candidate)
        else:
            LOGGER.warning("Aux list source missing: %s", candidate)

    if not csv_paths:
        raise FileNotFoundError(
            "No target definition CSVs found; unable to build all_targets.csv"
        )

    dataframes = [pd.read_csv(path) for path in csv_paths]

    common_columns = set(dataframes[0].columns)
    for frame in dataframes[1:]:
        common_columns &= set(frame.columns)

    if not common_columns:
        raise ValueError("Target lists share no common columns; cannot build aux list")

    primary_column_order = [
        col for col in dataframes[0].columns if col in common_columns
    ]

    trimmed = [frame.loc[:, primary_column_order] for frame in dataframes]
    combined = (
        pd.concat(trimmed, ignore_index=True).drop_duplicates().reset_index(drop=True)
    )

    output_path = data_dir / "all_targets.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, index=False)
    return output_path
