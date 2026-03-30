"""Build science calendar XML from observation schedules.

This module generates Pandora science calendar XML files from CSV schedules.
It handles:
- Creating visit and observation sequence XML structure
- Managing occultation and transit observations
- Splitting long observations into sequences
- Integrating target parameters from manifest files
- Outputting properly formatted XML calendars

This is the refactored version of xml_builder.py with clearer naming
and improved documentation.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime, timedelta
from numbers import Number
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from xml.dom import minidom

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from tqdm import tqdm

from pandorascheduler_rework import observation_utils
from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.utils.array_ops import (
    break_long_sequences,
)
from pandorascheduler_rework.utils.io import read_csv_cached, read_parquet_cached
from pandorascheduler_rework.xml import observation_sequence

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class ScienceCalendarInputs:
    """Filesystem pointers consumed by the XML builder."""

    schedule_csv: Path
    data_dir: Path
    schedule_row_start: Optional[int] = None
    schedule_row_end: Optional[int] = None


def generate_science_calendar(
    inputs: ScienceCalendarInputs,
    config: PandoraSchedulerConfig,
    output_path: Optional[Path] = None,
) -> Path:
    """Generate the science calendar XML, matching the legacy behaviour."""

    builder = _ScienceCalendarBuilder(inputs, config)
    calendar_element = builder.build_calendar()
    xml_string = _serialise_calendar(calendar_element)

    destination = output_path or (inputs.data_dir / "Pandora_science_calendar.xml")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(xml_string, encoding="utf-8")
    builder.write_sequence_provenance(
        destination.with_name(f"{destination.stem}_sequence_provenance.csv")
    )
    return destination


class _ScienceCalendarBuilder:
    """Encapsulates the translation from CSV schedules to XML."""

    def __init__(
        self, inputs: ScienceCalendarInputs, config: PandoraSchedulerConfig
    ) -> None:
        self.inputs = inputs
        self.config = config
        self.schedule = read_csv_cached(str(inputs.schedule_csv))
        if self.schedule is None:
            raise FileNotFoundError(f"Schedule CSV missing: {inputs.schedule_csv}")

        row_start = inputs.schedule_row_start
        row_end = inputs.schedule_row_end
        if row_start is not None or row_end is not None:
            start_idx = max(int(row_start or 1) - 1, 0)
            end_idx = None if row_end is None else max(int(row_end), 0)
            if end_idx is not None and end_idx < start_idx:
                raise ValueError(
                    "Invalid schedule row range: schedule_row_end must be >= schedule_row_start"
                )
            self.schedule = self.schedule.iloc[start_idx:end_idx].reset_index(drop=True)

        if self.schedule.empty:
            raise ValueError("Schedule CSV is empty; nothing to convert into XML")

        self.data_dir = inputs.data_dir
        self.target_catalog = _read_catalog(self.data_dir / "exoplanet_targets.csv")
        self.aux_catalog = _read_catalog(self.data_dir / "all_targets.csv")
        self.occ_catalog = _read_catalog(
            self.data_dir / "occultation-standard_targets.csv"
        )

        obs_minutes, occ_minutes = observation_utils.general_parameters(
            config.obs_sequence_duration_min,
            config.occ_sequence_limit_min,
        )
        self.sequence_duration = timedelta(minutes=obs_minutes)
        self.occultation_limit = timedelta(minutes=occ_minutes + 1)

        # Track cumulative observation time for each occultation target
        self.occultation_obs_time: Dict[str, timedelta] = {}
        self.sequence_provenance: list[dict[str, object]] = []

    def _record_sequence_provenance(
        self,
        visit_id: str,
        sequence_id: str,
        target: str,
        priority: str,
        start: datetime,
        stop: datetime,
        sequence_type: str,
        assignment_source: str,
        occultation_pass: str = "",
        visibility_fraction: float = 1.0,
        visible_minutes: float = float("nan"),
        non_visible_minutes: float = float("nan"),
    ) -> None:
        self.sequence_provenance.append(
            {
                "visit_id": visit_id,
                "sequence_id": sequence_id,
                "target": target,
                "priority": priority,
                "start_utc": start.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "stop_utc": stop.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "duration_minutes": (stop - start).total_seconds() / 60.0,
                "sequence_type": sequence_type,
                "assignment_source": assignment_source,
                "occultation_pass": occultation_pass,
                "visibility_fraction": visibility_fraction,
                "visible_minutes": visible_minutes,
                "non_visible_minutes": non_visible_minutes,
            }
        )

    def write_sequence_provenance(self, output_path: Path) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(self.sequence_provenance).to_csv(output_path, index=False)

    def _visibility_stats_in_df(
        self, vis_df: Optional[pd.DataFrame], seg_start: datetime, seg_stop: datetime
    ) -> tuple[float, float, float]:
        if vis_df is None or vis_df.empty:
            return float("nan"), float("nan"), float("nan")

        if "Time_UTC" in vis_df.columns and pd.api.types.is_datetime64_any_dtype(
            vis_df["Time_UTC"]
        ):
            times = pd.to_datetime(vis_df["Time_UTC"], errors="coerce")
        elif "Time(MJD_UTC)" in vis_df.columns:
            times = pd.to_datetime(
                Time(
                    vis_df["Time(MJD_UTC)"].to_numpy(dtype=float),
                    format="mjd",
                    scale="utc",
                ).to_datetime()
            )
        else:
            return float("nan"), float("nan"), float("nan")

        index = pd.DatetimeIndex(times)
        if getattr(index, "tz", None) is not None:
            index = index.tz_localize(None)
        index = index.round("min")

        prepared = pd.DataFrame(
            {"Visible": (vis_df["Visible"].to_numpy(dtype=float) > 0.5)},
            index=index,
        )
        prepared = prepared.groupby(level=0)["Visible"].max().to_frame()

        n_minutes = int(round((seg_stop - seg_start).total_seconds() / 60.0))
        if n_minutes <= 0:
            return float("nan"), float("nan"), float("nan")
        minute_index = pd.date_range(seg_start, periods=n_minutes, freq="min")
        aligned = prepared["Visible"].reindex(minute_index, fill_value=False)
        visible_minutes = float(aligned.sum())
        non_visible_minutes = float((~aligned).sum())
        return float(aligned.mean()), visible_minutes, non_visible_minutes

    def _visibility_fraction_in_df(
        self, vis_df: Optional[pd.DataFrame], seg_start: datetime, seg_stop: datetime
    ) -> float:
        fraction, _, _ = self._visibility_stats_in_df(vis_df, seg_start, seg_stop)
        return fraction

    def _occultation_visibility_stats(
        self, target_name: str, seg_start: datetime, seg_stop: datetime
    ) -> tuple[float, float, float]:
        vis_df = _read_visibility(self.data_dir / "aux_targets" / target_name, target_name)
        return self._visibility_stats_in_df(vis_df, seg_start, seg_stop)

    def _occultation_visibility_fraction(
        self, target_name: str, seg_start: datetime, seg_stop: datetime
    ) -> float:
        fraction, _, _ = self._occultation_visibility_stats(target_name, seg_start, seg_stop)
        return fraction

    def _science_min_duration(self) -> timedelta:
        """Resolved minimum standalone science fragment duration."""
        return timedelta(
            minutes=self.config.effective_min_science_sequence_minutes
        )

    def _occultation_min_duration(self) -> timedelta:
        """Resolved minimum standalone occultation fragment duration."""
        return timedelta(
            minutes=self.config.effective_min_occultation_sequence_minutes
        )

    def _get_occultation_time_limit(self, target_name: str) -> timedelta:
        """Get the time limit for an occultation target.

        Looks up 'Number of Hours Requested' from the occultation manifest.

        When ``requested_occ_time_override`` is False (default), raises
        ValueError if the catalog is missing, the target is not found, or
        the required column is missing.

        When ``requested_occ_time_override`` is True (override enabled), logs a
        warning and returns a very large effective limit so scheduling
        can continue.
        """
        enforce_requested_hours = not self.config.requested_occ_time_override
        _RELAXED_FALLBACK = timedelta(hours=1_000_000)

        if self.occ_catalog is None or self.occ_catalog.empty:
            msg = (
                f"Cannot get time limit for occultation target '{target_name}': "
                "occultation catalog is not loaded"
            )
            if enforce_requested_hours:
                raise ValueError(msg)
            LOGGER.warning("%s — using unlimited fallback", msg)
            return _RELAXED_FALLBACK

        match = self.occ_catalog[self.occ_catalog["Star Name"] == target_name]
        if match.empty:
            msg = f"Occultation target '{target_name}' not found in catalog"
            if enforce_requested_hours:
                raise ValueError(msg)
            LOGGER.warning("%s — using unlimited fallback", msg)
            return _RELAXED_FALLBACK

        if "Number of Hours Requested" not in match.columns:
            msg = (
                "Occultation catalog is missing required 'Number of Hours Requested' "
                "column"
            )
            if enforce_requested_hours:
                raise ValueError(msg)
            LOGGER.warning("%s — using unlimited fallback", msg)
            return _RELAXED_FALLBACK

        hours_req = match.iloc[0]["Number of Hours Requested"]
        if pd.isna(hours_req):
            msg = (
                f"Occultation target '{target_name}' has missing "
                "'Number of Hours Requested' value"
            )
            if enforce_requested_hours:
                raise ValueError(msg)
            LOGGER.warning("%s — using unlimited fallback", msg)
            return _RELAXED_FALLBACK

        return timedelta(hours=float(hours_req))

    def _next_chunk_end(
        self, current: datetime, step: timedelta, segment_stop: datetime
    ) -> datetime:
        """Compute the end of the next chunk, absorbing a short tail.

        If emitting a chunk of *step* would leave a remainder shorter than
        the resolved science threshold, extend this chunk to *segment_stop*
        so no short trailing science sequence is created.
        """
        candidate = min(current + step, segment_stop)
        if candidate >= segment_stop:
            return segment_stop
        remainder = segment_stop - candidate
        if remainder < self._science_min_duration():
            return segment_stop
        return candidate

    def _occ_chunk_end(
        self,
        current: datetime,
        segment_stop: datetime,
        occ_target: Optional[str] = None,
    ) -> datetime:
        """Occultation-aware chunk end, respecting break_occultation_sequences."""
        if self.config.break_occultation_sequences:
            candidate = min(current + self.occultation_limit, segment_stop)
            if candidate >= segment_stop:
                return segment_stop

            remainder = segment_stop - candidate
            if remainder >= self._occultation_min_duration():
                return candidate

            if occ_target:
                acceptable, _ = self._occ_visibility_score(
                    occ_target,
                    current,
                    segment_stop,
                )
                if acceptable:
                    return segment_stop
            return candidate
        return segment_stop

    @staticmethod
    def _iterate_segments(
        augmented_changes: List[int],
        visit_times: List[datetime],
        visibility_flags: List[int],
        start: datetime,
        final_time: datetime,
    ):
        """Yield ``(segment_start, segment_stop, is_visible)`` for each
        visibility-change segment within a visit."""
        last = len(augmented_changes) - 1
        for pos, change_idx in enumerate(augmented_changes):
            seg_start = (
                start if pos == 0
                else visit_times[augmented_changes[pos - 1] + 1]
            )
            seg_stop = (
                final_time if pos == last
                else visit_times[change_idx + 1]
            )
            yield seg_start, seg_stop, bool(visibility_flags[change_idx])

    def _raw_visit_segments(
        self,
        visit_times: Sequence[datetime],
        visibility_flags: Sequence[int],
        start: datetime,
        final_time: datetime,
    ) -> List[Tuple[datetime, datetime, bool]]:
        """Build raw visibility-change segments without short-fragment cleanup."""
        if not visit_times or not visibility_flags:
            return []

        changes = _visibility_change_indices(visibility_flags)
        segments: List[Tuple[datetime, datetime, bool]] = []
        seg_start = start
        for change_idx in changes:
            seg_stop = visit_times[change_idx + 1]
            if seg_stop > seg_start:
                segments.append(
                    (seg_start, seg_stop, bool(visibility_flags[change_idx]))
                )
            seg_start = seg_stop

        if final_time > seg_start:
            segments.append((seg_start, final_time, bool(visibility_flags[-1])))
        return segments

    @staticmethod
    def _coalesce_segments(
        segments: Sequence[Tuple[datetime, datetime, bool]]
    ) -> List[Tuple[datetime, datetime, bool]]:
        """Merge adjacent same-kind segments created by policy rewrites."""
        if not segments:
            return []

        merged: List[Tuple[datetime, datetime, bool]] = [segments[0]]
        tolerance = timedelta(seconds=1)
        for seg_start, seg_stop, is_visible in segments[1:]:
            prev_start, prev_stop, prev_visible = merged[-1]
            if is_visible == prev_visible and seg_start <= prev_stop + tolerance:
                merged[-1] = (prev_start, max(prev_stop, seg_stop), prev_visible)
            else:
                merged.append((seg_start, seg_stop, is_visible))
        return merged

    def _apply_science_fragment_policy(
        self,
        segments: Sequence[Tuple[datetime, datetime, bool]],
    ) -> List[Tuple[datetime, datetime, bool]]:
        """Handle short science-visible fragments before occultation fill.

        Policy:
        1. if a short science fragment can extend a contiguous preceding
           science fragment, merge it into that science chunk
        2. otherwise reclassify it as an occultation-fill interval
        """
        threshold = self._science_min_duration()
        if threshold <= timedelta(0):
            return self._coalesce_segments(list(segments))

        adjusted: List[Tuple[datetime, datetime, bool]] = []
        converted_short_science = 0
        tolerance = timedelta(seconds=1)

        for seg_start, seg_stop, is_visible in segments:
            duration = seg_stop - seg_start
            if is_visible and duration < threshold:
                can_absorb_into_prev_science = (
                    adjusted
                    and adjusted[-1][2]
                    and seg_start <= adjusted[-1][1] + tolerance
                )
                if can_absorb_into_prev_science:
                    prev_start, prev_stop, _ = adjusted[-1]
                    adjusted[-1] = (prev_start, max(prev_stop, seg_stop), True)
                else:
                    adjusted.append((seg_start, seg_stop, False))
                    converted_short_science += 1
                continue

            adjusted.append((seg_start, seg_stop, is_visible))

        if converted_short_science:
            LOGGER.debug(
                "Converted %d short science-visible fragment(s) shorter than %d "
                "min into occultation-fill intervals",
                converted_short_science,
                self.config.effective_min_science_sequence_minutes,
            )

        return self._coalesce_segments(adjusted)

    @staticmethod
    def _occultation_windows_from_segments(
        segments: Sequence[Tuple[datetime, datetime, bool]]
    ) -> tuple[List[datetime], List[datetime]]:
        """Extract occultation intervals from the policy-adjusted segment list."""
        starts: List[datetime] = []
        stops: List[datetime] = []
        for seg_start, seg_stop, is_visible in segments:
            if is_visible or seg_stop <= seg_start:
                continue
            starts.append(seg_start)
            stops.append(seg_stop)
        return starts, stops

    def _merged_segments(
        self,
        augmented_changes: List[int],
        visit_times: List[datetime],
        visibility_flags: List[int],
        start: datetime,
        final_time: datetime,
    ):
        """Yield ``(segment_start, segment_stop, is_visible)`` like
        :meth:`_iterate_segments`, but absorb any segment whose duration is
        below ``config.min_sequence_minutes`` into its following neighbour
        (or preceding neighbour when no following one exists).

        When ``min_sequence_minutes`` is 0 the output is identical to
        :meth:`_iterate_segments`."""
        threshold = timedelta(minutes=self.config.min_sequence_minutes)
        segments = list(
            self._iterate_segments(
                augmented_changes, visit_times, visibility_flags, start, final_time,
            )
        )
        if not threshold or len(segments) <= 1:
            yield from segments
            return

        # Iteratively absorb short segments until stable.
        changed = True
        while changed:
            changed = False
            merged: list = []
            i = 0
            while i < len(segments):
                seg_start, seg_stop, is_vis = segments[i]
                duration = seg_stop - seg_start
                if duration < threshold:
                    if i + 1 < len(segments):
                        # Absorb this short segment into the next one.
                        nxt_start, nxt_stop, nxt_vis = segments[i + 1]
                        segments[i + 1] = (seg_start, nxt_stop, nxt_vis)
                        changed = True
                    elif merged:
                        # No following segment: absorb into the previous one.
                        prev_start, prev_stop, prev_vis = merged[-1]
                        merged[-1] = (prev_start, seg_stop, prev_vis)
                        changed = True
                    # If there is neither a previous nor a next segment, we
                    # must have len(segments) == 1, which is handled by the
                    # early return above, so no action is needed here.
                else:
                    merged.append((seg_start, seg_stop, is_vis))
                i += 1
            segments = merged

        yield from segments

    def _emit_science_sequences(
        self,
        visit_element: ET.Element,
        visit_id: str,
        seq_counter: int,
        target_name: str,
        segment_start: datetime,
        segment_stop: datetime,
        ra_value: float,
        dec_value: float,
        target_info: Optional[pd.DataFrame],
        science_visibility_df: Optional[pd.DataFrame],
        priority_flag: bool,
        transit_start: Sequence[datetime],
        transit_stop: Sequence[datetime],
    ) -> int:
        """Emit chunked science observation sequences.  Returns updated
        *seq_counter*."""
        min_td = self._science_min_duration()
        if min_td and (segment_stop - segment_start) < min_td:
            LOGGER.debug(
                "Skipping short science segment for %s (%s < %d min)",
                target_name,
                segment_stop - segment_start,
                self.config.effective_min_science_sequence_minutes,
            )
            return seq_counter
        current = segment_start
        while current < segment_stop:
            next_value = self._next_chunk_end(
                current, self.sequence_duration, segment_stop
            )
            priority = _target_priority(
                priority_flag, transit_start, transit_stop, current, next_value,
            )
            sequence_id = f"{seq_counter:03d}"
            observation_sequence(
                visit_element,
                sequence_id,
                target_name,
                priority,
                current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                next_value.strftime("%Y-%m-%dT%H:%M:%SZ"),
                ra_value,
                dec_value,
                target_info if target_info is not None else pd.DataFrame(),
            )
            science_visibility_fraction, science_visible_minutes, science_non_visible_minutes = (
                self._visibility_stats_in_df(science_visibility_df, current, next_value)
            )
            self._record_sequence_provenance(
                visit_id,
                sequence_id,
                target_name,
                priority,
                current,
                next_value,
                "science",
                "science_schedule",
                visibility_fraction=science_visibility_fraction,
                visible_minutes=science_visible_minutes,
                non_visible_minutes=science_non_visible_minutes,
            )
            seq_counter += 1
            current = next_value
        return seq_counter

    def _emit_occultation_sequences(
        self,
        visit_element: ET.Element,
        visit_id: str,
        seq_counter: int,
        occ_target: str,
        segment_start: datetime,
        segment_stop: datetime,
        ra_occ: float,
        dec_occ: float,
        occ_info: Optional[pd.DataFrame],
        reference_ra: Optional[float] = None,
        reference_dec: Optional[float] = None,
        assignment_source: str = "catalog_fallback",
        occultation_pass: str = "",
    ) -> int:
        """Emit chunked occultation observation sequences.  Returns updated
        *seq_counter*."""
        current = segment_start
        while current < segment_stop:
            next_value = self._occ_chunk_end(current, segment_stop, occ_target)

            chunk_target = occ_target
            chunk_ra = ra_occ
            chunk_dec = dec_occ
            chunk_info = occ_info

            acceptable, _ = self._occ_visibility_score(
                chunk_target,
                current,
                next_value,
            )
            if (
                not acceptable
                and reference_ra is not None
                and reference_dec is not None
            ):
                fallback = self._select_fallback_occultation_target(
                    reference_ra,
                    reference_dec,
                    seg_start=current,
                    seg_stop=next_value,
                )
                if fallback is not None:
                    chunk_target, chunk_ra, chunk_dec, chunk_info = fallback
                    acceptable = True

            if not acceptable:
                LOGGER.warning(
                    "No visible occultation target for interval %s–%s",
                    current,
                    next_value,
                )
                current = next_value
                continue

            current_occ_time = self.occultation_obs_time.get(
                chunk_target, timedelta()
            )
            target_time_limit = self._get_occultation_time_limit(chunk_target)
            if current_occ_time >= target_time_limit:
                LOGGER.info(
                    "Skipping %s: exceeded occultation time limit (%.1f/%.1f hrs)",
                    chunk_target,
                    current_occ_time.total_seconds() / 3600,
                    target_time_limit.total_seconds() / 3600,
                )
                break
            sequence_id = f"{seq_counter:03d}"
            observation_sequence(
                visit_element,
                sequence_id,
                chunk_target,
                "0",
                current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                next_value.strftime("%Y-%m-%dT%H:%M:%SZ"),
                chunk_ra,
                chunk_dec,
                chunk_info if chunk_info is not None else pd.DataFrame(),
            )
            occ_visibility_fraction, occ_visible_minutes, occ_non_visible_minutes = (
                self._occultation_visibility_stats(chunk_target, current, next_value)
            )
            self._record_sequence_provenance(
                visit_id,
                sequence_id,
                chunk_target,
                "0",
                current,
                next_value,
                "occultation",
                assignment_source,
                occultation_pass,
                visibility_fraction=occ_visibility_fraction,
                visible_minutes=occ_visible_minutes,
                non_visible_minutes=occ_non_visible_minutes,
            )
            self.occultation_obs_time[chunk_target] = (
                self.occultation_obs_time.get(chunk_target, timedelta())
                + (next_value - current)
            )
            seq_counter += 1
            current = next_value
        return seq_counter

    def build_calendar(self) -> ET.Element:
        root = ET.Element("ScienceCalendar", xmlns="/pandora/calendar/")
        self._add_meta(root)

        visits = self.schedule
        if self.config.visit_limit is not None:
            visits = visits.head(self.config.visit_limit)

        iterator = tqdm(
            visits.iterrows(),
            total=len(visits),
            desc="Building science calendar",
            disable=not self.config.show_progress,
        )

        for visit_counter, (_, row) in enumerate(iterator, start=1):
            self._add_visit(root, visit_counter, row)

        return root

    def _add_meta(self, root: ET.Element) -> None:
        weights = ", ".join(
            f"{value:.1f}" for value in self.config.transit_scheduling_weights
        )
        keepout = ", ".join(
            f"{value:.1f}"
            for value in (
                self.config.sun_avoidance_deg,
                self.config.moon_avoidance_deg,
                self.config.earth_avoidance_deg,
            )
        )

        valid_from = str(self.schedule.iloc[0]["Observation Start"])
        expires = str(self.schedule.iloc[self.schedule.index[-1]]["Observation Stop"])

        raw_created = self.config.created_timestamp
        if isinstance(raw_created, str):
            created_value = raw_created
        else:
            timestamp = raw_created or datetime.now()
            # Round to nearest second
            created_value = str(
                (timestamp + timedelta(microseconds=500_000)).replace(microsecond=0)
            )

        attrs = {
            "Valid_From": valid_from,
            "Expires": expires,
            "Calendar_Weights": weights,
            "Keepout_Angles": keepout,
            "Observation_Sequence_Duration_hrs_max": str(self.sequence_duration),
            "Removed_Sequences_Shorter_Than_min": str(
                self.config.effective_min_science_sequence_minutes
            ),
            "Occultation_Sequences_Shorter_Than_min": str(
                self.config.effective_min_occultation_sequence_minutes
            ),
            "Created": created_value,
            "Delivery_Id": "",
        }
        if self.config.author:
            attrs["Author"] = self.config.author

        ET.SubElement(root, "Meta", attrib=attrs)

    def _add_visit(self, root: ET.Element, visit_counter: int, row: pd.Series) -> None:
        id_padding = 4 - len(str(visit_counter))
        target_label = str(row.get("Target", ""))

        if not target_label or target_label == "Free Time":
            return

        if target_label.startswith("WARNING"):
            LOGGER.warning("Need visible STD during %s", target_label)
            return

        target_name, star_name = _normalise_target_name(target_label)

        visit_element = ET.SubElement(root, "Visit")
        visit_id = f"{'0' * id_padding}{visit_counter}"
        ET.SubElement(visit_element, "ID").text = visit_id

        start = _parse_datetime(row.get("Observation Start"))
        stop = _parse_datetime(row.get("Observation Stop"))
        if start is None or stop is None:
            raise ValueError(f"Unable to parse observation window for {target_label}")

        planet_row = _lookup_planet_row(self.target_catalog, target_name)
        has_transit = _is_transit_entry(row)

        if planet_row is not None and has_transit:
            visibility_df = _read_visibility(
                self.data_dir / "targets" / star_name, star_name
            )
            transit_df = _read_planet_visibility(
                self.data_dir / "targets" / star_name / target_name, target_name
            )
            target_info = planet_row
            transit_windows = _transit_windows(transit_df)
            transit_start, transit_stop = (
                transit_windows if transit_windows else ([], [])
            )
            priority_flag = True
        else:
            visibility_df = _read_visibility(
                self.data_dir / "aux_targets" / star_name, target_name
            )
            target_info = _lookup_auxiliary_row(self.aux_catalog, target_name)
            transit_start, transit_stop = ([], [])
            priority_flag = False

        if visibility_df is None or visibility_df.empty:
            LOGGER.error(
                "No visibility data for %s. Aborting schedule build.", target_name
            )
            return

        try:
            ra_value = (
                float(target_info.iloc[0]["RA"])
                if target_info is not None
                else float("nan")
            )
            dec_value = (
                float(target_info.iloc[0]["DEC"])
                if target_info is not None
                else float("nan")
            )
        except (KeyError, ValueError, TypeError, AttributeError):
            ra_value, dec_value = _resolve_coordinates(star_name)

        visit_times, visibility_flags = _extract_visibility_segment(
            visibility_df,
            start,
            stop,
            self.config.min_sequence_minutes,
        )
        if not visit_times:
            if has_transit and transit_start and transit_stop:
                try:
                    visit_times = [transit_start[0], transit_stop[0]]
                    visibility_flags = [1, 1]
                except Exception:
                    source, time_min, time_max, in_window = _visibility_diagnostics(
                        visibility_df, start, stop
                    )
                    LOGGER.warning(
                        "No visibility samples within visit for %s (visit=%s..%s, vis_file=%s, vis_range=%s..%s, samples_in_window=%s)",
                        target_name,
                        start,
                        stop,
                        source,
                        time_min,
                        time_max,
                        in_window,
                    )
                    return
            else:
                source, time_min, time_max, in_window = _visibility_diagnostics(
                    visibility_df, start, stop
                )
                LOGGER.warning(
                    "No visibility samples within visit for %s (visit=%s..%s, vis_file=%s, vis_range=%s..%s, samples_in_window=%s)",
                    target_name,
                    start,
                    stop,
                    source,
                    time_min,
                    time_max,
                    in_window,
                )
                return

        # Use the CSV stop boundary (not visit_times[-1]) so the XML visit
        # spans the full scheduled window and no 1-minute inter-visit gap
        # is introduced by visibility-sample rounding.
        final_time = stop
        raw_segments = self._raw_visit_segments(
            visit_times,
            visibility_flags,
            start,
            final_time,
        )
        if not raw_segments:
            LOGGER.warning(
                "No usable visibility segments within visit for %s (%s..%s)",
                target_name,
                start,
                stop,
            )
            return
        seq_counter = 1

        # --- Path 1: occultation XML disabled — science-only ---------------
        if not self.config.enable_occultation_xml:
            for seg_start, seg_stop, is_visible in raw_segments:
                if is_visible:
                    seq_counter = self._emit_science_sequences(
                        visit_element, visit_id, seq_counter, target_name,
                        seg_start, seg_stop, ra_value, dec_value,
                        target_info, visibility_df, priority_flag, transit_start, transit_stop,
                    )
            return

        visit_segments = self._apply_science_fragment_policy(raw_segments)
        oc_starts, oc_stops = self._occultation_windows_from_segments(
            visit_segments
        )
        if not oc_starts:
            for seg_start, seg_stop, is_visible in visit_segments:
                if is_visible:
                    seq_counter = self._emit_science_sequences(
                        visit_element, visit_id, seq_counter, target_name,
                        seg_start, seg_stop, ra_value, dec_value,
                        target_info, visibility_df, priority_flag, transit_start, transit_stop,
                    )
            return

        # --- Resolve the occultation source for this visit ------------------
        occultation_info = self._find_occultation_target(
            oc_starts, oc_stops, start, final_time, ra_value, dec_value,
        )

        # Determine whether we have a scheduled occ_df or need a fallback.
        occ_df: Optional[pd.DataFrame] = None
        fallback_occultation: Optional[
            tuple[str, float, float, Optional[pd.DataFrame]]
        ] = None

        if occultation_info is not None:
            occ_df, scheduled = occultation_info
            if occ_df is not None and not occ_df.empty and "Target" in occ_df.columns:
                assigned_mask = occ_df["Target"].astype(str).str.strip().ne("")
                has_assigned_rows = bool(assigned_mask.any())
            else:
                has_assigned_rows = False

            if occ_df is None or occ_df.empty or not has_assigned_rows:
                occ_df = None
            elif not scheduled:
                LOGGER.info(
                    "Using partial scheduled occultation rows for %s between %s and %s; "
                    "unassigned occultation intervals will fall back per segment",
                    target_name,
                    start,
                    final_time,
                )

        if occ_df is None:
            LOGGER.warning(
                "Unable to schedule occultation target for %s between %s and %s",
                target_name, start, final_time,
            )

        # --- Path 2: catalog-fallback occultation (no occ_df) ---------------
        if occ_df is None:
            for seg_start, seg_stop, is_visible in visit_segments:
                if is_visible:
                    seq_counter = self._emit_science_sequences(
                        visit_element, visit_id, seq_counter, target_name,
                        seg_start, seg_stop, ra_value, dec_value,
                        target_info, visibility_df, priority_flag, transit_start, transit_stop,
                    )
                else:
                    # Select per-segment so the visibility check uses the
                    # actual occultation-gap times — a target visible in one gap
                    # may violate sun keepout in another.
                    fallback_occultation = self._select_fallback_occultation_target(
                        ra_value, dec_value,
                        seg_start=seg_start, seg_stop=seg_stop,
                    )
                    if fallback_occultation is not None:
                        occ_target, ra_occ, dec_occ, occ_info = fallback_occultation
                        seq_counter = self._emit_occultation_sequences(
                            visit_element, visit_id, seq_counter, occ_target,
                            seg_start, seg_stop, ra_occ, dec_occ, occ_info,
                            reference_ra=ra_value,
                            reference_dec=dec_value,
                            assignment_source="catalog_fallback",
                        )
            return

        # --- Path 3: scheduled occ_df available -----------------------------
        occ_time_index: Optional[pd.DataFrame] = None
        if {"start", "stop", "Target"}.issubset(set(occ_df.columns)):
            occ_time_index = occ_df.copy()
            try:
                occ_time_index["_start_dt"] = pd.to_datetime(
                    occ_time_index["start"], utc=True, format="ISO8601", errors="coerce"
                ).dt.tz_localize(None)
                occ_time_index["_stop_dt"] = pd.to_datetime(
                    occ_time_index["stop"], utc=True, format="ISO8601", errors="coerce"
                ).dt.tz_localize(None)
                occ_time_index = occ_time_index.dropna(
                    subset=["_start_dt", "_stop_dt", "Target"]
                )
            except Exception:
                occ_time_index = None

        oc_index = 0
        for seg_start, seg_stop, is_visible in visit_segments:
            if is_visible:
                seq_counter = self._emit_science_sequences(
                    visit_element, visit_id, seq_counter, target_name,
                    seg_start, seg_stop, ra_value, dec_value,
                    target_info, visibility_df, priority_flag, transit_start, transit_stop,
                )
                continue

            # Occultation segment — iterate using scheduled occ_df.
            current = seg_start
            while current < seg_stop:
                if oc_index >= len(occ_df):
                    # Pre-built schedule exhausted — fall back to catalog
                    # after the loop via the post-loop fallback block.
                    break

                lookup_stop = self._occ_chunk_end(current, seg_stop)

                # Prefer time-based lookup; fall back to positional index.
                occ_row = None
                used_fallback_row = False
                if occ_time_index is not None and not occ_time_index.empty:
                    exact_mask = (
                        (occ_time_index["_start_dt"] <= current)
                        & (occ_time_index["_stop_dt"] >= lookup_stop)
                    )
                    if exact_mask.any():
                        occ_row = occ_time_index.loc[exact_mask].iloc[0]
                    else:
                        overlap_mask = (
                            (occ_time_index["_start_dt"] < lookup_stop)
                            & (occ_time_index["_stop_dt"] > current)
                        )
                        if overlap_mask.any():
                            occ_row = occ_time_index.loc[overlap_mask].iloc[0]

                if occ_row is None:
                    occ_row = occ_df.iloc[oc_index]
                    used_fallback_row = True

                occ_target = str(occ_row["Target"]).strip()
                next_value = self._occ_chunk_end(
                    current,
                    seg_stop,
                    occ_target if occ_target and occ_target.lower() != "no target" else None,
                )

                if not occ_target or occ_target.lower() == "no target":
                    fallback = self._select_fallback_occultation_target(
                        ra_value,
                        dec_value,
                        seg_start=current,
                        seg_stop=next_value,
                    )
                    if fallback is not None:
                        fb_target, fb_ra, fb_dec, fb_info = fallback
                        seq_counter = self._emit_occultation_sequences(
                            visit_element,
                            visit_id,
                            seq_counter,
                            fb_target,
                            current,
                            next_value,
                            fb_ra,
                            fb_dec,
                            fb_info,
                            reference_ra=ra_value,
                            reference_dec=dec_value,
                            assignment_source="catalog_fallback",
                        )
                    else:
                        LOGGER.debug(
                            "No scheduled occultation target for %s–%s and no fallback "
                            "candidate was visible",
                            current,
                            next_value,
                        )
                    if used_fallback_row:
                        oc_index += 1
                    current = next_value
                    continue

                # Check if this occultation target has exceeded its time limit
                current_occ_time = self.occultation_obs_time.get(
                    occ_target, timedelta()
                )
                target_time_limit = self._get_occultation_time_limit(occ_target)
                if current_occ_time >= target_time_limit:
                    LOGGER.info(
                        "Skipping %s: exceeded occultation time limit "
                        "(%.1f/%.1f hrs)",
                        occ_target,
                        current_occ_time.total_seconds() / 3600,
                        target_time_limit.total_seconds() / 3600,
                    )
                    oc_index += 1
                    continue

                # Visibility gate: reject targets that violate keepout
                # constraints (e.g. sun keepout) at the scheduled time.
                # When allow_occ_startracker_violation is enabled, targets
                # with ST-only violations are still accepted.
                acceptable, _st_frac = self._occ_visibility_score(
                    occ_target, current, next_value,
                )
                if not acceptable:
                    LOGGER.debug(
                        "Skipping occultation target %s (%s–%s): "
                        "not visible (keepout constraint violated)",
                        occ_target, current, next_value,
                    )
                    oc_index += 1
                    continue

                occ_info = _lookup_occultation_info(
                    occ_target,
                    self.target_catalog,
                    self.aux_catalog,
                    self.occ_catalog,
                )
                ra_occ = _fallback_float(
                    occ_row.get("RA"), occ_info, "RA"
                )
                dec_occ = _fallback_float(
                    occ_row.get("DEC"), occ_info, "DEC"
                )

                sequence_id = f"{seq_counter:03d}"
                observation_sequence(
                    visit_element,
                    sequence_id,
                    occ_target,
                    "0",
                    current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    next_value.strftime("%Y-%m-%dT%H:%M:%SZ"),
                    ra_occ,
                    dec_occ,
                    occ_info if occ_info is not None else pd.DataFrame(),
                )
                occ_pass = str(occ_row.get("Occultation Pass", "") or "")
                occ_visibility_fraction, occ_visible_minutes, occ_non_visible_minutes = (
                    self._occultation_visibility_stats(occ_target, current, next_value)
                )
                self._record_sequence_provenance(
                    visit_id,
                    sequence_id,
                    occ_target,
                    "0",
                    current,
                    next_value,
                    "occultation",
                    "scheduled_occultation",
                    occ_pass,
                    visibility_fraction=occ_visibility_fraction,
                    visible_minutes=occ_visible_minutes,
                    non_visible_minutes=occ_non_visible_minutes,
                )

                sequence_duration = next_value - current
                self.occultation_obs_time[occ_target] = (
                    self.occultation_obs_time.get(occ_target, timedelta())
                    + sequence_duration
                )

                seq_counter += 1
                if used_fallback_row:
                    oc_index += 1
                current = next_value

            # After exhausting occ_df rows, try catalog fallback for any
            # remaining occultation time in this segment.
            if current < seg_stop:
                fallback = self._select_fallback_occultation_target(
                    ra_value, dec_value,
                    seg_start=current, seg_stop=seg_stop,
                )
                if fallback is not None:
                    fb_target, fb_ra, fb_dec, fb_info = fallback
                    seq_counter = self._emit_occultation_sequences(
                        visit_element, visit_id, seq_counter, fb_target,
                        current, seg_stop, fb_ra, fb_dec, fb_info,
                        reference_ra=ra_value,
                        reference_dec=dec_value,
                    )
                    current = seg_stop
                else:
                    gap_minutes = (seg_stop - current).total_seconds() / 60
                    LOGGER.warning(
                        "No visible occultation target for %s during "
                        "occultation gap %s–%s (%.0f min gap in XML)",
                        target_name, current, seg_stop, gap_minutes,
                    )

    def _emit_full_visibility(
        self,
        visit_element: ET.Element,
        target_name: str,
        start: datetime,
        stop: datetime,
        ra_value: float,
        dec_value: float,
        target_info: pd.DataFrame | None,
        priority_flag: bool,
        transit_start: Sequence[datetime],
        transit_stop: Sequence[datetime],
    ) -> None:
        segments = break_long_sequences(
            start, stop, self.sequence_duration,
            min_chunk=self._science_min_duration(),
        )
        seq_counter = 1
        for seg_start, seg_stop in segments:
            priority = _target_priority(
                priority_flag,
                transit_start,
                transit_stop,
                seg_start,
                seg_stop,
            )
            observation_sequence(
                visit_element,
                f"{seq_counter:03d}",
                target_name,
                priority,
                seg_start,
                seg_stop,
                ra_value,
                dec_value,
                target_info if target_info is not None else pd.DataFrame(),
            )
            seq_counter += 1

    def _is_target_visible_in_segment(
        self,
        star_name: str,
        seg_start: datetime,
        seg_stop: datetime,
    ) -> bool:
        """Return True if *star_name* has any visible minutes in [seg_start, seg_stop).

        Reads the target's ``aux_targets`` visibility parquet and checks the
        pre-computed ``Visible`` flag.  Returns False when the target is not
        visible or when visibility data cannot be loaded.
        """
        vis_df = _read_visibility(
            self.data_dir / "aux_targets" / star_name, star_name,
        )
        if vis_df is None or vis_df.empty:
            return False

        if "Time_UTC" in vis_df.columns and pd.api.types.is_datetime64_any_dtype(
            vis_df["Time_UTC"]
        ):
            times = vis_df["Time_UTC"].to_numpy(dtype="datetime64[ns]")
            mask = (
                (times >= np.datetime64(seg_start))
                & (times < np.datetime64(seg_stop))
            )
        elif "Time(MJD_UTC)" in vis_df.columns:
            start_mjd = Time(seg_start).mjd
            stop_mjd = Time(seg_stop).mjd
            mjd = vis_df["Time(MJD_UTC)"].to_numpy(dtype=float)
            mask = (mjd >= start_mjd) & (mjd < stop_mjd)
        else:
            return False

        segment_vis = vis_df.loc[mask, "Visible"]
        if segment_vis.empty:
            return False
        return bool((segment_vis > 0).any())

    def _occ_visibility_score(
        self,
        star_name: str,
        seg_start: datetime,
        seg_stop: datetime,
    ) -> tuple[bool, float]:
        """Score an occultation candidate's visibility in *[seg_start, seg_stop)*.

        Returns ``(acceptable, st_violation_frac)`` where:

        * *acceptable* — True if the target may be scheduled.  Always True when
          fully visible.  When ``allow_occ_startracker_violation`` is enabled,
          also True if the **only** keepout failure is star-tracker (boresight
          Sun/Moon/Earth constraints all pass).
        * *st_violation_frac* — fraction of segment minutes where the target
          is not visible (0.0 when fully visible, >0 for ST-only violations).
          Used to rank candidates: lower is better.
        """
        # Fast path: fully visible → always acceptable.
        if self._is_target_visible_in_segment(star_name, seg_start, seg_stop):
            return True, 0.0

        if not self.config.allow_occ_startracker_violation:
            return False, 1.0

        # Extended check: read boresight separation columns to determine
        # whether non-visibility is solely due to star-tracker keepout.
        vis_df = _read_visibility_extended(
            self.data_dir / "aux_targets" / star_name, star_name,
        )
        if vis_df is None or vis_df.empty:
            return False, 1.0

        # Build time mask for the segment.
        if "Time_UTC" in vis_df.columns and pd.api.types.is_datetime64_any_dtype(
            vis_df["Time_UTC"]
        ):
            times = vis_df["Time_UTC"].to_numpy(dtype="datetime64[ns]")
            mask = (
                (times >= np.datetime64(seg_start))
                & (times < np.datetime64(seg_stop))
            )
        elif "Time(MJD_UTC)" in vis_df.columns:
            start_mjd = Time(seg_start).mjd
            stop_mjd = Time(seg_stop).mjd
            mjd = vis_df["Time(MJD_UTC)"].to_numpy(dtype=float)
            mask = (mjd >= start_mjd) & (mjd < stop_mjd)
        else:
            return False, 1.0

        seg = vis_df.loc[mask]
        if seg.empty:
            return False, 1.0

        # Verify that the required separation columns exist.
        needed = {"Sun_Sep", "Moon_Sep", "Earth_Sep"}
        if not needed.issubset(seg.columns):
            return False, 1.0

        # Check boresight constraints on every minute in the segment.
        boresight_ok = (
            (seg["Sun_Sep"] >= self.config.sun_avoidance_deg)
            & (seg["Moon_Sep"] >= self.config.moon_avoidance_deg)
            & (seg["Earth_Sep"] >= self.config.earth_avoidance_deg)
        )

        if not boresight_ok.all():
            # At least one minute fails a boresight constraint — reject.
            return False, 1.0

        # Boresight passes everywhere → failure is star-tracker only.
        n_total = len(seg)
        n_not_visible = int((seg["Visible"] <= 0).sum())
        st_frac = n_not_visible / n_total if n_total else 1.0

        LOGGER.debug(
            "Occultation candidate %s: ST-only violation %.0f%% of segment "
            "(%s–%s)",
            star_name, st_frac * 100, seg_start, seg_stop,
        )
        return True, st_frac

    def _select_fallback_occultation_target(
        self,
        reference_ra: float,
        reference_dec: float,
        seg_start: Optional[datetime] = None,
        seg_stop: Optional[datetime] = None,
    ) -> Optional[tuple[str, float, float, Optional[pd.DataFrame]]]:
        if self.occ_catalog is None or self.occ_catalog.empty:
            return None
        if "Star Name" not in self.occ_catalog.columns:
            return None

        candidates = self.occ_catalog.copy()
        # Each accepted row carries its ST-violation fraction for ranking.
        scored_rows: list[tuple[pd.Series, float]] = []
        for _, row in candidates.iterrows():
            name = str(row.get("Star Name", "")).strip()
            if not name:
                continue
            current_occ_time = self.occultation_obs_time.get(name, timedelta())
            if current_occ_time >= self._get_occultation_time_limit(name):
                continue
            # When segment times are provided, filter by visibility score
            # (respects allow_occ_startracker_violation).
            if seg_start is not None and seg_stop is not None:
                acceptable, st_frac = self._occ_visibility_score(
                    name, seg_start, seg_stop,
                )
                if not acceptable:
                    LOGGER.debug(
                        "Skipping occultation candidate %s: not visible "
                        "between %s and %s",
                        name, seg_start, seg_stop,
                    )
                    continue
                scored_rows.append((row, st_frac))
            else:
                scored_rows.append((row, 0.0))

        if not scored_rows:
            return None

        # Sort: prefer fully-visible (st_frac == 0) first, then least ST
        # violation.  Within the same tier, preserve catalog/slew order.
        scored_rows.sort(key=lambda r: r[1])

        available_rows = [r[0] for r in scored_rows]
        candidates = pd.DataFrame(available_rows)
        if self.config.prioritise_occultations_by_slew:
            # Slew-priority sorting only among the best ST-violation tier.
            # Extract the best (lowest) fraction and filter to that tier,
            # then slew-sort within it.
            best_frac = scored_rows[0][1]
            tier_rows = [r for r, f in scored_rows if f <= best_frac + 1e-9]
            tier_df = pd.DataFrame(tier_rows)
            tier_df = _prioritise_occultation_targets(
                tier_df, reference_ra, reference_dec,
            )
            candidates = tier_df

        chosen = candidates.iloc[0]
        occ_target = str(chosen.get("Star Name", "")).strip()
        if not occ_target:
            return None

        occ_info = _lookup_occultation_info(
            occ_target,
            self.target_catalog,
            self.aux_catalog,
            self.occ_catalog,
        )
        ra_occ = _fallback_float(chosen.get("RA"), occ_info, "RA")
        dec_occ = _fallback_float(chosen.get("DEC"), occ_info, "DEC")
        return occ_target, ra_occ, dec_occ, occ_info

    def _find_occultation_target(
        self,
        starts: Sequence[datetime],
        stops: Sequence[datetime],
        visit_start: datetime,
        visit_stop: datetime,
        reference_ra: float,
        reference_dec: float,
    ) -> Optional[tuple[pd.DataFrame, bool]]:
        if not starts or not stops:
            return None

        if self.config.break_occultation_sequences:
            expanded_starts: list[datetime] = []
            expanded_stops: list[datetime] = []
            for start, stop in zip(starts, stops):
                segments = break_long_sequences(start, stop, self.occultation_limit)
                if not segments:
                    expanded_starts.append(start)
                    expanded_stops.append(stop)
                    continue

                for segment_start, segment_stop in segments:
                    expanded_starts.append(segment_start)
                    expanded_stops.append(segment_stop)
        else:
            expanded_starts = list(starts)
            expanded_stops = list(stops)

        if not expanded_starts:
            return None

        candidates: List[Tuple[Path, str, Path]] = [
            (
                self.data_dir / "occultation-standard_targets.csv",
                "occ list",
                self.data_dir / "aux_targets",
            ),
        ]
        if not self.config.use_target_list_for_occultations:
            candidates.reverse()

        # Build set of targets that have exceeded their time limit
        excluded_targets = {
            name
            for name, obs_time in self.occultation_obs_time.items()
            if obs_time >= self._get_occultation_time_limit(name)
        }

        def _try_candidates(excluded: Optional[set]) -> Optional[tuple[pd.DataFrame, bool]]:
            for csv_path, label, vis_root in candidates:
                result_df, flag = _build_occultation_schedule(
                    expanded_starts,
                    expanded_stops,
                    visit_start,
                    visit_stop,
                    csv_path,
                    vis_root,
                    label,
                    reference_ra,
                    reference_dec,
                    self.config.prioritise_occultations_by_slew,
                    excluded,
                    show_progress=self.config.show_progress,
                    use_pass1=self.config.enable_occultation_pass1,
                )
                if flag and result_df is not None:
                    return result_df, True
            return None

        result = _try_candidates(excluded_targets)
        if result is not None:
            return result

        # If strict limits are disabled, retry without exclusions.
        if (self.config.requested_occ_time_override) and excluded_targets:
            LOGGER.warning(
                "No occultation target assigned for %s..%s with %d targets excluded "
                "by time limits; retrying without exclusions",
                visit_start,
                visit_stop,
                len(excluded_targets),
            )
            result = _try_candidates(set())
            if result is not None:
                return result

        LOGGER.warning(
            "Occultation assignment failed for %s..%s (excluded_targets=%d, "
            "requested_hours_override=%s)",
            visit_start,
            visit_stop,
            len(excluded_targets),
            self.config.requested_occ_time_override,
        )
        return None


def _read_catalog(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required target catalog missing: {path}")
    df = read_csv_cached(str(path))
    if df is None:
        raise FileNotFoundError(f"Unable to read catalog: {path}")
    return df


def _normalise_target_name(target: str) -> tuple[str, str]:
    if target.endswith("STD"):
        stripped = target[:-4]
        return stripped, stripped
    # Star names whose trailing letter coincidentally matches a planet-suffix
    # character (b/c/d/e/f).  Add entries here when a new catalogue target
    # ends with one of those letters but is NOT a planet designation.
    _STAR_NAME_EXCEPTIONS = frozenset(
        {"EV_Lac", "AF_Psc", "AU_Mic", "55_Cnc", "AO_Cassiopeiae"}
    )
    if target.endswith(tuple("bcdef")) and target not in _STAR_NAME_EXCEPTIONS:
        return target, target[:-1].strip()
    return target, target


def _parse_datetime(value: object) -> Optional[datetime]:
    if isinstance(value, datetime):
        return value
    if isinstance(value, str):
        for pattern in (
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%dT%H:%M:%SZ",
            "%Y-%m-%d %H:%M:%S.%f",
            "%Y-%m-%d",  # Date only (time assumed to be 00:00:00)
        ):
            try:
                return datetime.strptime(value, pattern)
            except ValueError:
                continue
    return None


def _is_transit_entry(row: pd.Series) -> bool:
    value = row.get("Transit Coverage")
    if value is None:
        return False
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False
    return np.isfinite(numeric)


def _lookup_planet_row(
    catalog: pd.DataFrame, target_name: str
) -> Optional[pd.DataFrame]:
    if catalog.empty or "Planet Name" not in catalog.columns:
        return None
    match = catalog.loc[catalog["Planet Name"] == target_name]
    if match.empty:
        return None
    return match.head(1)


def _lookup_auxiliary_row(
    catalog: pd.DataFrame, target_name: str
) -> Optional[pd.DataFrame]:
    if catalog.empty or "Star Name" not in catalog.columns:
        return None
    match = catalog.loc[catalog["Star Name"] == target_name]
    if match.empty:
        return None
    # Filter out exoplanet rows (those with a Planet Name) so we return
    # the auxiliary row when both planet and auxiliary entries exist.
    if "Planet Name" in match.columns:
        aux_only = match.loc[match["Planet Name"].isna() | (match["Planet Name"] == "")]
        if not aux_only.empty:
            return aux_only.head(1)
    return match.head(1)


def _lookup_occultation_info(
    target_name: str,
    planet_catalog: pd.DataFrame,
    aux_catalog: pd.DataFrame,
    occ_catalog: pd.DataFrame,
) -> Optional[pd.DataFrame]:
    if "Star Name" in planet_catalog.columns:
        planet_match = planet_catalog.loc[planet_catalog["Star Name"] == target_name]
        if not planet_match.empty:
            return planet_match.head(1)
    if "Planet Name" in planet_catalog.columns:
        planet_match = planet_catalog.loc[planet_catalog["Planet Name"] == target_name]
        if not planet_match.empty:
            return planet_match.head(1)

    aux_match = aux_catalog.loc[
        aux_catalog.get("Star Name", pd.Series(dtype=object)) == target_name
    ]
    if not aux_match.empty:
        return aux_match.head(1)

    occ_match = occ_catalog.loc[
        occ_catalog.get("Star Name", pd.Series(dtype=object)) == target_name
    ]
    if not occ_match.empty:
        return occ_match.head(1)

    return None


def _read_visibility_extended(
    directory: Path, name: str,
) -> Optional[pd.DataFrame]:
    """Read visibility parquet with boresight separation columns.

    Returns a DataFrame with Time, Visible, Sun_Sep, Moon_Sep, Earth_Sep
    (when available).  Falls back gracefully when separation columns are
    absent (older parquet files).
    """
    path = directory / f"Visibility for {name}.parquet"
    cols = [
        "Time(MJD_UTC)", "Time_UTC", "Visible",
        "Sun_Sep", "Moon_Sep", "Earth_Sep",
    ]
    df = read_parquet_cached(str(path), columns=cols)
    if df is None:
        # Fall back to minimal columns (older parquet).
        df = read_parquet_cached(
            str(path), columns=["Time(MJD_UTC)", "Visible"],
        )
    return df


def _read_visibility(directory: Path, name: str) -> Optional[pd.DataFrame]:
    """Read star visibility file with caching."""
    path = directory / f"Visibility for {name}.parquet"
    df = read_parquet_cached(
        str(path),
        columns=["Time(MJD_UTC)", "Time_UTC", "Visible"],
    )
    if df is None:
        # Older visibility parquet fixtures (and some historical outputs) may not
        # include Time_UTC. Fall back to MJD-only visibility timeline.
        df = read_parquet_cached(
            str(path),
            columns=["Time(MJD_UTC)", "Visible"],
        )
    if df is None:
        LOGGER.debug("Visibility file missing for %s", name)
        return None
    # Keep the source path for later diagnostics/logging.
    df.attrs["_source_path"] = str(path)
    if df.empty:
        LOGGER.debug("DF is empty for %s", name)
    return df


def _visibility_diagnostics(
    visibility_df: pd.DataFrame,
    start: datetime,
    stop: datetime,
) -> tuple[str, str, str, str]:
    """Return (source_path, time_min, time_max, samples_in_window) for warnings."""
    source = str(visibility_df.attrs.get("_source_path", "<unknown>"))
    if "Time_UTC" in visibility_df.columns and pd.api.types.is_datetime64_any_dtype(
        visibility_df["Time_UTC"]
    ):
        times = visibility_df["Time_UTC"].to_numpy(dtype="datetime64[ns]")
        if times.size == 0:
            return source, "<empty>", "<empty>", "0"
        time_min = pd.to_datetime(times.min()).to_pydatetime().isoformat(sep=" ")
        time_max = pd.to_datetime(times.max()).to_pydatetime().isoformat(sep=" ")
        mask = (times >= np.datetime64(start)) & (times <= np.datetime64(stop))
        in_window = int(mask.sum())
        return source, time_min, time_max, str(in_window)

    if "Time(MJD_UTC)" in visibility_df.columns:
        mjd = visibility_df["Time(MJD_UTC)"].to_numpy(dtype=float)
        if mjd.size == 0:
            return source, "<empty>", "<empty>", "0"
        mjd_min = float(np.nanmin(mjd))
        mjd_max = float(np.nanmax(mjd))
        time_min = (
            Time(mjd_min, format="mjd", scale="utc").to_datetime().isoformat(sep=" ")
        )
        time_max = (
            Time(mjd_max, format="mjd", scale="utc").to_datetime().isoformat(sep=" ")
        )
        start_mjd = Time(start).mjd
        stop_mjd = Time(stop).mjd
        in_window = int(((mjd >= start_mjd) & (mjd <= stop_mjd)).sum())
        return source, time_min, time_max, str(in_window)

    return source, "<unknown>", "<unknown>", "<unknown>"


def _read_planet_visibility(directory: Path, name: str) -> Optional[pd.DataFrame]:
    """Read planet transit-visibility file with caching.

    Planet visibility parquet files contain transit windows (MJD start/stop), not a
    per-timestep visibility timeline. Therefore we load the transit window columns.
    """
    path = directory / f"Visibility for {name}.parquet"
    df = read_parquet_cached(
        str(path),
        columns=["Transit_Start", "Transit_Stop"],
    )
    if df is None:
        LOGGER.debug("Planet visibility file missing for %s", name)
    return df


def _extract_visibility_segment(
    visibility_df: pd.DataFrame,
    start: datetime,
    stop: datetime,
    min_sequence_minutes: int,
) -> tuple[List[datetime], List[int]]:
    if "Time_UTC" in visibility_df.columns and pd.api.types.is_datetime64_any_dtype(
        visibility_df["Time_UTC"]
    ):
        times_dt64 = visibility_df["Time_UTC"].to_numpy(dtype="datetime64[ns]")
        start_dt64 = np.datetime64(start)
        stop_dt64 = np.datetime64(stop)
        mask = (times_dt64 >= start_dt64) & (times_dt64 <= stop_dt64)
        if not bool(mask.any()):
            return [], []

        window_indices = np.flatnonzero(mask)

        # Round to nearest second with legacy half-up semantics.
        # (pandas .round('S') uses bankers rounding; we need >=0.5s to round up.)
        ns = times_dt64[window_indices].astype("datetime64[ns]").view("int64")
        rounded_ns = ((ns + 500_000_000) // 1_000_000_000) * 1_000_000_000
        rounded_dt64 = rounded_ns.view("datetime64[ns]")
        visit_times = pd.to_datetime(rounded_dt64).to_pydatetime().tolist()
    else:
        raw_times = Time(
            visibility_df["Time(MJD_UTC)"].to_numpy(),
            format="mjd",
            scale="utc",
        ).to_datetime()
        # Normalise astropy datetimes to naive Python datetimes (UTC) so comparisons
        # with schedule start/stop (which are naive datetimes) behave predictably.
        raw_times = [
            (rt.replace(tzinfo=None) if getattr(rt, "tzinfo", None) is not None else rt)
            for rt in raw_times
        ]

        mask = [start <= value <= stop for value in raw_times]
        if not any(mask):
            return [], []

        window_indices = [idx for idx, include in enumerate(mask) if include]
        filtered_times = [raw_times[idx] for idx in window_indices]
        # Round each time to nearest second
        visit_times = [
            (t + timedelta(microseconds=500_000)).replace(microsecond=0)
            for t in filtered_times
        ]

    flags = [int(float(visibility_df.iloc[int(idx)]["Visible"])) for idx in window_indices]
    return visit_times, flags


def _visibility_change_indices(flags: Sequence[int]) -> List[int]:
    return [idx for idx in range(len(flags) - 1) if flags[idx] != flags[idx + 1]]


def _merge_short_occultation_segments(
    starts: Sequence[datetime],
    stops: Sequence[datetime],
    min_sequence_minutes: int,
) -> tuple[List[datetime], List[datetime]]:
    """Merge occultation segments shorter than *min_sequence_minutes*.

    Merge policy per contiguous run:
    - short segment at run start -> merge forward
    - short segment at run end -> merge backward
    - isolated short segment -> drop
    """
    if not starts or not stops:
        return [], []
    if min_sequence_minutes <= 0:
        return list(starts), list(stops)

    threshold = timedelta(minutes=min_sequence_minutes)
    ordered = sorted(zip(starts, stops), key=lambda item: item[0])

    # Group segments into contiguous runs (allow tiny boundary jitter).
    runs: List[List[List[datetime]]] = []
    current_run: List[List[datetime]] = []
    adjacency_tolerance = timedelta(seconds=1)

    for start, stop in ordered:
        if stop <= start:
            continue
        if not current_run:
            current_run = [[start, stop]]
            continue
        if start <= current_run[-1][1] + adjacency_tolerance:
            current_run.append([start, stop])
        else:
            runs.append(current_run)
            current_run = [[start, stop]]
    if current_run:
        runs.append(current_run)

    merged: List[tuple[datetime, datetime]] = []
    dropped_short_isolated = 0
    for run in runs:
        if len(run) == 1:
            seg_start, seg_stop = run[0]
            if (seg_stop - seg_start) < threshold:
                dropped_short_isolated += 1
                continue

        # Iteratively merge short boundary segments into neighbours.
        changed = True
        while changed and len(run) > 1:
            changed = False
            for idx_seg, (seg_start, seg_stop) in enumerate(run):
                if (seg_stop - seg_start) >= threshold:
                    continue
                if idx_seg == 0:
                    run[1][0] = seg_start
                    del run[0]
                    changed = True
                    break
                if idx_seg == len(run) - 1:
                    run[idx_seg - 1][1] = seg_stop
                    del run[idx_seg]
                    changed = True
                    break
                run[idx_seg - 1][1] = seg_stop
                del run[idx_seg]
                changed = True
                break

        merged.extend((segment[0], segment[1]) for segment in run)

    if dropped_short_isolated > 0:
        LOGGER.info(
            "Dropped %d isolated occultation segment(s) shorter than %d min",
            dropped_short_isolated,
            min_sequence_minutes,
        )

    if not merged:
        return [], []
    return [item[0] for item in merged], [item[1] for item in merged]


def _occultation_windows(
    visit_times: Sequence[datetime],
    visibility_flags: Sequence[int],
    visibility_changes: Sequence[int],
) -> tuple[List[datetime], List[datetime], List[int]]:
    changes = list(visibility_changes)
    flags = list(visibility_flags)
    times = list(visit_times)

    if not flags:
        return [], [], []

    if flags[-1] == 0 and len(times) >= 2:
        changes.append(len(times) - 2)

    occ_starts: List[datetime] = []
    occ_stops: List[datetime] = []

    if flags[0] == 0 and changes:
        occ_starts.append(times[0])
        occ_stops.append(times[changes[0]])

    for idx in range(len(changes) - 1):
        change_idx = changes[idx]
        next_idx = changes[idx + 1]
        if flags[next_idx] == 0:
            occ_starts.append(times[change_idx + 1])
            occ_stops.append(times[next_idx])

    if flags[-1] == 1 and len(times) >= 2:
        changes.append(len(times) - 2)

    if not occ_starts:
        return [], [], changes

    if len(occ_starts) != len(occ_stops):
        raise ValueError("Occultation start/stop lists are mismatched")

    # Remove degenerate windows produced by boundary/rounding effects.
    filtered_pairs = [
        (start, stop) for start, stop in zip(occ_starts, occ_stops) if stop > start
    ]
    if len(filtered_pairs) != len(occ_starts):
        LOGGER.debug(
            "Dropped %d degenerate occultation window(s) at extraction",
            len(occ_starts) - len(filtered_pairs),
        )
    if not filtered_pairs:
        return [], [], changes
    return (
        [pair[0] for pair in filtered_pairs],
        [pair[1] for pair in filtered_pairs],
        changes,
    )


def _prioritise_occultation_targets(
    occ_list: pd.DataFrame,
    reference_ra: float,
    reference_dec: float,
) -> pd.DataFrame:
    if occ_list.empty:
        return occ_list

    if not (np.isfinite(reference_ra) and np.isfinite(reference_dec)):
        return occ_list

    if "RA" not in occ_list.columns or "DEC" not in occ_list.columns:
        return occ_list

    ra_values = pd.to_numeric(occ_list["RA"], errors="coerce")
    dec_values = pd.to_numeric(occ_list["DEC"], errors="coerce")
    if ra_values.isna().all() or dec_values.isna().all():
        return occ_list

    try:
        origin = SkyCoord(ra=reference_ra, dec=reference_dec, unit="deg")
        target_coords = SkyCoord(
            ra=ra_values.to_numpy(),
            dec=dec_values.to_numpy(),
            unit="deg",
        )
        separations = origin.separation(target_coords).deg
    except Exception as exc:  # SkyCoord failure should not abort scheduling
        LOGGER.debug("Unable to rank occultation targets by slew distance: %s", exc)
        return occ_list

    priorities = np.where(np.isfinite(separations), separations, np.inf)
    reordered = (
        occ_list.assign(_separation=priorities)
        .sort_values(by="_separation", kind="mergesort")
        .drop(columns="_separation")
        .reset_index(drop=True)
    )
    return reordered


def _build_occultation_schedule(
    starts: Sequence[datetime],
    stops: Sequence[datetime],
    visit_start: datetime,
    visit_stop: datetime,
    list_path: Path,
    vis_root: Path,
    label: str,
    reference_ra: float,
    reference_dec: float,
    prioritise_by_slew: bool,
    excluded_targets: Optional[set] = None,
    show_progress: bool = False,
    use_pass1: bool = True,
) -> tuple[Optional[pd.DataFrame], bool]:
    if not starts or not stops:
        return None, False

    # Guard against degenerate windows introduced by boundary rounding.
    filtered_pairs = [
        (start, stop) for start, stop in zip(starts, stops) if stop > start
    ]
    dropped_intervals = len(list(zip(starts, stops))) - len(filtered_pairs)
    if dropped_intervals > 0:
        LOGGER.info(
            "%s..%s: dropped %d degenerate occultation interval(s)",
            visit_start,
            visit_stop,
            dropped_intervals,
        )
    if not filtered_pairs:
        return None, False
    starts = [pair[0] for pair in filtered_pairs]
    stops = [pair[1] for pair in filtered_pairs]

    schedule_rows = [
        [
            "",
            start.strftime("%Y-%m-%dT%H:%M:%SZ"),
            stop.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "",
            "",
        ]
        for start, stop in zip(starts, stops)
    ]
    occ_df = pd.DataFrame(
        schedule_rows, columns=["Target", "start", "stop", "RA", "DEC"]
    )

    try:
        occ_list = read_csv_cached(str(list_path))
    except FileNotFoundError:
        LOGGER.warning("Occultation list missing: %s", list_path)
        return None, False

    # Filter out excluded targets (those that have exceeded their time limit)
    if excluded_targets and "Star Name" in occ_list.columns:
        before_count = len(occ_list)
        occ_list = occ_list[~occ_list["Star Name"].isin(excluded_targets)].reset_index(
            drop=True
        )
        if len(occ_list) < before_count:
            LOGGER.debug(
                "Excluded %d occultation targets that exceeded time limit",
                before_count - len(occ_list),
            )

    # If no targets remain after exclusion, fail early
    if occ_list.empty:
        LOGGER.warning("No occultation targets available after exclusion filter")
        return None, False

    if prioritise_by_slew:
        occ_list = _prioritise_occultation_targets(
            occ_list,
            reference_ra,
            reference_dec,
        )

    target_names = occ_list.get("Star Name", pd.Series(dtype=object)).to_numpy()
    starts_mjd = Time(list(starts), format="datetime", scale="utc").to_value("mjd")
    stops_mjd = Time(list(stops), format="datetime", scale="utc").to_value("mjd")

    occ_df, flag = observation_utils.schedule_occultation_targets(
        target_names,
        starts_mjd,
        stops_mjd,
        visit_start,
        visit_stop,
        str(vis_root),
        occ_df,
        occ_list,
        label,
        show_progress=show_progress,
        use_pass1=use_pass1,
    )
    return occ_df, flag


def _transit_windows(
    transit_df: Optional[pd.DataFrame],
) -> Optional[tuple[List[datetime], List[datetime]]]:
    if transit_df is None or transit_df.empty:
        return None

    start_times = Time(
        transit_df["Transit_Start"].to_numpy(),
        format="mjd",
        scale="utc",
    ).to_datetime()
    stop_times = Time(
        transit_df["Transit_Stop"].to_numpy(),
        format="mjd",
        scale="utc",
    ).to_datetime()

    # Round each time to nearest second to match legacy
    # `round_to_nearest_second` behaviour.
    start = [
        (t + timedelta(microseconds=500_000)).replace(microsecond=0)
        for t in start_times
    ]
    stop = [
        (t + timedelta(microseconds=500_000)).replace(microsecond=0) for t in stop_times
    ]
    return start, stop


def _target_priority(
    priority_flag: bool,
    transit_start: Sequence[datetime],
    transit_stop: Sequence[datetime],
    sequence_start: datetime,
    sequence_stop: datetime,
) -> str:
    if not priority_flag or not transit_start or not transit_stop:
        return "0"

    for start, stop in zip(transit_start, transit_stop):
        if start <= sequence_stop and stop >= sequence_start:
            return "2"
    return "1"


def _fallback_float(value: object, info: Optional[pd.DataFrame], column: str) -> float:
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            pass

    if info is not None and column in info.columns:
        try:
            candidate = info.iloc[0][column]
            if isinstance(candidate, (int, float, np.integer, np.floating)):
                return float(candidate)
            if isinstance(candidate, str):
                return float(candidate)
        except (KeyError, ValueError, TypeError):
            pass

    return float("nan")


def _resolve_coordinates(star_name: str) -> tuple[float, float]:
    """Raise error if coordinates not in catalog - no Simbad lookups allowed."""
    raise RuntimeError(f"No coordinates found in catalog for {star_name}")


def _serialise_calendar(root: ET.Element) -> str:
    _convert_numeric_content(root)
    xml_bytes = ET.tostring(root, encoding="utf-8", xml_declaration=True)
    xml_doc = minidom.parseString(xml_bytes)
    return xml_doc.toprettyxml(indent="\t")


def _convert_numeric_content(element: ET.Element) -> None:
    for child in element.iter():
        for key, value in list(child.attrib.items()):
            if isinstance(value, Number):
                child.set(key, str(value))
        if child.text is not None and isinstance(child.text, Number):
            child.text = str(child.text)
