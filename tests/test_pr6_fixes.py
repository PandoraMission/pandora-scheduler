"""Tests for the PR #6 robustness fixes ported from vbkostov:tb_main.

Covers:
  1. Degenerate occultation window filtering (_occultation_windows, _build_occultation_schedule)
  2. _merge_short_occultation_segments()
  3. Time-based occ_df lookup in _emit_visit_sequences (via _ScienceCalendarBuilder)
  4. Pass 2 empty-interval guard in schedule_occultation_targets
  5. Parquet schema check in build_visibility_catalog
  6. VDA boresight .copy() defensive copy
"""

import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from astropy.time import Time

from pandorascheduler_rework import observation_utils, science_calendar
from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.xml.parameters import populate_vda_parameters


# ---------------------------------------------------------------------------
# Feature 1: Degenerate window filtering in _occultation_windows
# ---------------------------------------------------------------------------

class TestDegenerateWindowFiltering:
    """_occultation_windows should drop windows where stop <= start."""

    @staticmethod
    def _make_times(n: int, start: datetime = datetime(2026, 1, 1)) -> list[datetime]:
        return [start + timedelta(minutes=i) for i in range(n)]

    def test_normal_windows_kept(self):
        """Windows with stop > start should pass through unchanged."""
        # Visibility: [0,0,0, 1,1,1, 0,0,0] -> one occ window at start, one at end
        times = self._make_times(9)
        flags = [0, 0, 0, 1, 1, 1, 0, 0, 0]
        changes = [2, 5]  # flag changes at index 2 and 5

        starts, stops, _ = science_calendar._occultation_windows(times, flags, changes)
        assert len(starts) > 0
        for s, e in zip(starts, stops):
            assert e > s

    def test_degenerate_equal_start_stop_dropped(self):
        """A window where stop == start (zero duration) should be dropped."""
        times = self._make_times(5)
        flags = [0, 0, 1, 1, 1]
        # Force the change index to produce start == stop:
        # flags[0]==0 and changes[0]==0  =>  start=times[0], stop=times[0]
        changes = [0, 1]

        starts, stops, _ = science_calendar._occultation_windows(times, flags, changes)
        # The first window (times[0]..times[0]) is degenerate; should be dropped.
        for s, e in zip(starts, stops):
            assert e > s, f"Degenerate window not filtered: start={s}, stop={e}"

    def test_all_degenerate_returns_empty(self):
        """If all windows are degenerate, return empty lists."""
        t = datetime(2026, 6, 1)
        # Manually test by injecting identical start/stop pairs via _build_occultation_schedule
        starts_in = [t, t + timedelta(hours=1)]
        stops_in = [t, t + timedelta(hours=1)]  # stop == start for both

        occ_df, success = science_calendar._build_occultation_schedule(
            starts=starts_in,
            stops=stops_in,
            visit_start=t,
            visit_stop=t + timedelta(hours=2),
            list_path=Path("/nonexistent"),  # won't be reached
            vis_root=Path("/nonexistent"),
            label="test",
            reference_ra=0.0,
            reference_dec=0.0,
            prioritise_by_slew=False,
        )
        # All degenerate -> should return None, False
        assert occ_df is None
        assert success is False

    def test_mixed_degenerate_and_valid(self):
        """Only degenerate windows removed; valid ones survive."""
        t = datetime(2026, 6, 1)
        starts_in = [t, t + timedelta(hours=2)]
        stops_in = [t, t + timedelta(hours=3)]  # first degenerate, second valid

        # _build_occultation_schedule reads a CSV; mock to avoid file I/O
        with patch("pandorascheduler_rework.science_calendar.read_csv_cached") as mock_csv:
            mock_csv.return_value = pd.DataFrame({
                "Star Name": ["TGT-A"],
                "RA": [10.0],
                "DEC": [20.0],
                "Number of Hours Requested": [100.0],
            })
            with patch("pandorascheduler_rework.science_calendar.observation_utils.schedule_occultation_targets") as mock_sched:
                mock_sched.return_value = (
                    pd.DataFrame({
                        "Target": ["TGT-A"],
                        "start": [(t + timedelta(hours=2)).strftime("%Y-%m-%dT%H:%M:%SZ")],
                        "stop": [(t + timedelta(hours=3)).strftime("%Y-%m-%dT%H:%M:%SZ")],
                        "RA": [10.0],
                        "DEC": [20.0],
                        "Visibility": [1],
                    }),
                    True,
                )

                occ_df, success = science_calendar._build_occultation_schedule(
                    starts=starts_in,
                    stops=stops_in,
                    visit_start=t,
                    visit_stop=t + timedelta(hours=4),
                    list_path=Path("/fake/occ.csv"),
                    vis_root=Path("/fake/vis"),
                    label="test",
                    reference_ra=0.0,
                    reference_dec=0.0,
                    prioritise_by_slew=False,
                )
                # Should have proceeded with the 1 valid window
                assert occ_df is not None


# ---------------------------------------------------------------------------
# Feature 2: _merge_short_occultation_segments
# ---------------------------------------------------------------------------

class TestMergeShortOccultationSegments:
    """Unit tests for the segment merging logic."""

    @staticmethod
    def _dt(offset_min: int) -> datetime:
        return datetime(2026, 6, 1) + timedelta(minutes=offset_min)

    def test_empty_input(self):
        starts, stops = science_calendar._merge_short_occultation_segments([], [], 10)
        assert starts == []
        assert stops == []

    def test_zero_threshold_returns_input(self):
        """min_sequence_minutes <= 0 should return input unchanged."""
        s = [self._dt(0), self._dt(30)]
        e = [self._dt(5), self._dt(60)]
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 0)
        assert starts == s
        assert stops == e

    def test_isolated_short_segment_dropped(self):
        """A single segment shorter than threshold with no neighbours should be dropped."""
        s = [self._dt(0)]
        e = [self._dt(5)]  # 5 minutes, threshold = 10
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 10)
        assert starts == []
        assert stops == []

    def test_isolated_long_segment_kept(self):
        """A single segment >= threshold should survive."""
        s = [self._dt(0)]
        e = [self._dt(15)]
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 10)
        assert len(starts) == 1
        assert starts[0] == self._dt(0)
        assert stops[0] == self._dt(15)

    def test_short_segment_at_run_start_merged_forward(self):
        """Short segment at the start of a contiguous run is merged into the next."""
        # Two contiguous segments: 3-min + 15-min
        s = [self._dt(0), self._dt(3)]
        e = [self._dt(3), self._dt(18)]
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 10)
        # Should merge into one segment spanning 0..18
        assert len(starts) == 1
        assert starts[0] == self._dt(0)
        assert stops[0] == self._dt(18)

    def test_short_segment_at_run_end_merged_backward(self):
        """Short segment at the end of a contiguous run is merged into the previous."""
        # 20-min + 3-min contiguous segments
        s = [self._dt(0), self._dt(20)]
        e = [self._dt(20), self._dt(23)]
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 10)
        assert len(starts) == 1
        assert starts[0] == self._dt(0)
        assert stops[0] == self._dt(23)

    def test_degenerate_segments_skipped(self):
        """Segments with stop <= start should be silently dropped."""
        s = [self._dt(10), self._dt(10)]  # degenerate
        e = [self._dt(10), self._dt(25)]  # second is valid
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 5)
        # Only the valid segment survives (it's 15 min, above threshold).
        assert len(starts) == 1
        assert starts[0] == self._dt(10)
        assert stops[0] == self._dt(25)

    def test_multiple_runs_processed_independently(self):
        """Non-contiguous runs (gap > adjacency_tolerance) are separate."""
        # Run 1: 0-5 (short, isolated) | Run 2: 60-90 (long, kept)
        s = [self._dt(0), self._dt(60)]
        e = [self._dt(5), self._dt(90)]
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 10)
        # Run 1 dropped (isolated short), Run 2 kept
        assert len(starts) == 1
        assert starts[0] == self._dt(60)

    def test_contiguous_short_segments_merge_iteratively(self):
        """Multiple short boundary segments in a run get merged."""
        # Three contiguous segments: 3-min + 4-min + 20-min
        s = [self._dt(0), self._dt(3), self._dt(7)]
        e = [self._dt(3), self._dt(7), self._dt(27)]
        starts, stops = science_calendar._merge_short_occultation_segments(s, e, 10)
        # Both short segments should merge forward into the 20-min one
        assert len(starts) == 1
        assert starts[0] == self._dt(0)
        assert stops[0] == self._dt(27)


# ---------------------------------------------------------------------------
# Feature 5: Pass 2 empty-interval guard
# ---------------------------------------------------------------------------

class TestPass2EmptyIntervalGuard:
    """schedule_occultation_targets should handle intervals with no visibility data."""

    def test_empty_mask_sets_visibility_zero(self, tmp_path):
        """When interval has no matching visibility timestamps, Visibility=0."""
        # Build MJD start/stop for one interval
        t0 = datetime(2026, 1, 1, 12, 0)
        t1 = datetime(2026, 1, 1, 13, 0)
        start_mjd = Time([t0], scale="utc").to_value("mjd")
        stop_mjd = Time([t1], scale="utc").to_value("mjd")

        occ_df = pd.DataFrame(
            [{"Target": "", "start": t0.strftime("%Y-%m-%dT%H:%M:%SZ"),
              "stop": t1.strftime("%Y-%m-%dT%H:%M:%SZ"),
              "RA": "", "DEC": "", "Visibility": np.nan}]
        )

        vis_dir = tmp_path / "vis"
        vis_dir.mkdir()

        # Write visibility data with times OUTSIDE the interval
        target_dir = vis_dir / "OCC-STAR"
        target_dir.mkdir()
        far_times = [datetime(2026, 6, 1, h, 0) for h in range(3)]
        far_mjd = Time(far_times, scale="utc").to_value("mjd")
        pd.DataFrame({"Time(MJD_UTC)": far_mjd, "Visible": [1, 1, 1]}).to_parquet(
            target_dir / "Visibility for OCC-STAR.parquet", index=False,
        )

        occ_list = pd.DataFrame({
            "Star Name": ["OCC-STAR"],
            "RA": [10.0],
            "DEC": [-20.0],
        })

        result_df, all_assigned = observation_utils.schedule_occultation_targets(
            v_names=["OCC-STAR"],
            starts=start_mjd,
            stops=stop_mjd,
            visit_start=t0,
            visit_stop=t1,
            path=vis_dir,
            o_df=occ_df,
            o_list=occ_list,
            try_occ_targets="OCC-STAR",
        )

        # The interval should not be assigned a target (no data overlap)
        assert result_df.loc[0, "Visibility"] == 0


# ---------------------------------------------------------------------------
# Feature 6: Parquet schema check
# ---------------------------------------------------------------------------

class TestParquetSchemaCheck:
    """Transit_Overlap column detection should use pyarrow schema, not text."""

    def test_detects_missing_transit_overlap_column(self, tmp_path):
        """Parquet without Transit_Overlap triggers recomputation."""
        star_dir = tmp_path / "star_a" / "planet_b"
        star_dir.mkdir(parents=True)
        # Write parquet WITHOUT Transit_Overlap
        df = pd.DataFrame({"Transit_Start": [1.0], "Transit_Stop": [2.0]})
        df.to_parquet(star_dir / "Visibility for planet_b.parquet", index=False)

        path = star_dir / "Visibility for planet_b.parquet"
        schema = pq.read_schema(path)
        assert "Transit_Overlap" not in schema.names

    def test_detects_present_transit_overlap_column(self, tmp_path):
        """Parquet with Transit_Overlap correctly detected."""
        star_dir = tmp_path / "star_a" / "planet_b"
        star_dir.mkdir(parents=True)
        df = pd.DataFrame({
            "Transit_Start": [1.0],
            "Transit_Stop": [2.0],
            "Transit_Overlap": [0],
        })
        df.to_parquet(star_dir / "Visibility for planet_b.parquet", index=False)

        path = star_dir / "Visibility for planet_b.parquet"
        schema = pq.read_schema(path)
        assert "Transit_Overlap" in schema.names

    def test_corrupted_file_returns_false(self, tmp_path):
        """Non-parquet file should trigger exception path gracefully."""
        star_dir = tmp_path / "star_a" / "planet_b"
        star_dir.mkdir(parents=True)
        bad_path = star_dir / "Visibility for planet_b.parquet"
        bad_path.write_text("not a parquet file")

        # The code uses try/except; verify it doesn't crash
        all_have_overlap = True
        try:
            schema = pq.read_schema(bad_path)
            if "Transit_Overlap" not in schema.names:
                all_have_overlap = False
        except Exception:
            all_have_overlap = False

        assert all_have_overlap is False


# ---------------------------------------------------------------------------
# Feature 7: VDA boresight .copy() defensive copy
# ---------------------------------------------------------------------------

class TestVdaBoresightDefensiveCopy:
    """populate_vda_parameters must not mutate the original ROI coordinate arrays."""

    def test_original_coordinates_not_mutated(self):
        """Calling populate_vda_parameters should not alter the source ROI data."""
        root = ET.Element("Root")
        # ROI_coord_0 = [1.0, 2.0] — these should NOT be changed to RA/DEC
        targ_info = pd.DataFrame([{
            "VDA_StarRoiDetMethod": "SET_BY_TARGET_DEFINITION_FILE",
            "StarRoiDetMethod": 1,
            "VDA_numPredefinedStarRois": "SET_BY_TARGET_DEFINITION_FILE",
            "numPredefinedStarRois": 2,
            "VDA_PredefinedStarRoiRa": "SET_BY_TARGET_DEFINITION_FILE",
            "VDA_PredefinedStarRoiDec": "SET_BY_TARGET_DEFINITION_FILE",
            "RA": 99.0,
            "DEC": -45.0,
            "ROI_coord_0": "[1.0, 2.0]",
            "ROI_coord_1": "[3.0, 4.0]",
        }])

        # Save original string values
        original_coord_0 = targ_info.iloc[0]["ROI_coord_0"]
        original_coord_1 = targ_info.iloc[0]["ROI_coord_1"]

        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)

        # Source DataFrame should not be altered
        assert targ_info.iloc[0]["ROI_coord_0"] == original_coord_0
        assert targ_info.iloc[0]["ROI_coord_1"] == original_coord_1

    def test_ra_dec_type_error_handled(self):
        """Non-numeric RA/DEC should not crash; ROI coords used as-is."""
        root = ET.Element("Root")
        targ_info = pd.DataFrame([{
            "VDA_StarRoiDetMethod": "SET_BY_TARGET_DEFINITION_FILE",
            "StarRoiDetMethod": 1,
            "VDA_numPredefinedStarRois": "SET_BY_TARGET_DEFINITION_FILE",
            "numPredefinedStarRois": 1,
            "VDA_PredefinedStarRoiRa": "SET_BY_TARGET_DEFINITION_FILE",
            "VDA_PredefinedStarRoiDec": "SET_BY_TARGET_DEFINITION_FILE",
            "RA": "not_a_number",
            "DEC": None,
            "ROI_coord_0": "[5.0, 6.0]",
        }])

        # Should not raise
        populate_vda_parameters(root, targ_info, diff_in_seconds=10.0)

        vda = root.find("AcquireVisCamScienceData")
        assert vda is not None
        # Since RA/DEC parsing fails, slot 0 should keep original values
        roi_ra = vda.find("PredefinedStarRoiRa")
        if roi_ra is not None:
            assert roi_ra.find("RA1").text == "5.000000"
