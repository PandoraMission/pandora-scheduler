"""Tests for science_calendar._ScienceCalendarBuilder internals."""

from datetime import datetime, timedelta
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.science_calendar import (
    ScienceCalendarInputs,
    _ScienceCalendarBuilder,
    _extract_visibility_segment,
    _is_transit_entry,
    _lookup_auxiliary_row,
    _lookup_planet_row,
    _normalise_target_name,
    _parse_datetime,
    _read_catalog,
    _split_at_priority_buffer_boundaries,
    _target_priority,
    _read_visibility_extended,
    _visibility_change_indices,
    generate_science_calendar,
)


def _seed_catalogs(data_dir: Path) -> None:
    """Create minimal catalog CSVs so _ScienceCalendarBuilder.__init__ succeeds."""
    data_dir.mkdir(exist_ok=True)
    pd.DataFrame({"Planet Name": [], "Star Name": []}).to_csv(
        data_dir / "exoplanet_targets.csv", index=False
    )
    pd.DataFrame({"Star Name": []}).to_csv(
        data_dir / "all_targets.csv", index=False
    )
    pd.DataFrame({"Star Name": []}).to_csv(
        data_dir / "occultation-standard_targets.csv", index=False
    )


# ---------------------------------------------------------------------------
# _ScienceCalendarBuilder.__init__ / build_calendar
# ---------------------------------------------------------------------------


class TestBuilderInit:
    """Tests for _ScienceCalendarBuilder initialization."""

    def _make_schedule(self, tmp_path, rows=None):
        csv = tmp_path / "sched.csv"
        if rows is None:
            rows = [
                {
                    "Target": "Free Time",
                    "Observation Start": "2026-03-01T00:00:00Z",
                    "Observation Stop": "2026-03-01T06:00:00Z",
                }
            ]
        pd.DataFrame(rows).to_csv(csv, index=False)
        return csv

    def test_raises_on_missing_csv(self, tmp_path):
        _seed_catalogs(tmp_path)
        inputs = ScienceCalendarInputs(
            schedule_csv=tmp_path / "nope.csv",
            data_dir=tmp_path,
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
        )
        with pytest.raises(FileNotFoundError, match="Schedule CSV missing"):
            _ScienceCalendarBuilder(inputs, config)

    def test_raises_on_empty_csv(self, tmp_path):
        _seed_catalogs(tmp_path)
        csv = tmp_path / "empty.csv"
        pd.DataFrame(columns=["Target"]).to_csv(csv, index=False)
        inputs = ScienceCalendarInputs(schedule_csv=csv, data_dir=tmp_path)
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
        )
        with pytest.raises(ValueError, match="empty"):
            _ScienceCalendarBuilder(inputs, config)

    def test_init_loads_schedule(self, tmp_path):
        _seed_catalogs(tmp_path)
        csv = self._make_schedule(tmp_path)
        inputs = ScienceCalendarInputs(schedule_csv=csv, data_dir=tmp_path)
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
        )
        builder = _ScienceCalendarBuilder(inputs, config)
        assert builder.schedule is not None
        assert len(builder.schedule) == 1


# ---------------------------------------------------------------------------
# _add_meta
# ---------------------------------------------------------------------------


class TestAddMeta:
    """Tests for XML Meta element construction."""

    def test_meta_contains_expected_attributes(self, tmp_path):
        _seed_catalogs(tmp_path)
        csv = tmp_path / "sched.csv"
        pd.DataFrame(
            [
                {
                    "Target": "Star b",
                    "Observation Start": "2026-03-01T00:00:00Z",
                    "Observation Stop": "2026-03-02T00:00:00Z",
                }
            ]
        ).to_csv(csv, index=False)

        inputs = ScienceCalendarInputs(schedule_csv=csv, data_dir=tmp_path)
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
            created_timestamp="2026-03-01T00:00:00Z",
            author="test",
        )
        builder = _ScienceCalendarBuilder(inputs, config)

        root = ET.Element("ScienceCalendar")
        builder._add_meta(root)

        meta = root.find("Meta")
        assert meta is not None
        assert meta.get("Valid_From") == "2026-03-01T00:00:00Z"
        assert meta.get("Created") == "2026-03-01T00:00:00Z"
        assert meta.get("Author") == "test"
        assert "Calendar_Weights" in meta.attrib


# ---------------------------------------------------------------------------
# _read_catalog
# ---------------------------------------------------------------------------


class TestReadCatalog:
    """Tests for catalog reading."""

    def test_raises_for_missing(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="Required target catalog missing"):
            _read_catalog(tmp_path / "nope.csv")

    def test_reads_existing_csv(self, tmp_path):
        path = tmp_path / "cat.csv"
        pd.DataFrame({"Star Name": ["Alpha"]}).to_csv(path, index=False)
        result = _read_catalog(path)
        assert result is not None
        assert len(result) == 1


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


class TestHelperFunctions:
    """Tests for module-level helper functions."""

    def test_normalise_target_name_planet(self):
        name, star = _normalise_target_name("TOI-700 b")
        assert name == "TOI-700 b"
        assert star == "TOI-700"

    def test_normalise_target_name_no_suffix(self):
        name, star = _normalise_target_name("Vega")
        assert name == "Vega"
        assert star == "Vega"

    def test_normalise_target_name_star_exceptions(self):
        """Star names whose trailing letter looks like a planet suffix."""
        for star in ("AU_Mic", "EV_Lac", "AF_Psc", "55_Cnc", "AO_Cassiopeiae"):
            name, sname = _normalise_target_name(star)
            assert name == star, f"{star}: target_name mangled"
            assert sname == star, f"{star}: star_name mangled to {sname!r}"

    def test_parse_datetime_iso(self):
        dt = _parse_datetime("2026-03-01T12:30:00Z")
        assert dt is not None
        assert dt.hour == 12

    def test_parse_datetime_none(self):
        assert _parse_datetime(None) is None

    def test_parse_datetime_already_datetime(self):
        dt = datetime(2026, 3, 1, 12, 30)
        assert _parse_datetime(dt) == dt

    def test_is_transit_entry_true(self):
        row = pd.Series({"Transit Coverage": 0.85})
        assert _is_transit_entry(row) == True

    def test_is_transit_entry_false(self):
        row = pd.Series({"Transit Coverage": None})
        assert _is_transit_entry(row) == False

    def test_is_transit_entry_numeric_1(self):
        row = pd.Series({"Transit Coverage": 1.0})
        assert _is_transit_entry(row) == True

    def test_target_priority_without_buffer_does_not_promote_adjacent_sequence(self):
        transit_start = [datetime(2026, 3, 1, 12, 0)]
        transit_stop = [datetime(2026, 3, 1, 13, 0)]
        priority = _target_priority(
            True,
            transit_start,
            transit_stop,
            datetime(2026, 3, 1, 10, 45),
            datetime(2026, 3, 1, 11, 30),
        )
        assert priority == "1"

    def test_target_priority_with_buffer_promotes_adjacent_sequence(self):
        transit_start = [datetime(2026, 3, 1, 12, 0)]
        transit_stop = [datetime(2026, 3, 1, 13, 0)]
        priority = _target_priority(
            True,
            transit_start,
            transit_stop,
            datetime(2026, 3, 1, 10, 45),
            datetime(2026, 3, 1, 11, 30),
            buffer_enabled=True,
            buffer_minutes=80,
        )
        assert priority == "2"

    def test_priority_buffer_boundaries_produce_exact_priority_regions(self):
        transit_start = [datetime(2026, 3, 1, 12, 0)]
        transit_stop = [datetime(2026, 3, 1, 13, 0)]
        regions = _split_at_priority_buffer_boundaries(
            datetime(2026, 3, 1, 10, 30),
            datetime(2026, 3, 1, 15, 0),
            True,
            transit_start,
            transit_stop,
            buffer_enabled=True,
            buffer_minutes=30,
        )
        assert regions == [
            (datetime(2026, 3, 1, 10, 30), datetime(2026, 3, 1, 11, 30)),
            (datetime(2026, 3, 1, 11, 30), datetime(2026, 3, 1, 13, 30)),
            (datetime(2026, 3, 1, 13, 30), datetime(2026, 3, 1, 15, 0)),
        ]

    def test_priority_buffer_disabled_preserves_original_segment(self):
        start = datetime(2026, 3, 1, 10, 30)
        stop = datetime(2026, 3, 1, 15, 0)
        regions = _split_at_priority_buffer_boundaries(
            start,
            stop,
            True,
            [datetime(2026, 3, 1, 12, 0)],
            [datetime(2026, 3, 1, 13, 0)],
            buffer_enabled=False,
            buffer_minutes=30,
        )
        assert regions == [(start, stop)]

    def test_target_priority_does_not_promote_touching_buffer_boundary(self):
        priority = _target_priority(
            True,
            [datetime(2026, 3, 1, 12, 0)],
            [datetime(2026, 3, 1, 13, 0)],
            datetime(2026, 3, 1, 10, 30),
            datetime(2026, 3, 1, 11, 30),
            buffer_enabled=True,
            buffer_minutes=30,
        )
        assert priority == "1"

    def test_adjacent_mode_selects_closest_sequence_on_each_side(self):
        builder = object.__new__(_ScienceCalendarBuilder)
        builder.config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 3, 2),
            obs_sequence_duration_min=90,
            priority_buffer=True,
            priority_buffer_mode="adjacent_sequences",
        )
        builder.sequence_duration = timedelta(minutes=90)
        segments = [
            (datetime(2026, 3, 1, 10, 0), datetime(2026, 3, 1, 11, 0), True),
            (datetime(2026, 3, 1, 11, 0), datetime(2026, 3, 1, 12, 0), False),
            (datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 13, 0), True),
            (datetime(2026, 3, 1, 13, 0), datetime(2026, 3, 1, 14, 0), False),
            (datetime(2026, 3, 1, 14, 0), datetime(2026, 3, 1, 15, 0), True),
        ]

        selected = builder._adjacent_priority_sequences(
            segments,
            True,
            [datetime(2026, 3, 1, 12, 15)],
            [datetime(2026, 3, 1, 12, 45)],
        )

        assert selected == {
            (datetime(2026, 3, 1, 10, 0), datetime(2026, 3, 1, 11, 0)),
            (datetime(2026, 3, 1, 14, 0), datetime(2026, 3, 1, 15, 0)),
        }

    def test_adjacent_mode_ignores_transits_outside_visit(self):
        builder = object.__new__(_ScienceCalendarBuilder)
        builder.config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 3, 2),
            priority_buffer=True,
            priority_buffer_mode="adjacent_sequences",
        )
        builder.sequence_duration = timedelta(minutes=90)
        segments = [
            (datetime(2026, 3, 1, 10, 0), datetime(2026, 3, 1, 11, 0), True),
        ]

        selected = builder._adjacent_priority_sequences(
            segments,
            True,
            [datetime(2026, 3, 2, 12, 0)],
            [datetime(2026, 3, 2, 13, 0)],
        )

        assert selected == set()

    def test_lookup_planet_row_found(self):
        cat = pd.DataFrame({"Planet Name": ["TOI-700 b"], "RA": [120.0]})
        result = _lookup_planet_row(cat, "TOI-700 b")
        assert result is not None
        assert result["RA"].iloc[0] == 120.0

    def test_lookup_planet_row_not_found(self):
        cat = pd.DataFrame({"Planet Name": ["TOI-700 b"], "RA": [120.0]})
        assert _lookup_planet_row(cat, "NOPE") is None

    def test_lookup_planet_row_empty_catalog(self):
        cat = pd.DataFrame({"Planet Name": [], "RA": []})
        assert _lookup_planet_row(cat, "anything") is None

    def test_lookup_auxiliary_row_found(self):
        cat = pd.DataFrame({"Star Name": ["Vega"], "DEC": [38.7]})
        result = _lookup_auxiliary_row(cat, "Vega")
        assert result is not None

    def test_lookup_auxiliary_row_not_found(self):
        cat = pd.DataFrame({"Star Name": ["Vega"]})
        assert _lookup_auxiliary_row(cat, "Sirius") is None


# ---------------------------------------------------------------------------
# Visibility segment extraction
# ---------------------------------------------------------------------------


class TestExtractVisibilitySegment:
    """Tests for _extract_visibility_segment."""

    def test_basic_split(self):
        times = pd.to_datetime(
            ["2026-03-01 00:00", "2026-03-01 00:01", "2026-03-01 00:02",
             "2026-03-01 00:03", "2026-03-01 00:04"]
        )
        vis = pd.DataFrame({
            "Time_UTC": times,
            "Time(MJD_UTC)": [61370.0 + i / 1440 for i in range(5)],
            "Visible": [1, 1, 0, 1, 1],
        })
        start = datetime(2026, 3, 1, 0, 0)
        stop = datetime(2026, 3, 1, 0, 4)
        visit_times, flags = _extract_visibility_segment(vis, start, stop, min_sequence_minutes=1)
        assert len(visit_times) >= 1
        assert len(visit_times) == len(flags)

    def test_empty_when_no_visible(self):
        times = pd.to_datetime(["2026-03-01 00:00", "2026-03-01 00:01"])
        vis = pd.DataFrame({
            "Time_UTC": times,
            "Time(MJD_UTC)": [61370.0, 61370.0 + 1 / 1440],
            "Visible": [0, 0],
        })
        visit_times, flags = _extract_visibility_segment(
            vis, datetime(2026, 3, 1), datetime(2026, 3, 1, 0, 1),
            min_sequence_minutes=1,
        )
        # Times are returned but all flags should be 0 (not visible)
        assert all(f == 0 for f in flags)


class TestVisibilityChangeIndices:
    """Tests for _visibility_change_indices."""

    def test_all_visible(self):
        indices = _visibility_change_indices(np.array([1, 1, 1, 1]))
        assert isinstance(indices, list)

    def test_alternating(self):
        indices = _visibility_change_indices(np.array([1, 0, 1, 0]))
        # Should detect each transition
        assert len(indices) >= 2


class TestVisibilityStats:
    """Tests for overlap-weighted provenance visibility statistics."""

    def test_second_aligned_fully_visible_sequence(self):
        index = pd.date_range("2026-08-12 05:29:00", periods=23, freq="min")
        prepared = pd.Series(True, index=index)

        fraction, visible, non_visible = (
            _ScienceCalendarBuilder._visibility_stats_in_series(
                prepared,
                datetime(2026, 8, 12, 5, 29, 43),
                datetime(2026, 8, 12, 5, 52, 0),
            )
        )

        assert fraction == pytest.approx(1.0)
        assert visible == pytest.approx(22 + 17 / 60)
        assert non_visible == pytest.approx(0.0)

    def test_second_aligned_sequence_weights_partial_minutes(self):
        index = pd.date_range("2026-08-12 05:29:00", periods=2, freq="min")
        prepared = pd.Series([False, True], index=index)

        fraction, visible, non_visible = (
            _ScienceCalendarBuilder._visibility_stats_in_series(
                prepared,
                datetime(2026, 8, 12, 5, 29, 30),
                datetime(2026, 8, 12, 5, 30, 30),
            )
        )

        assert fraction == pytest.approx(0.5)
        assert visible == pytest.approx(0.5)
        assert non_visible == pytest.approx(0.5)


# ---------------------------------------------------------------------------
# Full calendar generation (smoke test)
# ---------------------------------------------------------------------------


class TestGenerateScienceCalendar:
    """Smoke test: full generate_science_calendar with free-time-only schedule."""

    def test_free_time_only_schedule(self, tmp_path):
        """Schedule with only 'Free Time' produces a valid but visit-less XML."""
        csv = tmp_path / "sched.csv"
        pd.DataFrame(
            [
                {
                    "Target": "Free Time",
                    "Observation Start": "2026-03-01T00:00:00Z",
                    "Observation Stop": "2026-03-01T06:00:00Z",
                    "Primary Target": False,
                }
            ]
        ).to_csv(csv, index=False)

        data_dir = tmp_path / "data"
        _seed_catalogs(data_dir)

        inputs = ScienceCalendarInputs(schedule_csv=csv, data_dir=data_dir)
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
        )

        out = generate_science_calendar(inputs, config)
        assert out.exists()
        tree = ET.parse(out)
        root = tree.getroot()
        # The XML uses a namespace; use wildcard to find elements
        ns = {"ns": "/pandora/calendar/"}
        meta = root.find("ns:Meta", ns)
        assert meta is not None


# ---------------------------------------------------------------------------
# _merged_segments
# ---------------------------------------------------------------------------


class TestMergedSegments:
    """Tests for _ScienceCalendarBuilder._merged_segments."""

    def _make_builder(self, tmp_path, min_seq_min=10):
        _seed_catalogs(tmp_path)
        csv = tmp_path / "sched.csv"
        pd.DataFrame(
            [
                {
                    "Target": "Star b",
                    "Observation Start": "2026-03-01T00:00:00Z",
                    "Observation Stop": "2026-03-02T00:00:00Z",
                }
            ]
        ).to_csv(csv, index=False)
        inputs = ScienceCalendarInputs(
            schedule_csv=csv, data_dir=tmp_path
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
            min_sequence_minutes=min_seq_min,
        )
        return _ScienceCalendarBuilder(inputs, config)

    def test_short_segment_absorbed(self, tmp_path):
        """A 2-minute segment between two long segments is absorbed."""
        builder = self._make_builder(tmp_path, min_seq_min=10)
        t0 = datetime(2026, 3, 1, 0, 0)
        # 24 vis, 2 novis, 16 vis, 57 novis → 100 minute-resolution
        flags = [1] * 24 + [0] * 2 + [1] * 16 + [0] * 58
        times = [t0 + timedelta(minutes=i) for i in range(100)]
        changes = _visibility_change_indices(flags)

        raw = list(
            _ScienceCalendarBuilder._iterate_segments(
                changes, times, flags, t0, times[-1],
            )
        )
        merged = list(
            builder._merged_segments(
                changes, times, flags, t0, times[-1],
            )
        )
        # Raw has a short 2-min non-visible segment
        durations_raw = [
            (e - s).total_seconds() / 60 for s, e, v in raw
        ]
        assert any(d < 10 for d in durations_raw)
        # Merged has no segments shorter than threshold
        durations_merged = [
            (e - s).total_seconds() / 60 for s, e, v in merged
        ]
        assert all(d >= 10 for d in durations_merged)

    def test_no_merge_when_disabled(self, tmp_path):
        """min_sequence_minutes=0 disables merging."""
        builder = self._make_builder(tmp_path, min_seq_min=0)
        t0 = datetime(2026, 3, 1, 0, 0)
        flags = [1] * 5 + [0] * 2 + [1] * 5
        times = [t0 + timedelta(minutes=i) for i in range(12)]
        changes = _visibility_change_indices(flags)
        merged = list(
            builder._merged_segments(
                changes, times, flags, t0, times[-1],
            )
        )
        raw = list(
            _ScienceCalendarBuilder._iterate_segments(
                changes, times, flags, t0, times[-1],
            )
        )
        assert len(merged) == len(raw)

    def test_trailing_short_segment_absorbed(self, tmp_path):
        """A short segment at the end of the visit is absorbed backward."""
        builder = self._make_builder(tmp_path, min_seq_min=10)
        t0 = datetime(2026, 3, 1, 0, 0)
        # 30 vis, 50 novis, 3 vis, 1 novis  →  change indices at 29,79,82
        # The trailing 3-min visible segment is below the 10-min threshold.
        flags = [1] * 30 + [0] * 50 + [1] * 3 + [0] * 1
        times = [t0 + timedelta(minutes=i) for i in range(len(flags))]
        changes = _visibility_change_indices(flags)

        raw = list(
            _ScienceCalendarBuilder._iterate_segments(
                changes, times, flags, t0, times[-1],
            )
        )
        # Raw should have a trailing 3-min visible segment
        trailing_vis = [(s, e, v) for s, e, v in raw if v]
        assert any(
            (e - s).total_seconds() / 60 < 10 for s, e, v in trailing_vis
        ), f"Expected a short visible segment; got {trailing_vis}"

        merged = list(
            builder._merged_segments(
                changes, times, flags, t0, times[-1],
            )
        )
        # After merging, no segment should be shorter than threshold
        durations = [(e - s).total_seconds() / 60 for s, e, v in merged]
        assert all(d >= 10 for d in durations)


# ---------------------------------------------------------------------------
# _occ_visibility_score
# ---------------------------------------------------------------------------


def _write_visibility(
    data_dir: Path, star_name: str, df: pd.DataFrame,
) -> None:
    """Write a visibility parquet for *star_name* under ``aux_targets/``."""
    out_dir = data_dir / "aux_targets" / star_name
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / f"Visibility for {star_name}.parquet", index=False)


class TestOccVisibilityScore:
    """Tests for _occ_visibility_score and allow_occ_startracker_violation."""

    def _make_builder(self, tmp_path, *, allow_st=False):
        _seed_catalogs(tmp_path)
        csv = tmp_path / "sched.csv"
        pd.DataFrame(
            [
                {
                    "Target": "Star b",
                    "Observation Start": "2026-03-01T00:00:00Z",
                    "Observation Stop": "2026-03-02T00:00:00Z",
                }
            ]
        ).to_csv(csv, index=False)
        inputs = ScienceCalendarInputs(schedule_csv=csv, data_dir=tmp_path)
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 3, 1),
            window_end=datetime(2026, 4, 1),
            sun_avoidance_deg=91.0,
            moon_avoidance_deg=25.0,
            earth_avoidance_deg=110.0,
            allow_occ_startracker_violation=allow_st,
        )
        return _ScienceCalendarBuilder(inputs, config)

    def _make_vis_df(self, n_minutes, *, visible, sun_sep, moon_sep, earth_sep):
        """Create a visibility DataFrame for *n_minutes* starting 2026-03-01."""
        t0 = datetime(2026, 3, 1, 12, 0)
        from astropy.time import Time as AstropyTime
        times_utc = [t0 + timedelta(minutes=i) for i in range(n_minutes)]
        mjds = [AstropyTime(t).mjd for t in times_utc]
        return pd.DataFrame({
            "Time(MJD_UTC)": mjds,
            "Time_UTC": pd.to_datetime(times_utc),
            "Visible": [float(visible)] * n_minutes,
            "Sun_Sep": [float(sun_sep)] * n_minutes,
            "Moon_Sep": [float(moon_sep)] * n_minutes,
            "Earth_Sep": [float(earth_sep)] * n_minutes,
        })

    def test_fully_visible_returns_zero_frac(self, tmp_path):
        builder = self._make_builder(tmp_path)
        df = self._make_vis_df(10, visible=1, sun_sep=120, moon_sep=50, earth_sep=130)
        _write_visibility(tmp_path, "OccStar", df)
        ok, frac = builder._occ_visibility_score(
            "OccStar", datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 12, 10),
        )
        assert ok is True
        assert frac == 0.0

    def test_boresight_fail_rejected_even_with_st_relaxed(self, tmp_path):
        """Sun constraint violated → reject even with allow_occ_startracker_violation."""
        builder = self._make_builder(tmp_path, allow_st=True)
        # Sun_Sep < 91 → boresight fail
        df = self._make_vis_df(10, visible=0, sun_sep=50, moon_sep=50, earth_sep=130)
        _write_visibility(tmp_path, "OccStar", df)
        ok, frac = builder._occ_visibility_score(
            "OccStar", datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 12, 10),
        )
        assert ok is False
        assert frac == 1.0

    def test_st_only_fail_rejected_when_flag_off(self, tmp_path):
        """ST-only failure rejected when allow_occ_startracker_violation=False."""
        builder = self._make_builder(tmp_path, allow_st=False)
        # Boresight passes, but Visible=0 → star tracker killed it
        df = self._make_vis_df(10, visible=0, sun_sep=120, moon_sep=50, earth_sep=130)
        _write_visibility(tmp_path, "OccStar", df)
        ok, frac = builder._occ_visibility_score(
            "OccStar", datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 12, 10),
        )
        assert ok is False

    def test_st_only_fail_accepted_when_flag_on(self, tmp_path):
        """ST-only failure accepted when allow_occ_startracker_violation=True."""
        builder = self._make_builder(tmp_path, allow_st=True)
        # Boresight passes, but Visible=0 → star tracker only
        df = self._make_vis_df(10, visible=0, sun_sep=120, moon_sep=50, earth_sep=130)
        _write_visibility(tmp_path, "OccStar", df)
        ok, frac = builder._occ_visibility_score(
            "OccStar", datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 12, 10),
        )
        assert ok is True
        assert frac == 1.0  # 100% of segment is Visible=0

    def test_partial_visibility_st_violation_scored(self, tmp_path):
        """Some minutes not visible → falls through to extended check,
        computes real ST violation fraction (not the fast path)."""
        builder = self._make_builder(tmp_path, allow_st=True)
        t0 = datetime(2026, 3, 1, 12, 0)
        from astropy.time import Time as AstropyTime
        n = 10
        times_utc = [t0 + timedelta(minutes=i) for i in range(n)]
        mjds = [AstropyTime(t).mjd for t in times_utc]
        # First 6 visible, last 4 not (but boresight OK everywhere)
        df = pd.DataFrame({
            "Time(MJD_UTC)": mjds,
            "Time_UTC": pd.to_datetime(times_utc),
            "Visible": [1.0] * 6 + [0.0] * 4,
            "Sun_Sep": [120.0] * n,
            "Moon_Sep": [50.0] * n,
            "Earth_Sep": [130.0] * n,
        })
        _write_visibility(tmp_path, "OccStar", df)
        ok, frac = builder._occ_visibility_score(
            "OccStar", datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 12, 10),
        )
        assert ok is True
        assert frac == pytest.approx(0.4)  # 4/10 minutes not visible (ST-only)

    def test_missing_separation_columns_rejected(self, tmp_path):
        """If parquet lacks separation columns, treat as normal reject."""
        builder = self._make_builder(tmp_path, allow_st=True)
        # Write parquet with only the basic columns (no Sun_Sep etc.)
        t0 = datetime(2026, 3, 1, 12, 0)
        from astropy.time import Time as AstropyTime
        n = 10
        times_utc = [t0 + timedelta(minutes=i) for i in range(n)]
        mjds = [AstropyTime(t).mjd for t in times_utc]
        df = pd.DataFrame({
            "Time(MJD_UTC)": mjds,
            "Time_UTC": pd.to_datetime(times_utc),
            "Visible": [0.0] * n,
        })
        out_dir = tmp_path / "aux_targets" / "OccStar"
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_parquet(out_dir / "Visibility for OccStar.parquet", index=False)
        ok, frac = builder._occ_visibility_score(
            "OccStar", datetime(2026, 3, 1, 12, 0), datetime(2026, 3, 1, 12, 10),
        )
        assert ok is False
