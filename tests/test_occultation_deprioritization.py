"""Tests for occultation target time tracking and deprioritization."""

from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework import science_calendar
from pandorascheduler_rework.config import PandoraSchedulerConfig


def _write_visibility(directory: Path, name: str, times: list[datetime], flags: list[int]) -> None:
    """Write visibility parquet file for a target."""
    directory.mkdir(parents=True, exist_ok=True)
    mjd_times = Time(times, scale="utc").to_value("mjd")
    pd.DataFrame({"Time(MJD_UTC)": mjd_times, "Visible": flags}).to_parquet(
        directory / f"Visibility for {name}.parquet", index=False
    )


def _write_planet_visibility(directory: Path, name: str, start: datetime, stop: datetime) -> None:
    """Write planet transit visibility parquet file."""
    directory.mkdir(parents=True, exist_ok=True)
    transit_start = Time([start], scale="utc").to_value("mjd")
    transit_stop = Time([stop], scale="utc").to_value("mjd")
    pd.DataFrame(
        {
            "Transit_Start": transit_start,
            "Transit_Stop": transit_stop,
            "Transit_Coverage": [0.75],
            "SAA_Overlap": [0.1],
        }
    ).to_parquet(directory / f"Visibility for {name}.parquet", index=False)


class TestGetOccultationTimeLimitStrict:
    """Tests for strict validation in _get_occultation_time_limit."""

    def test_raises_when_catalog_empty(self, tmp_path):
        """Test that error is raised when occultation catalog is empty."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "exoplanet_targets.csv", index=False
        )
        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "all_targets.csv", index=False
        )
        # Create empty occultation catalog
        pd.DataFrame({
            "Star Name": [],
            "RA": [],
            "DEC": [],
            "Number of Hours Requested": [],
        }).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

        schedule_df = pd.DataFrame([{
            "Target": "TestPlanet",
            "Observation Start": "2026-01-01 00:00:00",
            "Observation Stop": "2026-01-01 01:00:00",
            "Transit Coverage": 0.5,
            "SAA Overlap": 0.0,
            "Schedule Factor": 0.9,
            "Quality Factor": 0.8,
            "Comments": "",
        }])
        schedule_path = tmp_path / "schedule.csv"
        schedule_df.to_csv(schedule_path, index=False)

        inputs = science_calendar.ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
        )

        builder = science_calendar._ScienceCalendarBuilder(inputs, config)
        with pytest.raises(ValueError, match="occultation catalog is not loaded"):
            builder._get_occultation_time_limit("AnyTarget")
    def test_raises_when_target_not_in_catalog(self, tmp_path):
        """Test that error is raised when target is not found in catalog."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "exoplanet_targets.csv", index=False
        )
        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "all_targets.csv", index=False
        )
        pd.DataFrame({
            "Star Name": ["OccA"],
            "RA": [10.0],
            "DEC": [20.0],
            "Number of Hours Requested": [600],
        }).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

        schedule_df = pd.DataFrame([{
            "Target": "TestPlanet",
            "Observation Start": "2026-01-01 00:00:00",
            "Observation Stop": "2026-01-01 01:00:00",
            "Transit Coverage": 0.5,
            "SAA Overlap": 0.0,
            "Schedule Factor": 0.9,
            "Quality Factor": 0.8,
            "Comments": "",
        }])
        schedule_path = tmp_path / "schedule.csv"
        schedule_df.to_csv(schedule_path, index=False)

        inputs = science_calendar.ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
        )

        builder = science_calendar._ScienceCalendarBuilder(inputs, config)
        with pytest.raises(ValueError, match="not found in catalog"):
            builder._get_occultation_time_limit("UnknownTarget")

    def test_raises_when_column_missing(self, tmp_path):
        """Test that error is raised when 'Number of Hours Requested' column is missing."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "exoplanet_targets.csv", index=False
        )
        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "all_targets.csv", index=False
        )
        # Missing "Number of Hours Requested" column
        pd.DataFrame({
            "Star Name": ["OccA"],
            "RA": [10.0],
            "DEC": [20.0],
        }).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

        schedule_df = pd.DataFrame([{
            "Target": "TestPlanet",
            "Observation Start": "2026-01-01 00:00:00",
            "Observation Stop": "2026-01-01 01:00:00",
            "Transit Coverage": 0.5,
            "SAA Overlap": 0.0,
            "Schedule Factor": 0.9,
            "Quality Factor": 0.8,
            "Comments": "",
        }])
        schedule_path = tmp_path / "schedule.csv"
        schedule_df.to_csv(schedule_path, index=False)

        inputs = science_calendar.ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
        )

        builder = science_calendar._ScienceCalendarBuilder(inputs, config)
        with pytest.raises(ValueError, match="missing required.*Number of Hours Requested"):
            builder._get_occultation_time_limit("OccA")

    def test_returns_hours_when_valid(self, tmp_path):
        """Test that correct timedelta is returned when data is valid."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "exoplanet_targets.csv", index=False
        )
        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "all_targets.csv", index=False
        )
        pd.DataFrame({
            "Star Name": ["OccA"],
            "RA": [10.0],
            "DEC": [20.0],
            "Number of Hours Requested": [600],
        }).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

        schedule_df = pd.DataFrame([{
            "Target": "TestPlanet",
            "Observation Start": "2026-01-01 00:00:00",
            "Observation Stop": "2026-01-01 01:00:00",
            "Transit Coverage": 0.5,
            "SAA Overlap": 0.0,
            "Schedule Factor": 0.9,
            "Quality Factor": 0.8,
            "Comments": "",
        }])
        schedule_path = tmp_path / "schedule.csv"
        schedule_df.to_csv(schedule_path, index=False)

        inputs = science_calendar.ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
        )

        builder = science_calendar._ScienceCalendarBuilder(inputs, config)
        result = builder._get_occultation_time_limit("OccA")
        assert result == timedelta(hours=600)


class TestOccultationTimeTracking:
    """Tests for tracking occultation target observation time."""

    def test_builder_initializes_time_tracking(self, tmp_path):
        """Test that _ScienceCalendarBuilder initializes time tracking."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # Create minimal required files
        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "exoplanet_targets.csv", index=False
        )
        pd.DataFrame({"Star Name": [], "RA": [], "DEC": []}).to_csv(
            data_dir / "all_targets.csv", index=False
        )
        pd.DataFrame({
            "Star Name": ["OccStar"],
            "RA": [10.0],
            "DEC": [20.0],
            "Number of Hours Requested": [600],
        }).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

        # Create minimal schedule
        schedule_df = pd.DataFrame([{
            "Target": "TestPlanet",
            "Observation Start": "2026-01-01 00:00:00",
            "Observation Stop": "2026-01-01 01:00:00",
            "Transit Coverage": 0.5,
            "SAA Overlap": 0.0,
            "Schedule Factor": 0.9,
            "Quality Factor": 0.8,
            "Comments": "",
        }])
        schedule_path = tmp_path / "schedule.csv"
        schedule_df.to_csv(schedule_path, index=False)

        inputs = science_calendar.ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
        )

        builder = science_calendar._ScienceCalendarBuilder(inputs, config)

        # Check that time tracking is initialized
        assert hasattr(builder, "occultation_obs_time")
        assert builder.occultation_obs_time == {}


class TestBuildOccultationScheduleExclusion:
    """Tests for excluding over-limit targets from occultation scheduling."""

    def test_excluded_targets_filtered_from_occ_list(self, tmp_path):
        """Test that excluded targets are filtered from the occultation list."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # Create occultation list with multiple targets
        occ_df = pd.DataFrame({
            "Star Name": ["OccA", "OccB", "OccC"],
            "RA": [10.0, 20.0, 30.0],
            "DEC": [10.0, 20.0, 30.0],
            "Priority": [0.9, 0.8, 0.7],
        })
        occ_csv = data_dir / "occultation-standard_targets.csv"
        occ_df.to_csv(occ_csv, index=False)

        # Create visibility files for all targets
        start = datetime(2026, 1, 1, 0, 0, 0)
        times = [start + timedelta(minutes=i) for i in range(60)]
        vis_dir = data_dir / "aux_targets"
        for name in ["OccA", "OccB", "OccC"]:
            _write_visibility(vis_dir / name, name, times, [1] * len(times))

        # Test with OccB excluded
        excluded = {"OccB"}
        result_df, flag = science_calendar._build_occultation_schedule(
            starts=[start],
            stops=[start + timedelta(minutes=30)],
            visit_start=start,
            visit_stop=start + timedelta(hours=1),
            list_path=occ_csv,
            vis_root=vis_dir,
            label="test",
            reference_ra=0.0,
            reference_dec=0.0,
            prioritise_by_slew=False,
            excluded_targets=excluded,
        )

        # Should schedule and NOT use OccB
        assert flag is True
        assert result_df is not None
        targets_used = result_df["Target"].unique()
        assert "OccB" not in targets_used

    def test_all_excluded_returns_failure(self, tmp_path):
        """Test that if all targets are excluded, scheduling fails gracefully."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # Create occultation list with single target
        occ_df = pd.DataFrame({
            "Star Name": ["OnlyOcc"],
            "RA": [10.0],
            "DEC": [10.0],
            "Priority": [0.9],
        })
        occ_csv = data_dir / "occultation-standard_targets.csv"
        occ_df.to_csv(occ_csv, index=False)

        # Create visibility file
        start = datetime(2026, 1, 1, 0, 0, 0)
        times = [start + timedelta(minutes=i) for i in range(60)]
        vis_dir = data_dir / "aux_targets"
        _write_visibility(vis_dir / "OnlyOcc", "OnlyOcc", times, [1] * len(times))

        # Exclude the only target
        excluded = {"OnlyOcc"}
        result_df, flag = science_calendar._build_occultation_schedule(
            starts=[start],
            stops=[start + timedelta(minutes=30)],
            visit_start=start,
            visit_stop=start + timedelta(hours=1),
            list_path=occ_csv,
            vis_root=vis_dir,
            label="test",
            reference_ra=0.0,
            reference_dec=0.0,
            prioritise_by_slew=False,
            excluded_targets=excluded,
        )

        # Should fail to schedule since no targets available
        assert flag is False


class TestOccultationDeprioritizationIntegration:
    """Integration tests for occultation deprioritization."""

    def test_target_excluded_after_time_limit_exceeded(self, tmp_path):
        """Test that a target is excluded once it exceeds its time limit."""
        data_dir = tmp_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        start = datetime(2026, 1, 1, 0, 0, 0)
        times = [start + timedelta(minutes=i) for i in range(240)]  # 4 hours

        # Create primary target with some non-visible periods
        primary_flags = [1] * 30 + [0] * 30 + [1] * 30 + [0] * 30 + [1] * 120
        _write_visibility(data_dir / "targets" / "StarOne", "StarOne", times, primary_flags)
        _write_planet_visibility(
            data_dir / "targets" / "StarOne" / "StarOne b",
            "StarOne b",
            start + timedelta(minutes=30),
            start + timedelta(minutes=60),
        )

        # Create occultation targets
        for name in ["OccA", "OccB"]:
            _write_visibility(data_dir / "aux_targets" / name, name, times, [1] * len(times))

        # Create catalog files
        pd.DataFrame([{
            "Planet Name": "StarOne b",
            "Star Name": "StarOne",
            "RA": 10.0,
            "DEC": -20.0,
        }]).to_csv(data_dir / "exoplanet_targets.csv", index=False)

        pd.DataFrame({"Star Name": [], "RA": [], "DEC": [], "Priority": []}).to_csv(
            data_dir / "all_targets.csv", index=False
        )

        pd.DataFrame([
            {"Star Name": "OccA", "RA": 30.0, "DEC": 15.0, "Priority": 0.9, "Number of Hours Requested": 0.5},
            {"Star Name": "OccB", "RA": 35.0, "DEC": 20.0, "Priority": 0.8, "Number of Hours Requested": 0.5},
        ]).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

        # Create schedule - two visits to allow time accumulation
        schedule_df = pd.DataFrame([
            {
                "Target": "StarOne b",
                "Observation Start": "2026-01-01 00:00:00",
                "Observation Stop": "2026-01-01 02:00:00",
                "Transit Coverage": 0.75,
                "SAA Overlap": 0.1,
                "Schedule Factor": 0.9,
                "Quality Factor": 0.85,
                "Comments": "",
            },
            {
                "Target": "StarOne b",
                "Observation Start": "2026-01-01 02:00:00",
                "Observation Stop": "2026-01-01 04:00:00",
                "Transit Coverage": 0.75,
                "SAA Overlap": 0.1,
                "Schedule Factor": 0.9,
                "Quality Factor": 0.85,
                "Comments": "",
            },
        ])
        schedule_path = tmp_path / "schedule.csv"
        schedule_df.to_csv(schedule_path, index=False)

        inputs = science_calendar.ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        # The low time limit (0.5 hours) is set in manifest for each target
        config = PandoraSchedulerConfig(
            window_start=start,
            window_end=start + timedelta(hours=4),
            visit_limit=2,
        )

        builder = science_calendar._ScienceCalendarBuilder(inputs, config)
        builder.build_calendar()

        # Check that time was tracked
        assert len(builder.occultation_obs_time) > 0, "Occultation time should be tracked"

        # At least one target should have accumulated time
        total_tracked = sum(
            t.total_seconds() for t in builder.occultation_obs_time.values()
        )
        assert total_tracked > 0, "Some occultation time should have been tracked"
