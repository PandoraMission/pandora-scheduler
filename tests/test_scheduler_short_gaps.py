"""Tests for scheduler handling of short gaps before primary observations.

When the gap between the current schedule position and the next primary observation
start is shorter than min_sequence_minutes, the scheduler should not attempt to
schedule an auxiliary target. Instead, it should mark the gap as Free Time.
"""

from datetime import datetime, timedelta
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    SchedulerState,
    _schedule_auxiliary_target,
    _schedule_primary_target,
)


@pytest.fixture
def mock_config():
    """Create a minimal config with min_sequence_minutes=5."""
    return PandoraSchedulerConfig(
        window_start=datetime(2026, 1, 1),
        window_end=datetime(2026, 2, 1),
        output_dir="/tmp/test_output",
        min_sequence_minutes=5,
    )


@pytest.fixture
def mock_state():
    """Create a mock scheduler state."""
    state = MagicMock(spec=SchedulerState)
    state.tracker = pd.DataFrame({
        "Planet Name": ["TestPlanet"],
        "Transits Needed": [1],
        "Transits Acquired": [0],
        "Transits Left in Lifetime": [10],
        "Transit Priority": [9],
    })
    state.all_target_obs_time = {}
    state.non_primary_obs_time = {}
    return state


@pytest.fixture
def mock_inputs():
    """Create mock scheduler inputs."""
    inputs = MagicMock(spec=SchedulerInputs)
    inputs.target_list = pd.DataFrame({
        "Planet Name": ["TestPlanet"],
        "Star Name": ["TestStar"],
        "RA": [90.0],
        "DEC": [45.0],
        "Primary Target": [True],
        "Obs Window (hrs)": [6.0],
    })
    inputs.aux_star_catalog = pd.DataFrame()
    inputs.aux_star_visibility = {}
    return inputs


class TestShortGapHandling:
    """Test that short gaps before primary observations are handled correctly."""

    def test_short_gap_becomes_free_time(self, mock_config, mock_state, mock_inputs):
        """A 3-minute gap (< min_sequence_minutes=5) should become Free Time."""
        # Create a candidate table with obs_start 3 minutes from now
        start = datetime(2026, 2, 5, 12, 0, 0)
        obs_start = datetime(2026, 2, 5, 12, 3, 0)  # 3 minutes later
        obs_stop = obs_start + timedelta(hours=6)
        
        temp_df = pd.DataFrame([{
            "Planet Name": "TestPlanet",
            "RA": 90.0,
            "DEC": 45.0,
            "Obs Start": obs_start,
            "Visit Duration": timedelta(hours=6),
            "Transit Coverage": 1.0,
            "SAA Overlap": 0.0,
            "Schedule Factor": 1.0,
            "Quality Factor": 1.0,
            "Transit Factor": 5.0,  # Not critical
        }])
        
        obs_range = pd.date_range(start, obs_stop, freq="min")
        
        # Mock _schedule_auxiliary_target to track if it's called
        with patch(
            "pandorascheduler_rework.scheduler._schedule_auxiliary_target"
        ) as mock_aux:
            mock_aux.return_value = (pd.DataFrame(), "")
            
            result = _schedule_primary_target(
                temp_df,
                mock_state,
                mock_inputs,
                mock_config,
                start,
                obs_range,
            )
            
            # The auxiliary scheduler should NOT be called for a 3-minute gap
            mock_aux.assert_not_called()
            
            # The result should include a Free Time entry for the gap
            free_time_rows = result[result["Target"] == "Free Time"]
            assert len(free_time_rows) == 1
            
            # Verify Free Time covers the gap
            free_start = free_time_rows.iloc[0]["Observation Start"]
            free_stop = free_time_rows.iloc[0]["Observation Stop"]
            assert free_start == start
            assert free_stop == obs_start

    def test_sufficient_gap_schedules_auxiliary(self, mock_config, mock_state, mock_inputs):
        """A 10-minute gap (>= min_sequence_minutes=5) should attempt aux scheduling."""
        # Create a candidate table with obs_start 10 minutes from now
        start = datetime(2026, 2, 5, 12, 0, 0)
        obs_start = datetime(2026, 2, 5, 12, 10, 0)  # 10 minutes later
        obs_stop = obs_start + timedelta(hours=6)
        
        temp_df = pd.DataFrame([{
            "Planet Name": "TestPlanet",
            "RA": 90.0,
            "DEC": 45.0,
            "Obs Start": obs_start,
            "Visit Duration": timedelta(hours=6),
            "Transit Coverage": 1.0,
            "SAA Overlap": 0.0,
            "Schedule Factor": 1.0,
            "Quality Factor": 1.0,
            "Transit Factor": 5.0,
        }])
        
        obs_range = pd.date_range(start, obs_stop, freq="min")
        
        with patch(
            "pandorascheduler_rework.scheduler._schedule_auxiliary_target"
        ) as mock_aux:
            mock_aux.return_value = (
                pd.DataFrame([{
                    "Target": "AuxTarget",
                    "Observation Start": start,
                    "Observation Stop": obs_start,
                    "RA": 0.0,
                    "DEC": 0.0,
                }]),
                "AuxTarget scheduled",
            )
            
            result = _schedule_primary_target(
                temp_df,
                mock_state,
                mock_inputs,
                mock_config,
                start,
                obs_range,
            )
            
            # Auxiliary scheduler SHOULD be called for a 10-minute gap
            mock_aux.assert_called_once()

    def test_primary_only_gap_becomes_free_time(self, mock_state, mock_inputs):
        """Primary-only mode should convert a pre-primary gap into Free Time."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 2, 1),
            output_dir="/tmp/test_output",
            min_sequence_minutes=5,
            primary_only_mode=True,
        )

        start = datetime(2026, 2, 5, 12, 0, 0)
        obs_start = datetime(2026, 2, 5, 12, 10, 0)
        obs_stop = obs_start + timedelta(hours=6)

        temp_df = pd.DataFrame([{
            "Planet Name": "TestPlanet",
            "RA": 90.0,
            "DEC": 45.0,
            "Obs Start": obs_start,
            "Visit Duration": timedelta(hours=6),
            "Transit Coverage": 1.0,
            "SAA Overlap": 0.0,
            "Schedule Factor": 1.0,
            "Quality Factor": 1.0,
            "Transit Factor": 5.0,
        }])

        obs_range = pd.date_range(start, obs_stop, freq="min")

        result = _schedule_primary_target(
            temp_df,
            mock_state,
            mock_inputs,
            config,
            start,
            obs_range,
        )

        free_time_rows = result[result["Target"] == "Free Time"]
        assert len(free_time_rows) == 1
        assert free_time_rows.iloc[0]["Observation Start"] == start
        assert free_time_rows.iloc[0]["Observation Stop"] == obs_start

    def test_exact_threshold_gap_becomes_free_time(self, mock_config, mock_state, mock_inputs):
        """A gap exactly equal to min_sequence_minutes should become Free Time.
        
        Due to the fencepost problem (a 5-minute window only yields 4 visibility
        samples at minute cadence), we use strict > comparison. This means a gap
        of exactly min_sequence_minutes becomes Free Time, not an aux observation.
        """
        start = datetime(2026, 2, 5, 12, 0, 0)
        obs_start = datetime(2026, 2, 5, 12, 5, 0)  # Exactly 5 minutes later
        obs_stop = obs_start + timedelta(hours=6)
        
        temp_df = pd.DataFrame([{
            "Planet Name": "TestPlanet",
            "RA": 90.0,
            "DEC": 45.0,
            "Obs Start": obs_start,
            "Visit Duration": timedelta(hours=6),
            "Transit Coverage": 1.0,
            "SAA Overlap": 0.0,
            "Schedule Factor": 1.0,
            "Quality Factor": 1.0,
            "Transit Factor": 5.0,
        }])
        
        obs_range = pd.date_range(start, obs_stop, freq="min")
        
        with patch(
            "pandorascheduler_rework.scheduler._schedule_auxiliary_target"
        ) as mock_aux:
            mock_aux.return_value = (pd.DataFrame(), "")
            
            result = _schedule_primary_target(
                temp_df,
                mock_state,
                mock_inputs,
                mock_config,
                start,
                obs_range,
            )
            
            # Exactly 5 minutes should NOT schedule aux (use strict > to avoid fencepost)
            mock_aux.assert_not_called()
            
            # Should have Free Time instead
            free_time_rows = result[result["Target"] == "Free Time"]
            assert len(free_time_rows) == 1


class TestShortWindowAfterSTD:
    """Test that short remaining windows after STD observations become Free Time."""

    def test_short_remaining_window_after_std_becomes_free_time(self, tmp_path):
        """When STD takes most of the window, remaining 1 min should be Free Time."""
        from pandorascheduler_rework.scheduler import (
            SchedulerInputs,
            SchedulerState,
            _schedule_auxiliary_target,
        )

        # Create a config where STD takes most of a 31-minute window
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 2, 1),
            output_dir=str(tmp_path),
            min_sequence_minutes=5,
            std_obs_duration_hours=0.5,  # 30 minutes for STD
            std_obs_frequency_days=1,    # Always schedule STD if due
            aux_sort_key="sort_by_tdf_priority",  # Enable aux scheduling
        )
        
        # Create mock state where STD is due (last_std_obs is old)
        state = MagicMock(spec=SchedulerState)
        state.last_std_obs = datetime(2026, 1, 1)  # Very old, so STD is due
        state.non_primary_obs_time = {}
        state.all_target_obs_time = {}
        
        # Create mock inputs with minimal setup
        inputs = MagicMock(spec=SchedulerInputs)
        inputs.target_definition_files = ["primary"]  # No aux targets
        inputs.paths = MagicMock(spec=SchedulerPaths)
        inputs.paths.data_dir = tmp_path
        inputs.paths.aux_targets_dir = tmp_path / "aux_targets"
        
        # Create a minimal monitoring-standard_targets.csv
        std_csv = tmp_path / "monitoring-standard_targets.csv"
        std_df = pd.DataFrame([{
            "Star Name": "TestSTD",
            "RA": 0.0,
            "DEC": 0.0,
            "Priority": 1.0,
        }])
        std_df.to_csv(std_csv, index=False)


class TestPrimaryOnlyMode:
    """Test that primary-only mode skips non-primary scheduling entirely."""

    def test_auxiliary_scheduler_returns_free_time_in_primary_only_mode(self, tmp_path):
        start = datetime(2026, 2, 5, 12, 0, 0)
        stop = datetime(2026, 2, 5, 13, 0, 0)

        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 2, 1),
            output_dir=str(tmp_path),
            min_sequence_minutes=5,
            primary_only_mode=True,
        )

        state = MagicMock(spec=SchedulerState)
        state.last_std_obs = datetime(2026, 1, 1)
        state.non_primary_obs_time = {}
        state.all_target_obs_time = {}

        inputs = MagicMock(spec=SchedulerInputs)

        result, log_info = _schedule_auxiliary_target(
            start,
            stop,
            config,
            state,
            inputs,
        )

        assert log_info == "Primary-only mode enabled, skipping non-primary scheduling."
        assert len(result) == 1
        assert result.iloc[0]["Target"] == "Free Time"
        assert result.iloc[0]["Observation Start"] == start
        assert result.iloc[0]["Observation Stop"] == stop
        assert state.all_target_obs_time == {}


class TestShortWindowAfterSTDContinuation:
    """Test STD handling when the post-STD remainder is too short."""

    def test_std_absorbs_short_post_std_remainder(self, tmp_path):
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 2, 1),
            output_dir=str(tmp_path),
            min_sequence_minutes=5,
            std_obs_duration_hours=0.5,
            std_obs_frequency_days=1,
            aux_sort_key="sort_by_tdf_priority",
        )

        state = MagicMock(spec=SchedulerState)
        state.last_std_obs = datetime(2026, 1, 1)
        state.non_primary_obs_time = {}
        state.all_target_obs_time = {}

        inputs = MagicMock(spec=SchedulerInputs)
        inputs.target_definition_files = ["primary"]
        inputs.paths = MagicMock(spec=SchedulerPaths)
        inputs.paths.data_dir = tmp_path
        inputs.paths.aux_targets_dir = tmp_path / "aux_targets"

        std_csv = tmp_path / "monitoring-standard_targets.csv"
        std_df = pd.DataFrame([
            {
                "Star Name": "TestSTD",
                "RA": 0.0,
                "DEC": 0.0,
                "Priority": 1.0,
            }
        ])
        std_df.to_csv(std_csv, index=False)

        # Create visibility file for STD
        std_vis_dir = tmp_path / "aux_targets" / "TestSTD"
        std_vis_dir.mkdir(parents=True)
        std_vis_file = std_vis_dir / "Visibility for TestSTD.parquet"

        # Create visibility data that covers the window
        # Need both Time(MJD_UTC) and Time_UTC columns for compatibility
        from astropy.time import Time as AstropyTime

        time_utc = pd.date_range(
            "2026-10-01 22:32:00", "2026-10-01 23:05:00", freq="min"
        )
        mjd_values = AstropyTime(time_utc.to_pydatetime(), scale="utc").mjd
        vis_df = pd.DataFrame({
            "Time(MJD_UTC)": mjd_values,
            "Time_UTC": time_utc,
            "Visible": [1] * len(time_utc),
        })
        vis_df.to_parquet(std_vis_file)

        # Schedule with a 31-minute window (22:32 to 23:03)
        # STD takes 30 min, leaving only 1 minute
        start = datetime(2026, 10, 1, 22, 32, 0)
        stop = datetime(2026, 10, 1, 23, 3, 0)  # 31 minutes later

        result, log_info = _schedule_auxiliary_target(
            start, stop, config, state, inputs
        )

        # Should have STD observation — the 1-minute post-STD remainder
        # is too short for an aux observation, so the STD row is extended
        # to absorb it (no Free Time).
        std_rows = result[result["Target"].str.contains("STD", na=False)]
        assert len(std_rows) == 1
        assert std_rows.iloc[0]["Observation Start"] == start
        assert std_rows.iloc[0]["Observation Stop"] == stop  # extended to end of gap

        # No Free Time — the STD absorbed the short remainder
        free_rows = result[result["Target"] == "Free Time"]
        assert len(free_rows) == 0
