"""Error handling tests for scheduler module.

These tests verify that the scheduler handles error conditions gracefully
and provides helpful error messages when things go wrong.
"""

from datetime import datetime

import pandas as pd
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    _load_too_table,
    run_scheduler,
)


class TestSchedulerErrorHandling:
    """Test scheduler error handling and edge cases."""

    def test_load_too_table_prefers_output_dir_over_data_dir(self, tmp_path):
        output_too = tmp_path / "ToO_list.csv"
        output_too.write_text(
            "Target,Obs Window Start,Obs Window Stop\n"
            "RootTarget,2026-06-02 00:00:00,2026-06-02 02:00:00\n",
            encoding="utf-8",
        )

        data_dir = tmp_path / "data_91_20_111"
        data_dir.mkdir()
        data_too = data_dir / "ToO_list.csv"
        data_too.write_text(
            "Target,Obs Window Start,Obs Window Stop\n"
            "StaleTarget,2026-06-01 00:00:00,2026-06-01 02:00:00\n",
            encoding="utf-8",
        )

        targets, starts, stops = _load_too_table(data_dir, tmp_path)

        assert targets == ["RootTarget"]
        assert starts == [datetime(2026, 6, 2, 0, 0, 0)]
        assert stops == [datetime(2026, 6, 2, 2, 0, 0)]
    
    def test_scheduler_handles_missing_visibility_file(self, tmp_path):
        """Verify graceful error when visibility file is missing."""
        # Create minimal directory structure
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        targets_dir = data_dir / "targets"
        targets_dir.mkdir()
        
        # Create target CSV
        target_df = pd.DataFrame([{
            "Planet Name": "TestPlanet",
            "Star Name": "TestStar",
            "RA": 0.0,
            "DEC": 0.0,
            "Primary Target": True,
            "Number of Transits to Capture": 1,
            "Transit Duration (hrs)": 2.0,
            "Period (days)": 1.0,
            "Transit Epoch (BJD_TDB-2400000.5)": 60000.0,
        }])
        primary_csv = tmp_path / "primary.csv"
        target_df.to_csv(primary_csv, index=False)
        
        # Create config
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            output_dir=tmp_path,
        )
        
        # Create paths
        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=targets_dir,
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )
        
        # Create inputs WITHOUT visibility files
        inputs = SchedulerInputs(
            pandora_start=config.window_start,
            pandora_stop=config.window_end,
            sched_start=config.window_start,
            sched_stop=config.window_end,
            target_list=target_df,
            paths=paths,
            target_definition_files=[],
            primary_target_csv=primary_csv,
            auxiliary_target_csv=tmp_path / "aux.csv",
            occultation_target_csv=tmp_path / "occ.csv",
            output_dir=tmp_path,
        )
        
        # Create empty auxiliary/occultation CSVs
        pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(inputs.auxiliary_target_csv, index=False)
        pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(inputs.occultation_target_csv, index=False)
        
        # Should raise FileNotFoundError when trying to find visibility
        with pytest.raises(FileNotFoundError):
            run_scheduler(inputs, config)
    
    def test_scheduler_handles_corrupt_tracker_pickle(self, tmp_path):
        """Document behavior with corrupt tracker pickle.
        
        Currently the scheduler may or may not validate the pickle before using it.
        This test documents that corrupt pickles should be handled gracefully.
        """
        # Create corrupt pickle file
        tracker_path = tmp_path / "tracker.pkl"
        tracker_path.write_bytes(b"corrupt data not a pickle")
        
        # When we have full error handling, this test should verify
        # that the scheduler either:
        # 1. Detects the corrupt pickle and falls back to new tracker, OR
        # 2. Raises a clear error message
        
        # For now, just verify the file exists
        assert tracker_path.exists()
        
        # TODO: Implement proper error handling in scheduler
        # and update this test to verify the behavior
    
    def test_empty_tracker_initialization(self, tmp_path):
        """Verify tracker correctly initialized when no previous data."""
        # Create minimal valid setup with NO tracker pickle
        
        # Empty target list
        target_df = pd.DataFrame(columns=[
            "Planet Name", "Star Name", "RA", "DEC", "Primary Target",
            "Number of Transits to Capture", "Transit Duration (hrs)",
            "Period (days)", "Transit Epoch (BJD_TDB-2400000.5)"
        ])
        target_csv = tmp_path / "targets.csv"
        target_df.to_csv(target_csv, index=False)
        
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            output_dir=tmp_path,
        )
        
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        
        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=data_dir / "targets",
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )
        
        inputs = SchedulerInputs(
            pandora_start=config.window_start,
            pandora_stop=config.window_end,
            sched_start=config.window_start,
            sched_stop=config.window_end,
            target_list=target_df,
            paths=paths,
            target_definition_files=[],
            primary_target_csv=target_csv,
            auxiliary_target_csv=tmp_path / "aux.csv",
            occultation_target_csv=tmp_path / "occ.csv",
            output_dir=tmp_path,
            tracker_pickle_path=None,  # No previous tracker
        )
        
        # Create empty auxiliary/occultation
        pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(inputs.auxiliary_target_csv, index=False)
        pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(inputs.occultation_target_csv, index=False)
        
        # Should complete successfully with empty schedule
        outputs = run_scheduler(inputs, config)
        
        # Verify tracker was initialized
        assert outputs.tracker is not None
        assert isinstance(outputs.tracker, pd.DataFrame)
        
        # Verify schedule is valid (may be empty)
        assert outputs.schedule is not None
        assert isinstance(outputs.schedule, pd.DataFrame)
    
    def test_invalid_observation_window(self):
        """Test that invalid observation window is handled."""
        # Create config with window_end before window_start
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 2),
            window_end=datetime(2026, 1, 1),  # Before start!
        )
        
        # Currently no validation, but documents the behavior
        assert config.window_end < config.window_start
        
        # TODO: Consider adding validation to config.__post_init__
        # This would prevent invalid configs from being created
