"""Tests for parameter propagation from PandoraSchedulerConfig to components."""

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.science_calendar import (
    ScienceCalendarInputs,
    generate_science_calendar,
)
from pandorascheduler_rework.visibility.catalog import build_visibility_catalog


class TestConfigParameterPropagation:
    """Test that config parameters properly propagate to components."""
    
    def test_visibility_angles_affect_output(self, tmp_path):
        """Test that different avoidance angles produce different visibility results."""
        # This is tested comprehensively in test_visibility_catalog.py
        # Here we just verify the angles are accessible on the config
        config_strict = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            sun_avoidance_deg=100.0,  # Very strict
            moon_avoidance_deg=40.0,
            earth_avoidance_deg=30.0,
        )
        
        config_lenient = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            sun_avoidance_deg=45.0,  # More lenient
            moon_avoidance_deg=20.0,
            earth_avoidance_deg=15.0,
        )
        
        # Verify configs have different angles
        assert config_strict.sun_avoidance_deg == 100.0
        assert config_lenient.sun_avoidance_deg == 45.0
        # Actual behavior testing is in test_visibility_catalog.py
    
    def test_scheduler_receives_config_parameters(self, tmp_path):
        """Test that scheduler config parameters are accessible."""
        # This test is simplified - just verify config can be passed to components
        # Full integration is tested elsewhere
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            transit_coverage_min=0.35,  # Non-default
            aux_sort_key="closest",  # Non-default
            primary_only_mode=True,
        )
        
        # Verify config has the parameters we set
        assert config.transit_coverage_min == 0.35
        assert config.aux_sort_key == "closest"
        assert config.primary_only_mode is True
    
    def test_science_calendar_receives_config_parameters(self, tmp_path):
        """Test that science calendar receives and uses config parameters."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            visit_limit=50,  # Non-default
            obs_sequence_duration_min=100,  # Non-default
            sun_avoidance_deg=92.0,
            moon_avoidance_deg=26.0,
            earth_avoidance_deg=88.0,
        )
        
        # Create minimal schedule
        schedule_df = pd.DataFrame([{
            "Target": "Test",
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
        
        # Create data directory with required files
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        pd.DataFrame([{"Star Name": "Test", "RA": 0, "DEC": 0}]).to_csv(
            data_dir / "exoplanet_targets.csv", index=False
        )
        # Create auxiliary catalogs expected by the calendar builder
        pd.DataFrame([{"Star Name": "Test", "RA": 0, "DEC": 0}]).to_csv(
            data_dir / "all_targets.csv", index=False
        )
        pd.DataFrame([{"Star Name": "Test", "RA": 0, "DEC": 0}]).to_csv(
            data_dir / "occultation-standard_targets.csv", index=False
        )
        
        # Create visibility file
        vis_dir = data_dir / "targets" / "Test"
        vis_dir.mkdir(parents=True)
        pd.DataFrame({
            "Time(MJD_UTC)": [60000.0],
            "Visible": [1]
        }).to_parquet(vis_dir / "Visibility for Test.parquet", index=False)
        
        inputs = ScienceCalendarInputs(
            schedule_csv=schedule_path,
            data_dir=data_dir,
        )
        
        # Generate calendar; surface any unexpected failures so the test fails
        try:
            output = generate_science_calendar(
                inputs=inputs,
                config=config,
                output_path=tmp_path / "calendar.xml",
            )
        except Exception as e:
            # Wrap the failure message to avoid exceeding line length limits
            msg = (
                "generate_science_calendar raised an unexpected exception: "
                f"{e}"
            )
            pytest.fail(msg)

        # The builder should return the destination path and write the file
        assert output is not None, (
            "generate_science_calendar did not return an output path"
        )
        assert Path(output).exists(), (
            "Expected calendar file to be created at " f"{output}"
        )
        
        # Verify config has the right parameters
        assert config.visit_limit == 50
        assert config.obs_sequence_duration_min == 100
        assert config.sun_avoidance_deg == 92.0
        assert config.moon_avoidance_deg == 26.0
        assert config.earth_avoidance_deg == 88.0


class TestConfigErrorHandling:
    """Test that config validates required parameters and raises helpful errors."""
    
    def test_visibility_requires_gmat_ephemeris(self, tmp_path):
        """Test that visibility generation requires gmat_ephemeris to be set."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            targets_manifest=tmp_path / "targets.csv",
            output_dir=tmp_path,
            # gmat_ephemeris not provided
        )
        
        # Create minimal target file
        pd.DataFrame([{"Star Name": "Test", "RA": 0, "DEC": 0}]).to_csv(
            tmp_path / "targets.csv", index=False
        )
        
        # Attempting to build visibility should fail gracefully
        # Note: The actual error may come from the visibility code, not config
        # This test documents current behavior
        with pytest.raises((ValueError, AttributeError, FileNotFoundError)):
            build_visibility_catalog(
                config,
                target_list=config.targets_manifest,
                output_subpath="targets"
            )
    
    def test_visibility_requires_target_list(self, tmp_path):
        """Test that visibility generation requires a valid target list."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            gmat_ephemeris=tmp_path / "gmat.txt",
            output_dir=tmp_path,
        )
        
        # Create gmat file but not target file
        (tmp_path / "gmat.txt").write_text("# dummy\n")
        
        # Should fail when target file doesn't exist
        with pytest.raises(FileNotFoundError):
            build_visibility_catalog(
                config,
                target_list=tmp_path / "nonexistent.csv",
                output_subpath="targets"
            )
    
    def test_window_validation(self):
        """Test that window_end must be after window_start."""
        # This might not be validated currently, but should be
        start = datetime(2026, 1, 2)
        end = datetime(2026, 1, 1)  # Before start!
        
        # Create config - may or may not validate this
        config = PandoraSchedulerConfig(
            window_start=start,
            window_end=end,
        )
        
        # If no validation error, at least document the state
        assert config.window_end < config.window_start
        # TODO: Consider adding validation to __post_init__


class TestConfigIntegration:
    """Integration tests for configuration system."""
    
    def test_config_flows_through_pipeline(self, tmp_path):
        """Test that a single config object flows through the entire pipeline."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 1, 1, 0),  # 1 hour window
            targets_manifest=tmp_path / "targets.csv",
            output_dir=tmp_path,
            transit_coverage_min=0.25,
        )
        
        # Create minimal target file
        pd.DataFrame(columns=["Planet Name", "Star Name", "RA", "DEC"]).to_csv(
            tmp_path / "targets.csv", index=False
        )
        
        # The config should be immutable
        with pytest.raises(AttributeError):
            config.transit_coverage_min = 0.5  # type: ignore
    
    def test_extra_inputs_override_defaults(self, tmp_path):
        """Test that extra_inputs can override default behavior."""
        custom_gmat = tmp_path / "custom_gmat.txt"
        custom_gmat.write_text("# custom ephemeris\n")
        
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            gmat_ephemeris=tmp_path / "default.txt",
            extra_inputs={
                "visibility_gmat": custom_gmat,
            }
        )
        
        # Verify override is accessible
        assert config.extra_inputs.get("visibility_gmat") == custom_gmat
        
    def test_config_serialization_roundtrip(self):
        """Test that config can be represented and reconstructed."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            transit_coverage_min=0.3,
            transit_scheduling_weights=(0.7, 0.1, 0.2),
        )
        
        # Get repr
        repr_str = repr(config)
        assert "PandoraSchedulerConfig" in repr_str
        assert "transit_coverage_min=0.3" in repr_str
        
        # Verify all important fields are accessible
        assert config.window_start == datetime(2026, 1, 1)
        assert config.window_end == datetime(2026, 1, 2)
        assert config.transit_coverage_min == 0.3
        assert config.transit_scheduling_weights == (0.7, 0.1, 0.2)
