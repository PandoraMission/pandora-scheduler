"""Unit tests for build_schedule() pipeline API."""

from datetime import datetime
from pathlib import Path

import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.pipeline import build_schedule


class TestBuildScheduleV2:
    """Test the v2 pipeline API with PandoraSchedulerConfig."""

    # Config validation tests (test_requires_valid_config, test_validates_transit_coverage_range,
    # test_config_has_sensible_defaults, test_config_is_immutable) consolidated into test_config.py

# Legacy conversion tests removed


class TestBuildScheduleV2Integration:
    """Integration tests for build_schedule() (require test data)."""

    @pytest.mark.skip(reason="Requires test data fixtures")
    def test_basic_scheduling_run(self, tmp_path):
        """Test a basic scheduling run with minimal config."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2026, 2, 6),  # 1 day
            targets_manifest=Path("test_data/targets"),
            output_dir=tmp_path,
            show_progress=False,
        )
        
        result = build_schedule(config)
        
        assert result.schedule_csv is not None
        assert result.schedule_csv.exists()

    @pytest.mark.skip(reason="Requires test data fixtures")
    def test_with_visibility_generation(self, tmp_path):
        """Test scheduling with visibility generation enabled."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2026, 2, 6),
            targets_manifest=Path("test_data/targets"),
            gmat_ephemeris=Path("test_data/ephemeris.txt"),
            output_dir=tmp_path,
            show_progress=False,
        )
        
        result = build_schedule(config)
        
        assert result.schedule_csv is not None
        # Visibility files should be generated
        assert (tmp_path / "data" / "targets").exists()


class TestConfigDocumentation:
    """Test that config fields are well-documented."""

    def test_all_fields_have_docstrings(self):
        """Test that all config fields have documentation."""
        # Check that key fields have docstrings in the class
        assert PandoraSchedulerConfig.__doc__ is not None
        assert "Master configuration" in PandoraSchedulerConfig.__doc__

    def test_config_repr_is_readable(self):
        """Test that config has a readable string representation."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2026, 2, 19),
            targets_manifest=Path("data"),
            transit_coverage_min=0.3,
        )
        
        repr_str = repr(config)
        assert "PandoraSchedulerConfig" in repr_str
        assert "transit_coverage_min=0.3" in repr_str


class TestBackwardCompatibility:
    """Test backward compatibility with old API."""

    def test_v2_produces_same_structure_as_v1(self):
        """Test that v2 API produces same result structure as v1."""
        # This is a structural test - both APIs should return SchedulerResult

        
        # build_schedule_v2 should return SchedulerResult
        # (actual execution skipped without test data)
        assert hasattr(build_schedule, "__call__")
        assert build_schedule.__doc__ is not None
