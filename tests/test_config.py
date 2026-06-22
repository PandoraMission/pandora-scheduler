"""Unit tests for unified configuration system."""

from datetime import datetime, timedelta
from pathlib import Path

import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig


class TestPandoraSchedulerConfig:
    """Tests for unified configuration class."""
    
    def test_minimal_config(self):
        """Test creating config with minimal required parameters."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
        )
        
        assert config.window_start == datetime(2026, 2, 5)
        assert config.window_end == datetime(2027, 2, 5)
        assert config.schedule_step == timedelta(hours=24)  # Default
        
    def test_full_config(self):
        """Test creating config with all parameters."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            targets_manifest=Path("targets"),
            gmat_ephemeris=Path("ephemeris.csv"),
            output_dir=Path("output"),
            transit_coverage_min=0.3,
            transit_scheduling_weights=(0.6, 0.2, 0.2),
            show_progress=True,
        )
        
        assert config.transit_coverage_min == 0.3
        assert config.transit_scheduling_weights == (0.6, 0.2, 0.2)
        assert config.show_progress is True
        
    def test_transit_scheduling_weights_validation_pass(self):
        """Test that transit_scheduling_weights summing to 1.0 validates."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            transit_scheduling_weights=(0.7, 0.1, 0.2),  # Sums to 1.0
        )
        assert sum(config.transit_scheduling_weights) == pytest.approx(1.0)
        
    def test_transit_scheduling_weights_validation_fail(self):
        """Test that transit_scheduling_weights not summing to 1.0 raises error."""
        with pytest.raises(ValueError, match="transit_scheduling_weights must sum to 1.0"):
            PandoraSchedulerConfig(
                window_start=datetime(2026, 2, 5),
                window_end=datetime(2027, 2, 5),
                transit_scheduling_weights=(1.0, 0.5, 0.5),  # Sums to 2.0!
            )
            
    def test_transit_coverage_min_validation_pass(self):
        """Test that transit_coverage_min in [0, 1] validates."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            transit_coverage_min=0.5,
        )
        assert config.transit_coverage_min == 0.5
        
    def test_transit_coverage_min_validation_fail_negative(self):
        """Test that negative transit_coverage_min raises error."""
        with pytest.raises(ValueError, match="transit_coverage_min must be in"):
            PandoraSchedulerConfig(
                window_start=datetime(2026, 2, 5),
                window_end=datetime(2027, 2, 5),
                transit_coverage_min=-0.1,
            )
            
    def test_transit_coverage_min_validation_fail_too_large(self):
        """Test that transit_coverage_min > 1.0 raises error."""
        with pytest.raises(ValueError, match="transit_coverage_min must be in"):
            PandoraSchedulerConfig(
                window_start=datetime(2026, 2, 5),
                window_end=datetime(2027, 2, 5),
                transit_coverage_min=1.5,
            )
            
# Legacy conversion tests removed
            
    def test_defaults_are_sensible(self):
        """Test that all default values are sensible."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
        )
        
        # Defaults should be the same as in the classes they replace
        assert config.transit_coverage_min == 0.2
        assert config.transit_scheduling_weights == (0.8, 0.0, 0.2)
        assert config.sun_avoidance_deg == 91.0
        assert config.moon_avoidance_deg == 25.0
        assert config.earth_avoidance_deg == 110.0
        assert config.earth_avoidance_day_deg is None
        assert config.earth_avoidance_night_deg is None
        assert config.st_required == 1
        assert config.roll_step_deg == 2.0
        assert config.min_power_frac == 0.7
        assert config.obs_sequence_duration_min == 90
        assert config.priority_buffer is False
        assert config.priority_buffer_minutes == 80
        assert config.show_progress is False

    def test_priority_buffer_minutes_validation_fail_negative(self):
        """Test that negative priority_buffer_minutes raises error."""
        with pytest.raises(ValueError, match="priority_buffer_minutes must be >= 0"):
            PandoraSchedulerConfig(
                window_start=datetime(2026, 2, 5),
                window_end=datetime(2027, 2, 5),
                priority_buffer_minutes=-1,
            )
        
    def test_frozen_immutable(self):
        """Test that config is immutable (frozen dataclass)."""
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
        )
        
        with pytest.raises(AttributeError):
            config.transit_coverage_min = 0.5  # Can't modify frozen dataclass
