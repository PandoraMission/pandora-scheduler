"""Tests for per-target visit duration feature."""

from datetime import datetime, timedelta

import pandas as pd
import pytest

from pandorascheduler_rework.observation_utils import (
    MissingObsWindowError,
    TransitUnschedulableError,
    compute_edge_buffer,
    compute_transit_start_bounds,
    get_target_visit_duration,
    validate_transit_schedulable,
)


class TestGetTargetVisitDuration:
    """Tests for get_target_visit_duration function."""

    def test_returns_value_from_column(self):
        """Should return duration from Obs Window (hrs) column."""
        target_list = pd.DataFrame({
            "Planet Name": ["Planet-A", "Planet-B"],
            "Obs Window (hrs)": [6.0, 12.0],
        })

        result = get_target_visit_duration("Planet-A", target_list)
        assert result == timedelta(hours=6)

        result = get_target_visit_duration("Planet-B", target_list)
        assert result == timedelta(hours=12)

    def test_raises_when_column_missing(self):
        """Should raise when Obs Window (hrs) column doesn't exist."""
        target_list = pd.DataFrame({
            "Planet Name": ["Planet-A"],
        })

        with pytest.raises(MissingObsWindowError, match="missing required column"):
            get_target_visit_duration("Planet-A", target_list)

    def test_raises_when_planet_not_found(self):
        """Should raise when planet is not in target list."""
        target_list = pd.DataFrame({
            "Planet Name": ["Planet-A"],
            "Obs Window (hrs)": [6.0],
        })

        with pytest.raises(MissingObsWindowError, match="not present"):
            get_target_visit_duration("Planet-X", target_list)

    def test_raises_for_nan_value(self):
        """Should raise when value is NaN."""
        target_list = pd.DataFrame({
            "Planet Name": ["Planet-A"],
            "Obs Window (hrs)": [float("nan")],
        })

        with pytest.raises(MissingObsWindowError, match="missing"):
            get_target_visit_duration("Planet-A", target_list)

    def test_raises_for_invalid_value(self):
        """Should raise for non-positive values."""
        target_list = pd.DataFrame({
            "Planet Name": ["Planet-A", "Planet-B"],
            "Obs Window (hrs)": [0.0, -5.0],
        })

        with pytest.raises(MissingObsWindowError, match="invalid"):
            get_target_visit_duration("Planet-A", target_list)

        with pytest.raises(MissingObsWindowError, match="invalid"):
            get_target_visit_duration("Planet-B", target_list)


class TestComputeEdgeBuffer:
    """Tests for compute_edge_buffer function."""

    def test_short_visit_uses_short_buffer(self):
        """Visits < threshold should use short buffer."""
        result = compute_edge_buffer(
            timedelta(hours=6),
            short_visit_threshold_hours=12.0,
            short_visit_edge_buffer_hours=1.5,
            long_visit_edge_buffer_hours=4.0,
        )
        assert result == timedelta(hours=1.5)

    def test_long_visit_uses_long_buffer(self):
        """Visits >= threshold should use long buffer."""
        result = compute_edge_buffer(
            timedelta(hours=24),
            short_visit_threshold_hours=12.0,
            short_visit_edge_buffer_hours=1.5,
            long_visit_edge_buffer_hours=4.0,
        )
        assert result == timedelta(hours=4)

    def test_threshold_boundary_uses_long_buffer(self):
        """Visits exactly at threshold should use long buffer."""
        result = compute_edge_buffer(
            timedelta(hours=12),
            short_visit_threshold_hours=12.0,
            short_visit_edge_buffer_hours=1.5,
            long_visit_edge_buffer_hours=4.0,
        )
        assert result == timedelta(hours=4)

    def test_custom_parameters(self):
        """Should respect custom parameter values."""
        result = compute_edge_buffer(
            timedelta(hours=8),
            short_visit_threshold_hours=10.0,
            short_visit_edge_buffer_hours=2.0,
            long_visit_edge_buffer_hours=5.0,
        )
        assert result == timedelta(hours=2)


class TestValidateTransitSchedulable:
    """Tests for validate_transit_schedulable function."""

    def test_valid_transit_passes(self):
        """Transit that fits within visit should not raise."""
        # 2h transit + 2 * 1.5h buffer = 5h, fits in 6h visit
        validate_transit_schedulable(
            "Planet-A",
            transit_duration=timedelta(hours=2),
            visit_duration=timedelta(hours=6),
            edge_buffer=timedelta(hours=1.5),
        )

    def test_transit_too_long_raises(self):
        """Transit that doesn't fit should raise TransitUnschedulableError."""
        # 4h transit + 2 * 1.5h buffer = 7h, doesn't fit in 6h visit
        with pytest.raises(TransitUnschedulableError) as exc_info:
            validate_transit_schedulable(
                "Planet-A",
                transit_duration=timedelta(hours=4),
                visit_duration=timedelta(hours=6),
                edge_buffer=timedelta(hours=1.5),
            )
        assert "Planet-A" in str(exc_info.value)
        assert "cannot be scheduled" in str(exc_info.value)

    def test_exact_fit_passes(self):
        """Transit that exactly fits should pass."""
        # 3h transit + 2 * 1.5h buffer = 6h, exactly fits in 6h visit
        validate_transit_schedulable(
            "Planet-A",
            transit_duration=timedelta(hours=3),
            visit_duration=timedelta(hours=6),
            edge_buffer=timedelta(hours=1.5),
        )


class TestComputeTransitStartBounds:
    """Tests for compute_transit_start_bounds function."""

    def test_24h_visit_matches_legacy(self):
        """24h visit with 4h buffer should match legacy -20h/-4h logic."""
        transit_start = datetime(2026, 3, 1, 12, 0, 0)
        transit_stop = datetime(2026, 3, 1, 14, 0, 0)  # 2h transit

        earliest, latest = compute_transit_start_bounds(
            transit_start,
            transit_stop,
            visit_duration=timedelta(hours=24),
            edge_buffer=timedelta(hours=4),
        )

        # Legacy: early_start = transit_stop - 20h = 2026-03-01 14:00 - 20h = 2026-02-28 18:00
        # Legacy: late_start = transit_start - 4h = 2026-03-01 12:00 - 4h = 2026-03-01 08:00
        expected_earliest = transit_stop - timedelta(hours=20)
        expected_latest = transit_start - timedelta(hours=4)

        assert earliest == expected_earliest
        assert latest == expected_latest

    def test_6h_visit_with_short_buffer(self):
        """6h visit with 1.5h buffer should compute correct bounds."""
        transit_start = datetime(2026, 3, 1, 12, 0, 0)
        transit_stop = datetime(2026, 3, 1, 14, 0, 0)  # 2h transit

        earliest, latest = compute_transit_start_bounds(
            transit_start,
            transit_stop,
            visit_duration=timedelta(hours=6),
            edge_buffer=timedelta(hours=1.5),
        )

        # earliest = transit_stop + edge_buffer - visit_duration
        #          = 14:00 + 1.5h - 6h = 09:30
        # latest = transit_start - edge_buffer
        #        = 12:00 - 1.5h = 10:30
        expected_earliest = datetime(2026, 3, 1, 9, 30, 0)
        expected_latest = datetime(2026, 3, 1, 10, 30, 0)

        assert earliest == expected_earliest
        assert latest == expected_latest

    def test_no_valid_window_returns_inverted_bounds(self):
        """When no valid window exists, earliest > latest."""
        transit_start = datetime(2026, 3, 1, 12, 0, 0)
        transit_stop = datetime(2026, 3, 1, 15, 0, 0)  # 3h transit

        # 3h transit + 2 * 2h buffer = 7h, but visit is only 6h
        # This should result in earliest > latest
        earliest, latest = compute_transit_start_bounds(
            transit_start,
            transit_stop,
            visit_duration=timedelta(hours=6),
            edge_buffer=timedelta(hours=2),
        )

        # earliest = 15:00 + 2h - 6h = 11:00
        # latest = 12:00 - 2h = 10:00
        # earliest (11:00) > latest (10:00) -> no valid window
        assert earliest > latest


class TestIntegration:
    """Integration tests for per-target visit duration."""

    def test_6h_target_produces_valid_bounds(self):
        """A 6h target with typical transit should have valid start bounds."""
        # Simulate TOI-700b: 6h visit, ~2.17h transit
        transit_duration = timedelta(hours=2.17)
        visit_duration = timedelta(hours=6)
        edge_buffer = compute_edge_buffer(
            visit_duration,
            short_visit_threshold_hours=12.0,
            short_visit_edge_buffer_hours=1.5,
            long_visit_edge_buffer_hours=4.0,
        )

        # Should not raise
        validate_transit_schedulable(
            "TOI-700b",
            transit_duration,
            visit_duration,
            edge_buffer,
        )

        # Start bounds should be valid (earliest <= latest)
        transit_start = datetime(2026, 3, 1, 12, 0, 0)
        transit_stop = transit_start + transit_duration

        earliest, latest = compute_transit_start_bounds(
            transit_start, transit_stop, visit_duration, edge_buffer
        )

        assert earliest <= latest, f"No valid window: {earliest} > {latest}"

    def test_target_list_lookup_integration(self):
        """Full integration of target list lookup with edge buffer computation."""
        target_list = pd.DataFrame({
            "Planet Name": ["WASP-107b", "TOI-700b"],
            "Obs Window (hrs)": [24.0, 6.0],
        })

        # WASP-107b: 24h visit -> 4h buffer
        duration = get_target_visit_duration("WASP-107b", target_list)
        buffer = compute_edge_buffer(duration)
        assert duration == timedelta(hours=24)
        assert buffer == timedelta(hours=4)

        # TOI-700b: 6h visit -> 1.5h buffer
        duration = get_target_visit_duration("TOI-700b", target_list)
        buffer = compute_edge_buffer(duration)
        assert duration == timedelta(hours=6)
        assert buffer == timedelta(hours=1.5)
