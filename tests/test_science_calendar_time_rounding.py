"""Unit tests for time rounding in science_calendar module.

These tests ensure that all time rounding operations use consistent
rounding-to-nearest-second logic, matching the legacy behavior from
the removed time_utils.round_to_nearest_second() function.
"""

from datetime import datetime, timedelta

import pandas as pd
from astropy.time import Time

from pandorascheduler_rework.science_calendar import (
    _extract_visibility_segment,
    _parse_datetime,
    _transit_windows,
)


class TestTimeRounding:
    """Test that time rounding is consistent and matches legacy behavior."""

    def test_round_to_nearest_second_logic(self):
        """Test the rounding logic used throughout science_calendar.
        
        This is the expected behavior from the legacy round_to_nearest_second():
        - Times with microseconds < 500,000 should round down
        - Times with microseconds >= 500,000 should round up
        """
        # Test case: exactly 499,999 microseconds (should round down)
        dt_down = datetime(2026, 1, 1, 12, 0, 0, 499_999)
        rounded_down = (dt_down + timedelta(microseconds=500_000)).replace(microsecond=0)
        assert rounded_down == datetime(2026, 1, 1, 12, 0, 0)
        
        # Test case: exactly 500,000 microseconds (should round up)
        dt_up = datetime(2026, 1, 1, 12, 0, 0, 500_000)
        rounded_up = (dt_up + timedelta(microseconds=500_000)).replace(microsecond=0)
        assert rounded_up == datetime(2026, 1, 1, 12, 0, 1)
        
        # Test case: 750,000 microseconds (should round up)
        dt_up2 = datetime(2026, 1, 1, 12, 0, 0, 750_000)
        rounded_up2 = (dt_up2 + timedelta(microseconds=500_000)).replace(microsecond=0)
        assert rounded_up2 == datetime(2026, 1, 1, 12, 0, 1)
        
        # Test case: 1 microsecond (should round down)
        dt_min = datetime(2026, 1, 1, 12, 0, 0, 1)
        rounded_min = (dt_min + timedelta(microseconds=500_000)).replace(microsecond=0)
        assert rounded_min == datetime(2026, 1, 1, 12, 0, 0)

    def test_transit_windows_rounds_to_nearest_second(self):
        """Test that _transit_windows rounds transit times to nearest second.
        
        This test would have caught the bug where transit times were being
        truncated instead of rounded.
        """
        # Create a transit dataframe with sub-second precision times
        # MJD with fractional seconds that when converted to datetime
        # will have microseconds >= 500,000 (should round up)
        
        # 2026-02-05 00:00:00.750 UTC (should round to 00:00:01)
        mjd_start = Time("2026-02-05T00:00:00.750", format="isot", scale="utc").mjd
        # 2026-02-05 01:00:00.250 UTC (should round to 01:00:00)
        mjd_stop = Time("2026-02-05T01:00:00.250", format="isot", scale="utc").mjd
        
        transit_df = pd.DataFrame({
            "Transit_Start": [mjd_start],
            "Transit_Stop": [mjd_stop],
        })
        
        result = _transit_windows(transit_df)
        assert result is not None
        
        start_times, stop_times = result
        assert len(start_times) == 1
        assert len(stop_times) == 1
        
        # Start time with 750ms should round UP to next second
        expected_start = datetime(2026, 2, 5, 0, 0, 1)
        assert start_times[0] == expected_start, (
            f"Expected start time to round up to {expected_start}, "
            f"but got {start_times[0]}"
        )
        
        # Stop time with 250ms should round DOWN (stay at same second)
        expected_stop = datetime(2026, 2, 5, 1, 0, 0)
        assert stop_times[0] == expected_stop, (
            f"Expected stop time to round down to {expected_stop}, "
            f"but got {stop_times[0]}"
        )

    def test_transit_windows_no_microseconds(self):
        """Test that _transit_windows works correctly with already-rounded times."""
        # Create times that are already on exact seconds
        mjd_start = Time("2026-02-05T00:00:00.000", format="isot", scale="utc").mjd
        mjd_stop = Time("2026-02-05T01:00:00.000", format="isot", scale="utc").mjd
        
        transit_df = pd.DataFrame({
            "Transit_Start": [mjd_start],
            "Transit_Stop": [mjd_stop],
        })
        
        result = _transit_windows(transit_df)
        assert result is not None
        
        start_times, stop_times = result
        
        # Should remain unchanged
        assert start_times[0] == datetime(2026, 2, 5, 0, 0, 0)
        assert stop_times[0] == datetime(2026, 2, 5, 1, 0, 0)

    def test_transit_windows_edge_case_exactly_half_second(self):
        """Test the exact boundary case: 0.5 seconds should round UP."""
        mjd_start = Time("2026-02-05T12:30:45.500", format="isot", scale="utc").mjd
        mjd_stop = Time("2026-02-05T13:15:30.500", format="isot", scale="utc").mjd
        
        transit_df = pd.DataFrame({
            "Transit_Start": [mjd_start],
            "Transit_Stop": [mjd_stop],
        })
        
        result = _transit_windows(transit_df)
        assert result is not None
        
        start_times, stop_times = result
        
        # Exactly 0.5 seconds should round UP
        assert start_times[0] == datetime(2026, 2, 5, 12, 30, 46)
        assert stop_times[0] == datetime(2026, 2, 5, 13, 15, 31)

    def test_extract_visibility_segment_rounds_consistently(self):
        """Test that visibility times are also rounded to nearest second."""
        # Create a visibility dataframe with sub-second times
        mjd_times = [
            Time("2026-02-05T00:00:00.250", format="isot", scale="utc").mjd,
            Time("2026-02-05T00:01:00.750", format="isot", scale="utc").mjd,
            Time("2026-02-05T00:02:00.499", format="isot", scale="utc").mjd,
            Time("2026-02-05T00:03:00.500", format="isot", scale="utc").mjd,
        ]
        
        visibility_df = pd.DataFrame({
            "Time(MJD_UTC)": mjd_times,
            "Visible": [1, 1, 0, 1],
        })
        
        start = datetime(2026, 2, 5, 0, 0, 0)
        stop = datetime(2026, 2, 5, 0, 4, 0)
        
        times, flags = _extract_visibility_segment(visibility_df, start, stop, min_sequence_minutes=1)
        
        assert len(times) == 4
        
        # Check that times are properly rounded
        assert times[0] == datetime(2026, 2, 5, 0, 0, 0)  # 0.25s rounds down
        assert times[1] == datetime(2026, 2, 5, 0, 1, 1)  # 0.75s rounds up
        assert times[2] == datetime(2026, 2, 5, 0, 2, 0)  # 0.499s rounds down
        assert times[3] == datetime(2026, 2, 5, 0, 3, 1)  # 0.500s rounds up

    def test_parse_datetime_handles_various_formats(self):
        """Test that datetime parsing is robust across different formats."""
        # ISO format with Z
        dt1 = _parse_datetime("2026-02-05T12:30:45Z")
        assert dt1 == datetime(2026, 2, 5, 12, 30, 45)
        
        # Space-separated format
        dt2 = _parse_datetime("2026-02-05 12:30:45")
        assert dt2 == datetime(2026, 2, 5, 12, 30, 45)
        
        # With microseconds
        dt3 = _parse_datetime("2026-02-05 12:30:45.123456")
        assert dt3 == datetime(2026, 2, 5, 12, 30, 45, 123456)
        
        # Already a datetime object
        dt_obj = datetime(2026, 2, 5, 12, 30, 45)
        dt4 = _parse_datetime(dt_obj)
        assert dt4 == dt_obj
        
        # Invalid format returns None
        dt5 = _parse_datetime("invalid")
        assert dt5 is None


class TestTransitWindowsRegression:
    """Specific regression tests for the transit time rounding bug."""
    
    def test_transit_times_not_truncated(self):
        """Regression test: ensure transit times are NOT truncated to whole seconds.
        
        The bug was that transit times used .replace(second=0, microsecond=0)
        which truncates the seconds field instead of rounding microseconds.
        This test would have caught that error.
        """
        # Create a time at 12:30:45.750 - the seconds should be preserved!
        mjd_time = Time("2026-02-05T12:30:45.750", format="isot", scale="utc").mjd
        
        transit_df = pd.DataFrame({
            "Transit_Start": [mjd_time],
            "Transit_Stop": [mjd_time],
        })
        
        result = _transit_windows(transit_df)
        start_times, stop_times = result
        
        # The bug would have set seconds=0, giving 12:30:00
        # The correct behavior rounds to 12:30:46
        buggy_result = datetime(2026, 2, 5, 12, 30, 0)  # What the bug would produce
        correct_result = datetime(2026, 2, 5, 12, 30, 46)  # What we should get
        
        assert start_times[0] != buggy_result, (
            "Transit time appears to be truncated (seconds set to 0). "
            "This suggests the old buggy behavior."
        )
        assert start_times[0] == correct_result, (
            f"Transit time should round to {correct_result}, got {start_times[0]}"
        )

    def test_multiple_transits_all_rounded(self):
        """Test that all transit times in a batch are properly rounded."""
        mjd_times_start = [
            Time("2026-02-05T00:00:00.100", format="isot", scale="utc").mjd,
            Time("2026-02-05T01:00:00.600", format="isot", scale="utc").mjd,
            Time("2026-02-05T02:00:00.499", format="isot", scale="utc").mjd,
            Time("2026-02-05T03:00:00.999", format="isot", scale="utc").mjd,
        ]
        
        mjd_times_stop = [
            Time("2026-02-05T00:30:00.200", format="isot", scale="utc").mjd,
            Time("2026-02-05T01:30:00.800", format="isot", scale="utc").mjd,
            Time("2026-02-05T02:30:00.400", format="isot", scale="utc").mjd,
            Time("2026-02-05T03:30:00.501", format="isot", scale="utc").mjd,
        ]
        
        transit_df = pd.DataFrame({
            "Transit_Start": mjd_times_start,
            "Transit_Stop": mjd_times_stop,
        })
        
        result = _transit_windows(transit_df)
        start_times, stop_times = result
        
        # Expected rounded times
        expected_starts = [
            datetime(2026, 2, 5, 0, 0, 0),  # 100ms rounds down
            datetime(2026, 2, 5, 1, 0, 1),  # 600ms rounds up
            datetime(2026, 2, 5, 2, 0, 0),  # 499ms rounds down
            datetime(2026, 2, 5, 3, 0, 1),  # 999ms rounds up
        ]
        
        expected_stops = [
            datetime(2026, 2, 5, 0, 30, 0),  # 200ms rounds down
            datetime(2026, 2, 5, 1, 30, 1),  # 800ms rounds up
            datetime(2026, 2, 5, 2, 30, 0),  # 400ms rounds down
            datetime(2026, 2, 5, 3, 30, 1),  # 501ms rounds up
        ]
        
        assert start_times == expected_starts
        assert stop_times == expected_stops
