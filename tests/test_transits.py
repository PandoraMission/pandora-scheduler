"""Regression tests for the legacy transit helpers."""

from datetime import datetime

import pytest

from pandorascheduler.transits import (
    _merged_interval_overlap_fraction,
    _round_datetime_to_second,
)


def test_overlap_union_does_not_double_count_overlapping_intervals():
    target_start = datetime(2026, 1, 1, 0, 0)
    target_stop = datetime(2026, 1, 1, 1, 40)
    partner_intervals = [
        (datetime(2026, 1, 1, 0, 10), datetime(2026, 1, 1, 0, 50)),
        (datetime(2026, 1, 1, 0, 30), datetime(2026, 1, 1, 1, 10)),
        (datetime(2026, 1, 1, 1, 20), datetime(2026, 1, 1, 1, 30)),
    ]

    # The union covers 70 minutes of a 100-minute target interval.
    assert _merged_interval_overlap_fraction(
        target_start, target_stop, partner_intervals
    ) == pytest.approx(0.7)


def test_overlap_union_combines_disjoint_intervals_from_multiple_partners():
    target_start = datetime(2026, 1, 1, 0, 0)
    target_stop = datetime(2026, 1, 1, 1, 0)
    partner_intervals = [
        (datetime(2025, 12, 31, 23, 50), datetime(2026, 1, 1, 0, 15)),
        (datetime(2026, 1, 1, 0, 45), datetime(2026, 1, 1, 1, 10)),
    ]

    assert _merged_interval_overlap_fraction(
        target_start, target_stop, partner_intervals
    ) == pytest.approx(0.5)


def test_overlap_union_handles_no_overlap_and_full_coverage():
    target_start = datetime(2026, 1, 1, 0, 0)
    target_stop = datetime(2026, 1, 1, 1, 0)

    assert _merged_interval_overlap_fraction(
        target_start,
        target_stop,
        [(datetime(2026, 1, 1, 2, 0), datetime(2026, 1, 1, 3, 0))],
    ) == 0.0
    assert _merged_interval_overlap_fraction(
        target_start,
        target_stop,
        [(datetime(2025, 12, 31, 23, 0), datetime(2026, 1, 1, 2, 0))],
    ) == 1.0


def test_round_datetime_to_second_uses_nearest_second():
    assert _round_datetime_to_second(
        datetime(2026, 1, 1, 0, 0, 0, 499999)
    ) == datetime(2026, 1, 1, 0, 0, 0)
    assert _round_datetime_to_second(
        datetime(2026, 1, 1, 0, 0, 0, 500000)
    ) == datetime(2026, 1, 1, 0, 0, 1)
