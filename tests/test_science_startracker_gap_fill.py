"""Tests for science-calendar star-tracker-only gap filling."""

from datetime import datetime, timedelta
from types import MethodType

import numpy as np

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.science_calendar import _ScienceCalendarBuilder


def _builder(*, enabled=True, max_gap_minutes=5):
    builder = object.__new__(_ScienceCalendarBuilder)
    builder.config = PandoraSchedulerConfig(
        window_start=datetime(2026, 3, 1),
        window_end=datetime(2026, 3, 2),
        allow_science_startracker_gap_fill=enabled,
        science_startracker_gap_max_minutes=max_gap_minutes,
        min_science_sequence_minutes=8,
    )
    return builder


def _constant_visibility(builder, visible):
    def interval_visibility(self, ra, dec, start, stop, visibility_config):
        minutes = int((stop - start).total_seconds() / 60)
        return np.full(minutes, visible, dtype=bool)

    builder._science_interval_visibility = MethodType(
        interval_visibility, builder
    )


def test_fills_short_boresight_visible_gap():
    builder = _builder()
    _constant_visibility(builder, True)
    t0 = datetime(2026, 3, 1, 0, 0)
    segments = [
        (t0, t0 + timedelta(minutes=20), True),
        (t0 + timedelta(minutes=20), t0 + timedelta(minutes=24), False),
        (t0 + timedelta(minutes=24), t0 + timedelta(minutes=40), True),
    ]

    adjusted, filled = builder._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    )

    assert adjusted == [(t0, t0 + timedelta(minutes=40), True)]
    assert filled == [
        (t0 + timedelta(minutes=20), t0 + timedelta(minutes=24))
    ]


def test_does_not_fill_gap_when_boresight_is_not_visible():
    builder = _builder()
    _constant_visibility(builder, False)
    t0 = datetime(2026, 3, 1, 0, 0)
    gap = (t0 + timedelta(minutes=20), t0 + timedelta(minutes=24), False)
    segments = [
        (t0, t0 + timedelta(minutes=20), True),
        gap,
        (t0 + timedelta(minutes=24), t0 + timedelta(minutes=40), True),
    ]

    adjusted, filled = builder._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    )

    assert adjusted == segments
    assert filled == []


def test_does_not_fill_edge_gap():
    builder = _builder()
    _constant_visibility(builder, True)
    t0 = datetime(2026, 3, 1, 0, 0)
    segments = [
        (t0, t0 + timedelta(minutes=4), False),
        (t0 + timedelta(minutes=4), t0 + timedelta(minutes=20), True),
    ]

    adjusted, filled = builder._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    )

    assert adjusted == segments
    assert filled == []


def test_does_not_fill_gap_after_subminimum_visible_fragment():
    builder = _builder()
    _constant_visibility(builder, True)
    t0 = datetime(2026, 3, 1, 0, 0)
    segments = [
        (t0, t0 + timedelta(minutes=3), True),
        (t0 + timedelta(minutes=3), t0 + timedelta(minutes=15), False),
        (t0 + timedelta(minutes=15), t0 + timedelta(minutes=35), True),
    ]

    adjusted, filled = builder._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    )

    assert adjusted == segments
    assert filled == []


def test_fills_gap_before_subminimum_trailing_fragment():
    builder = _builder(max_gap_minutes=15)
    _constant_visibility(builder, True)
    t0 = datetime(2026, 3, 1, 0, 0)
    segments = [
        (t0, t0 + timedelta(minutes=20), True),
        (t0 + timedelta(minutes=20), t0 + timedelta(minutes=32), False),
        (t0 + timedelta(minutes=32), t0 + timedelta(minutes=35), True),
    ]

    adjusted, filled = builder._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    )

    assert adjusted == [(t0, t0 + timedelta(minutes=35), True)]
    assert filled == [(t0 + timedelta(minutes=20), t0 + timedelta(minutes=32))]


def test_does_not_fill_isolated_gap():
    builder = _builder()
    _constant_visibility(builder, True)
    t0 = datetime(2026, 3, 1, 0, 0)
    segments = [
        (t0, t0 + timedelta(minutes=20), True),
        (t0 + timedelta(minutes=20), t0 + timedelta(minutes=24), False),
        (t0 + timedelta(minutes=24), t0 + timedelta(minutes=28), False),
        (t0 + timedelta(minutes=28), t0 + timedelta(minutes=40), True),
    ]

    adjusted, filled = builder._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    )

    assert adjusted == [
        (t0, t0 + timedelta(minutes=20), True),
        (t0 + timedelta(minutes=20), t0 + timedelta(minutes=28), False),
        (t0 + timedelta(minutes=28), t0 + timedelta(minutes=40), True),
    ]
    assert filled == []


def test_disabled_or_long_gap_is_unchanged():
    t0 = datetime(2026, 3, 1, 0, 0)
    segments = [
        (t0, t0 + timedelta(minutes=20), True),
        (t0 + timedelta(minutes=20), t0 + timedelta(minutes=26), False),
        (t0 + timedelta(minutes=26), t0 + timedelta(minutes=40), True),
    ]

    disabled = _builder(enabled=False)
    _constant_visibility(disabled, True)
    assert disabled._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    ) == (segments, [])

    too_long = _builder(max_gap_minutes=5)
    _constant_visibility(too_long, True)
    assert too_long._fill_science_segments_with_startracker_gaps(
        segments, 10.0, 20.0
    ) == (segments, [])
