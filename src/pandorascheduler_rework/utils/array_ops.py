"""Array processing utilities for observation scheduling."""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import List, Optional, Tuple

import numpy as np


def remove_short_sequences(
    array: np.ndarray, sequence_too_short: int
) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
    """Remove visibility sequences shorter than threshold.

    This performs two passes:
    1. Zero out short visible runs (short 1-runs → 0) so that
       briefly-visible windows are not scheduled.
    2. Fill short non-visible gaps (short 0-runs between 1-runs → 1)
       so that brief occultation dips do not fragment the schedule
       into unusably short sequences.

    Args:
        array: Binary visibility array (1=visible, 0=not visible)
        sequence_too_short: Minimum sequence length to keep

    Returns:
        Tuple of (cleaned array, list of removed visible spans)
    """
    cleaned = np.asarray(array, dtype=float).copy()
    start_index = None
    spans: List[Tuple[int, int]] = []

    # Pass 1: remove short visible runs (1-runs shorter than threshold).
    for idx, value in enumerate(cleaned):
        if value == 1 and start_index is None:
            start_index = idx
            continue
        if value == 0 and start_index is not None:
            if idx - start_index < sequence_too_short:
                spans.append((start_index, idx - 1))
            start_index = None

    if start_index is not None and len(cleaned) - start_index < sequence_too_short:
        spans.append((start_index, len(cleaned) - 1))

    for start_idx, stop_idx in spans:
        cleaned[start_idx : stop_idx + 1] = 0.0

    # Pass 2: fill short non-visible gaps (0-runs shorter than threshold
    # that sit between two visible regions).  Trailing gaps at the end of
    # the array are left alone.
    gap_start: Optional[int] = None
    for idx, value in enumerate(cleaned):
        if value == 0 and gap_start is None:
            gap_start = idx
        elif value != 0 and gap_start is not None:
            if idx - gap_start < sequence_too_short and gap_start > 0:
                cleaned[gap_start:idx] = 1.0
            gap_start = None

    return cleaned, spans


def break_long_sequences(
    start: datetime,
    end: datetime,
    step: timedelta,
    min_chunk: Optional[timedelta] = None,
) -> List[Tuple[datetime, datetime]]:
    """Break long time range into smaller chunks.

    Args:
        start: Start time
        end: End time
        step: Maximum chunk duration
        min_chunk: If the last chunk would be shorter than this,
            merge it into the previous chunk instead of emitting it
            as a standalone segment.  ``None`` disables merging.

    Returns:
        List of (start, end) tuples for each chunk
    """
    ranges: List[Tuple[datetime, datetime]] = []
    current = start
    while current < end:
        next_val = min(current + step, end)
        ranges.append((current, next_val))
        current = next_val

    # Absorb a short trailing chunk into the previous one.
    if (
        min_chunk is not None
        and len(ranges) >= 2
        and (ranges[-1][1] - ranges[-1][0]) < min_chunk
    ):
        ranges[-2] = (ranges[-2][0], ranges[-1][1])
        ranges.pop()

    return ranges
