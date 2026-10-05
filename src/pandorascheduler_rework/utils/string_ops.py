"""String processing utilities for target names."""

from __future__ import annotations

import re

import pandas as pd

# Pre-compile regex pattern for performance
_PLANET_SUFFIX_PATTERN = re.compile(r"\s+[a-z]$", flags=re.ASCII)


def remove_suffix(value: str) -> str:
    """Return value with trailing planet suffix removed.

    The legacy helper accepted identifiers like "WASP-107 b" and stripped the
    final " b" so that callers could recover the stellar host name. The rework
    pipeline relies on the same behaviour when looking up per-star visibility
    files.

    Args:
        value: Target identifier (e.g., "WASP-107 b")

    Returns:
        Identifier with suffix removed (e.g., "WASP-107")

    Examples:
        >>> remove_suffix("WASP-107 b")
        "WASP-107"
        >>> remove_suffix("HD-189733")
        "HD-189733"
    """
    return _PLANET_SUFFIX_PATTERN.sub("", value)


def target_identifier(row: pd.Series) -> str:
    """Extract target identifier from manifest row.

    Prefers Planet Name if present, otherwise uses Star Name.
    Removes trailing single-letter suffix from planet names.

    Args:
        row: DataFrame row with target information

    Returns:
        Target identifier string
    """
    planet = row.get("Planet Name") if isinstance(row, pd.Series) else None
    if planet is not None and pd.notna(planet):
        return re.sub(r"\s+([A-Za-z])$", r"\1", str(planet))

    star = row.get("Star Name") if isinstance(row, pd.Series) else None
    if star is not None and pd.notna(star):
        return str(star)

    return ""
