"""Visibility catalog generation for the Pandora scheduler rework.

This module reimplements the legacy `pandorascheduler.transits` functionality in a
structured, testable fashion while preserving output parity. The public API will
be fleshed out incrementally; keep an eye on `docs/visibility-plan.md` for the
latest roadmap.
"""

from .catalog import build_visibility_catalog  # noqa: F401
