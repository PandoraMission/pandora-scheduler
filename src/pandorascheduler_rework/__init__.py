"""Modernized scheduler implementation living alongside the legacy pipeline.

The original entry points remain under :mod:`pandorascheduler`.  Modules here
provide a clean-room reimplementation that we can validate against the existing
code before switching over.
"""

from .pipeline import SchedulerResult, build_schedule

__all__ = [
    "SchedulerResult",
    "build_schedule",
]
