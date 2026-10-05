from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    SchedulerState,
    _schedule_primary_target,
)


def _minimal_inputs(tmp_path: Path) -> SchedulerInputs:
    paths = SchedulerPaths.from_package_root(tmp_path)
    dummy = tmp_path / "dummy.csv"
    return SchedulerInputs(
        pandora_start=datetime(2026, 1, 1),
        pandora_stop=datetime(2026, 1, 2),
        sched_start=datetime(2026, 1, 1),
        sched_stop=datetime(2026, 1, 2),
        target_list=pd.DataFrame(),
        paths=paths,
        target_definition_files=[],
        primary_target_csv=dummy,
        auxiliary_target_csv=dummy,
        occultation_target_csv=dummy,
        output_dir=tmp_path,
        tracker_pickle_path=None,
    )


def test_full_coverage_beats_higher_quality_when_not_critical(tmp_path: Path):
    start = datetime(2026, 1, 1, 0, 0, 0)
    obs_range = pd.date_range(start, start + timedelta(hours=1), freq="min")

    # Both candidates are non-critical (Transit Factor > 2), but one is full coverage.
    temp_df = pd.DataFrame(
        [
            {
                "Planet Name": "Full",
                "RA": 0.0,
                "DEC": 0.0,
                "Obs Start": start,
                "Obs Gap Time": timedelta(),
                "Visit Duration": timedelta(hours=10),
                "Transit Coverage": 1.0,
                "SAA Overlap": 0.0,
                "Schedule Factor": 0.0,
                "Transit Factor": 3.0,
                "Quality Factor": 0.1,
                "Comments": pd.NA,
            },
            {
                "Planet Name": "Partial",
                "RA": 0.0,
                "DEC": 0.0,
                "Obs Start": start,
                "Obs Gap Time": timedelta(),
                "Visit Duration": timedelta(hours=10),
                "Transit Coverage": 0.9,
                "SAA Overlap": 0.0,
                "Schedule Factor": 0.0,
                "Transit Factor": 3.0,
                "Quality Factor": 0.99,
                "Comments": pd.NA,
            },
        ]
    )

    state = SchedulerState(
        tracker=pd.DataFrame(
            {
                "Planet Name": ["Full", "Partial"],
                "Transits Needed": [1, 1],
                "Transits Acquired": [0, 0],
                "Transits Left in Lifetime": [10, 10],
                "Transit Priority": [0, 0],
            }
        )
    )

    inputs = _minimal_inputs(tmp_path)
    config = PandoraSchedulerConfig(window_start=start, window_end=start + timedelta(days=1))

    scheduled = _schedule_primary_target(temp_df, state, inputs, config, start, obs_range)

    assert scheduled.iloc[0]["Target"] == "Full"


def test_full_coverage_tiebreaks_when_critical(tmp_path: Path):
    start = datetime(2026, 1, 1, 0, 0, 0)
    obs_range = pd.date_range(start, start + timedelta(hours=1), freq="min")

    # Critical case triggered (Transit Factor <= 2 exists). Both have same urgency.
    temp_df = pd.DataFrame(
        [
            {
                "Planet Name": "Full",
                "RA": 0.0,
                "DEC": 0.0,
                "Obs Start": start,
                "Obs Gap Time": timedelta(),
                "Visit Duration": timedelta(hours=10),
                "Transit Coverage": 1.0,
                "SAA Overlap": 0.0,
                "Schedule Factor": 0.0,
                "Transit Factor": 1.5,
                "Quality Factor": 0.1,
                "Comments": pd.NA,
            },
            {
                "Planet Name": "Partial",
                "RA": 0.0,
                "DEC": 0.0,
                "Obs Start": start,
                "Obs Gap Time": timedelta(),
                "Visit Duration": timedelta(hours=10),
                "Transit Coverage": 0.8,
                "SAA Overlap": 0.0,
                "Schedule Factor": 0.0,
                "Transit Factor": 1.5,
                "Quality Factor": 0.99,
                "Comments": pd.NA,
            },
        ]
    )

    state = SchedulerState(
        tracker=pd.DataFrame(
            {
                "Planet Name": ["Full", "Partial"],
                "Transits Needed": [1, 1],
                "Transits Acquired": [0, 0],
                "Transits Left in Lifetime": [1, 1],
                "Transit Priority": [0, 0],
            }
        )
    )

    inputs = _minimal_inputs(tmp_path)
    config = PandoraSchedulerConfig(window_start=start, window_end=start + timedelta(days=1))

    scheduled = _schedule_primary_target(temp_df, state, inputs, config, start, obs_range)

    assert scheduled.iloc[0]["Target"] == "Full"
