import logging
from datetime import datetime, timedelta

import pandas as pd
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    run_scheduler,
)
from pandorascheduler_rework.utils.io import build_star_visibility_path


def test_aux_falls_back_to_over_requested_with_warning(tmp_path, caplog):
    caplog.set_level(logging.WARNING)

    window_start = datetime(2026, 1, 1, 0, 0, 0)
    window_end = datetime(2026, 1, 1, 6, 0, 0)

    # Arrange minimal scheduler inputs
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    aux_targets_dir = data_dir / "aux_targets"
    aux_targets_dir.mkdir()

    paths = SchedulerPaths(
        package_dir=tmp_path,
        data_dir=data_dir,
        targets_dir=data_dir / "targets",
        aux_targets_dir=aux_targets_dir,
        baseline_dir=data_dir / "baseline",
    )
    paths.targets_dir.mkdir()

    # Empty primary target list so scheduler always goes to aux scheduling
    target_df = pd.DataFrame(
        columns=[
            "Planet Name",
            "Star Name",
            "RA",
            "DEC",
            "Primary Target",
            "Number of Transits to Capture",
            "Transit Duration (hrs)",
            "Period (days)",
            "Transit Epoch (BJD_TDB-2400000.5)",
        ]
    )
    primary_csv = tmp_path / "primary.csv"
    target_df.to_csv(primary_csv, index=False)

    # Create a single auxiliary target with a small requested-hour budget.
    # The scheduling window is 2 hours, so after the first visit the target will
    # have met/exceeded its requested hours.
    aux_df = pd.DataFrame(
        [
            {
                "Star Name": "AuxStar",
                "RA": 0.0,
                "DEC": 0.0,
                "Priority": 1.0,
                "Number of Hours Requested": 1.0,
            }
        ]
    )
    (data_dir / "auxiliary-standard_targets.csv").write_text(aux_df.to_csv(index=False))

    # Make the target fully visible for the whole window
    vis_path = build_star_visibility_path(aux_targets_dir, "AuxStar")
    vis_path.parent.mkdir(parents=True, exist_ok=True)
    start_mjd = Time(window_start).mjd
    two_hours_days = 2.0 / 24.0
    pd.DataFrame(
        {
            "Time(MJD_UTC)": [
                start_mjd,
                start_mjd + two_hours_days,
                start_mjd + 2 * two_hours_days,
                start_mjd + 3 * two_hours_days,
            ],
            "Visible": [1, 1, 1, 1],
        }
    ).to_parquet(vis_path, index=False)

    config = PandoraSchedulerConfig(
        window_start=window_start,
        window_end=window_end,
        schedule_step=timedelta(hours=2),
        std_obs_frequency_days=999999.0,  # avoid STD scheduling in this test
    )

    inputs = SchedulerInputs(
        pandora_start=config.window_start,
        pandora_stop=config.window_end,
        sched_start=config.window_start,
        sched_stop=config.window_end,
        target_list=target_df,
        paths=paths,
        target_definition_files=["exoplanet", "auxiliary-standard"],
        primary_target_csv=primary_csv,
        auxiliary_target_csv=tmp_path / "aux.csv",
        occultation_target_csv=tmp_path / "occ.csv",
        output_dir=tmp_path,
        tracker_pickle_path=None,
    )

    # These CSVs act as requested-hours catalogs for the observation report
    pd.DataFrame(
        [{"Star Name": "AuxStar", "Number of Hours Requested": 1.0}]
    ).to_csv(inputs.auxiliary_target_csv, index=False)
    pd.DataFrame(columns=["Star Name", "Number of Hours Requested"]).to_csv(
        inputs.occultation_target_csv, index=False
    )

    # Act
    outputs = run_scheduler(inputs, config)

    # Assert: we scheduled AuxStar at least once
    assert (outputs.schedule["Target"] == "AuxStar").any()

    # Assert: warning emitted for over-requested fallback
    assert any(
        "has met requested observation time" in rec.getMessage() for rec in caplog.records
    )
