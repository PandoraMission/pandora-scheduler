from __future__ import annotations

import importlib.util
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.pipeline import SchedulerResult
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    run_scheduler,
)


def _load_run_scheduler_module():
    module_path = Path(__file__).resolve().parents[1] / "run_scheduler.py"
    spec = importlib.util.spec_from_file_location("run_scheduler", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load run_scheduler module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


scheduler_driver = _load_run_scheduler_module()


def test_require_manifest_exists_raises_when_manifest_missing(tmp_path):
    missing = tmp_path / "run_config_manifest.json"

    with pytest.raises(FileNotFoundError):
        scheduler_driver._require_manifest_exists(missing, tmp_path)


def test_run_scheduler_uses_free_time_when_primary_only_and_no_primary_targets(tmp_path):
    window_start = datetime(2026, 1, 1, 0, 0, 0)
    window_end = datetime(2026, 1, 1, 6, 0, 0)

    data_dir = tmp_path / "data"
    data_dir.mkdir()
    aux_targets_dir = data_dir / "aux_targets"
    aux_targets_dir.mkdir()
    targets_dir = data_dir / "targets"
    targets_dir.mkdir()

    paths = SchedulerPaths(
        package_dir=tmp_path,
        data_dir=data_dir,
        targets_dir=targets_dir,
        aux_targets_dir=aux_targets_dir,
        baseline_dir=data_dir / "baseline",
    )

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

    pd.DataFrame(
        [
            {
                "Star Name": "AuxStar",
                "RA": 0.0,
                "DEC": 0.0,
                "Priority": 1.0,
                "Number of Hours Requested": 1.0,
            }
        ]
    ).to_csv(data_dir / "auxiliary-standard_targets.csv", index=False)
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        tmp_path / "occ.csv", index=False
    )

    config = PandoraSchedulerConfig(
        window_start=window_start,
        window_end=window_end,
        schedule_step=timedelta(hours=2),
        primary_only_mode=True,
        std_obs_frequency_days=999999.0,
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
        auxiliary_target_csv=data_dir / "auxiliary-standard_targets.csv",
        occultation_target_csv=tmp_path / "occ.csv",
        output_dir=tmp_path,
        tracker_pickle_path=None,
    )

    outputs = run_scheduler(inputs, config)

    assert not outputs.schedule.empty
    assert set(outputs.schedule["Target"].astype(str)) == {"Free Time"}


def test_main_propagates_primary_only_cli_flag(monkeypatch, tmp_path):
    captured: dict[str, PandoraSchedulerConfig] = {}

    def fake_build_schedule(config: PandoraSchedulerConfig) -> SchedulerResult:
        captured["config"] = config
        return SchedulerResult(schedule_csv=None, reports={}, diagnostics={})

    monkeypatch.setattr(scheduler_driver, "build_schedule", fake_build_schedule)
    monkeypatch.setattr(
        scheduler_driver,
        "print_summary",
        lambda result, xml_path: None,
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_scheduler.py",
            "--start",
            "2026-02-05",
            "--end",
            "2026-02-06",
            "--output",
            str(tmp_path),
            "--primary-only",
        ],
    )

    exit_code = scheduler_driver.main()

    assert exit_code == 0
    assert captured["config"].primary_only_mode is True
