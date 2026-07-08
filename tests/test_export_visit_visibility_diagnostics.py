from __future__ import annotations

import argparse
import importlib.util
from datetime import datetime
from pathlib import Path

import numpy as np


def _load_script_module():
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "export_visit_visibility_diagnostics.py"
    )
    spec = importlib.util.spec_from_file_location(
        "export_visit_visibility_diagnostics_test_module",
        script_path,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_resolve_visit_window_from_provenance_uses_visit_and_sequence_ids(tmp_path):
    module = _load_script_module()
    provenance_csv = tmp_path / "Pandora_science_calendar_sequence_provenance.csv"
    provenance_csv.write_text(
        "visit_id,sequence_id,target,start_utc,stop_utc\n"
        "12,3,WASP-18b,2026-07-09T01:03:59Z,2026-07-09T03:16:35Z\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        visit_id="12",
        sequence_id="3",
        provenance_csv=provenance_csv,
        output=None,
        config=tmp_path / "scheduler_config.json",
        run_dir=None,
    )

    target_name, start, stop = module._resolve_visit_window_from_provenance(args)

    assert target_name == "WASP-18b"
    assert start.strftime("%Y-%m-%d %H:%M:%S") == "2026-07-09 01:03:59"
    assert stop.strftime("%Y-%m-%d %H:%M:%S") == "2026-07-09 03:16:35"


def test_resolve_visit_window_from_provenance_accepts_zero_padded_ids(tmp_path):
    module = _load_script_module()
    provenance_csv = tmp_path / "Pandora_science_calendar_sequence_provenance.csv"
    provenance_csv.write_text(
        "visit_id,sequence_id,target,start_utc,stop_utc\n"
        "0001,003,WASP-18b,2026-07-13T00:29:00Z,2026-07-13T00:45:00Z\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        visit_id="1",
        sequence_id="3",
        provenance_csv=provenance_csv,
        output=None,
        config=tmp_path / "scheduler_config.json",
        run_dir=None,
    )

    target_name, start, stop = module._resolve_visit_window_from_provenance(args)

    assert target_name == "WASP-18b"
    assert start.strftime("%Y-%m-%d %H:%M:%S") == "2026-07-13 00:29:00"
    assert stop.strftime("%Y-%m-%d %H:%M:%S") == "2026-07-13 00:45:00"


def test_default_output_uses_visit_and_sequence_ids(tmp_path):
    module = _load_script_module()
    args = argparse.Namespace(
        visit_id="12",
        sequence_id="3",
        config=tmp_path / "scheduler_config.json",
        run_dir=None,
        target_name=None,
        start=None,
    )

    result = module._default_output(args)

    assert result == tmp_path / "visit12_seq3_visibility.csv"


def test_default_output_prefers_run_dir_for_visit_sequence_mode(tmp_path):
    module = _load_script_module()
    run_dir = tmp_path / "output_20260622_20260629_EarthDay111"
    args = argparse.Namespace(
        visit_id="12",
        sequence_id="3",
        config=tmp_path / "scheduler_config.json",
        run_dir=run_dir,
        target_name=None,
        start=None,
    )

    result = module._default_output(args)

    assert result == run_dir / "visit12_seq3_visibility.csv"


def test_default_provenance_prefers_run_dir(tmp_path):
    module = _load_script_module()
    run_dir = tmp_path / "output_20260622_20260629_EarthDay111"
    args = argparse.Namespace(
        run_dir=run_dir,
        output=None,
        config=tmp_path / "scheduler_config.json",
    )

    result = module._default_provenance_csv(args)

    assert result == run_dir / "Pandora_science_calendar_sequence_provenance.csv"


def test_lookup_radec_from_run_manifests_finds_non_exoplanet_target(tmp_path):
    module = _load_script_module()
    run_dir = tmp_path / "output_20260713_20260720_EarthDay111_test"
    data_dir = run_dir / "data_91_20_sc111_86_occ121_96"
    data_dir.mkdir(parents=True)
    (data_dir / "occultation-standard_targets.csv").write_text(
        "Star Name,RA,DEC\n"
        "G1203023711162117120,123.45,-54.321\n",
        encoding="utf-8",
    )

    result = module._lookup_radec_from_run_manifests(
        run_dir,
        "G1203023711162117120",
    )

    assert result == (123.45, -54.321)


def test_build_config_for_target_uses_science_keepouts_for_exoplanets(tmp_path):
    module = _load_script_module()
    run_dir = tmp_path / "output"
    data_dir = run_dir / "data_91_20_sc111_86_occ121_96"
    data_dir.mkdir(parents=True)
    (data_dir / "exoplanet_targets.csv").write_text(
        "Planet Name,Star Name,RA,DEC\n"
        "WASP-18b,WASP-18,24.3,-45.6\n",
        encoding="utf-8",
    )
    cfg = {
        "earth_keepouts": "different",
        "earth_avoidance_day_deg_science": 111.0,
        "earth_avoidance_night_deg_science": 86.0,
        "earth_avoidance_day_deg_occultation": 121.0,
        "earth_avoidance_night_deg_occultation": 96.0,
    }

    result = module._build_config_for_target(
        cfg,
        start=datetime(2026, 7, 13, 0, 0, 0),
        stop=datetime(2026, 7, 13, 1, 0, 0),
        target_name="WASP-18b",
        run_dir=run_dir,
    )

    assert result.earth_avoidance_day_deg == 111.0
    assert result.earth_avoidance_night_deg == 86.0


def test_build_config_for_target_uses_occultation_keepouts_for_non_exoplanets(tmp_path):
    module = _load_script_module()
    run_dir = tmp_path / "output"
    data_dir = run_dir / "data_91_20_sc111_86_occ121_96"
    data_dir.mkdir(parents=True)
    (data_dir / "occultation-standard_targets.csv").write_text(
        "Star Name,RA,DEC\n"
        "G1203023711162117120,123.45,-54.321\n",
        encoding="utf-8",
    )
    cfg = {
        "earth_keepouts": "different",
        "earth_avoidance_day_deg_science": 111.0,
        "earth_avoidance_night_deg_science": 86.0,
        "earth_avoidance_day_deg_occultation": 121.0,
        "earth_avoidance_night_deg_occultation": 96.0,
    }

    result = module._build_config_for_target(
        cfg,
        start=datetime(2026, 7, 13, 0, 0, 0),
        stop=datetime(2026, 7, 13, 1, 0, 0),
        target_name="G1203023711162117120",
        run_dir=run_dir,
    )

    assert result.earth_avoidance_day_deg == 121.0
    assert result.earth_avoidance_night_deg == 96.0


def test_earth_threshold_branch_labels_for_exoplanet_different_mode():
    module = _load_script_module()

    result = module._earth_threshold_branch_labels(
        target_category="exoplanet",
        earth_keepouts="different",
        sunlit_subsatellite=np.array([True, False]),
    )

    assert result.tolist() == ["science_day", "science_night"]


def test_earth_threshold_branch_labels_for_occultation_different_mode():
    module = _load_script_module()

    result = module._earth_threshold_branch_labels(
        target_category="occultation-standard",
        earth_keepouts="different",
        sunlit_subsatellite=np.array([True, False]),
    )

    assert result.tolist() == ["occultation_day", "occultation_night"]
