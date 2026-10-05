"""Phase 4 tests: debug_occultation_visit script helpers."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest

# Import helpers from the scripts directory.
import importlib.util, sys

_spec = importlib.util.spec_from_file_location(
    "debug_occultation_visit",
    Path(__file__).resolve().parents[1] / "scripts" / "debug_occultation_visit.py",
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _mod
_spec.loader.exec_module(_mod)

_load_schedule = _mod._load_schedule
_find_occultation_candidates = _mod._find_occultation_candidates
_infer_data_dir = _mod._infer_data_dir
debug_visit = _mod.debug_visit


class TestLoadSchedule:
    def test_loads_csv_and_parses_datetimes(self, tmp_path):
        csv = tmp_path / "sched.csv"
        csv.write_text(
            "Target,Observation Start,Observation Stop,RA,DEC\n"
            "Star-A,2026-03-12 00:00:00,2026-03-12 06:00:00,120.0,-30.0\n"
        )
        df = _load_schedule(csv)
        assert len(df) == 1
        assert df.iloc[0]["Target"] == "Star-A"
        assert isinstance(df.iloc[0]["Observation Start"], pd.Timestamp)


class TestFindOccultationCandidates:
    def test_finds_stars_with_parquet(self, tmp_path):
        data_dir = tmp_path / "data"
        aux = data_dir / "aux_targets"
        for name in ("Star-A", "Star-B"):
            d = aux / name
            d.mkdir(parents=True)
            (d / f"Visibility for {name}.parquet").write_bytes(b"dummy")
        result = _find_occultation_candidates(data_dir)
        assert set(result) == {"Star-A", "Star-B"}

    def test_empty_when_no_aux_dir(self, tmp_path):
        assert _find_occultation_candidates(tmp_path) == []


class TestInferDataDir:
    def test_prefers_legacy_data_dir(self, tmp_path):
        sched = tmp_path / "sched.csv"
        sched.write_text("Target,Observation Start,Observation Stop\n")
        (tmp_path / "data").mkdir()
        (tmp_path / "data_91_25_110").mkdir()
        assert _infer_data_dir(sched) == (tmp_path / "data")

    def test_uses_single_data_star_candidate(self, tmp_path):
        sched = tmp_path / "sched.csv"
        sched.write_text("Target,Observation Start,Observation Stop\n")
        candidate = tmp_path / "data_91_25_110"
        candidate.mkdir()
        assert _infer_data_dir(sched) == candidate

    def test_errors_when_no_candidates(self, tmp_path):
        sched = tmp_path / "sched.csv"
        sched.write_text("Target,Observation Start,Observation Stop\n")
        with pytest.raises(FileNotFoundError, match="Could not infer --data-dir"):
            _infer_data_dir(sched)

    def test_errors_when_multiple_candidates(self, tmp_path):
        sched = tmp_path / "sched.csv"
        sched.write_text("Target,Observation Start,Observation Stop\n")
        (tmp_path / "data_91_25_110").mkdir()
        (tmp_path / "data_100_30_90").mkdir()
        with pytest.raises(ValueError, match="multiple candidate directories"):
            _infer_data_dir(sched)


class TestDebugVisitOutput:
    def test_prints_visit_info(self, tmp_path, capsys):
        """debug_visit prints structured output even with no visibility data."""
        sched = pd.DataFrame(
            {
                "Target": ["Test-Star"],
                "Observation Start": [pd.Timestamp("2026-03-12 00:00:00")],
                "Observation Stop": [pd.Timestamp("2026-03-12 01:00:00")],
                "RA": [120.0],
                "DEC": [-30.0],
            }
        )
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        debug_visit(sched, 0, data_dir)
        out = capsys.readouterr().out
        assert "Visit 0: Test-Star" in out
        assert "Duration: 60.0 min" in out
        assert "No auxiliary target" in out
