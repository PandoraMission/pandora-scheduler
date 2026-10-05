from __future__ import annotations

import types
import warnings
from importlib import util as importlib_util
from pathlib import Path

import pandas as pd
import pytest
from erfa import ErfaWarning

# Import the rework manifest builder (place imports at top to satisfy linters)
from pandorascheduler_rework.targets.manifest import build_target_manifest

warnings.filterwarnings("ignore", category=ErfaWarning)

pytestmark = pytest.mark.filterwarnings(
    r"ignore:ERFA function .*:erfa.ErfaWarning"
)


def _find_target_definition_dir() -> Path:
    """Find target definition files in order of preference."""
    # Option 1: Limited fixture set for CI/quick testing
    limited_dir = (
        Path(__file__).resolve().parents[1]
        / "comparison_outputs"
        / "target_definition_files_limited"
    )
    if limited_dir.is_dir():
        return limited_dir
    
    # Option 2: Full target list (typical development setup)
    full_dir = Path(__file__).resolve().parents[2] / "PandoraTargetList" / "target_definition_files"
    if full_dir.is_dir():
        return full_dir
    
    # Return limited dir path even if it doesn't exist (will be caught by skip check)
    return limited_dir


BASE_DIR = _find_target_definition_dir()


def _load_legacy_helper_codes() -> types.ModuleType:
    module_name = "_legacy_helper_codes_for_manifest_test"
    legacy_path = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "pandorascheduler"
        / "helper_codes.py"
    )
    spec = importlib_util.spec_from_file_location(module_name, legacy_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load legacy helper_codes module for tests")
    module = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


legacy_helper = _load_legacy_helper_codes()


@pytest.mark.parametrize(
    "category",
    [
        "exoplanet",
        "primary-exoplanet",
        "secondary-exoplanet",
        "auxiliary-exoplanet",
        "auxiliary-standard",
        "monitoring-standard",
        "occultation-standard",
    ],
)
def test_manifest_matches_legacy(category: str, monkeypatch: pytest.MonkeyPatch) -> None:
    if not BASE_DIR.is_dir():
        pytest.skip("Target definition fixtures are not available")

    monkeypatch.setattr(legacy_helper, "TARGET_DEF_BASE", str(BASE_DIR))
    monkeypatch.setattr(legacy_helper, "tqdm", lambda iterable, **_: iterable)

    legacy_df = legacy_helper.process_target_files(category)
    rework_df = build_target_manifest(category, BASE_DIR)

    assert isinstance(legacy_df, pd.DataFrame)
    assert isinstance(rework_df, pd.DataFrame)

    # For occultation-standard, rework intentionally adds "Number of Hours Requested"
    # which was missing in legacy code (bug fix - the priority table has hours_req)
    if category == "occultation-standard":
        expected_extra = {"Number of Hours Requested"}
        assert set(rework_df.columns) - set(legacy_df.columns) == expected_extra
        # Compare only common columns
        common_cols = list(legacy_df.columns)
        rework_df = rework_df[common_cols]
    else:
        assert set(legacy_df.columns) == set(rework_df.columns)
        rework_df = rework_df[legacy_df.columns]

    sort_key = "Original Filename" if "Original Filename" in legacy_df.columns else "Star Name"
    legacy_sorted = legacy_df.sort_values(by=sort_key).reset_index(drop=True)
    rework_sorted = rework_df.sort_values(by=sort_key).reset_index(drop=True)

    pd.testing.assert_frame_equal(rework_sorted, legacy_sorted)


class TestManifestStrictValidation:
    """Tests for strict validation of required fields in manifests."""

    def test_missing_hours_req_column_raises_error(self):
        """Test that missing 'hours_req' column in priority table raises an error."""
        from pandorascheduler_rework.targets.manifest import (
            TargetDefinitionError,
            _apply_priority,
        )

        category = "auxiliary-standard"
        
        # Create priority table WITHOUT hours_req column
        priority_df = pd.DataFrame({
            "rank": [1],
            "target": ["test_star"],
            "priority": [0.9],
            # Intentionally missing: "hours_req"
        })
        
        row = {"Original Filename": "test_star", "Star Name": "TestStar"}

        with pytest.raises(TargetDefinitionError, match="missing required.*hours_req"):
            _apply_priority(row, category, priority_df)

    def test_missing_hours_req_value_raises_error(self):
        """Test that missing 'hours_req' value for a target raises an error."""
        import numpy as np

        from pandorascheduler_rework.targets.manifest import (
            TargetDefinitionError,
            _apply_priority,
        )

        category = "auxiliary-standard"
        
        # Create priority table with NaN hours_req value
        priority_df = pd.DataFrame({
            "rank": [1],
            "target": ["test_star"],
            "priority": [0.9],
            "hours_req": [np.nan],  # Missing value
        })
        
        row = {"Original Filename": "test_star", "Star Name": "TestStar"}

        with pytest.raises(TargetDefinitionError, match="missing.*hours_req.*value"):
            _apply_priority(row, category, priority_df)

    def test_valid_hours_req_works(self):
        """Test that valid 'hours_req' value is read correctly."""
        from pandorascheduler_rework.targets.manifest import _apply_priority

        category = "auxiliary-standard"
        
        # Create priority table with valid hours_req
        priority_df = pd.DataFrame({
            "rank": [1],
            "target": ["test_star"],
            "priority": [0.9],
            "hours_req": [100],
        })
        
        row = {"Original Filename": "test_star", "Star Name": "TestStar"}
        _apply_priority(row, category, priority_df)
        
        assert "Number of Hours Requested" in row
        assert row["Number of Hours Requested"] == 100

