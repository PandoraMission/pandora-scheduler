"""Tests for pipeline.py helper functions (coverage gap)."""

from datetime import datetime
from pathlib import Path

import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.pipeline import (
    _as_bool,
    _coerce_optional_path,
    _coerce_path,
    _resolve_target_definition_files,
    _stage_too_list,
    _target_definition_from_csv,
    _visibility_config_for_target_definition,
)


class TestCoercePath:
    def test_none_returns_default(self):
        default = Path("/fallback")
        assert _coerce_path(None, default) == default

    def test_string_converted(self):
        result = _coerce_path("some/dir", Path("/default"))
        assert isinstance(result, Path)
        assert result.is_absolute()

    def test_path_passthrough(self):
        result = _coerce_path(Path("/absolute/path"), Path("/default"))
        assert result == Path("/absolute/path")


class TestCoerceOptionalPath:
    def test_none_returns_none(self):
        assert _coerce_optional_path(None) is None

    def test_string_converted(self):
        result = _coerce_optional_path("some/path")
        assert isinstance(result, Path)
        assert result.is_absolute()


class TestResolveTargetDefinitionFiles:
    def test_none_uses_fallback(self):
        result = _resolve_target_definition_files(None, ["exoplanet", "auxiliary-standard"])
        assert result == ["exoplanet", "auxiliary-standard"]

    def test_list_passthrough(self):
        result = _resolve_target_definition_files(["a", "b"], [])
        assert result == ["a", "b"]

    def test_tuple_converted(self):
        result = _resolve_target_definition_files(("x", "y"), [])
        assert result == ["x", "y"]

    def test_comma_separated_string(self):
        result = _resolve_target_definition_files("exoplanet, auxiliary-standard", [])
        assert result == ["exoplanet", "auxiliary-standard"]

    def test_single_string(self):
        result = _resolve_target_definition_files("exoplanet", [])
        assert result == ["exoplanet"]

    def test_invalid_type_raises(self):
        with pytest.raises(TypeError):
            _resolve_target_definition_files(42, [])


class TestTargetDefinitionFromCsv:
    def test_strips_targets_suffix(self):
        assert _target_definition_from_csv(Path("exoplanet_targets.csv")) == "exoplanet"

    def test_no_suffix(self):
        assert _target_definition_from_csv(Path("custom.csv")) == "custom"

    def test_nested_path(self):
        assert _target_definition_from_csv(Path("/data/output/auxiliary-standard_targets.csv")) == "auxiliary-standard"


class TestAsBool:
    def test_none_returns_default(self):
        assert _as_bool(None, True) is True
        assert _as_bool(None, False) is False

    def test_bool_passthrough(self):
        assert _as_bool(True, False) is True
        assert _as_bool(False, True) is False

    def test_string_true_values(self):
        for val in ("1", "true", "True", "yes", "y", "on", "YES"):
            assert _as_bool(val, False) is True

    def test_string_false_values(self):
        for val in ("0", "false", "no", "off"):
            assert _as_bool(val, True) is False

    def test_empty_string_returns_default(self):
        assert _as_bool("", True) is True
        assert _as_bool("  ", False) is False

    def test_int_values(self):
        assert _as_bool(1, False) is True
        assert _as_bool(0, True) is False

    def test_float_values(self):
        assert _as_bool(1.0, False) is True
        assert _as_bool(0.0, True) is False

    def test_other_type_returns_default(self):
        assert _as_bool([], True) is True


class TestStageTooList:
    def test_copies_configured_too_list_to_output_root(self, tmp_path):
        source = tmp_path / "source_too.csv"
        source.write_text(
            "Target,Obs Window Start,Obs Window Stop\n"
            "TOI-674b,2026-06-04 08:25:29,2026-06-04 09:32:07\n",
            encoding="utf-8",
        )

        output_dir = tmp_path / "output"
        out_data = output_dir / "data_91_20_111"
        out_data.mkdir(parents=True)

        resolved = _stage_too_list(
            {"too_list_csv": str(source)},
            output_dir,
            out_data,
        )

        assert resolved == output_dir / "ToO_list.csv"
        assert resolved.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")
        assert not (out_data / "ToO_list.csv").exists()

    def test_prefers_existing_output_root_too_list_when_unconfigured(self, tmp_path):
        output_dir = tmp_path / "output"
        out_data = output_dir / "data_91_20_111"
        out_data.mkdir(parents=True)

        root_too = output_dir / "ToO_list.csv"
        root_too.write_text(
            "Target,Obs Window Start,Obs Window Stop\n"
            "HATS-72b,2026-06-03 03:11:15,2026-06-03 06:16:20\n",
            encoding="utf-8",
        )
        (out_data / "ToO_list.csv").write_text(
            "Target,Obs Window Start,Obs Window Stop\n"
            "StaleTarget,2026-06-01 00:00:00,2026-06-01 01:00:00\n",
            encoding="utf-8",
        )

        resolved = _stage_too_list({}, output_dir, out_data)

        assert resolved == root_too


class TestVisibilityConfigForTargetDefinition:
    def test_same_mode_returns_original_config(self):
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            earth_keepouts="same",
            earth_avoidance_day_deg=111.0,
            earth_avoidance_night_deg=86.0,
        )

        result = _visibility_config_for_target_definition(config, "exoplanet")

        assert result is config

    def test_same_mode_uses_science_keepouts_without_legacy_fields(self):
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            earth_keepouts="same",
            earth_avoidance_day_deg_science=111.0,
            earth_avoidance_night_deg_science=86.0,
        )

        result = _visibility_config_for_target_definition(config, "exoplanet")

        assert result is config
        assert result.earth_avoidance_day_deg == 111.0
        assert result.earth_avoidance_night_deg == 86.0

    def test_different_mode_uses_science_keepouts_for_exoplanets(self):
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            earth_keepouts="different",
            earth_avoidance_day_deg_science=112.0,
            earth_avoidance_night_deg_science=87.0,
            earth_avoidance_day_deg_occultation=121.0,
            earth_avoidance_night_deg_occultation=96.0,
        )

        result = _visibility_config_for_target_definition(config, "exoplanet")

        assert result is config
        assert result.earth_avoidance_day_deg == 112.0
        assert result.earth_avoidance_night_deg == 87.0

    def test_different_mode_uses_occultation_keepouts_for_non_exoplanets(self):
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 1, 2),
            earth_keepouts="different",
            earth_avoidance_day_deg_science=111.0,
            earth_avoidance_night_deg_science=86.0,
            earth_avoidance_day_deg_occultation=121.0,
            earth_avoidance_night_deg_occultation=96.0,
        )

        result = _visibility_config_for_target_definition(
            config, "occultation-standard"
        )

        assert result is not config
        assert result.earth_avoidance_day_deg == 121.0
        assert result.earth_avoidance_night_deg == 96.0
