"""Phase 3 tests: visibility parquet metadata and generate_visibility three-way logic."""

from __future__ import annotations

import io
from datetime import datetime
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.pipeline import _maybe_generate_visibility, SchedulerPaths
from pandorascheduler_rework.visibility.catalog import (
    _write_visibility_parquet,
    _apply_transit_overlaps,
)


# ---------------------------------------------------------------------------
# _write_visibility_parquet metadata tests
# ---------------------------------------------------------------------------


class TestWriteVisibilityParquet:
    @pytest.fixture()
    def config(self) -> PandoraSchedulerConfig:
        return PandoraSchedulerConfig(
            targets_manifest=Path("/dev/null"),
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            sun_avoidance_deg=93.0,
            moon_avoidance_deg=45.0,
            earth_avoidance_deg=86.0,
        )

    @pytest.fixture()
    def sample_df(self) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "Time(MJD_UTC)": [60000.0, 60000.001],
                "Visible": [True, False],
            }
        )

    def test_file_contains_keepout_metadata(self, tmp_path, config, sample_df):
        out = tmp_path / "test.parquet"
        _write_visibility_parquet(sample_df, out, config)

        meta = pq.read_schema(out).metadata
        assert meta[b"pandora.visibility_sun_deg"] == b"93.0"
        assert meta[b"pandora.visibility_moon_deg"] == b"45.0"
        assert meta[b"pandora.visibility_earth_deg"] == b"86.0"

    def test_bytesio_contains_keepout_metadata(self, config, sample_df):
        buf = io.BytesIO()
        _write_visibility_parquet(sample_df, buf, config)
        buf.seek(0)

        meta = pq.read_schema(buf).metadata
        assert meta[b"pandora.visibility_sun_deg"] == b"93.0"
        assert meta[b"pandora.visibility_moon_deg"] == b"45.0"
        assert meta[b"pandora.visibility_earth_deg"] == b"86.0"

    def test_data_round_trips(self, tmp_path, config, sample_df):
        out = tmp_path / "test.parquet"
        _write_visibility_parquet(sample_df, out, config)

        result = pd.read_parquet(out)
        pd.testing.assert_frame_equal(result, sample_df)

    def test_preserves_existing_pandas_metadata(self, tmp_path, config, sample_df):
        """pandas embeds its own schema metadata; ours should not clobber it."""
        out = tmp_path / "test.parquet"
        _write_visibility_parquet(sample_df, out, config)

        meta = pq.read_schema(out).metadata
        # pandas always writes b"pandas" key
        assert b"pandas" in meta
        # and our keys coexist
        assert b"pandora.visibility_sun_deg" in meta


# ---------------------------------------------------------------------------
# _apply_transit_overlaps with config parameter
# ---------------------------------------------------------------------------


class TestApplyTransitOverlapsConfig:
    def test_with_config_embeds_metadata(self, tmp_path):
        """When config is provided, rewritten parquet files get metadata."""
        config = PandoraSchedulerConfig(
            targets_manifest=Path("/dev/null"),
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            sun_avoidance_deg=93.0,
            moon_avoidance_deg=45.0,
            earth_avoidance_deg=86.0,
        )

        star = "HD-1234"
        p1, p2 = "Planet-A", "Planet-B"
        for planet in (p1, p2):
            d = tmp_path / star / planet
            d.mkdir(parents=True)
            df = pd.DataFrame(
                {
                    "Transit_Start": [60000.0],
                    "Transit_Stop": [60000.01],
                    "Transit_Start_UTC": ["2026-02-05T00:00:00"],
                    "Transit_Stop_UTC": ["2026-02-05T00:15:00"],
                    "Transit_Coverage": [1.0],
                    "SAA_Overlap": [0.0],
                }
            )
            _write_visibility_parquet(df, d / f"Visibility for {planet}.parquet", config)

        _apply_transit_overlaps([(star, p1), (star, p2)], tmp_path, config)

        for planet in (p1, p2):
            path = tmp_path / star / planet / f"Visibility for {planet}.parquet"
            meta = pq.read_schema(path).metadata
            assert meta[b"pandora.visibility_sun_deg"] == b"93.0"

    def test_without_config_still_works(self, tmp_path):
        """Backward compat: config=None falls back to plain to_parquet."""
        star = "HD-1234"
        p1, p2 = "Planet-A", "Planet-B"
        for planet in (p1, p2):
            d = tmp_path / star / planet
            d.mkdir(parents=True)
            df = pd.DataFrame(
                {
                    "Transit_Start": [60000.0],
                    "Transit_Stop": [60000.01],
                    "Transit_Start_UTC": ["2026-02-05T00:00:00"],
                    "Transit_Stop_UTC": ["2026-02-05T00:15:00"],
                    "Transit_Coverage": [1.0],
                    "SAA_Overlap": [0.0],
                }
            )
            df.to_parquet(
                d / f"Visibility for {planet}.parquet",
                index=False,
                engine="pyarrow",
            )

        _apply_transit_overlaps([(star, p1), (star, p2)], tmp_path)

        for planet in (p1, p2):
            path = tmp_path / star / planet / f"Visibility for {planet}.parquet"
            result = pd.read_parquet(path)
            assert "Transit_Overlap" in result.columns


# ---------------------------------------------------------------------------
# generate_visibility three-way logic tests
# ---------------------------------------------------------------------------


class TestGenerateVisibilityThreeWay:
    """Test the three-way generate_visibility precedence in _maybe_generate_visibility."""

    @pytest.fixture()
    def _mock_builder(self, monkeypatch):
        """Patch out the actual visibility builder to just count calls."""
        self.call_count = 0

        def fake_builder(cfg, target_list, partner_list=None, output_subpath="targets"):
            self.call_count += 1

        monkeypatch.setattr(
            "pandorascheduler_rework.pipeline.build_visibility_catalog",
            fake_builder,
        )

    @pytest.fixture()
    def _setup_paths(self, tmp_path):
        pkg = tmp_path / "package"
        data = pkg / "data"
        data.mkdir(parents=True)
        for name in (
            "exoplanet_targets.csv",
            "auxiliary-standard_targets.csv",
            "monitoring-standard_targets.csv",
            "occultation-standard_targets.csv",
        ):
            (data / name).write_text("Star Name,RA,DEC\n")
        self.paths = SchedulerPaths.from_package_root(pkg)
        self.csv_files = {
            "primary": data / "exoplanet_targets.csv",
            "auxiliary": data / "auxiliary-standard_targets.csv",
            "monitoring": data / "monitoring-standard_targets.csv",
            "occultation": data / "occultation-standard_targets.csv",
        }
        self.tmp_path = tmp_path

    def _call(self, config):
        _maybe_generate_visibility(
            config,
            self.paths,
            config.window_start,
            config.window_end,
            self.csv_files["primary"],
            self.csv_files["auxiliary"],
            self.csv_files["monitoring"],
            self.csv_files["occultation"],
        )

    @pytest.mark.usefixtures("_mock_builder", "_setup_paths")
    def test_explicit_true_without_gmat_generates(self):
        """generate_visibility='true' forces generation even without GMAT."""
        config = PandoraSchedulerConfig(
            targets_manifest=self.csv_files["primary"],
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            output_dir=self.tmp_path,
            extra_inputs={"generate_visibility": "true"},
        )
        self._call(config)
        assert self.call_count == 4

    @pytest.mark.usefixtures("_mock_builder", "_setup_paths")
    def test_explicit_false_with_gmat_skips(self):
        """generate_visibility='false' disables generation even with GMAT present."""
        gmat = self.tmp_path / "gmat.txt"
        gmat.write_text("# dummy\n")
        config = PandoraSchedulerConfig(
            targets_manifest=self.csv_files["primary"],
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            output_dir=self.tmp_path,
            gmat_ephemeris=gmat,
            extra_inputs={"generate_visibility": "false"},
        )
        self._call(config)
        assert self.call_count == 0

    @pytest.mark.usefixtures("_mock_builder", "_setup_paths")
    def test_unset_with_gmat_generates(self):
        """When generate_visibility is unset, presence of GMAT triggers generation."""
        gmat = self.tmp_path / "gmat.txt"
        gmat.write_text("# dummy\n")
        config = PandoraSchedulerConfig(
            targets_manifest=self.csv_files["primary"],
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            output_dir=self.tmp_path,
            gmat_ephemeris=gmat,
        )
        self._call(config)
        assert self.call_count == 4

    @pytest.mark.usefixtures("_mock_builder", "_setup_paths")
    def test_unset_without_gmat_skips(self):
        """When generate_visibility is unset and no GMAT, skip generation."""
        config = PandoraSchedulerConfig(
            targets_manifest=self.csv_files["primary"],
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            output_dir=self.tmp_path,
        )
        self._call(config)
        assert self.call_count == 0

    @pytest.mark.usefixtures("_mock_builder", "_setup_paths")
    def test_explicit_no_skips(self):
        """generate_visibility='no' is treated as explicit false."""
        gmat = self.tmp_path / "gmat.txt"
        gmat.write_text("# dummy\n")
        config = PandoraSchedulerConfig(
            targets_manifest=self.csv_files["primary"],
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            output_dir=self.tmp_path,
            gmat_ephemeris=gmat,
            extra_inputs={"generate_visibility": "no"},
        )
        self._call(config)
        assert self.call_count == 0
