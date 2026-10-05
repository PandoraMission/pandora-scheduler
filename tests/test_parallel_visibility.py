"""Tests that parallel visibility generation produces identical output to serial."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.catalog import build_visibility_catalog


def _to_gmat_mod_julian(times: Time) -> np.ndarray:
    mjd = np.asarray(times.to_value("mjd"), dtype=float)
    return mjd - 29999.5


def _make_gmat_csv(tmp_path: Path, window_start: datetime, window_end: datetime) -> Path:
    """Build a synthetic GMAT ephemeris CSV spanning the given window."""
    gmat_samples = Time(
        [
            window_start - timedelta(minutes=10),
            window_start,
            window_start + timedelta(hours=1.5),
            window_start + timedelta(hours=3),
            window_end,
            window_end + timedelta(minutes=10),
        ],
        format="datetime",
        scale="utc",
    )
    n = gmat_samples.size

    # Put spacecraft at origin, bodies at various fixed positions
    # so that Sun, Moon, Earth are at well-separated angles from boresight.
    gmat_df = pd.DataFrame(
        {
            "Earth.UTCModJulian": _to_gmat_mod_julian(gmat_samples),
            "Earth.EarthMJ2000Eq.X": np.full(n, -7000.0),
            "Earth.EarthMJ2000Eq.Y": np.zeros(n),
            "Earth.EarthMJ2000Eq.Z": np.zeros(n),
            "Pandora.EarthMJ2000Eq.X": np.zeros(n),
            "Pandora.EarthMJ2000Eq.Y": np.zeros(n),
            "Pandora.EarthMJ2000Eq.Z": np.zeros(n),
            "Sun.EarthMJ2000Eq.X": np.zeros(n),
            "Sun.EarthMJ2000Eq.Y": np.full(n, 7000.0),
            "Sun.EarthMJ2000Eq.Z": np.zeros(n),
            "Luna.EarthMJ2000Eq.X": np.zeros(n),
            "Luna.EarthMJ2000Eq.Y": np.zeros(n),
            "Luna.EarthMJ2000Eq.Z": np.full(n, 7000.0),
            "Pandora.Earth.Latitude": np.full(n, -20.0),
            "Pandora.Earth.Longitude": np.full(n, -50.0),
        }
    )
    path = tmp_path / "gmat.csv"
    gmat_df.to_csv(path, index=False)
    return path


def _make_target_manifest(tmp_path: Path, star_names: list[str]) -> Path:
    """Make a CSV manifest with *star_names* at distinct RA/DEC."""
    rows = []
    period_days = 30.0 / (24.0 * 60.0)
    epoch_time = Time(datetime(2025, 2, 4), scale="tdb", format="datetime")
    epoch_bjd_tdb = float(epoch_time.jd) - 2400000.5
    for i, name in enumerate(star_names):
        rows.append(
            {
                "Star Name": name,
                "Star Simbad Name": name,
                "Planet Name": f"{name}b",
                "Planet Simbad Name": f"{name}b",
                "RA": 20.0 * i,  # distinct RA so each star has unique visibility
                "DEC": 10.0 + 5.0 * i,
                "Transit Duration (hrs)": 1.0,
                "Period (days)": period_days,
                "Transit Epoch (BJD_TDB-2400000.5)": epoch_bjd_tdb,
            }
        )
    df = pd.DataFrame(rows)
    path = tmp_path / "exoplanet_targets.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture
def window():
    start = datetime(2025, 2, 5)
    end = start + timedelta(hours=6)
    return start, end


def _run_catalog(
    tmp_path: Path,
    window: tuple[datetime, datetime],
    parallel_workers: int,
) -> dict[str, pd.DataFrame]:
    """Run build_visibility_catalog and return {star_name: DataFrame}."""
    start, end = window
    out_dir = tmp_path / f"out_{parallel_workers}"
    gmat_path = _make_gmat_csv(tmp_path, start, end)
    manifest = _make_target_manifest(tmp_path, ["StarA", "StarB", "StarC"])

    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=end,
        gmat_ephemeris=gmat_path,
        targets_manifest=manifest,
        output_dir=out_dir,
        extra_inputs={"data_subdir": "data"},
        force_regenerate=True,
        sun_avoidance_deg=45.0,
        moon_avoidance_deg=30.0,
        earth_avoidance_deg=20.0,
        parallel_workers=parallel_workers,
    )
    build_visibility_catalog(
        config,
        target_list=manifest,
        output_subpath="targets",
    )
    result = {}
    data_root = out_dir / "data" / "targets"
    for name in ["StarA", "StarB", "StarC"]:
        parquet = data_root / name / f"Visibility for {name}.parquet"
        assert parquet.exists(), f"Missing {parquet}"
        result[name] = pd.read_parquet(parquet)
    return result


def test_parallel_matches_serial(tmp_path, window):
    """Parallel (2 workers) must produce bit-identical parquet output to serial."""
    serial = _run_catalog(tmp_path, window, parallel_workers=1)
    parallel = _run_catalog(tmp_path, window, parallel_workers=2)

    for star in serial:
        pd.testing.assert_frame_equal(
            serial[star],
            parallel[star],
            check_exact=True,
            obj=f"Visibility for {star}",
        )


def test_parallel_workers_config_default():
    """parallel_workers defaults to 0 (auto)."""
    cfg = PandoraSchedulerConfig(
        window_start=datetime(2025, 1, 1),
        window_end=datetime(2025, 1, 2),
    )
    assert cfg.parallel_workers == 0


def test_parallel_workers_config_set():
    """parallel_workers can be set explicitly."""
    cfg = PandoraSchedulerConfig(
        window_start=datetime(2025, 1, 1),
        window_end=datetime(2025, 1, 2),
        parallel_workers=4,
    )
    assert cfg.parallel_workers == 4


def test_parallel_workers_config_rejects_negative():
    """parallel_workers must be non-negative."""
    with pytest.raises(ValueError, match="parallel_workers"):
        PandoraSchedulerConfig(
            window_start=datetime(2025, 1, 1),
            window_end=datetime(2025, 1, 2),
            parallel_workers=-1,
        )
