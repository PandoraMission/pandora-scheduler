from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.catalog import build_visibility_catalog
from pandorascheduler_rework.visibility.geometry import (
    build_minute_cadence,
    compute_saa_crossings,
    interpolate_gmat_ephemeris,
)

MJD_EPOCH = datetime(1858, 11, 17)


def mjd_to_datetime(mjd: float) -> datetime:
    return MJD_EPOCH + timedelta(days=mjd)


def _to_gmat_mod_julian(times: Time) -> np.ndarray:
    mjd = np.asarray(times.to_value("mjd"), dtype=float)
    return mjd - 29999.5


@pytest.fixture
def mock_config(tmp_path):
    return PandoraSchedulerConfig(
        window_start=datetime(2025, 2, 5),
        window_end=datetime(2025, 2, 5, 6, 0, 0),
        gmat_ephemeris=tmp_path / "gmat.csv",
        targets_manifest=tmp_path / "targets.csv",
        output_dir=tmp_path / "output",
        extra_inputs={"data_subdir": "data"},
        force_regenerate=True,
        sun_avoidance_deg=45.0,
        moon_avoidance_deg=30.0,
        earth_avoidance_deg=20.0,
    )


def test_build_visibility_catalog_generates_star_and_planet_outputs(tmp_path, mock_config):
    window_start = datetime(2025, 2, 5, 0, 0, 0)
    window_end = window_start + timedelta(hours=6)

    cadence = build_minute_cadence(window_start, window_end)

    gmat_samples = Time(
        [
            window_start - timedelta(minutes=10),
            window_start,
            window_start + timedelta(hours=3),
            window_end,
            window_end + timedelta(minutes=10),
        ],
        format="datetime",
        scale="utc",
    )

    gmat_df = pd.DataFrame(
        {
            "Earth.UTCModJulian": _to_gmat_mod_julian(gmat_samples),
            "Earth.EarthMJ2000Eq.X": np.full(gmat_samples.size, -7000.0),
            "Earth.EarthMJ2000Eq.Y": np.zeros(gmat_samples.size),
            "Earth.EarthMJ2000Eq.Z": np.zeros(gmat_samples.size),
            "Pandora.EarthMJ2000Eq.X": np.zeros(gmat_samples.size),
            "Pandora.EarthMJ2000Eq.Y": np.zeros(gmat_samples.size),
            "Pandora.EarthMJ2000Eq.Z": np.zeros(gmat_samples.size),
            "Sun.EarthMJ2000Eq.X": np.zeros(gmat_samples.size),
            "Sun.EarthMJ2000Eq.Y": np.full(gmat_samples.size, 7000.0),
            "Sun.EarthMJ2000Eq.Z": np.zeros(gmat_samples.size),
            "Luna.EarthMJ2000Eq.X": np.zeros(gmat_samples.size),
            "Luna.EarthMJ2000Eq.Y": np.zeros(gmat_samples.size),
            "Luna.EarthMJ2000Eq.Z": np.full(gmat_samples.size, 7000.0),
            "Pandora.Earth.Latitude": np.full(gmat_samples.size, -20.0),
            "Pandora.Earth.Longitude": np.full(gmat_samples.size, -50.0),
        }
    )

    gmat_path = tmp_path / "gmat.csv"
    gmat_df.to_csv(gmat_path, index=False)

    period_minutes = 30.0
    period_days = period_minutes / (24.0 * 60.0)
    epoch_time = Time(window_start - timedelta(minutes=150), scale="tdb", format="datetime")
    epoch_bjd_tdb = float(epoch_time.jd)

    target_df = pd.DataFrame(
        {
            "Star Name": ["TestStar"],
            "Star Simbad Name": ["TestStar"],
            "Planet Name": ["TestPlanet"],
            "Planet Simbad Name": ["TestPlanet"],
            "RA": [0.0],
            "DEC": [0.0],
            "Transit Duration (hrs)": [1.0],
            "Period (days)": [period_days],
            "Transit Epoch (BJD_TDB-2400000.5)": [epoch_bjd_tdb - 2400000.5],
        }
    )
    target_path = mock_config.targets_manifest
    target_df.to_csv(target_path, index=False)

    build_visibility_catalog(
        mock_config,
        target_list=mock_config.targets_manifest,
        output_subpath="targets"
    )

    output_root = mock_config.output_dir / "data" / "targets"
    star_output = output_root / "TestStar" / "Visibility for TestStar.parquet"
    planet_output = output_root / "TestStar" / "TestPlanet" / "Visibility for TestPlanet.parquet"

    assert star_output.exists()
    assert planet_output.exists()

    # Verify star visibility Parquet was created with correct structure
    star_df = pd.read_parquet(star_output)
    assert not star_df.empty
    expected_cols = {"Time(MJD_UTC)", "Time_UTC", "SAA_Crossing", "Visible", "Earth_Sep", "Earth_Threshold", "Solar_Power_Frac", "Moon_Sep", "Sun_Sep", "Roll_Deg", "N_ST_Pass"}
    assert set(star_df.columns) == expected_cols

    star_visibility = pd.read_parquet(star_output)
    planet_visibility = pd.read_parquet(planet_output)

    ephemeris = interpolate_gmat_ephemeris(gmat_path, cadence)

    star_unit_vector = np.array([1.0, 0.0, 0.0])

    def _separation_deg(vectors: np.ndarray) -> np.ndarray:
        norms = np.linalg.norm(vectors, axis=1)
        dot = vectors @ star_unit_vector
        cos_theta = np.clip(dot / norms, -1.0, 1.0)
        return np.round(np.degrees(np.arccos(cos_theta)), 3)

    earth_sep = _separation_deg(ephemeris.earth_pc)
    sun_sep = _separation_deg(ephemeris.sun_pc)
    moon_sep = _separation_deg(ephemeris.moon_pc)

    visible = (sun_sep > 45.0) & (moon_sep > 30.0) & (earth_sep > 20.0)
    visible = np.round(visible.astype(float), 1)

    saa_crossing = np.round(
        compute_saa_crossings(
            ephemeris.spacecraft_lat_deg, ephemeris.spacecraft_lon_deg
        ),
        1,
    )

    expected_star = pd.DataFrame(
        {
            "Time(MJD_UTC)": np.round(cadence.mjd_utc, 6),
            "Time_UTC": Time(cadence.mjd_utc, format="mjd", scale="utc").to_datetime(),
            "SAA_Crossing": saa_crossing,
            "Visible": visible,
            "Earth_Sep": earth_sep,
            "Moon_Sep": moon_sep,
            "Sun_Sep": sun_sep,
        }
    )

    # Drop diagnostic columns that are not part of the legacy geometry check;
    # datetime conversion may also have microsecond differences.
    drop_cols = ["Time_UTC", "Earth_Threshold", "Solar_Power_Frac", "Roll_Deg", "N_ST_Pass"]
    star_visibility_compare = star_visibility.drop(columns=[c for c in drop_cols if c in star_visibility.columns])
    expected_star_compare = expected_star.drop(columns=[c for c in drop_cols if c in expected_star.columns])
    pd.testing.assert_frame_equal(star_visibility_compare, expected_star_compare)

    assert not planet_visibility.empty
    np.testing.assert_array_equal(
        planet_visibility["Transits"].to_numpy(),
        np.arange(len(planet_visibility)),
    )

    star_time_mjd = star_visibility["Time(MJD_UTC)"].to_numpy(dtype=float)
    star_visible = star_visibility["Visible"].to_numpy()
    star_saa = star_visibility["SAA_Crossing"].to_numpy()

    star_time = Time(star_time_mjd, format="mjd", scale="utc")
    star_time_iso = Time(star_time.iso, format="iso", scale="utc")
    star_datetimes = star_time_iso.to_value("datetime")

    visible_times = {
        dt for dt, flag in zip(star_datetimes, star_visible) if flag == 1.0
    }
    saa_times = {
        dt for dt, flag in zip(star_datetimes, star_saa) if flag == 1.0
    }

    assert "SAA_Overlap" in planet_visibility.columns
    assert "Transit_Overlap" not in planet_visibility.columns

    for _, row in planet_visibility.iterrows():
        start_mjd = float(row["Transit_Start"])
        end_mjd = float(row["Transit_Stop"])
        start_dt = mjd_to_datetime(start_mjd).replace(second=0, microsecond=0)
        end_dt = mjd_to_datetime(end_mjd).replace(second=0, microsecond=0)
        # Use half-open [start, end) range to match catalog.py denominator fix
        n_minutes = int((end_dt - start_dt).total_seconds() // 60)
        minute_range = pd.date_range(start_dt, periods=n_minutes, freq="min").to_pydatetime()
        if len(minute_range) == 0:
            assert row["Transit_Coverage"] == 0.0
            continue
        overlap = visible_times.intersection(minute_range)
        expected_cov = len(overlap) / len(minute_range)
        assert np.isclose(row["Transit_Coverage"], expected_cov)

        saa_overlap = saa_times.intersection(minute_range)
        expected_saa = len(saa_overlap) / len(minute_range)
        assert np.isclose(row["SAA_Overlap"], expected_saa)


def test_resolve_star_coord_uses_catalog_coordinates():
    """Test that _resolve_star_coord uses catalog coordinates (no Simbad)."""
    from pandorascheduler_rework.visibility import catalog

    row = pd.Series(
        {
            "Star Name": "Foo",
            "RA": 123.4,
            "DEC": 56.7,
        }
    )
    metadata = {"Foo": (123.4, 56.7)}

    coord = catalog._resolve_star_coord(row, metadata)

    assert np.isclose(coord.ra.deg, 123.4)
    assert np.isclose(coord.dec.deg, 56.7)



def test_resolve_star_coord_uses_metadata_fallback():
    """Test that _resolve_star_coord uses metadata as fallback."""
    from pandorascheduler_rework.visibility import catalog

    row = pd.Series(
        {
            "Star Name": "Bar",
            "RA": np.nan,  # Missing in row
            "DEC": np.nan,  # Missing in row
        }
    )
    metadata = {"Bar": (10.0, -11.0)}

    coord = catalog._resolve_star_coord(row, metadata)

    assert coord.ra.deg == 10.0
    assert coord.dec.deg == -11.0



def test_resolve_star_coord_raises_error_when_missing():
    """Test that _resolve_star_coord raises error when coordinates missing."""
    import pytest

    from pandorascheduler_rework.visibility import catalog

    row = pd.Series(
        {
            "Star Name": "Baz",
            "RA": np.nan,
            "DEC": np.nan,
        }
    )
    metadata = {}  # No fallback

    with pytest.raises(RuntimeError, match="No coordinates found in catalog for Baz"):
        catalog._resolve_star_coord(row, metadata)


def _write_csv(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)

