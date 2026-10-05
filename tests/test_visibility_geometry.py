from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework.visibility.geometry import (
    build_minute_cadence,
    compute_saa_crossings,
    interpolate_gmat_ephemeris,
)


def test_build_minute_cadence_matches_legacy_conversion():
    start = datetime(2025, 2, 5, 0, 0, 0)
    end = start + timedelta(minutes=2)

    cadence = build_minute_cadence(start, end)

    dt_iso = pd.date_range(start, end, freq="min")
    expected_jd = Time(dt_iso.to_julian_date(), format="jd", scale="utc")
    expected_mjd = Time(expected_jd.value - 2400000.5, format="mjd", scale="utc").value

    assert np.allclose(cadence.timestamps.jd, expected_jd.jd)
    assert np.allclose(cadence.mjd_utc, expected_mjd)



def _expected_longitudes(source_times, source_lon, target_times):
    lon_cont = np.copy(source_lon)
    for idx in range(lon_cont.size - 1):
        if abs(lon_cont[idx + 1] - lon_cont[idx]) > 100:
            lon_cont[idx + 1 :] = lon_cont[idx + 1 :] - 360

    interpolated = np.interp(target_times, source_times, lon_cont)
    adjusted = np.copy(interpolated)
    for idx in range(adjusted.size):
        while adjusted[0] < -180:
            adjusted += 360
        if adjusted[idx] < -180:
            adjusted[idx:] += 360
    return adjusted


def test_interpolate_gmat_ephemeris_respects_minute_grid(tmp_path):
    start = datetime(2025, 5, 1, 0, 0, 0)
    end = start + timedelta(minutes=3)
    cadence = build_minute_cadence(start, end)

    base_mod = cadence.timestamps.jd[0] - 2430000.0
    raw_times = base_mod + np.array(
        [-0.0100, -0.0010, -0.0005, 0.0, 0.0005, 0.0010, 0.0015, 0.0100]
    )
    values = np.linspace(0.0, 7.0, raw_times.size)

    gmat_df = pd.DataFrame(
        {
            "Earth.UTCModJulian": raw_times,
            "Earth.EarthMJ2000Eq.X": values,
            "Earth.EarthMJ2000Eq.Y": values + 1,
            "Earth.EarthMJ2000Eq.Z": values + 2,
            "Pandora.EarthMJ2000Eq.X": values + 10,
            "Pandora.EarthMJ2000Eq.Y": values + 11,
            "Pandora.EarthMJ2000Eq.Z": values + 12,
            "Sun.EarthMJ2000Eq.X": values + 20,
            "Sun.EarthMJ2000Eq.Y": values + 21,
            "Sun.EarthMJ2000Eq.Z": values + 22,
            "Luna.EarthMJ2000Eq.X": values + 30,
            "Luna.EarthMJ2000Eq.Y": values + 31,
            "Luna.EarthMJ2000Eq.Z": values + 32,
            "Pandora.Earth.Latitude": np.linspace(-50, 10, raw_times.size),
            "Pandora.Earth.Longitude": np.array(
                [140.0, 150.0, 170.0, 175.0, 179.0, -179.0, -175.0, -160.0]
            ),
        }
    )

    gmat_path = tmp_path / "sample_gmat.csv"
    gmat_df.to_csv(gmat_path, index=False)

    result = interpolate_gmat_ephemeris(gmat_path, cadence)

    trim_mask = (raw_times >= cadence.timestamps.jd[0] - 2430000.0 - 0.0007) & (
        raw_times <= cadence.timestamps.jd[-1] - 2430000.0 + 0.0007
    )
    trimmed_times = raw_times[trim_mask] + 2430000.0 - 2400000.5
    trimmed_values = values[trim_mask]

    assert np.allclose(
        result.earth_ec[:, 0], np.interp(cadence.mjd_utc, trimmed_times, trimmed_values)
    )
    assert np.allclose(
        result.spacecraft_ec[:, 2],
        np.interp(cadence.mjd_utc, trimmed_times, trimmed_values + 12),
    )
    assert np.allclose(
        result.sun_ec[:, 1], np.interp(cadence.mjd_utc, trimmed_times, trimmed_values + 21)
    )

    trimmed_lat = np.linspace(-50, 10, raw_times.size)[trim_mask]
    trimmed_lon = np.array(
        [140.0, 150.0, 170.0, 175.0, 179.0, -179.0, -175.0, -160.0]
    )[trim_mask]
    assert np.allclose(
        result.spacecraft_lat_deg,
        np.interp(cadence.mjd_utc, trimmed_times, trimmed_lat),
    )
    expected_lon = _expected_longitudes(trimmed_times, trimmed_lon, cadence.mjd_utc)
    assert np.allclose(result.spacecraft_lon_deg, expected_lon)


@pytest.mark.parametrize(
    "lat, lon, expected",
    [
        (np.array([-20.0, 5.0]), np.array([-80.0, 0.0]), np.array([1.0, 0.0])),
        (np.array([-39.9, -41.0]), np.array([29.9, -100.0]), np.array([1.0, 0.0])),
    ],
)
def test_compute_saa_crossings(lat, lon, expected):
    mask = compute_saa_crossings(lat, lon)
    assert np.array_equal(mask, expected)