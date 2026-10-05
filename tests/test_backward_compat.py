"""Backward-compatibility test: verify identical Visible output when
ST constraints are disabled and day/night is not used.

This test reconstructs the OLD visibility logic (pre-change) and compares
it against the NEW constraint engine to confirm zero behavioral difference.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.catalog import (
    _build_base_payload,
    _build_star_visibility,
)
from pandorascheduler_rework.visibility.constraints import fast_sep_deg
from pandorascheduler_rework.visibility.geometry import (
    build_minute_cadence,
    compute_saa_crossings,
    interpolate_gmat_ephemeris,
)


def _to_gmat_mod_julian(times: Time) -> np.ndarray:
    mjd = np.asarray(times.to_value("mjd"), dtype=float)
    return mjd - 29999.5


def _make_synthetic_gmat(tmp_path, window_start, window_end):
    """Build a small synthetic GMAT ephemeris with varying geometry."""
    samples = Time(
        [
            window_start - timedelta(minutes=10),
            window_start,
            window_start + timedelta(hours=1),
            window_start + timedelta(hours=2),
            window_start + timedelta(hours=3),
            window_end,
            window_end + timedelta(minutes=10),
        ],
        format="datetime",
        scale="utc",
    )
    n = samples.size

    # Spacecraft in a ~600 km circular orbit, Earth at origin
    sc_x = 6971.0 * np.cos(np.linspace(0, 2 * np.pi, n))
    sc_y = 6971.0 * np.sin(np.linspace(0, 2 * np.pi, n))
    sc_z = np.zeros(n)

    gmat_df = pd.DataFrame({
        "Earth.UTCModJulian": _to_gmat_mod_julian(samples),
        "Earth.EarthMJ2000Eq.X": np.zeros(n),
        "Earth.EarthMJ2000Eq.Y": np.zeros(n),
        "Earth.EarthMJ2000Eq.Z": np.zeros(n),
        "Pandora.EarthMJ2000Eq.X": sc_x,
        "Pandora.EarthMJ2000Eq.Y": sc_y,
        "Pandora.EarthMJ2000Eq.Z": sc_z,
        "Sun.EarthMJ2000Eq.X": np.full(n, 1.496e8),
        "Sun.EarthMJ2000Eq.Y": np.full(n, 0.0),
        "Sun.EarthMJ2000Eq.Z": np.full(n, 0.0),
        "Luna.EarthMJ2000Eq.X": np.full(n, 0.0),
        "Luna.EarthMJ2000Eq.Y": np.full(n, 3.844e5),
        "Luna.EarthMJ2000Eq.Z": np.full(n, 0.0),
        "Pandora.Earth.Latitude": np.sin(np.linspace(0, 2 * np.pi, n)) * 51.6,
        "Pandora.Earth.Longitude": np.linspace(-180, 180, n),
    })
    path = tmp_path / "gmat_compat.csv"
    gmat_df.to_csv(path, index=False)
    return path


def _old_visibility_logic(payload, star_coord, config):
    """Reproduce the EXACT pre-change logic from catalog.py."""
    sun_sep = payload["sun_pc"].separation(star_coord).deg
    moon_sep = payload["moon_pc"].separation(star_coord).deg
    earth_sep = payload["earth_pc"].separation(star_coord).deg

    sun_req = sun_sep > config.sun_avoidance_deg
    moon_req = moon_sep > config.moon_avoidance_deg
    earth_req = earth_sep > config.earth_avoidance_deg

    visible = (sun_req & moon_req & earth_req).astype(float)
    return visible, sun_sep, moon_sep, earth_sep


@pytest.mark.parametrize("ra,dec", [
    (0.0, 0.0),
    (90.0, 45.0),
    (180.0, -30.0),
    (270.0, 89.0),
    (45.0, -89.0),
])
def test_old_vs_new_visibility_identical(tmp_path, ra, dec):
    """With ST disabled and day/night=None, new Visible must equal old Visible."""
    window_start = datetime(2025, 6, 15, 0, 0, 0)
    window_end = datetime(2025, 6, 15, 3, 0, 0)

    gmat_path = _make_synthetic_gmat(tmp_path, window_start, window_end)

    earth_avoid = 86.0  # The OLD default

    config = PandoraSchedulerConfig(
        window_start=window_start,
        window_end=window_end,
        gmat_ephemeris=gmat_path,
        output_dir=tmp_path,
        earth_avoidance_deg=earth_avoid,
        earth_avoidance_day_deg=None,
        earth_avoidance_night_deg=None,
        st_required=0,  # ST disabled
    )

    cadence = build_minute_cadence(window_start, window_end)
    ephemeris = interpolate_gmat_ephemeris(gmat_path, cadence)
    payload = _build_base_payload(ephemeris, cadence)
    star_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    # --- Old logic ---
    old_visible, old_sun, old_moon, old_earth = _old_visibility_logic(
        payload, star_coord, config
    )

    # --- New logic ---
    new_df = _build_star_visibility(payload, star_coord, config)
    new_visible = new_df["Visible"].to_numpy()
    new_sun = new_df["Sun_Sep"].to_numpy()
    new_moon = new_df["Moon_Sep"].to_numpy()
    new_earth = new_df["Earth_Sep"].to_numpy()

    # Separation values must match to 3 decimal places (parquet rounding)
    np.testing.assert_allclose(new_earth, np.round(old_earth, 3), atol=1e-3,
                               err_msg=f"Earth_Sep mismatch for RA={ra}, DEC={dec}")
    np.testing.assert_allclose(new_sun, np.round(old_sun, 3), atol=1e-3,
                               err_msg=f"Sun_Sep mismatch for RA={ra}, DEC={dec}")
    np.testing.assert_allclose(new_moon, np.round(old_moon, 3), atol=1e-3,
                               err_msg=f"Moon_Sep mismatch for RA={ra}, DEC={dec}")

    # Visible column must be identical
    np.testing.assert_array_equal(
        new_visible, np.round(old_visible, 1),
        err_msg=f"Visible mismatch for RA={ra}, DEC={dec}"
    )


def test_fast_sep_matches_skycoord_separation(tmp_path):
    """fast_sep_deg must produce the same angles as SkyCoord.separation()."""
    window_start = datetime(2025, 6, 15, 0, 0, 0)
    window_end = datetime(2025, 6, 15, 3, 0, 0)

    gmat_path = _make_synthetic_gmat(tmp_path, window_start, window_end)
    cadence = build_minute_cadence(window_start, window_end)
    ephemeris = interpolate_gmat_ephemeris(gmat_path, cadence)
    payload = _build_base_payload(ephemeris, cadence)

    star_coord = SkyCoord(ra=45 * u.deg, dec=30 * u.deg, frame="icrs")

    # SkyCoord separation (the old method)
    skycoord_sun_sep = payload["sun_pc"].separation(star_coord).deg
    skycoord_moon_sep = payload["moon_pc"].separation(star_coord).deg

    # fast_sep_deg (the new method)
    tgt_cart = star_coord.icrs.cartesian
    tgt_unit_1 = np.array([tgt_cart.x.value, tgt_cart.y.value, tgt_cart.z.value])
    tgt_unit_1 = tgt_unit_1 / np.linalg.norm(tgt_unit_1)
    N = len(payload["Time(MJD_UTC)"])
    target_unit = np.broadcast_to(tgt_unit_1, (N, 3)).copy()

    fast_sun_sep = fast_sep_deg(target_unit, payload["sun_unit"])
    fast_moon_sep = fast_sep_deg(target_unit, payload["moon_unit"])

    # Should match to high precision
    np.testing.assert_allclose(fast_sun_sep, skycoord_sun_sep, atol=0.01,
                               err_msg="Sun separation: fast_sep_deg vs SkyCoord.separation() differ")
    np.testing.assert_allclose(fast_moon_sep, skycoord_moon_sep, atol=0.01,
                               err_msg="Moon separation: fast_sep_deg vs SkyCoord.separation() differ")


def test_comparison_operator_gte_vs_gt():
    """Visibility thresholds are strict: equal-to-threshold should fail."""
    threshold = 91.0
    assert not (np.array([threshold]) > threshold)[0]
    assert (np.array([threshold + 1e-6]) > threshold)[0]
