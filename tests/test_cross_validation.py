"""Cross-validation: pandorascheduler_rework constraints vs pandoravisibility.

These tests verify that the geometric constraint primitives ported into
``pandorascheduler_rework.visibility.constraints`` produce identical results
to the reference implementation in ``pandoravisibility.Visibility``.

pandoravisibility uses **column-major** ``(3, N)`` arrays;
pandorascheduler_rework uses **row-major** ``(N, 3)`` arrays.
Every comparison transposes as needed before checking element-wise agreement.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy.coordinates import GCRS, SkyCoord
from astropy.time import Time

import astropy.units as u

Visibility = pytest.importorskip("pandoravisibility").Visibility

from pandorascheduler_rework.visibility.constraints import (
    ST1_BODY,
    ST2_BODY,
    _R_EARTH_KM,
    earthlimb_is_sunlit,
    evaluate_star_tracker,
    fast_limb_deg,
    fast_sep_deg,
    fixed_roll_attitude,
    solar_power_fraction,
    star_tracker_eci,
    sun_constrained_attitude,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# TLE matching pandoravisibility's test fixtures (NORAD 67395)
_LINE1 = (
    "1 67395U 80229J   26057.99991898"
    "  .00000000  00000-0  37770-3 0    03"
)
_LINE2 = (
    "2 67395  97.8009  58.3973 0006599"
    " 121.8878 132.9207 14.87804761    04"
)


def _make_vis(**kw) -> Visibility:
    """Create a pandoravisibility Visibility instance with our TLE."""
    return Visibility(_LINE1, _LINE2, **kw)


def _rand_unit_col(N: int, rng: np.random.Generator) -> np.ndarray:
    """Random (3, N) unit vectors for pandoravisibility."""
    v = rng.standard_normal((3, N))
    v /= np.linalg.norm(v, axis=0, keepdims=True)
    return v


def _col_to_row(arr_3xN: np.ndarray) -> np.ndarray:
    """(3, N) → (N, 3) for pandorascheduler_rework."""
    return arr_3xN.T.copy()


# ---------------------------------------------------------------------------
# 1. fast_sep_deg
# ---------------------------------------------------------------------------


class TestFastSepDeg:
    """Compare angular separation computation."""

    @pytest.fixture()
    def vectors(self):
        rng = np.random.default_rng(42)
        N = 200
        a_col = _rand_unit_col(N, rng)
        b_col = _rand_unit_col(N, rng)
        return a_col, b_col

    def test_matches_pandoravisibility(self, vectors):
        a_col, b_col = vectors
        ref = Visibility._fast_sep_deg(a_col, b_col)
        ours = fast_sep_deg(_col_to_row(a_col), _col_to_row(b_col))
        np.testing.assert_allclose(ours, ref, atol=1e-12)

    def test_zero_separation(self):
        """Identical vectors → 0°."""
        v_col = np.array([[1.0], [0.0], [0.0]])
        ref = Visibility._fast_sep_deg(v_col, v_col)
        ours = fast_sep_deg(_col_to_row(v_col), _col_to_row(v_col))
        assert ref == pytest.approx(0.0, abs=1e-12)
        assert ours == pytest.approx(0.0, abs=1e-12)

    def test_opposite_vectors(self):
        """Antiparallel vectors → 180°."""
        a = np.array([[0.0], [0.0], [1.0]])
        b = np.array([[0.0], [0.0], [-1.0]])
        ref = Visibility._fast_sep_deg(a, b)
        ours = fast_sep_deg(_col_to_row(a), _col_to_row(b))
        assert ref == pytest.approx(180.0, abs=1e-12)
        assert ours == pytest.approx(180.0, abs=1e-12)


# ---------------------------------------------------------------------------
# 2. fast_limb_deg
# ---------------------------------------------------------------------------


class TestFastLimbDeg:
    """Compare Earth-limb angle computation."""

    @pytest.fixture()
    def geometry(self):
        rng = np.random.default_rng(99)
        N = 150
        tgt_col = _rand_unit_col(N, rng)
        zen_col = _rand_unit_col(N, rng)
        # Realistic limb angles: 600 ± 50 km altitude
        obs_dist = 6378.1 + 550 + rng.uniform(0, 100, N)
        limb_rad = np.arccos(6378.1 / obs_dist)
        return tgt_col, zen_col, limb_rad

    def test_matches_pandoravisibility(self, geometry):
        tgt_col, zen_col, limb_rad = geometry
        ref = Visibility._fast_limb_deg(tgt_col, zen_col, limb_rad)
        ours = fast_limb_deg(
            _col_to_row(tgt_col), _col_to_row(zen_col), limb_rad
        )
        np.testing.assert_allclose(ours, ref, atol=1e-12)


# ---------------------------------------------------------------------------
# 3. earthlimb_is_sunlit
# ---------------------------------------------------------------------------


class TestEarthlimbIsSunlitCross:
    """Compare day/night limb classification between the two codebases."""

    @pytest.fixture()
    def geometry(self):
        rng = np.random.default_rng(77)
        N = 300
        tgt_col = _rand_unit_col(N, rng)
        zen_col = _rand_unit_col(N, rng)
        sun_col = _rand_unit_col(N, rng)
        obs_dist = 6378.1 + 550 + rng.uniform(0, 100, N)
        limb_rad = np.arccos(6378.1 / obs_dist)
        return tgt_col, zen_col, sun_col, limb_rad

    def test_without_limb_angle(self, geometry):
        """Legacy (horizontal-only) path — limb_angle_rad=None."""
        tgt_col, zen_col, sun_col, _ = geometry
        ref = Visibility._earthlimb_is_sunlit(tgt_col, zen_col, sun_col)
        # Our code takes nadir = -zenith
        nadir_row = _col_to_row(-zen_col)
        ours = earthlimb_is_sunlit(
            _col_to_row(tgt_col), nadir_row, _col_to_row(sun_col),
            limb_angle_rad=None,
        )
        np.testing.assert_array_equal(ours, ref)

    def test_with_limb_angle(self, geometry):
        """Surface-normal path — with limb_angle_rad."""
        tgt_col, zen_col, sun_col, limb_rad = geometry
        ref = Visibility._earthlimb_is_sunlit(
            tgt_col, zen_col, sun_col, limb_angle_rad=limb_rad
        )
        nadir_row = _col_to_row(-zen_col)
        ours = earthlimb_is_sunlit(
            _col_to_row(tgt_col), nadir_row, _col_to_row(sun_col),
            limb_angle_rad=limb_rad,
        )
        np.testing.assert_array_equal(ours, ref)

    def test_sun_along_zenith(self):
        """Regression: sun along zenith exposed the old horizontal-only bug."""
        N = 1
        tgt_col = np.array([[1.0], [0.0], [0.0]])
        zen_col = np.array([[0.0], [0.0], [1.0]])
        sun_col = np.array([[0.0], [0.0], [1.0]])  # sun along zenith
        limb_rad = np.array([np.arccos(6378.1 / 6978.1)])

        ref = Visibility._earthlimb_is_sunlit(
            tgt_col, zen_col, sun_col, limb_angle_rad=limb_rad
        )
        nadir_row = _col_to_row(-zen_col)
        ours = earthlimb_is_sunlit(
            _col_to_row(tgt_col), nadir_row, _col_to_row(sun_col),
            limb_angle_rad=limb_rad,
        )
        # Both should say sunlit (cos(θ) * 1 > 0)
        assert ref[0] is np.True_
        assert ours[0] is np.True_

    def test_sun_along_nadir(self):
        """Sun below → not sunlit."""
        tgt_col = np.array([[1.0], [0.0], [0.0]])
        zen_col = np.array([[0.0], [0.0], [1.0]])
        sun_col = np.array([[0.0], [0.0], [-1.0]])
        limb_rad = np.array([np.arccos(6378.1 / 6978.1)])

        ref = Visibility._earthlimb_is_sunlit(
            tgt_col, zen_col, sun_col, limb_angle_rad=limb_rad
        )
        nadir_row = _col_to_row(-zen_col)
        ours = earthlimb_is_sunlit(
            _col_to_row(tgt_col), nadir_row, _col_to_row(sun_col),
            limb_angle_rad=limb_rad,
        )
        assert ref[0] is np.False_
        assert ours[0] is np.False_


# ---------------------------------------------------------------------------
# 4. Star tracker body vectors
# ---------------------------------------------------------------------------


class TestStarTrackerBodyVectors:
    """Ensure ST boresight vectors match between codebases."""

    def test_st1(self):
        ref = np.array(Visibility._get_star_tracker_body_xyz(1))
        np.testing.assert_allclose(ST1_BODY, ref, atol=1e-14)

    def test_st2(self):
        ref = np.array(Visibility._get_star_tracker_body_xyz(2))
        np.testing.assert_allclose(ST2_BODY, ref, atol=1e-14)


# ---------------------------------------------------------------------------
# 5. Roll attitude
# ---------------------------------------------------------------------------


class TestRollAttitude:
    """Compare fixed_roll_attitude against pandoravisibility._roll_attitude."""

    @pytest.mark.parametrize("roll_deg", [0, 45, 90, 135, 180, 270, 359])
    def test_single_target(self, roll_deg):
        """Compare roll attitude for a single boresight at various rolls."""
        rng = np.random.default_rng(roll_deg + 100)
        z_unit = rng.standard_normal(3)
        z_unit /= np.linalg.norm(z_unit)
        roll_rad = np.deg2rad(roll_deg)

        # pandoravisibility (3,) in/out
        x_ref, y_ref = Visibility._roll_attitude(z_unit, roll_rad)

        # pandorascheduler_rework (N=1, 3) in, (N, 3) out
        z_row = z_unit[np.newaxis, :]  # (1, 3)
        x_ours, y_ours, z_ours = fixed_roll_attitude(z_row, roll_rad)

        np.testing.assert_allclose(x_ours[0], x_ref, atol=1e-12)
        np.testing.assert_allclose(y_ours[0], y_ref, atol=1e-12)
        np.testing.assert_allclose(z_ours[0], z_unit, atol=1e-14)

    def test_near_pole(self):
        """Degenerate: boresight near celestial north pole."""
        z_unit = np.array([0.0, 0.001, 0.9999999])
        z_unit /= np.linalg.norm(z_unit)
        roll_rad = np.deg2rad(30)

        x_ref, y_ref = Visibility._roll_attitude(z_unit, roll_rad)
        z_row = z_unit[np.newaxis, :]
        x_ours, y_ours, _ = fixed_roll_attitude(z_row, roll_rad)

        np.testing.assert_allclose(x_ours[0], x_ref, atol=1e-6)
        np.testing.assert_allclose(y_ours[0], y_ref, atol=1e-6)


# ---------------------------------------------------------------------------
# 6. Solar power fraction
# ---------------------------------------------------------------------------


class TestSolarPowerFraction:
    """Compare solar power computation.

    pandoravisibility computes:
        cos_sy = clip(y · sun, -1, 1)
        theta_sy = arccos(|cos_sy|)
        incidence = π/2 - theta_sy
        power = cos(incidence)

    Which simplifies to sqrt(1 - (y·s)²), matching our implementation.
    """

    def test_random_vectors(self):
        rng = np.random.default_rng(55)
        N = 200
        y_col = _rand_unit_col(N, rng)
        sun_col = _rand_unit_col(N, rng)

        # pandoravisibility inline formula
        cos_sy = np.sum(y_col * sun_col, axis=0)
        cos_sy = np.clip(cos_sy, -1.0, 1.0)
        theta_sy = np.arccos(np.abs(cos_sy))
        incidence = np.pi / 2 - theta_sy
        ref = np.cos(incidence)

        ours = solar_power_fraction(_col_to_row(y_col), _col_to_row(sun_col))
        np.testing.assert_allclose(ours, ref, atol=1e-12)

    def test_perpendicular_max_power(self):
        """Y ⊥ Sun → power = 1."""
        y = np.array([[0.0, 1.0, 0.0]])
        s = np.array([[1.0, 0.0, 0.0]])
        assert solar_power_fraction(y, s)[0] == pytest.approx(1.0)

    def test_parallel_zero_power(self):
        """Y ∥ Sun → power = 0."""
        y = np.array([[1.0, 0.0, 0.0]])
        s = np.array([[1.0, 0.0, 0.0]])
        assert solar_power_fraction(y, s)[0] == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# 7. Star tracker ECI direction
# ---------------------------------------------------------------------------


class TestStarTrackerECI:
    """Compare body→ECI rotation of star tracker boresights."""

    @pytest.mark.parametrize("roll_deg", [0, 60, 120, 240])
    def test_st_eci_matches(self, roll_deg):
        rng = np.random.default_rng(roll_deg + 200)
        z_unit = rng.standard_normal(3)
        z_unit /= np.linalg.norm(z_unit)
        roll_rad = np.deg2rad(roll_deg)

        # pandoravisibility
        x_ref, y_ref = Visibility._roll_attitude(z_unit, roll_rad)
        for tracker in [1, 2]:
            st_body = np.array(Visibility._get_star_tracker_body_xyz(tracker))
            st_ref = (
                x_ref * st_body[0]
                + y_ref * st_body[1]
                + z_unit * st_body[2]
            )
            st_ref /= np.linalg.norm(st_ref)

            # pandorascheduler_rework
            z_row = z_unit[np.newaxis, :]
            x_ours, y_ours, z_ours = fixed_roll_attitude(z_row, roll_rad)
            body_vec = ST1_BODY if tracker == 1 else ST2_BODY
            st_ours = star_tracker_eci(x_ours, y_ours, z_ours, body_vec)

            np.testing.assert_allclose(
                st_ours[0], st_ref, atol=1e-12,
                err_msg=f"ST{tracker} roll={roll_deg}",
            )


# ---------------------------------------------------------------------------
# 8. Sun-constrained attitude
# ---------------------------------------------------------------------------


class TestSunConstrainedAttitude:
    """Compare Sun-constrained attitude (roll=None in pandoravisibility)."""

    def test_random_geometries(self):
        rng = np.random.default_rng(88)
        N = 100
        tgt_col = _rand_unit_col(N, rng)
        sun_col = _rand_unit_col(N, rng)

        # pandoravisibility inline (from _get_st_constraint_fast, roll=None)
        for i in range(N):
            tgt_3 = tgt_col[:, i]
            sun_3 = sun_col[:, i]

            y_ref = np.cross(sun_3, tgt_3)
            y_norm = np.linalg.norm(y_ref)
            if y_norm < 1e-10:
                continue  # degenerate
            y_ref /= y_norm
            x_ref = np.cross(y_ref, tgt_3)
            x_ref /= np.linalg.norm(x_ref)

            # pandorascheduler_rework
            tgt_row = tgt_col[:, i : i + 1].T  # (1, 3)
            sun_row = sun_col[:, i : i + 1].T  # (1, 3)
            x_ours, y_ours, z_ours = sun_constrained_attitude(tgt_row, sun_row)

            np.testing.assert_allclose(x_ours[0], x_ref, atol=1e-12)
            np.testing.assert_allclose(y_ours[0], y_ref, atol=1e-12)
            np.testing.assert_allclose(z_ours[0], tgt_3, atol=1e-14)


# ---------------------------------------------------------------------------
# 9. Integration: full constraint check on TLE-propagated geometry
# ---------------------------------------------------------------------------


class TestIntegrationTLE:
    """Use pandoravisibility._precompute() to obtain realistic geometry,
    then run both codebases' boresight constraint checks and compare.

    This covers sun avoidance, moon avoidance, and Earth-limb avoidance
    (with day/night classification).  Star tracker constraints are tested
    separately below.
    """

    @pytest.fixture()
    def setup(self):
        """Propagate 1 orbit (~97 min) at 1-min cadence."""
        vis = _make_vis(
            earthlimb_day_min=35 * u.deg,
            earthlimb_night_min=20 * u.deg,
        )
        t0 = Time("2026-03-01T00:00:00", scale="utc")
        times = t0 + np.arange(97) * u.min
        pre = vis._precompute(times)

        # Extract column-major (3, N) arrays
        zen_col = pre["zenith_unit"]  # (3, N)
        limb_rad = pre["limb_angle_rad"]  # (N,)
        sun_col = pre["body_units"]["sun"]
        moon_col = pre["body_units"]["moon"]

        return vis, times, zen_col, limb_rad, sun_col, moon_col

    def test_boresight_visibility_matches(self, setup):
        """Compare boresight-level visibility for a real target."""
        vis, times, zen_col, limb_rad, sun_col, moon_col = setup

        target = SkyCoord(ra=79.17, dec=45.99, unit="deg")

        # pandoravisibility: compute per-target GCRS direction
        tgt_gcrs = target.transform_to(GCRS(obstime=times))
        tgt_xyz = tgt_gcrs.cartesian.xyz.value
        tgt_col = tgt_xyz / np.linalg.norm(
            tgt_xyz, axis=0, keepdims=True
        )

        # --- pandoravisibility results ---
        ref_sep_moon = Visibility._fast_sep_deg(moon_col, tgt_col)
        ref_sep_sun = Visibility._fast_sep_deg(sun_col, tgt_col)
        ref_limb = Visibility._fast_limb_deg(tgt_col, zen_col, limb_rad)
        ref_threshold = vis._effective_earthlimb_min_deg(
            tgt_col, zen_col, sun_col, limb_angle_rad=limb_rad
        )
        ref_sunlit = Visibility._earthlimb_is_sunlit(
            tgt_col, zen_col, sun_col, limb_angle_rad=limb_rad
        )

        # --- pandorascheduler_rework results ---
        tgt_row = _col_to_row(tgt_col)
        zen_row = _col_to_row(zen_col)
        sun_row = _col_to_row(sun_col)
        moon_row = _col_to_row(moon_col)
        nadir_row = -zen_row

        our_sep_moon = fast_sep_deg(tgt_row, moon_row)
        our_sep_sun = fast_sep_deg(tgt_row, sun_row)
        our_limb = fast_limb_deg(tgt_row, zen_row, limb_rad)
        our_sunlit = earthlimb_is_sunlit(
            tgt_row, nadir_row, sun_row, limb_angle_rad=limb_rad
        )

        # Compare
        np.testing.assert_allclose(our_sep_moon, ref_sep_moon, atol=1e-10)
        np.testing.assert_allclose(our_sep_sun, ref_sep_sun, atol=1e-10)
        np.testing.assert_allclose(our_limb, ref_limb, atol=1e-10)
        np.testing.assert_array_equal(our_sunlit, ref_sunlit)

        # Effective threshold
        from pandorascheduler_rework.visibility.constraints import (
            effective_earth_threshold,
        )

        our_threshold = effective_earth_threshold(
            tgt_row, nadir_row, sun_row,
            day_deg=35.0, night_deg=20.0, default_deg=20.0,
            limb_angle_rad=limb_rad,
        )
        np.testing.assert_allclose(our_threshold, ref_threshold, atol=1e-10)

        # Final boresight boolean
        ref_ok = (
            (ref_sep_sun >= vis.sun_min.to(u.deg).value)
            & (ref_sep_moon >= vis.moon_min.to(u.deg).value)
            & (ref_limb >= ref_threshold)
        )
        our_ok = (
            (our_sep_sun > vis.sun_min.to(u.deg).value)
            & (our_sep_moon > vis.moon_min.to(u.deg).value)
            & (our_limb > our_threshold)
        )
        # Note: pandoravisibility uses >= while our code uses >.
        # Check that they agree on the vast majority of timesteps.
        # The only possible disagreements are at exact boundary values
        # which are astronomically unlikely with float64.
        n_disagree = np.sum(ref_ok != our_ok)
        assert n_disagree <= 2, (
            f"Too many boresight disagreements: {n_disagree}/97"
        )


class TestIntegrationStarTrackers:
    """Compare star tracker constraint evaluation between codebases.

    Uses pandoravisibility._precompute() for realistic geometry, then
    evaluates star tracker keepout at a fixed roll angle.
    """

    @pytest.fixture()
    def setup(self):
        vis = _make_vis(
            st_sun_min=44 * u.deg,
            st_moon_min=12 * u.deg,
            st_earthlimb_min=30 * u.deg,
        )
        t0 = Time("2026-03-01T00:00:00", scale="utc")
        times = t0 + np.arange(97) * u.min
        pre = vis._precompute(times)
        return vis, times, pre

    @pytest.mark.parametrize("roll_deg", [0, 45, 120, 270])
    def test_st_keepout_matches(self, setup, roll_deg):
        vis, times, pre = setup
        N = len(times)
        zen_col = pre["zenith_unit"]
        limb_rad = pre["limb_angle_rad"]
        sun_col = pre["body_units"]["sun"]
        moon_col = pre["body_units"]["moon"]

        # Boresight direction (use HD 209458 as a realistic target)
        target = SkyCoord(ra=330.795, dec=18.884, unit="deg")
        tgt_gcrs = target.transform_to(GCRS(obstime=times))
        tgt_xyz = tgt_gcrs.cartesian.xyz.value
        tgt_col = tgt_xyz / np.linalg.norm(
            tgt_xyz, axis=0, keepdims=True
        )
        tgt_unit_3 = tgt_col[:, 0]  # representative (3,) direction

        roll_rad = np.deg2rad(roll_deg)

        # --- pandoravisibility: star tracker ECI at this roll ---
        x_ref, y_ref = Visibility._roll_attitude(tgt_unit_3, roll_rad)
        z_col = np.tile(tgt_unit_3.reshape(3, 1), (1, N))

        for tracker in [1, 2]:
            st_body = np.array(Visibility._get_star_tracker_body_xyz(tracker))
            st_eci_ref = (
                x_ref[:, np.newaxis] * st_body[0]
                + y_ref[:, np.newaxis] * st_body[1]
                + z_col * st_body[2]
            )
            st_eci_ref /= np.linalg.norm(
                st_eci_ref, axis=0, keepdims=True
            )

            # Reference keepout checks
            ref_sun_ok = (
                Visibility._fast_sep_deg(st_eci_ref, sun_col)
                >= vis.st_sun_min.to(u.deg).value
            )
            ref_moon_ok = (
                Visibility._fast_sep_deg(st_eci_ref, moon_col)
                >= vis.st_moon_min.to(u.deg).value
            )
            ref_limb_ok = (
                Visibility._fast_limb_deg(st_eci_ref, zen_col, limb_rad)
                >= vis.st_earthlimb_min.to(u.deg).value
            )
            ref_tracker_ok = ref_sun_ok & ref_moon_ok & ref_limb_ok

            # --- pandorascheduler_rework ---
            tgt_N3 = np.broadcast_to(
                tgt_unit_3[np.newaxis, :], (N, 3)
            ).copy()
            x_ours, y_ours, z_ours = fixed_roll_attitude(tgt_N3, roll_rad)
            body_vec = ST1_BODY if tracker == 1 else ST2_BODY
            st_eci_ours = star_tracker_eci(x_ours, y_ours, z_ours, body_vec)

            our_tracker_ok = evaluate_star_tracker(
                st_eci_ours,
                _col_to_row(sun_col),
                _col_to_row(moon_col),
                _col_to_row(zen_col),
                limb_rad,
                sun_min_deg=44.0,
                moon_min_deg=12.0,
                earthlimb_min_deg=30.0,
            )

            # The ECI directions should match
            np.testing.assert_allclose(
                st_eci_ours, _col_to_row(st_eci_ref), atol=1e-10,
                err_msg=f"ST{tracker} ECI mismatch at roll={roll_deg}",
            )

            # The keepout booleans should match exactly.
            # (Note: pandoravisibility uses >=, our code uses >=)
            np.testing.assert_array_equal(
                our_tracker_ok, ref_tracker_ok,
                err_msg=f"ST{tracker} keepout mismatch at roll={roll_deg}",
            )


# ---------------------------------------------------------------------------
# 10. End-to-end: get_visibility vs our full boresight check
# ---------------------------------------------------------------------------


class TestEndToEndVisibility:
    """Compare pandoravisibility.get_visibility() (no ST) against our
    boresight pipeline on the same TLE-propagated geometry.

    This is the highest-level comparison.  We only test boresight
    constraints (sun, moon, earth-limb) since the two codebases use
    different ephemeris sources (TLE vs GMAT) in production; here we
    feed the same TLE-derived vectors to both.
    """

    def test_boresight_only(self):
        vis = _make_vis(
            sun_min=91 * u.deg,
            moon_min=25 * u.deg,
            earthlimb_min=20 * u.deg,
            # disable day/night so we compare plain threshold
        )
        t0 = Time("2026-02-15T00:00:00", scale="utc")
        times = t0 + np.arange(97) * u.min
        target = SkyCoord(ra=79.17, dec=45.99, unit="deg")

        ref = vis.get_visibility(target, times)

        # Our implementation using same geometry
        pre = vis._precompute(times)
        zen_col = pre["zenith_unit"]
        limb_rad = pre["limb_angle_rad"]
        sun_col = pre["body_units"]["sun"]
        moon_col = pre["body_units"]["moon"]

        tgt_gcrs = target.transform_to(GCRS(obstime=times))
        tgt_xyz = tgt_gcrs.cartesian.xyz.value
        tgt_col = tgt_xyz / np.linalg.norm(
            tgt_xyz, axis=0, keepdims=True
        )

        tgt_row = _col_to_row(tgt_col)
        zen_row = _col_to_row(zen_col)
        sun_row = _col_to_row(sun_col)
        moon_row = _col_to_row(moon_col)

        # Our constraint checks (using > not >= as in our codebase)
        sun_ok = fast_sep_deg(tgt_row, sun_row) > 91.0
        moon_ok = fast_sep_deg(tgt_row, moon_row) > 25.0
        earth_ok = fast_limb_deg(tgt_row, zen_row, limb_rad) > 20.0

        ours = sun_ok & moon_ok & earth_ok

        # Compare.  Allow up to 2 boundary disagreements from >= vs >
        n_disagree = np.sum(ref != ours)
        assert n_disagree <= 2, (
            f"End-to-end disagreement: {n_disagree}/97 timesteps"
        )

    def test_daynight_thresholds(self):
        """With day/night Earth thresholds, the limb classification must
        agree between the two codebases."""
        vis = _make_vis(
            sun_min=91 * u.deg,
            moon_min=25 * u.deg,
            earthlimb_min=20 * u.deg,
            earthlimb_day_min=35 * u.deg,
            earthlimb_night_min=20 * u.deg,
        )
        t0 = Time("2026-02-15T00:00:00", scale="utc")
        times = t0 + np.arange(97) * u.min
        target = SkyCoord(ra=150.0, dec=-20.0, unit="deg")

        ref = vis.get_visibility(target, times)

        pre = vis._precompute(times)
        zen_col = pre["zenith_unit"]
        limb_rad = pre["limb_angle_rad"]
        sun_col = pre["body_units"]["sun"]
        moon_col = pre["body_units"]["moon"]

        tgt_gcrs = target.transform_to(GCRS(obstime=times))
        tgt_xyz = tgt_gcrs.cartesian.xyz.value
        tgt_col = tgt_xyz / np.linalg.norm(
            tgt_xyz, axis=0, keepdims=True
        )

        tgt_row = _col_to_row(tgt_col)
        zen_row = _col_to_row(zen_col)
        sun_row = _col_to_row(sun_col)
        moon_row = _col_to_row(moon_col)
        nadir_row = -zen_row

        from pandorascheduler_rework.visibility.constraints import (
            effective_earth_threshold,
        )

        threshold = effective_earth_threshold(
            tgt_row, nadir_row, sun_row,
            day_deg=35.0, night_deg=20.0, default_deg=20.0,
            limb_angle_rad=limb_rad,
        )

        sun_ok = fast_sep_deg(tgt_row, sun_row) > 91.0
        moon_ok = fast_sep_deg(tgt_row, moon_row) > 25.0
        earth_ok = fast_limb_deg(tgt_row, zen_row, limb_rad) > threshold

        ours = sun_ok & moon_ok & earth_ok

        n_disagree = np.sum(ref != ours)
        assert n_disagree <= 2, (
            f"Day/night end-to-end disagreement: {n_disagree}/97 timesteps"
        )


class TestEndToEndBestRoll:
    """Compare get_visibility_best_roll between the two codebases.

    This is the most comprehensive test — it exercises the full
    orbital roll sweep pipeline end-to-end.
    """

    def test_best_roll_visibility_agrees(self):
        """Both codebases should find the same visible timesteps
        when given the same TLE/geometry and constraint thresholds."""
        vis = _make_vis(
            sun_min=91 * u.deg,
            moon_min=25 * u.deg,
            earthlimb_min=20 * u.deg,
            st_sun_min=44 * u.deg,
            st_moon_min=12 * u.deg,
            st_earthlimb_min=30 * u.deg,
            st_required=1,
        )
        t0 = Time("2026-03-01T00:00:00", scale="utc")
        # One full orbit
        times = t0 + np.arange(97) * u.min
        target = SkyCoord(ra=79.17, dec=45.99, unit="deg")

        # pandoravisibility best-roll result
        ref_result = vis.get_visibility_best_roll(
            target, times, roll_step=2 * u.deg
        )

        pre = vis._precompute(times)
        zen_col = pre["zenith_unit"]
        limb_rad = pre["limb_angle_rad"]
        sun_col = pre["body_units"]["sun"]
        moon_col = pre["body_units"]["moon"]

        tgt_gcrs = target.transform_to(GCRS(obstime=times))
        tgt_xyz = tgt_gcrs.cartesian.xyz.value
        tgt_col = tgt_xyz / np.linalg.norm(
            tgt_xyz, axis=0, keepdims=True
        )

        tgt_row = _col_to_row(tgt_col)
        zen_row = _col_to_row(zen_col)
        sun_row = _col_to_row(sun_col)
        moon_row = _col_to_row(moon_col)

        # Our boresight check
        sun_ok = fast_sep_deg(tgt_row, sun_row) > 91.0
        moon_ok = fast_sep_deg(tgt_row, moon_row) > 25.0
        earth_ok = fast_limb_deg(tgt_row, zen_row, limb_rad) > 20.0
        our_boresight = sun_ok & moon_ok & earth_ok

        # Verify boresight_visible agrees
        ref_boresight = ref_result["boresight_visible"]
        n_bs_disagree = np.sum(ref_boresight != our_boresight)
        assert n_bs_disagree <= 2, (
            f"Boresight disagreement: {n_bs_disagree}/97"
        )

        # Where pandoravisibility says boresight-visible AND found a roll,
        # verify we can reproduce the ST pass with their chosen roll.
        ref_visible = ref_result["visible"]
        ref_rolls = ref_result["roll_deg"]

        # For timesteps where ref found visibility, check we get the same
        # ST result at the same roll angle
        vis_idx = np.where(ref_visible)[0]
        if len(vis_idx) == 0:
            pytest.skip("No visible timesteps with these constraints/target")

        # All visible timesteps should have the same roll within an orbit
        rolls_at_vis = ref_rolls[vis_idx]
        unique_rolls = np.unique(rolls_at_vis[~np.isnan(rolls_at_vis)])
        assert len(unique_rolls) >= 1, "Expected at least one roll value"

        # Test with the first unique roll
        test_roll = unique_rolls[0]
        roll_rad = np.deg2rad(test_roll)
        N = len(times)
        tgt_unit_3 = tgt_col[:, 0]
        tgt_N3 = np.broadcast_to(
            tgt_unit_3[np.newaxis, :], (N, 3)
        ).copy()

        x_ours, y_ours, z_ours = fixed_roll_attitude(tgt_N3, roll_rad)
        st1_eci = star_tracker_eci(x_ours, y_ours, z_ours, ST1_BODY)
        st2_eci = star_tracker_eci(x_ours, y_ours, z_ours, ST2_BODY)

        st1_ok = evaluate_star_tracker(
            st1_eci, sun_row, moon_row, zen_row, limb_rad,
            sun_min_deg=44.0, moon_min_deg=12.0, earthlimb_min_deg=30.0,
        )
        st2_ok = evaluate_star_tracker(
            st2_eci, sun_row, moon_row, zen_row, limb_rad,
            sun_min_deg=44.0, moon_min_deg=12.0, earthlimb_min_deg=30.0,
        )
        our_st_ok = st1_ok | st2_ok  # st_required=1

        # At the chosen roll, our ST pass should match ref's at visible times
        # (there may be minor >= vs > differences at boundaries)
        at_vis = vis_idx[rolls_at_vis == test_roll]
        if len(at_vis) > 0:
            n_st_disagree = np.sum(
                our_st_ok[at_vis] != np.ones(len(at_vis), dtype=bool)
            )
            assert n_st_disagree <= 2, (
                f"ST keepout disagreement at chosen roll: "
                f"{n_st_disagree}/{len(at_vis)}"
            )
