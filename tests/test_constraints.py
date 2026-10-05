"""Unit tests for visibility constraint primitives.

Tests cover the geometric functions ported from pandoravisibility:
day/night Earth detection, star tracker keepout, roll attitude, and
per-orbit roll sweep.
"""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.constraints import (
    ST1_BODY,
    ST2_BODY,
    _normalise,
    detect_orbit_boundaries,
    earthlimb_is_sunlit,
    effective_earth_threshold,
    evaluate_star_tracker,
    fast_limb_deg,
    fast_sep_deg,
    find_best_roll_per_orbit,
    fixed_roll_attitude,
    orbit_slices_from_boundaries,
    solar_power_fraction,
    star_tracker_eci,
    sun_constrained_attitude,
    compute_visibility_with_constraints,
    subsatellite_is_sunlit,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_config(**overrides) -> PandoraSchedulerConfig:
    """Create a test config with sensible defaults; override any field."""
    defaults = dict(
        window_start=datetime(2026, 2, 5),
        window_end=datetime(2027, 2, 5),
    )
    defaults.update(overrides)
    return PandoraSchedulerConfig(**defaults)


def _unit(v):
    """Create a unit vector from v."""
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def _tile(v, n):
    """Tile a 1-D vector into (n, 3)."""
    return np.broadcast_to(np.asarray(v, dtype=float), (n, 3)).copy()


# ============================================================================
# fast_sep_deg
# ============================================================================


class TestFastSepDeg:
    def test_zero_separation(self):
        a = _tile([1, 0, 0], 3)
        b = _tile([1, 0, 0], 3)
        assert np.allclose(fast_sep_deg(a, b), 0.0)

    def test_ninety_degrees(self):
        a = _tile([1, 0, 0], 2)
        b = _tile([0, 1, 0], 2)
        result = fast_sep_deg(a, b)
        assert np.allclose(result, 90.0)

    def test_opposite(self):
        a = _tile([1, 0, 0], 1)
        b = _tile([-1, 0, 0], 1)
        assert np.allclose(fast_sep_deg(a, b), 180.0)

    def test_arbitrary_angle(self):
        # 60° between two vectors
        a = _tile([1, 0, 0], 1)
        b = _tile([0.5, np.sqrt(3) / 2, 0], 1)
        assert np.allclose(fast_sep_deg(a, b), 60.0, atol=1e-10)


# ============================================================================
# fast_limb_deg
# ============================================================================


class TestFastLimbDeg:
    def test_target_at_zenith(self):
        """Target directly above observer → max limb separation."""
        zenith = _tile([0, 0, 1], 1)
        target = _tile([0, 0, 1], 1)
        limb_rad = np.array([np.deg2rad(60)])  # ~60° Earth half-angle
        result = fast_limb_deg(target, zenith, limb_rad)
        # elev = arcsin(1) = 90°, result = 90 + 60 = 150
        assert np.allclose(result, 150.0)

    def test_target_at_horizon(self):
        """Target perpendicular to zenith → elev = 0."""
        zenith = _tile([0, 0, 1], 1)
        target = _tile([1, 0, 0], 1)
        limb_rad = np.array([np.deg2rad(60)])
        result = fast_limb_deg(target, zenith, limb_rad)
        # elev = arcsin(0) = 0°, result = 0 + 60 = 60
        assert np.allclose(result, 60.0)

    def test_target_below_horizon(self):
        """Target toward nadir → negative elevation, limb angle < limb_rad."""
        zenith = _tile([0, 0, 1], 1)
        target = _tile([0, 0, -1], 1)
        limb_rad = np.array([np.deg2rad(60)])
        result = fast_limb_deg(target, zenith, limb_rad)
        # elev = -90°, result = -90 + 60 = -30
        assert result[0] < 0


# ============================================================================
# earthlimb_is_sunlit
# ============================================================================


class TestEarthlimbIsSunlit:
    """Test day/night limb classification."""

    def test_sunlit_limb(self):
        """Target and Sun on the same side → nearest limb is sunlit."""
        nadir = _tile([0, 0, -1], 1)   # observer above Earth
        target = _tile([1, 0, 0], 1)   # target toward +X
        sun = _tile([1, 0, 0], 1)      # Sun same direction
        assert earthlimb_is_sunlit(target, nadir, sun)[0]

    def test_dark_limb(self):
        """Target and Sun on opposite sides → nearest limb is dark."""
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        sun = _tile([-1, 0, 0], 1)
        assert not earthlimb_is_sunlit(target, nadir, sun)[0]

    def test_array(self):
        """Vectorised over multiple timesteps."""
        nadir = _tile([0, 0, -1], 3)
        target = _tile([1, 0, 0], 3)
        sun = np.array([
            [1, 0, 0],     # sunlit
            [-1, 0, 0],    # dark
            [0, 1, 0],     # perpendicular — limb is at +X, Sun at +Y, dot > 0? No.
        ], dtype=float)
        sun = _normalise(sun)
        result = earthlimb_is_sunlit(target, nadir, sun)
        assert result[0] is np.True_
        assert result[1] is np.False_

    def test_legacy_none_limb_angle(self):
        """Without limb_angle_rad, falls back to horizontal projection."""
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        sun = _tile([1, 0, 0], 1)
        result_legacy = earthlimb_is_sunlit(target, nadir, sun, limb_angle_rad=None)
        result_no_arg = earthlimb_is_sunlit(target, nadir, sun)
        assert result_legacy[0] == result_no_arg[0]

    def test_surface_normal_geometry(self):
        """With limb_angle_rad, the zenith component of the surface normal matters.

        At 600 km altitude the limb half-angle θ ≈ 24° so cos(θ) ≈ 0.91.
        The surface normal n = cos(θ)·zenith + sin(θ)·limb_dir.  If the Sun
        is purely along −zenith (i.e. directly below the observer), the old
        horizontal-only check would return sunlit=False (limb_dir·sun = 0),
        but actually n·sun = cos(θ)·(zenith·sun) < 0 so it's still dark.

        More importantly: Sun along +zenith (above observer, opposite nadir)
        means zenith·sun = 1, so n·sun = cos(θ) > 0 → sunlit, regardless of
        limb_dir·sun.  With horizontal-only check, limb_dir·sun = 0 gives
        NOT sunlit.  The surface-normal check correctly identifies this as
        sunlit because the zenith component dominates.
        """
        nadir = _tile([0, 0, -1], 1)     # Earth below
        target = _tile([1, 0, 0], 1)     # target along +X
        sun = _tile([0, 0, 1], 1)        # Sun directly above (along zenith)

        # 600 km altitude: limb angle ≈ arccos(6371/6971)
        R_EARTH = 6371.0
        alt_km = 600.0
        la_rad = np.array([np.arccos(R_EARTH / (R_EARTH + alt_km))])

        # Without limb_angle_rad: horizontal projection only
        # limb_dir = [1,0,0], sun = [0,0,1] → dot = 0 → NOT sunlit
        legacy = earthlimb_is_sunlit(target, nadir, sun, limb_angle_rad=None)
        assert not legacy[0], "Legacy horizontal check should give dark"

        # With limb_angle_rad: surface normal has dominant zenith component
        # n·sun = cos(θ)·1 + sin(θ)·0 = cos(θ) > 0 → sunlit
        fixed = earthlimb_is_sunlit(target, nadir, sun, limb_angle_rad=la_rad)
        assert fixed[0], "Surface-normal check should detect sunlit"

    def test_surface_normal_sun_below(self):
        """Sun directly below (along nadir) → limb is NOT sunlit."""
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        sun = _tile([0, 0, -1], 1)   # Sun along nadir

        R_EARTH = 6371.0
        la_rad = np.array([np.arccos(R_EARTH / (R_EARTH + 600.0))])

        result = earthlimb_is_sunlit(target, nadir, sun, limb_angle_rad=la_rad)
        # n·sun = cos(θ)·(zenith·nadir-dir) + sin(θ)·(limb·nadir-dir)
        # zenith = [0,0,1], sun = [0,0,-1] → zenith·sun = -1
        # n·sun = cos(θ)·(-1) + sin(θ)·0 = -cos(θ) < 0 → NOT sunlit
        assert not result[0]

    # ---- twilight_margin_deg tests ----

    def test_twilight_margin_default_zero(self):
        """Default margin produces same result as explicit margin=0."""
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        sun = _tile([1, 0, 0], 1)
        default = earthlimb_is_sunlit(target, nadir, sun)
        explicit = earthlimb_is_sunlit(target, nadir, sun, twilight_margin_deg=0.0)
        np.testing.assert_array_equal(default, explicit)

    def test_twilight_margin_broadens_sunlit(self):
        """A positive margin classifies marginal cases as sunlit.

        Construct geometry where dot(n, sun) is slightly negative (just past
        terminator).  With margin=0 → dark; with a large enough margin → sunlit.
        """
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        # Sun barely below the horizon of the limb normal:
        # dot(limb_dir, sun) = cos(91°) ≈ -0.0175
        sun = _tile(_unit([np.cos(np.deg2rad(91)), np.sin(np.deg2rad(91)), 0]), 1)

        # Without limb_angle_rad (legacy path), limb_dir·sun ≈ -0.017
        assert not earthlimb_is_sunlit(target, nadir, sun, twilight_margin_deg=0.0)[0]
        # Margin of 5° → threshold = -sin(5°) ≈ -0.087, so -0.017 > -0.087 → sunlit
        assert earthlimb_is_sunlit(target, nadir, sun, twilight_margin_deg=5.0)[0]

    def test_twilight_margin_with_limb_angle(self):
        """Margin works with the full surface-normal path (limb_angle_rad given)."""
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        R_EARTH = 6371.0
        la_rad = np.array([np.arccos(R_EARTH / (R_EARTH + 600.0))])

        # Surface normal: n = [sin(la), 0, cos(la)]  (in x-z plane)
        # Choose sun so that n·sun = cos(angle_between) = -0.02
        # i.e. angle_between = arccos(-0.02) ≈ 91.1° from n
        # Parameterise sun = [sin(β), 0, cos(β)] in x-z plane:
        #   n·sun = cos(β - la) = -0.02  →  β = la + arccos(-0.02)
        beta = la_rad[0] + np.arccos(-0.02)
        sun = _tile([np.sin(beta), 0.0, np.cos(beta)], 1)

        dot_n_sun = (
            np.sin(la_rad[0]) * sun[0, 0] + np.cos(la_rad[0]) * sun[0, 2]
        )
        assert dot_n_sun < 0, "Sanity: n·sun should be slightly negative"

        assert not earthlimb_is_sunlit(target, nadir, sun, la_rad, twilight_margin_deg=0.0)[0]
        # Margin=10° → threshold = -sin(10°) ≈ -0.174, -0.02 > -0.174 → sunlit
        assert earthlimb_is_sunlit(target, nadir, sun, la_rad, twilight_margin_deg=10.0)[0]

    def test_twilight_margin_monotonic(self):
        """Increasing margin never reduces the number of sunlit timesteps."""
        nadir = _tile([0, 0, -1], 10)
        target = _tile([1, 0, 0], 10)
        # Sun sweeps from fully sunlit to fully dark
        angles = np.linspace(0, np.pi, 10)
        sun = np.column_stack([np.cos(angles), np.sin(angles), np.zeros(10)])

        count_0 = earthlimb_is_sunlit(target, nadir, sun, twilight_margin_deg=0.0).sum()
        count_5 = earthlimb_is_sunlit(target, nadir, sun, twilight_margin_deg=5.0).sum()
        count_18 = earthlimb_is_sunlit(target, nadir, sun, twilight_margin_deg=18.0).sum()
        assert count_0 <= count_5 <= count_18


# ============================================================================
# effective_earth_threshold
# ============================================================================


class TestEffectiveEarthThreshold:
    def test_uniform_when_none(self):
        """Both day/night None → returns scalar default."""
        target = _tile([1, 0, 0], 2)
        nadir = _tile([0, 0, -1], 2)
        sun = _tile([1, 0, 0], 2)
        result = effective_earth_threshold(target, nadir, sun, None, None, 86.0)
        assert result == 86.0

    def test_day_night_split(self):
        """Explicit day/night values → per-timestep array."""
        nadir = _tile([0, 0, -1], 2)
        target = _tile([1, 0, 0], 2)
        # zenith = [0,0,1]; sun must have +Z component to be sunlit
        sun = np.array([
            _unit([1, 0, 1]),    # dot(zenith, sun) > 0 → sunlit → day
            _unit([-1, 0, -1]),  # dot(zenith, sun) < 0 → dark   → night
        ])
        result = effective_earth_threshold(target, nadir, sun, 110.0, 80.0, 86.0)
        assert np.isclose(result[0], 110.0)
        assert np.isclose(result[1], 80.0)

    def test_twilight_margin_shifts_classification(self):
        """Margin moves a marginal timestep from night→day threshold.

        Sun at 91° from limb_dir → dot ≈ -0.017 < 0.  With margin=0 it is
        dark (→ night_deg); with margin=5° → threshold ≈ -0.087 so it becomes
        sunlit (→ day_deg).
        """
        nadir = _tile([0, 0, -1], 1)
        target = _tile([1, 0, 0], 1)
        sun = _tile(_unit([np.cos(np.deg2rad(91)), np.sin(np.deg2rad(91)), 0]), 1)

        res0 = effective_earth_threshold(
            target, nadir, sun, 110.0, 80.0, 86.0, twilight_margin_deg=0.0,
        )
        assert np.isclose(res0[0], 80.0), "Margin 0: marginal point → night"

        res5 = effective_earth_threshold(
            target, nadir, sun, 110.0, 80.0, 86.0, twilight_margin_deg=5.0,
        )
        assert np.isclose(res5[0], 110.0), "Margin 5°: marginal point → day"


# ============================================================================
# Config default
# ============================================================================


class TestTwilightMarginConfig:
    def test_default_zero(self):
        """Config defaults to 0.0 twilight margin."""
        cfg = _make_config()
        assert cfg.twilight_margin_deg == 0.0

    def test_custom_margin(self):
        """Config accepts a positive twilight margin."""
        cfg = _make_config(twilight_margin_deg=18.0)
        assert cfg.twilight_margin_deg == 18.0


# ============================================================================
# Attitude construction
# ============================================================================


class TestSunConstrainedAttitude:
    def test_orthonormal(self):
        """Payload axes must be orthonormal."""
        target = _tile(_unit([1, 1, 0]), 5)
        sun = _tile(_unit([0, 0, 1]), 5)
        x, y, z = sun_constrained_attitude(target, sun)
        for i in range(5):
            frame = np.stack([x[i], y[i], z[i]])
            assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-12)

    def test_boresight_is_target(self):
        target = _tile(_unit([0, 0, 1]), 1)
        sun = _tile(_unit([1, 0, 0]), 1)
        _, _, z = sun_constrained_attitude(target, sun)
        assert np.allclose(z[0], _unit([0, 0, 1]))


class TestFixedRollAttitude:
    def test_orthonormal_at_various_rolls(self):
        target = _tile(_unit([1, 0, 0]), 5)
        for roll_deg in [0, 45, 90, 180, 270]:
            x, y, z = fixed_roll_attitude(target, np.deg2rad(roll_deg))
            for i in range(5):
                frame = np.stack([x[i], y[i], z[i]])
                assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-12)

    def test_zero_roll_aligns_with_north(self):
        """At roll = 0, X_payload should be the projection of +Z onto the
        plane perpendicular to the boresight."""
        target = _tile(_unit([1, 0, 0]), 1)  # boresight along +X
        x, y, z = fixed_roll_attitude(target, 0.0)
        # North = [0,0,1], projection onto plane perp to [1,0,0] = [0,0,1]
        assert np.allclose(x[0], [0, 0, 1], atol=1e-12)


# ============================================================================
# Star tracker ECI transform
# ============================================================================


class TestStarTrackerECI:
    def test_identity_frame(self):
        """When payload frame = ECI axes, ST body vector comes through unchanged."""
        x = _tile([1, 0, 0], 1)
        y = _tile([0, 1, 0], 1)
        z = _tile([0, 0, 1], 1)
        st_body = np.array([1.0, 0.0, 0.0])
        result = star_tracker_eci(x, y, z, st_body)
        assert np.allclose(result[0], [1, 0, 0], atol=1e-12)

    def test_st1_st2_symmetric(self):
        """ST1 and ST2 should be symmetric about the XZ plane."""
        assert np.allclose(ST1_BODY[0], ST2_BODY[0])
        assert np.allclose(ST1_BODY[1], -ST2_BODY[1])
        assert np.allclose(ST1_BODY[2], ST2_BODY[2])


# ============================================================================
# evaluate_star_tracker
# ============================================================================


class TestEvaluateStarTracker:
    def test_all_disabled(self):
        """All thresholds = 0 → always passes."""
        st = _tile(_unit([1, 0, 0]), 3)
        sun = _tile(_unit([0, 1, 0]), 3)
        moon = _tile(_unit([0, 0, 1]), 3)
        zenith = _tile(_unit([0, 0, 1]), 3)
        limb = np.full(3, np.deg2rad(60))
        ok = evaluate_star_tracker(st, sun, moon, zenith, limb, 0.0, 0.0, 0.0)
        assert np.all(ok)

    def test_sun_constraint_fail(self):
        """ST pointing near Sun should fail."""
        st = _tile(_unit([1, 0, 0]), 1)
        sun = _tile(_unit([1, 0.01, 0]), 1)  # ~0.5° from ST
        moon = _tile(_unit([0, 0, 1]), 1)
        zenith = _tile(_unit([0, 0, 1]), 1)
        limb = np.full(1, np.deg2rad(60))
        ok = evaluate_star_tracker(st, sun, moon, zenith, limb, 45.0, 0.0, 0.0)
        assert not ok[0]


# ============================================================================
# solar_power_fraction
# ============================================================================


class TestSolarPowerFraction:
    def test_max_power(self):
        """Y axis perpendicular to Sun → panels face Sun → power = 1."""
        y = _tile([0, 1, 0], 1)
        sun = _tile([1, 0, 0], 1)
        assert np.allclose(solar_power_fraction(y, sun), 1.0)

    def test_zero_power(self):
        """Y axis aligned with Sun → panels edge-on → power = 0."""
        y = _tile([0, 1, 0], 1)
        sun = _tile([0, 1, 0], 1)
        assert np.allclose(solar_power_fraction(y, sun), 0.0, atol=1e-12)


# ============================================================================
# orbit boundaries
# ============================================================================


class TestOrbitBoundaries:
    def test_simple_sinusoid(self):
        """Ascending nodes of a sinusoidal latitude."""
        t = np.linspace(0, 10 * np.pi, 1000)
        lat = np.sin(t) * 50  # 5 full cycles
        boundaries = detect_orbit_boundaries(lat)
        # 5 ascending crossings + index 0
        assert boundaries[0] == 0
        assert len(boundaries) >= 5

    def test_single_orbit(self):
        lat = np.sin(np.linspace(0, 2 * np.pi, 100)) * 50
        boundaries = detect_orbit_boundaries(lat)
        slices = orbit_slices_from_boundaries(boundaries, 100)
        # At least covers all data
        total = sum(s.stop - s.start for s in slices)
        assert total == 100


# ============================================================================
# find_best_roll_per_orbit
# ============================================================================


class TestFindBestRollPerOrbit:
    def test_st_required_zero_skips(self):
        """When st_required=0, no roll sweep should happen."""
        config = _make_config(st_required=0)
        N = 10
        target = _tile(_unit([1, 0, 0]), N)
        zenith = _tile(_unit([0, 0, 1]), N)
        sun = _tile(_unit([0, 1, 0]), N)
        moon = _tile(_unit([0, 0, 1]), N)
        limb = np.full(N, np.deg2rad(60))
        bvis = np.ones(N, dtype=bool)

        # The function should still be callable; results don't matter when
        # the caller gates on _st_thresholds_active.
        roll, st_vis, pwr = find_best_roll_per_orbit(
            target, zenith, sun, moon, limb,
            [slice(0, N)], bvis, config,
        )
        assert roll.shape == (N,)

    def test_finds_a_roll(self):
        """With mild ST constraints, the sweep should find a valid roll."""
        config = _make_config(
            st_sun_min_deg=30.0,
            st_moon_min_deg=0.0,
            st_earthlimb_min_deg=0.0,
            st_required=1,
            roll_step_deg=10.0,
            min_power_frac=0.0,
        )
        N = 20
        # Target along +X, Sun along +Y (90° separation — comfortable)
        target = _tile(_unit([1, 0, 0]), N)
        zenith = _tile(_unit([0, 0, 1]), N)
        sun = _tile(_unit([0, 1, 0]), N)
        moon = _tile(_unit([0, 0, -1]), N)
        limb = np.full(N, np.deg2rad(60))
        bvis = np.ones(N, dtype=bool)

        roll, st_vis, pwr = find_best_roll_per_orbit(
            target, zenith, sun, moon, limb,
            [slice(0, N)], bvis, config,
        )
        # Should find at least some minutes with ST pass
        assert np.any(st_vis)
        assert np.any(~np.isnan(roll))


# ============================================================================
# compute_visibility_with_constraints (integration)
# ============================================================================


class TestComputeVisibilityIntegration:
    def test_uniform_earth_avoidance(self):
        """With day/night=None, should replicate simple threshold."""
        config = _make_config(
            earth_avoidance_deg=90.0,
            earth_avoidance_day_deg=None,
            earth_avoidance_night_deg=None,
            st_required=0,
        )
        N = 5
        target = _tile(_unit([1, 0, 0]), N)
        nadir = _tile(_unit([0, 0, -1]), N)
        sun = _tile(_unit([0, 1, 0]), N)
        moon = _tile(_unit([0, 0, 1]), N)
        zenith = -nadir
        dist_km = np.full(N, 6971.0)
        limb = np.arccos(6371.0 / dist_km)
        # Earth centre separation: angle between nadir and target
        # nadir=[0,0,-1], target=[1,0,0] → 90°. With threshold 90° → NOT visible (must be >)
        earth_sep = fast_sep_deg(nadir, target)
        assert np.allclose(earth_sep, 90.0)

        result = compute_visibility_with_constraints(
            target, nadir, sun, moon, dist_km, zenith, limb,
            [slice(0, N)], earth_sep, config,
        )
        assert not np.any(result["visible"])

    def test_daynight_makes_difference(self):
        """Day = 110, night = 80 should produce different visibility."""
        config = _make_config(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            st_required=0,
        )
        N = 2
        # Target at 95° from Earth centre — passes night (80), fails day (110)
        # Set up: nadir along -Z, target slightly off +X
        nadir = _tile(_unit([0, 0, -1]), N)
        zenith = -nadir
        # Make target vector such that angle from nadir ≈ 95°
        # nadir = [0,0,-1], we want dot(nadir, target) = cos(95°) ≈ -0.087
        theta = np.deg2rad(95)
        target = _tile([np.sin(theta), 0, -np.cos(theta)], N)
        target = _normalise(target)

        # Sun on same side as target (sunlit) for first timestep,
        # opposite for second (dark)
        sun = np.array([
            _unit([1, 0, 0]),    # sunlit limb → day threshold 110 → FAIL (95 < 110)
            _unit([-1, 0, 0]),   # dark limb   → night threshold 80 → PASS (95 > 80)
        ])
        moon = _tile(_unit([0, 1, 0]), N)
        dist_km = np.full(N, 6971.0)
        limb = np.arccos(6371.0 / dist_km)
        earth_sep = fast_sep_deg(nadir, target)

        result = compute_visibility_with_constraints(
            target, nadir, sun, moon, dist_km, zenith, limb,
            [slice(0, N)], earth_sep, config,
        )
        # Daytime (sunlit limb) should fail, nighttime should pass
        assert not result["visible"][0]
        assert result["visible"][1]


# ============================================================================
# Config validation
# ============================================================================


class TestConfigValidation:
    def test_st_required_invalid(self):
        with pytest.raises(ValueError, match="st_required"):
            _make_config(st_required=3)

    def test_roll_step_zero(self):
        with pytest.raises(ValueError, match="roll_step_deg"):
            _make_config(roll_step_deg=0.0)

    def test_min_power_frac_range(self):
        with pytest.raises(ValueError, match="min_power_frac"):
            _make_config(min_power_frac=1.5)

    def test_daynight_mode_invalid(self):
        with pytest.raises(ValueError, match="daynight_mode"):
            _make_config(daynight_mode="bogus")


# ============================================================================
# subsatellite_is_sunlit
# ============================================================================


class TestSubsatelliteIsSunlit:
    def test_dayside(self):
        """Spacecraft over dayside: zenith and sun on same hemisphere."""
        nadir = _tile([0, 0, -1], 1)  # zenith = +Z
        sun = _tile([0, 0, 1], 1)     # sun above → dayside
        result = subsatellite_is_sunlit(nadir, sun)
        assert result[0] is True or bool(result[0]) is True

    def test_nightside(self):
        """Spacecraft over nightside: zenith and sun on opposite hemispheres."""
        nadir = _tile([0, 0, -1], 1)
        sun = _tile([0, 0, -1], 1)    # sun below → nightside
        result = subsatellite_is_sunlit(nadir, sun)
        assert bool(result[0]) is False

    def test_terminator_sharp(self):
        """At exact terminator with margin=0 → nightside (not strictly dayside)."""
        nadir = _tile([0, 0, -1], 1)
        sun = _tile([1, 0, 0], 1)     # sun perpendicular → terminator
        result = subsatellite_is_sunlit(nadir, sun, twilight_margin_deg=0.0)
        assert bool(result[0]) is False

    def test_terminator_with_margin(self):
        """At terminator with margin → classified as dayside."""
        nadir = _tile([0, 0, -1], 1)
        sun = _tile([1, 0, 0], 1)
        result = subsatellite_is_sunlit(nadir, sun, twilight_margin_deg=10.0)
        assert bool(result[0]) is True

    def test_array(self):
        """Array input with mixed day/night timesteps."""
        nadir = np.array([
            [0, 0, -1],
            [0, 0, -1],
        ], dtype=float)
        sun = np.array([
            [0, 0, 1],     # day
            [0, 0, -1],    # night
        ], dtype=float)
        result = subsatellite_is_sunlit(nadir, sun)
        assert bool(result[0]) is True
        assert bool(result[1]) is False


# ============================================================================
# effective_earth_threshold with daynight_mode
# ============================================================================


class TestEffectiveEarthThresholdDaynightMode:
    def test_limb_mode_default(self):
        """Default mode='limb' uses earthlimb_is_sunlit."""
        N = 2
        nadir = _tile(_unit([0, 0, -1]), N)
        target = _tile(_unit([1, 0, 0]), N)
        # Sun same direction as target → limb is sunlit → day threshold
        sun = _tile(_unit([1, 0, 0]), N)
        limb = np.full(N, np.arccos(6371.0 / 6971.0))

        result = effective_earth_threshold(
            target, nadir, sun,
            day_deg=110.0, night_deg=80.0, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="limb",
        )
        # With limb mode, sun aligned with target → sunlit → day=110
        np.testing.assert_allclose(result, 110.0)

    def test_subsatellite_mode(self):
        """subsatellite mode uses zenith·sun instead of limb geometry."""
        N = 2
        nadir = _tile(_unit([0, 0, -1]), N)
        target = _tile(_unit([1, 0, 0]), N)
        limb = np.full(N, np.arccos(6371.0 / 6971.0))

        # Sun overhead (+Z): subsatellite is dayside regardless of target dir
        sun_day = _tile(_unit([0, 0, 1]), N)
        result_day = effective_earth_threshold(
            target, nadir, sun_day,
            day_deg=110.0, night_deg=80.0, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="subsatellite",
        )
        np.testing.assert_allclose(result_day, 110.0)

        # Sun below (-Z): subsatellite is nightside
        sun_night = _tile(_unit([0, 0, -1]), N)
        result_night = effective_earth_threshold(
            target, nadir, sun_night,
            day_deg=110.0, night_deg=80.0, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="subsatellite",
        )
        np.testing.assert_allclose(result_night, 80.0)

    def test_modes_can_disagree(self):
        """Show a geometry where limb and subsatellite give different answers.

        Target in +X, zenith in +Z, sun in +X:
        - limb mode: nearest limb to target is in +X; limb normal has +X
          component → sunlit → day threshold
        - subsatellite mode: dot(zenith, sun) = dot(+Z, +X) = 0
          → terminator → NOT sunlit → night threshold
        """
        N = 1
        nadir = _tile(_unit([0, 0, -1]), N)
        target = _tile(_unit([1, 0, 0]), N)
        sun = _tile(_unit([1, 0, 0]), N)
        limb = np.full(N, np.arccos(6371.0 / 6971.0))

        result_limb = effective_earth_threshold(
            target, nadir, sun,
            day_deg=110.0, night_deg=80.0, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="limb",
        )
        result_subsat = effective_earth_threshold(
            target, nadir, sun,
            day_deg=110.0, night_deg=80.0, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="subsatellite",
        )
        # limb: sunlit → 110; subsatellite: terminator → 80
        np.testing.assert_allclose(result_limb, 110.0)
        np.testing.assert_allclose(result_subsat, 80.0)

    def test_no_effect_when_both_none(self):
        """When day/night both None, mode doesn't matter."""
        N = 2
        nadir = _tile(_unit([0, 0, -1]), N)
        target = _tile(_unit([1, 0, 0]), N)
        sun = _tile(_unit([0, 0, 1]), N)
        limb = np.full(N, np.arccos(6371.0 / 6971.0))

        r_limb = effective_earth_threshold(
            target, nadir, sun,
            day_deg=None, night_deg=None, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="limb",
        )
        r_subsat = effective_earth_threshold(
            target, nadir, sun,
            day_deg=None, night_deg=None, default_deg=95.0,
            limb_angle_rad=limb, daynight_mode="subsatellite",
        )
        assert r_limb == 95.0
        assert r_subsat == 95.0

    def test_invalid_daynight_mode_raises(self):
        """Invalid daynight_mode raises ValueError when day/night branching is active."""
        N = 2
        nadir = _tile(_unit([0, 0, -1]), N)
        target = _tile(_unit([1, 0, 0]), N)
        sun = _tile(_unit([0, 0, 1]), N)
        limb = np.full(N, np.arccos(6371.0 / 6971.0))

        with pytest.raises(ValueError, match="daynight_mode"):
            effective_earth_threshold(
                target, nadir, sun,
                day_deg=110.0, night_deg=80.0, default_deg=95.0,
                limb_angle_rad=limb, daynight_mode="bogus",
            )


# ============================================================================
# compute_visibility_with_constraints + daynight_mode
# ============================================================================


class TestComputeVisibilitySubsatellite:
    def test_subsatellite_mode_flips_visibility(self):
        """A geometry that passes with subsatellite night but fails with limb day.

        Target at 95° from Earth centre:
        - day threshold 110 → FAIL  (95 < 110)
        - night threshold 80 → PASS (95 > 80)

        Sun in +Y (perpendicular to target), zenith in +Z:
        - limb mode: nearest limb in target direction (+X component)
          has surface normal with sin(la)*X + cos(la)*Z.  Sun in +Y
          → dot(n, sun) = 0, threshold = 0 → NOT sunlit → night → PASS
        We use sun in -Y for limb to be sunlit:
        Actually, to make limb sunlit we need sun aligned with the
        limb direction.

        Simpler approach: put sun overhead (+Z), well separated from target.
        - limb mode: zenith·sun = cos(la)>0, limb·sun has a component →
          sunlit → day threshold 110 → FAIL (95<110)
        - subsatellite: dot(zenith, sun) = 1 > 0 → dayside → day 110 →
          FAIL too.

        Better: use a geometry where the two modes disagree.
        Sun in +X but slightly below horizon: zenith in +Z, sun in [1,0,-0.1].
        - subsatellite: dot(+Z, normalize([1,0,-0.1])) ≈ -0.0995 < 0 → NIGHT
        - limb: target in +X direction, limb point in +X, surface normal
          has sin(la)*X component, dot(n, sun) has +X contribution → SUNLIT

        We must also ensure sun separation from target is >= 91°.
        Target at 95° from nadir [-Z]: target ≈ [sin95, 0, cos95] ≈ [0.996, 0, 0.087]
        Sun: [1, 0, -0.1] normalised → would be too close to target.

        Use different axes: target along +Y offset from -Z, sun in +X.
        """
        N = 1
        # nadir = -Z, zenith = +Z
        nadir = _tile(_unit([0, 0, -1]), N)
        zenith = -nadir

        # Target at 95° from nadir (Earth centre), in the Y-Z plane
        theta = np.deg2rad(95)
        target = _tile([0, np.sin(theta), -np.cos(theta)], N)
        target = _normalise(target)
        # target ≈ [0, 0.996, 0.087]

        # Sun in +X: sep from target ≈ arccos(0) = 90°
        # Need > 91°, so tilt sun a bit: [-0.1, 0, 0.995]... no.
        # Actually put sun at [-1, 0, 0]:
        # sep from target ≈ arccos(0) = 90°. Not enough.
        # Use sun at [-0.1, -1, 0] normalised:
        # dot(target, sun_unit) = 0 * (-0.0995) + 0.996 * (-0.995) + 0.087 * 0
        #                       ≈ -0.99 → sep ≈ 171° > 91 ✓
        sun = _tile(_unit([-0.1, -1, 0]), N)

        # Moon far away
        moon = _tile(_unit([0.5, 0.5, 0.5]), N)
        # moon-target sep: dot = 0*0.577+0.996*0.577+0.087*0.577 ≈ 0.625 → 51° > 25 ✓

        dist_km = np.full(N, 6971.0)
        limb = np.arccos(6371.0 / dist_km)
        earth_sep = fast_sep_deg(nadir, target)
        # earth_sep = angle between nadir and target
        # nadir = [0,0,-1], target ≈ [0, 0.996, 0.087]
        # dot = 0 + 0 + (-1)*0.087... wait, target has -cos(95)
        # cos(95°) = -0.087, so target ≈ [0, 0.996, 0.087]
        # nadir·target = 0*0 + 0*0.996 + (-1)*0.087 = -0.087
        # sep = arccos(-0.087) ≈ 95° ✓

        # Now check limb-is-sunlit for limb mode:
        # zenith = +Z, target ≈ [0, 0.996, 0.087]
        # dot(target, zenith) = 0.087
        # proj = target - zenith * 0.087 = [0, 0.996, 0]
        # limb_dir = [0, 1, 0]
        # sun_unit ≈ [-0.0995, -0.995, 0]
        # n = cos(la)*zenith + sin(la)*limb_dir
        #   = cos(la)*[0,0,1] + sin(la)*[0,1,0]
        # n·sun = cos(la)*0 + sin(la)*(-0.995) < 0 → NOT sunlit → NIGHT

        # For subsatellite mode:
        # zenith = [0,0,1], sun ≈ [-0.0995, -0.995, 0]
        # dot(zenith, sun) = 0 ≈ 0 → at terminator → NOT sunlit → NIGHT

        # Both modes give NIGHT here. We need them to disagree.
        # Let's use sun = [0.1, -1, 0.2] normalised:
        sun = _tile(_unit([0.1, -1, 0.2]), N)
        # zenith·sun_unit = 0*(...) + 0*(...) + 1*0.195 ≈ 0.195 > 0 → DAYSIDE
        # (for subsatellite)
        # limb n·sun: limb_dir = [0,1,0]
        #   n = cos(la)*[0,0,1] + sin(la)*[0,1,0]
        #   n·sun = cos(la)*0.195 + sin(la)*(-0.976) ≈ 0.91*0.195 + 0.41*(-0.976)
        #         ≈ 0.177 - 0.400 = -0.223 < 0 → NOT sunlit → NIGHT

        # So: subsatellite → DAYSIDE (110 threshold) → 95<110 → FAIL
        #     limb → NIGHT (80 threshold) → 95>80 → PASS

        # Check sun-target separation:
        # target ≈ [0, 0.996, 0.087], sun_unit ≈ [0.098, -0.976, 0.195]
        # dot = 0*0.098 + 0.996*(-0.976) + 0.087*0.195 ≈ -0.972 + 0.017 = -0.955
        # sep ≈ 163° > 91 ✓

        config_limb = _make_config(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            daynight_mode="limb",
            st_required=0,
        )
        config_subsat = _make_config(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            daynight_mode="subsatellite",
            st_required=0,
        )

        r_limb = compute_visibility_with_constraints(
            target, nadir, sun, moon, dist_km, zenith, limb,
            [slice(0, N)], earth_sep, config_limb,
        )
        r_subsat = compute_visibility_with_constraints(
            target, nadir, sun, moon, dist_km, zenith, limb,
            [slice(0, N)], earth_sep, config_subsat,
        )

        # limb mode: nearest limb NOT sunlit → night threshold 80 → PASS (95>80)
        assert r_limb["visible"][0], "limb mode: nightside limb → threshold 80 → should pass"
        # subsatellite mode: spacecraft on dayside → day threshold 110 → FAIL (95<110)
        assert not r_subsat["visible"][0], "subsatellite mode: dayside → threshold 110 → should fail"
