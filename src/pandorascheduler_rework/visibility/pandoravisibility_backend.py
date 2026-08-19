"""Adapter for the optional :mod:`pandoravisibility` visibility engine."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig


def earth_center_to_limb_threshold_deg(center_deg, limb_angle_rad):
    """Convert an Earth-center keepout to package apparent-limb clearance."""
    return np.asarray(center_deg) - 90.0 + np.rad2deg(limb_angle_rad)


def _read_tle(path: Path) -> tuple[str, str]:
    """Read the first valid line-1/line-2 pair from a TLE text file."""
    try:
        lines = [line.strip() for line in path.expanduser().read_text().splitlines()]
    except OSError as exc:
        raise ValueError(f"Could not read visibility TLE file {path}: {exc}") from exc

    line1 = next((line for line in lines if line.startswith("1 ")), None)
    line2 = next((line for line in lines if line.startswith("2 ")), None)
    if line1 is None or line2 is None:
        raise ValueError(
            f"Visibility TLE file {path} must contain lines beginning with "
            "'1 ' and '2 '"
        )
    return line1, line2


def resolve_tle(config: PandoraSchedulerConfig) -> tuple[str, str]:
    """Resolve inline or file-backed TLE configuration."""
    if config.visibility_tle_line1 and config.visibility_tle_line2:
        return config.visibility_tle_line1.strip(), config.visibility_tle_line2.strip()
    assert config.visibility_tle_file is not None
    return _read_tle(config.visibility_tle_file)


class PandoraVisibilityBackend:
    """One reusable package engine and time grid per scheduler worker."""

    def __init__(self, config: PandoraSchedulerConfig, payload: dict) -> None:
        try:
            from pandoravisibility import Visibility
        except ImportError as exc:
            raise RuntimeError(
                "visibility_backend='pandoravisibility' requires the optional "
                "pandoravisibility package. Install the project's crossval dependency "
                "group (for example: poetry install --with crossval)."
            ) from exc

        class EarthCenterVisibility(Visibility):
            """Use scheduler Earth-center limits on package-propagated geometry."""

            def __init__(
                self,
                tle_line1: str,
                tle_line2: str,
                *,
                earth_center_day_deg: float,
                earth_center_night_deg: float,
                **custom_limits,
            ) -> None:
                self.earth_center_day_deg = earth_center_day_deg
                self.earth_center_night_deg = earth_center_night_deg
                # The public package uses apparent-limb clearance. Its Earth
                # check calls the override below, which converts our center
                # limits at every propagated spacecraft position.
                custom_limits["earthlimb_min"] = -90 * u.deg
                custom_limits.pop("earthlimb_day_min", None)
                custom_limits.pop("earthlimb_night_min", None)
                super().__init__(tle_line1, tle_line2, **custom_limits)

            def _effective_earthlimb_min_deg(
                self,
                target_unit,
                zenith_unit,
                sun_unit,
                limb_angle_rad=None,
            ):
                if limb_angle_rad is None:
                    raise ValueError(
                        "Earth-center keepouts require the propagated Earth limb angle"
                    )
                sunlit = self._daynight_is_sunlit(
                    target_unit,
                    zenith_unit,
                    sun_unit,
                    limb_angle_rad=limb_angle_rad,
                )
                center_limit = np.where(
                    sunlit,
                    self.earth_center_day_deg,
                    self.earth_center_night_deg,
                )
                # package limb clearance = center separation - 90 + limb angle
                return earth_center_to_limb_threshold_deg(
                    center_limit, limb_angle_rad
                )

        line1, line2 = resolve_tle(config)
        earth_day = (
            config.earth_avoidance_day_deg
            if config.earth_avoidance_day_deg is not None
            else config.earth_avoidance_deg
        )
        earth_night = (
            config.earth_avoidance_night_deg
            if config.earth_avoidance_night_deg is not None
            else config.earth_avoidance_deg
        )
        limits = {
            "sun_min": config.sun_avoidance_deg * u.deg,
            "moon_min": config.moon_avoidance_deg * u.deg,
            "twilight_margin": config.twilight_margin_deg * u.deg,
            "daynight_mode": config.daynight_mode,
            "st_sun_min": config.st_sun_min_deg * u.deg,
            "st_moon_min": config.st_moon_min_deg * u.deg,
            "st_earthlimb_min": config.st_earthlimb_min_deg * u.deg,
            "st_required": config.st_required,
            "ephemeris_step": 1 * u.min,
        }
        if config.st1_earthlimb_min_deg is not None:
            limits["st1_earthlimb_min"] = config.st1_earthlimb_min_deg * u.deg
        if config.st2_earthlimb_min_deg is not None:
            limits["st2_earthlimb_min"] = config.st2_earthlimb_min_deg * u.deg

        self.engine = EarthCenterVisibility(
            line1,
            line2,
            earth_center_day_deg=earth_day,
            earth_center_night_deg=earth_night,
            **limits,
        )
        self.times = Time(
            np.asarray(payload["Time(MJD_UTC)"], dtype=float),
            format="mjd",
            scale="utc",
        )

    @property
    def package_version(self) -> str:
        try:
            return version("pandoravisibility")
        except PackageNotFoundError:
            return "unknown"

    def build_star(
        self,
        payload: dict,
        star_coord: SkyCoord,
        config: PandoraSchedulerConfig,
    ) -> pd.DataFrame:
        """Compute a scheduler-schema visibility frame with package geometry."""
        result = self.engine.get_visibility_best_roll(
            star_coord,
            self.times,
            roll_step=config.roll_step_deg * u.deg,
            orbit_time_step=1 * u.min,
        )
        separations = self.engine.get_separations(star_coord, self.times)

        visible = np.asarray(result["visible"], dtype=bool)
        power = np.asarray(result["solar_power_frac"], dtype=float)
        # pandoravisibility uses power as a roll tie-breaker. The scheduler's
        # stronger minimum-power contract is retained as a final constraint.
        finite_power = np.isfinite(power)
        visible &= ~finite_power | (power >= config.min_power_frac)

        earth_center_sep = payload["earth_pc"].separation(star_coord).deg
        return pd.DataFrame(
            {
                "Time(MJD_UTC)": payload["Time(MJD_UTC)"],
                "Time_UTC": payload["Time_UTC"],
                "SAA_Crossing": payload["SAA_Crossing"],
                "Visible": np.round(visible.astype(float), 1),
                "Earth_Sep": np.round(earth_center_sep, 3),
                "Moon_Sep": np.round(separations["moon"].to_value(u.deg), 3),
                "Sun_Sep": np.round(separations["sun"].to_value(u.deg), 3),
                "Roll_Deg": np.round(np.asarray(result["roll_deg"], dtype=float), 2),
                "N_ST_Pass": np.asarray(result["n_st_pass"], dtype=int),
            }
        )
