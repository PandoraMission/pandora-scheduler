"""Geometry utilities for the visibility pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from astropy.time import Time


@dataclass
class MinuteCadence:
    """Helper describing the interpolated minute grid for an observing window."""

    timestamps: Time
    mjd_utc: np.ndarray


@dataclass
class InterpolatedEphemeris:
    """Interpolated GMAT state vectors on the requested cadence."""

    cadence: MinuteCadence
    earth_ec: np.ndarray
    spacecraft_ec: np.ndarray
    sun_ec: np.ndarray
    moon_ec: np.ndarray
    spacecraft_lat_deg: np.ndarray
    spacecraft_lon_deg: np.ndarray

    @property
    def earth_pc(self) -> np.ndarray:
        """Earth position in the spacecraft-centric frame."""

        return self.earth_ec - self.spacecraft_ec

    @property
    def sun_pc(self) -> np.ndarray:
        """Sun position in the spacecraft-centric frame."""

        return self.sun_ec - self.spacecraft_ec

    @property
    def moon_pc(self) -> np.ndarray:
        """Moon position in the spacecraft-centric frame."""

        return self.moon_ec - self.spacecraft_ec


def build_minute_cadence(window_start: datetime, window_end: datetime) -> MinuteCadence:
    """Create a minute-spaced `astropy.time.Time` grid matching legacy behaviour."""

    if window_end < window_start:
        raise ValueError("window_end must not be earlier than window_start")

    dt_iso_utc = pd.date_range(window_start, window_end, freq="min")
    jd_utc = Time(dt_iso_utc.to_julian_date(), format="jd", scale="utc")
    mjd_utc = Time(jd_utc.value - 2400000.5, format="mjd", scale="utc").value
    return MinuteCadence(timestamps=jd_utc, mjd_utc=np.asarray(mjd_utc, dtype=float))


def interpolate_gmat_ephemeris(
    gmat_path: Path,
    cadence: MinuteCadence,
    spacecraft_name: str | None = None,
) -> InterpolatedEphemeris:
    """Load GMAT output and interpolate state vectors onto the cadence grid."""

    gmat_df = pd.read_csv(gmat_path, sep=",", engine="python")

    time_column = _detect_time_column(gmat_df.columns)
    spacecraft = spacecraft_name or _detect_spacecraft_name(gmat_df.columns)

    cadence_jd = cadence.timestamps.jd
    lower_bound = cadence_jd[0] - 2430000.0 - 0.0007
    upper_bound = cadence_jd[-1] - 2430000.0 + 0.0007

    gmat_df = gmat_df[
        (gmat_df[time_column] >= lower_bound) & (gmat_df[time_column] <= upper_bound)
    ].reset_index(drop=True)

    gmat_mjd_utc = gmat_df[time_column].to_numpy(dtype=float) + 2430000.0 - 2400000.5

    earth_ec = _interp_vectors(
        cadence.mjd_utc, gmat_mjd_utc, _extract_vectors(gmat_df, "Earth")
    )
    spacecraft_ec = _interp_vectors(
        cadence.mjd_utc,
        gmat_mjd_utc,
        _extract_vectors(gmat_df, spacecraft),
    )
    sun_ec = _interp_vectors(
        cadence.mjd_utc, gmat_mjd_utc, _extract_vectors(gmat_df, "Sun")
    )
    moon_ec = _interp_vectors(
        cadence.mjd_utc, gmat_mjd_utc, _extract_vectors(gmat_df, "Luna")
    )

    lat_series = gmat_df[f"{spacecraft}.Earth.Latitude"].to_numpy(dtype=float)
    lon_series = gmat_df[f"{spacecraft}.Earth.Longitude"].to_numpy(dtype=float)

    lat_interp = np.interp(cadence.mjd_utc, gmat_mjd_utc, lat_series)
    lon_interp = _interpolate_longitudes(cadence.mjd_utc, gmat_mjd_utc, lon_series)

    return InterpolatedEphemeris(
        cadence=cadence,
        earth_ec=earth_ec,
        spacecraft_ec=spacecraft_ec,
        sun_ec=sun_ec,
        moon_ec=moon_ec,
        spacecraft_lat_deg=lat_interp,
        spacecraft_lon_deg=lon_interp,
    )


def compute_saa_crossings(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    """Return the legacy SAA mask (1.0 inside, 0.0 outside)."""

    saa_lat_max = 0.0
    saa_lat_min = -40.0
    saa_lon_max = 30.0
    saa_lon_min = -90.0

    mask = (
        (saa_lat_min <= lat_deg)
        & (lat_deg <= saa_lat_max)
        & (saa_lon_min <= lon_deg)
        & (lon_deg <= saa_lon_max)
    )
    return mask.astype(float)


def _detect_time_column(columns: Iterable[str]) -> str:
    for candidate in ("Earth.UTCModJulian", "Earth.A1ModJulian"):
        if candidate in columns:
            return candidate
    raise ValueError("Unable to locate a GMAT time column")


def _detect_spacecraft_name(columns: Iterable[str]) -> str:
    excluded = {"Earth", "Sun", "Luna"}
    for column in columns:
        if column.endswith(".EarthMJ2000Eq.X"):
            name = column.split(".", 1)[0]
            if name not in excluded:
                return name
    raise ValueError("Unable to infer spacecraft name from GMAT columns")


def _extract_vectors(frame: pd.DataFrame, prefix: str) -> np.ndarray:
    cols = [
        f"{prefix}.EarthMJ2000Eq.X",
        f"{prefix}.EarthMJ2000Eq.Y",
        f"{prefix}.EarthMJ2000Eq.Z",
    ]
    missing = [col for col in cols if col not in frame]
    if missing:
        raise KeyError(f"Missing GMAT columns for {prefix}: {missing}")
    return frame[cols].to_numpy(dtype=float)


def _interp_vectors(
    target_times: np.ndarray,
    source_times: np.ndarray,
    source_vectors: np.ndarray,
) -> np.ndarray:
    result = np.empty((target_times.size, source_vectors.shape[1]), dtype=float)
    for axis in range(source_vectors.shape[1]):
        result[:, axis] = np.interp(target_times, source_times, source_vectors[:, axis])
    return result


def _interpolate_longitudes(
    target_times: np.ndarray,
    source_times: np.ndarray,
    longitudes: np.ndarray,
) -> np.ndarray:
    lon_cont = np.copy(longitudes)
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
