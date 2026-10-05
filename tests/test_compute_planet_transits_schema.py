"""Tests that _compute_planet_transits always returns a consistent column schema.

Every code path—including the four early-return branches for degenerate inputs—must
produce a DataFrame with the full set of columns so that downstream parquet readers
(observation_utils, _apply_transit_overlaps) never encounter a missing-column error.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time

from pandorascheduler_rework.visibility.catalog import (
    _compute_planet_transits,
    _write_visibility_parquet,
)
from pandorascheduler_rework.config import PandoraSchedulerConfig

# Canonical columns that every planet visibility DataFrame must contain.
REQUIRED_COLUMNS = {
    "Transits",
    "Transit_Start",
    "Transit_Stop",
    "Transit_Start_UTC",
    "Transit_Stop_UTC",
    "Transit_Coverage",
    "SAA_Overlap",
}


def _write_star_parquet(path: Path, t_mjd: np.ndarray) -> None:
    """Write a minimal star-visibility parquet file."""
    n = len(t_mjd)
    df = pd.DataFrame(
        {
            "Time(MJD_UTC)": t_mjd,
            "Visible": np.ones(n, dtype=float),
            "SAA_Crossing": np.zeros(n, dtype=float),
        }
    )
    df.to_parquet(path, index=False)


def _make_planet_row(
    *,
    star_name: str = "TestStar",
    planet_name: str = "TestPlanetb",
    ra: float = 10.0,
    dec: float = 20.0,
    transit_duration_hrs: float = 1.0,
    period_days: float = 3.0,
    epoch_bjd_tdb_minus: float | None = None,
) -> pd.Series:
    """Build a minimal planet row Series."""
    if epoch_bjd_tdb_minus is None:
        # Put the epoch well before a typical window so transits occur inside it.
        epoch_bjd_tdb_minus = Time(
            datetime(2025, 1, 1), scale="tdb", format="datetime"
        ).jd - 2400000.5
    return pd.Series(
        {
            "Star Name": star_name,
            "Planet Name": planet_name,
            "RA": ra,
            "DEC": dec,
            "Transit Duration (hrs)": transit_duration_hrs,
            "Period (days)": period_days,
            "Transit Epoch (BJD_TDB-2400000.5)": epoch_bjd_tdb_minus,
        }
    )


OBSERVER = EarthLocation.from_geodetic(0.0, 0.0, 0.0)
STAR_METADATA: dict[str, tuple[float, float]] = {"TestStar": (10.0, 20.0)}


def _assert_schema(df: pd.DataFrame) -> None:
    missing = REQUIRED_COLUMNS - set(df.columns)
    assert not missing, f"Missing columns in result: {missing}"


# --------------------------------------------------------------------------
# Early-return path 1: empty time grid (t_mjd.size == 0)
# This is actually guarded by the star_visibility.empty check which raises
# FileNotFoundError first. We verify it raises rather than returning a
# broken schema.
# --------------------------------------------------------------------------
class TestEmptyTimeGrid:
    def test_raises_for_empty_star_visibility(self, tmp_path):
        star_path = tmp_path / "star_vis.parquet"
        _write_star_parquet(star_path, np.array([], dtype=float))

        with pytest.raises(FileNotFoundError):
            _compute_planet_transits(
                star_path,
                _make_planet_row(),
                STAR_METADATA,
                OBSERVER,
            )


# --------------------------------------------------------------------------
# Early-return path 2: NaN ephemeris values
# --------------------------------------------------------------------------
class TestNaNEphemeris:
    @pytest.mark.parametrize(
        "field",
        ["Transit Duration (hrs)", "Period (days)", "Transit Epoch (BJD_TDB-2400000.5)"],
    )
    def test_nan_field_returns_full_schema(self, tmp_path, field):
        window_start = datetime(2025, 2, 5)
        cadence = pd.date_range(window_start, periods=360, freq="min")
        t_mjd = Time(list(cadence.to_pydatetime()), scale="utc").mjd

        star_path = tmp_path / "star_vis.parquet"
        _write_star_parquet(star_path, t_mjd)

        row = _make_planet_row()
        row[field] = float("nan")

        result = _compute_planet_transits(
            star_path, row, STAR_METADATA, OBSERVER
        )
        _assert_schema(result)
        assert result.empty


# --------------------------------------------------------------------------
# Early-return path 3: non-positive period
# --------------------------------------------------------------------------
class TestNonPositivePeriod:
    @pytest.mark.parametrize("bad_period", [0.0, -1.0])
    def test_returns_full_schema(self, tmp_path, bad_period):
        window_start = datetime(2025, 2, 5)
        cadence = pd.date_range(window_start, periods=360, freq="min")
        t_mjd = Time(list(cadence.to_pydatetime()), scale="utc").mjd

        star_path = tmp_path / "star_vis.parquet"
        _write_star_parquet(star_path, t_mjd)

        row = _make_planet_row(period_days=bad_period)

        result = _compute_planet_transits(
            star_path, row, STAR_METADATA, OBSERVER
        )
        _assert_schema(result)
        assert result.empty


# --------------------------------------------------------------------------
# Early-return path 4: no transits in window
# --------------------------------------------------------------------------
class TestNoTransitsInWindow:
    def test_returns_full_schema(self, tmp_path):
        # Very short window (10 min) with a long period so no transit falls inside
        window_start = datetime(2025, 2, 5)
        cadence = pd.date_range(window_start, periods=10, freq="min")
        t_mjd = Time(list(cadence.to_pydatetime()), scale="utc").mjd

        star_path = tmp_path / "star_vis.parquet"
        _write_star_parquet(star_path, t_mjd)

        # Period of 365 days, epoch far in the past → no transit in a 10-min window
        row = _make_planet_row(period_days=365.0)

        result = _compute_planet_transits(
            star_path, row, STAR_METADATA, OBSERVER
        )
        _assert_schema(result)
        assert result.empty


# --------------------------------------------------------------------------
# Happy path: transits found — schema must also be correct
# --------------------------------------------------------------------------
class TestHappyPath:
    def test_returns_full_schema_with_data(self, tmp_path):
        window_start = datetime(2025, 2, 5)
        cadence = pd.date_range(window_start, periods=1440 * 5, freq="min")  # 5 days
        t_mjd = Time(list(cadence.to_pydatetime()), scale="utc").mjd

        star_path = tmp_path / "star_vis.parquet"
        _write_star_parquet(star_path, t_mjd)

        # Short period so at least one transit lands in the 5-day window
        row = _make_planet_row(period_days=1.5, transit_duration_hrs=1.0)

        result = _compute_planet_transits(
            star_path, row, STAR_METADATA, OBSERVER
        )
        _assert_schema(result)
        assert not result.empty
        assert (result["SAA_Overlap"] >= 0).all()
        assert (result["Transit_Coverage"] >= 0).all()


# --------------------------------------------------------------------------
# Round-trip: write → read with column selection (what downstream code does)
# --------------------------------------------------------------------------
class TestParquetRoundTrip:
    def test_written_parquet_readable_with_saa_column(self, tmp_path):
        """Parquet written by _write_visibility_parquet must be readable with SAA_Overlap column."""
        window_start = datetime(2025, 2, 5)
        cadence = pd.date_range(window_start, periods=1440 * 5, freq="min")
        t_mjd = Time(list(cadence.to_pydatetime()), scale="utc").mjd

        star_path = tmp_path / "star_vis.parquet"
        _write_star_parquet(star_path, t_mjd)

        row = _make_planet_row(period_days=1.5, transit_duration_hrs=1.0)

        result = _compute_planet_transits(
            star_path, row, STAR_METADATA, OBSERVER
        )
        assert not result.empty

        config = PandoraSchedulerConfig(
            window_start=window_start,
            window_end=window_start + timedelta(days=5),
            output_dir=tmp_path,
            sun_avoidance_deg=45.0,
            moon_avoidance_deg=30.0,
            earth_avoidance_deg=20.0,
        )

        planet_parquet = tmp_path / "planet_visibility.parquet"
        _write_visibility_parquet(result, planet_parquet, config)

        # This is exactly what observation_utils and _apply_transit_overlaps do:
        read_back = pd.read_parquet(
            planet_parquet,
            columns=["Transit_Start", "Transit_Stop", "Transit_Coverage", "SAA_Overlap"],
        )
        assert set(read_back.columns) == {
            "Transit_Start",
            "Transit_Stop",
            "Transit_Coverage",
            "SAA_Overlap",
        }
        assert len(read_back) == len(result)

    def test_empty_df_parquet_readable_with_saa_column(self, tmp_path):
        """Even empty DataFrames must produce parquet files with the full schema."""
        # Build an empty DataFrame with the correct schema directly
        result = pd.DataFrame(
            {
                "Transits": np.array([], dtype=float),
                "Transit_Start": np.array([], dtype=float),
                "Transit_Stop": np.array([], dtype=float),
                "Transit_Start_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Stop_UTC": pd.Series([], dtype="datetime64[ns]"),
                "Transit_Coverage": np.array([], dtype=float),
                "SAA_Overlap": np.array([], dtype=float),
            }
        )
        assert result.empty

        config = PandoraSchedulerConfig(
            window_start=datetime(2025, 2, 5),
            window_end=datetime(2025, 2, 10),
            output_dir=tmp_path,
            sun_avoidance_deg=45.0,
            moon_avoidance_deg=30.0,
            earth_avoidance_deg=20.0,
        )

        planet_parquet = tmp_path / "planet_visibility.parquet"
        _write_visibility_parquet(result, planet_parquet, config)

        # Must not raise even for empty files
        read_back = pd.read_parquet(
            planet_parquet,
            columns=["Transit_Start", "Transit_Stop", "Transit_Coverage", "SAA_Overlap"],
        )
        assert read_back.empty
        assert "SAA_Overlap" in read_back.columns
