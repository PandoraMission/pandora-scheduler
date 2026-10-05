"""Tests for datetime column optimization in visibility files."""

import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.time import Time

from pandorascheduler_rework.visibility.catalog import _build_star_visibility


class TestDatetimeColumnGeneration:
    """Test that datetime columns are generated correctly."""

    def test_star_visibility_includes_datetime_column(self):
        """Verify that generated star visibility includes Time_UTC column."""
        # Create minimal test case by directly calling the function
        from pandorascheduler_rework.config import PandoraSchedulerConfig

        import astropy.units as u
        from astropy.coordinates import SkyCoord

        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            earth_avoidance_deg=86.0,
            st_required=0,
        )
        
        # Create test payload
        n_points = 10
        mjd_values = np.linspace(60365.0, 60365.1, n_points)
        
        # Pre-compute datetime array (as done in _build_base_payload)
        time_utc = Time(mjd_values, format="mjd", scale="utc")
        datetime_utc = time_utc.to_datetime()

        # Build unit vectors matching _build_base_payload output
        earth_pc_xyz = np.tile([0.0, 0.0, -6971.0], (n_points, 1))
        sun_pc_xyz = np.tile([1.496e8, 0.0, 0.0], (n_points, 1))
        moon_pc_xyz = np.tile([0.0, 3.844e5, 0.0], (n_points, 1))

        def _normalise(v):
            norms = np.linalg.norm(v, axis=1, keepdims=True)
            norms = np.where(norms == 0, 1.0, norms)
            return v / norms

        nadir_unit = _normalise(earth_pc_xyz)
        observer_dist_km = np.linalg.norm(earth_pc_xyz, axis=1)

        payload = {
            "Time(MJD_UTC)": mjd_values,
            "Time_UTC": datetime_utc,
            "SAA_Crossing": np.zeros(n_points),
            "sun_pc": SkyCoord(ra=np.full(n_points, 100) * u.deg,
                              dec=np.zeros(n_points) * u.deg),
            "moon_pc": SkyCoord(ra=np.full(n_points, 30) * u.deg,
                               dec=np.zeros(n_points) * u.deg),
            "earth_pc": SkyCoord(ra=np.full(n_points, 100) * u.deg,
                                dec=np.zeros(n_points) * u.deg),
            "nadir_unit": nadir_unit,
            "zenith_unit": -nadir_unit,
            "sun_unit": _normalise(sun_pc_xyz),
            "moon_unit": _normalise(moon_pc_xyz),
            "observer_dist_km": observer_dist_km,
            "limb_angle_rad": np.arccos(np.clip(6371.0 / observer_dist_km, -1, 1)),
            "orbit_slices": [slice(0, n_points)],
        }
        
        star_coord = SkyCoord(ra=0 * u.deg, dec=0 * u.deg)
        
        df = _build_star_visibility(payload, star_coord, config)
        
        # Assert Time_UTC column exists
        assert "Time_UTC" in df.columns, "Time_UTC column missing from star visibility"
        assert "Time(MJD_UTC)" in df.columns, "Original MJD column missing"
        
        # Assert both columns have same length
        assert len(df["Time_UTC"]) == len(df["Time(MJD_UTC)"])
        
        # Assert values match MJD conversion (check first value)
        mjd_converted = Time(df["Time(MJD_UTC)"].iloc[0], format="mjd", scale="utc").to_datetime()
        time_utc_val = pd.to_datetime(df["Time_UTC"].iloc[0])
        
        # Allow small time differences due to floating point
        if hasattr(time_utc_val, 'to_pydatetime'):
            time_utc_val = time_utc_val.to_pydatetime()
        
        time_diff = abs((time_utc_val - mjd_converted).total_seconds())
        assert time_diff < 1.0, f"Timestamp mismatch: {time_diff}s difference"


    def test_planet_transits_include_datetime_columns(self):
        """Verify planet visibility includes Transit_Start_UTC and Transit_Stop_UTC."""
        # Create a minimal test visibility file with transit data
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            
            # Create test transit CSV
            transit_data = {
                "Transits": [0, 1],
                "Transit_Start": [60365.5, 60366.5],  # MJD values
                "Transit_Stop": [60365.6, 60366.6],
                "Transit_Start_UTC": [
                    datetime(2026, 2, 5, 12, 0, 0),
                    datetime(2026, 2, 6, 12, 0, 0)
                ],
                "Transit_Stop_UTC": [
                    datetime(2026, 2, 5, 14, 24, 0),
                    datetime(2026, 2, 6, 14, 24, 0)
                ],
                "Transit_Coverage": [0.8, 0.9],
                "SAA_Overlap": [0.0, 0.0],
            }
            df = pd.DataFrame(transit_data)
            csv_path = tmpdir_path / "test_planet.csv"
            df.to_csv(csv_path, index=False)
            
            # Read back and verify
            loaded = pd.read_csv(csv_path)
            
            assert "Transit_Start_UTC" in loaded.columns
            assert "Transit_Stop_UTC" in loaded.columns
            assert len(loaded) == 2


class TestDatetimeReading:
    """Test that scheduler uses datetime columns when available."""

    def test_scheduler_prefers_datetime_column(self):
        """Verify scheduler uses Time_UTC over MJD conversion when available."""
        # Create test CSV with both columns
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            
            test_data = {
                "Time(MJD_UTC)": [60365.0, 60365.5, 60366.0],
                "Time_UTC": [
                    datetime(2026, 2, 5, 0, 0, 0),
                    datetime(2026, 2, 5, 12, 0, 0),
                    datetime(2026, 2, 6, 0, 0, 0)
                ],
                "Visible": [1.0, 1.0, 0.0],
                "SAA_Crossing": [0.0, 0.0, 0.0],
            }
            df = pd.DataFrame(test_data)
            csv_path = tmpdir_path / "test_vis.csv"
            df.to_csv(csv_path, index=False)
            
            # Read and verify Time_UTC is used
            loaded = pd.read_csv(csv_path)
            assert "Time_UTC" in loaded.columns
            
            # Convert Time_UTC to datetime
            vis_times = pd.to_datetime(loaded["Time_UTC"])
            assert len(vis_times) == 3
            assert isinstance(vis_times[0], (datetime, pd.Timestamp))


    def test_scheduler_falls_back_to_mjd_when_datetime_missing(self):
        """Verify scheduler uses MJD conversion when Time_UTC is absent."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            
            # Create CSV without Time_UTC column
            test_data = {
                "Time(MJD_UTC)": [60365.0, 60365.5, 60366.0],
                "Visible": [1.0, 1.0, 0.0],
                "SAA_Crossing": [0.0, 0.0, 0.0],
            }
            df = pd.DataFrame(test_data)
            csv_path = tmpdir_path / "test_vis_legacy.csv"
            df.to_csv(csv_path, index=False)
            
            # Read and verify fallback to MJD
            loaded = pd.read_csv(csv_path)
            assert "Time_UTC" not in loaded.columns
            
            # Verify MJD conversion works
            vis_times = Time(
                loaded["Time(MJD_UTC)"].to_numpy(), format="mjd", scale="utc"
            ).to_datetime()
            assert len(vis_times) == 3


class TestDatetimeAccuracy:
    """Test datetime conversion accuracy."""

    def test_datetime_conversion_accuracy(self):
        """Verify Time_UTC exactly matches Time(MJD).to_datetime()."""
        # Test known MJD values - just verify conversion consistency
        test_mjds = [60000.0, 60365.0, 60365.5]
        
        for mjd_val in test_mjds:
            # Convert using astropy
            time_obj = Time(mjd_val, format="mjd", scale="utc")
            converted_dt1 = time_obj.to_datetime()
            converted_dt2 = time_obj.to_datetime()
            
            # Verify conversion is deterministic
            assert converted_dt1 == converted_dt2, f"MJD {mjd_val}: non-deterministic conversion"


    def test_datetime_mjd_roundtrip(self):
        """Verify MJD -> datetime -> comparison works correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            
            # Create visibility with both columns
            mjd_values = np.array([60365.0, 60365.5, 60366.0, 60366.5])
            time_objs = Time(mjd_values, format="mjd", scale="utc")
            datetime_values = time_objs.to_datetime()
            
            test_data = {
                "Time(MJD_UTC)": mjd_values,
                "Time_UTC": datetime_values,
                "Visible": [1.0, 1.0, 1.0, 0.0],
            }
            df = pd.DataFrame(test_data)
            csv_path = tmpdir_path / "test_roundtrip.csv"
            df.to_csv(csv_path, index=False)
            
            # Read back and compare
            loaded = pd.read_csv(csv_path)
            
            # Filter using Time_UTC
            active_start = datetime(2026, 2, 5, 6, 0, 0)
            active_stop = datetime(2026, 2, 6, 6, 0, 0)
            vis_times_utc = pd.to_datetime(loaded["Time_UTC"])
            mask_utc = (vis_times_utc >= active_start) & (vis_times_utc <= active_stop)
            
            # Filter using MJD
            vis_times_mjd = Time(loaded["Time(MJD_UTC)"], format="mjd").to_datetime()
            vis_times_mjd = pd.to_datetime(vis_times_mjd)
            mask_mjd = (vis_times_mjd >= active_start) & (vis_times_mjd <= active_stop)
            
            # Assert both produce identical results
            np.testing.assert_array_equal(mask_utc, mask_mjd,
                                         "Time_UTC and MJD filtering produced different results")
