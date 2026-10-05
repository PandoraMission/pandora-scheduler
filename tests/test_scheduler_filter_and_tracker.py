"""Tests for scheduler._filter_visibility_by_time and _initialize_tracker."""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    _filter_visibility_by_time,
    _initialize_tracker,
)


# ---------------------------------------------------------------------------
# _filter_visibility_by_time
# ---------------------------------------------------------------------------


class TestFilterVisibilityByTime:
    """Tests for both legacy and modern filtering paths."""

    @staticmethod
    def _make_vis(n=5, start_mjd=61079.0, step=1 / 1440):
        """Build a small visibility DataFrame with naive-datetime Time_UTC."""
        mjds = [start_mjd + i * step for i in range(n)]
        times_utc = [
            Time(m, format="mjd", scale="utc").to_datetime().replace(tzinfo=None)
            for m in mjds
        ]
        return pd.DataFrame(
            {
                "Time(MJD_UTC)": mjds,
                # Use format without trailing Z to avoid pd.to_datetime producing tz-aware
                "Time_UTC": [t.strftime("%Y-%m-%d %H:%M:%S") for t in times_utc],
                "Visible": [1] * n,
            }
        )

    def test_legacy_mode_uses_mjd(self):
        vis = self._make_vis(n=10)
        start = Time(vis["Time(MJD_UTC)"].iloc[2], format="mjd").to_datetime()
        stop = Time(vis["Time(MJD_UTC)"].iloc[7], format="mjd").to_datetime()
        result = _filter_visibility_by_time(vis, start, stop, use_legacy_mode=True)
        assert len(result) <= 6  # at most indices 2..7

    def test_modern_mode_uses_datetime_column(self):
        vis = self._make_vis(n=10)
        start_mjd = vis["Time(MJD_UTC)"].iloc[2]
        stop_mjd = vis["Time(MJD_UTC)"].iloc[7]
        start = Time(start_mjd, format="mjd").to_datetime().replace(tzinfo=None)
        stop = Time(stop_mjd, format="mjd").to_datetime().replace(tzinfo=None)
        result = _filter_visibility_by_time(vis.copy(), start, stop, use_legacy_mode=False)
        assert len(result) >= 5  # inclusive

    def test_modern_mode_fallback_without_time_utc(self):
        """When Time_UTC column is absent, falls back to MJD conversion."""
        vis = self._make_vis(n=10).drop(columns=["Time_UTC"])
        start = Time(vis["Time(MJD_UTC)"].iloc[0], format="mjd").to_datetime().replace(tzinfo=None)
        stop = Time(vis["Time(MJD_UTC)"].iloc[9], format="mjd").to_datetime().replace(tzinfo=None)
        result = _filter_visibility_by_time(vis.copy(), start, stop, use_legacy_mode=False)
        assert len(result) == 10

    def test_modern_mode_with_datetime64_column(self):
        """Pre-parsed datetime64 column skips reparsing."""
        vis = self._make_vis(n=5)
        vis["Time_UTC"] = pd.to_datetime(vis["Time_UTC"])
        start = vis["Time_UTC"].iloc[1]
        stop = vis["Time_UTC"].iloc[3]
        result = _filter_visibility_by_time(vis, start, stop, use_legacy_mode=False)
        assert len(result) == 3

    def test_empty_result_when_window_has_no_data(self):
        vis = self._make_vis(n=5)
        future = datetime(2030, 1, 1)
        result = _filter_visibility_by_time(
            vis, future, future + timedelta(hours=1), use_legacy_mode=True
        )
        assert len(result) == 0


# ---------------------------------------------------------------------------
# _initialize_tracker
# ---------------------------------------------------------------------------


class TestInitializeTracker:
    """Tests for tracker initialization from target list and visibility files."""

    @staticmethod
    def _make_planet_parquet(path, n_transits=3, coverage=0.5, start_mjd=61079.0):
        """Write a planet visibility parquet with Transit_Start/Stop/Coverage."""
        period = 3.0  # days
        rows = []
        for i in range(n_transits):
            t_start = start_mjd + i * period
            t_stop = t_start + 0.05
            rows.append(
                {
                    "Transit_Start": t_start,
                    "Transit_Stop": t_stop,
                    "Transit_Coverage": coverage,
                    "Transit_Start_UTC": Time(t_start, format="mjd").to_datetime().strftime(
                        "%Y-%m-%dT%H:%M:%SZ"
                    ),
                    "Transit_Stop_UTC": Time(t_stop, format="mjd").to_datetime().strftime(
                        "%Y-%m-%dT%H:%M:%SZ"
                    ),
                }
            )
        pd.DataFrame(rows).to_parquet(path, index=False)

    def test_tracker_shape_matches_target_list(self, tmp_path):
        """Tracker has one row per target in target_list."""
        data_dir = tmp_path / "data"
        targets_dir = data_dir / "targets"

        # Create target list
        target_list = pd.DataFrame(
            {
                "Planet Name": ["Planet A b", "Planet B b"],
                "Star Name": ["Planet A", "Planet B"],
                "Primary Target": [True, False],
                "RA": [120.0, 180.0],
                "DEC": [30.0, -15.0],
                "Number of Transits to Capture": [3, 2],
            }
        )

        # Write planet visibility parquets
        for _, row in target_list.iterrows():
            planet_dir = targets_dir / row["Star Name"] / row["Planet Name"]
            planet_dir.mkdir(parents=True)
            self._make_planet_parquet(
                planet_dir / f"Visibility for {row['Planet Name']}.parquet",
                n_transits=5,
                coverage=0.5,
            )

        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=targets_dir,
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )
        inputs = SchedulerInputs(
            pandora_start=datetime(2026, 2, 5),
            pandora_stop=datetime(2027, 2, 5),
            sched_start=datetime(2026, 2, 5),
            sched_stop=datetime(2027, 2, 5),
            target_list=target_list,
            paths=paths,
            target_definition_files=["primary"],
            primary_target_csv=tmp_path / "prim.csv",
            auxiliary_target_csv=tmp_path / "aux.csv",
            occultation_target_csv=tmp_path / "occ.csv",
            output_dir=tmp_path / "output",
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            transit_coverage_min=0.2,
        )

        tracker = _initialize_tracker(
            inputs, config,
            datetime(2026, 2, 5), datetime(2027, 2, 5),
            datetime(2026, 2, 5), datetime(2027, 2, 5),
        )

        assert len(tracker) == 2
        assert "Planet Name" in tracker.columns
        assert "Transit Priority" in tracker.columns
        assert "Transits Left in Lifetime" in tracker.columns

    def test_tracker_filters_low_coverage_transits(self, tmp_path):
        """Transits below transit_coverage_min are dropped."""
        data_dir = tmp_path / "data"
        targets_dir = data_dir / "targets"

        target_list = pd.DataFrame(
            {
                "Planet Name": ["Lowcov b"],
                "Star Name": ["Lowcov"],
                "Primary Target": [True],
                "RA": [100.0],
                "DEC": [10.0],
                "Number of Transits to Capture": [5],
            }
        )

        planet_dir = targets_dir / "Lowcov" / "Lowcov b"
        planet_dir.mkdir(parents=True)
        self._make_planet_parquet(
            planet_dir / "Visibility for Lowcov b.parquet",
            n_transits=5,
            coverage=0.1,  # all below threshold
        )

        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=targets_dir,
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )
        inputs = SchedulerInputs(
            pandora_start=datetime(2026, 2, 5),
            pandora_stop=datetime(2027, 2, 5),
            sched_start=datetime(2026, 2, 5),
            sched_stop=datetime(2027, 2, 5),
            target_list=target_list,
            paths=paths,
            target_definition_files=["primary"],
            primary_target_csv=tmp_path / "prim.csv",
            auxiliary_target_csv=tmp_path / "aux.csv",
            occultation_target_csv=tmp_path / "occ.csv",
            output_dir=tmp_path / "output",
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            transit_coverage_min=0.5,
        )

        tracker = _initialize_tracker(
            inputs, config,
            datetime(2026, 2, 5), datetime(2027, 2, 5),
            datetime(2026, 2, 5), datetime(2027, 2, 5),
        )

        # All transits filtered out → 0 left
        assert tracker["Transits Left in Lifetime"].iloc[0] == 0

    def test_tracker_missing_visibility_raises(self, tmp_path):
        """Missing visibility file raises FileNotFoundError."""
        data_dir = tmp_path / "data"
        targets_dir = data_dir / "targets"
        targets_dir.mkdir(parents=True)

        target_list = pd.DataFrame(
            {
                "Planet Name": ["Nope b"],
                "Star Name": ["Nope"],
                "Primary Target": [True],
                "RA": [0.0],
                "DEC": [0.0],
                "Number of Transits to Capture": [1],
            }
        )

        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=targets_dir,
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )
        inputs = SchedulerInputs(
            pandora_start=datetime(2026, 2, 5),
            pandora_stop=datetime(2027, 2, 5),
            sched_start=datetime(2026, 2, 5),
            sched_stop=datetime(2027, 2, 5),
            target_list=target_list,
            paths=paths,
            target_definition_files=["primary"],
            primary_target_csv=tmp_path / "p.csv",
            auxiliary_target_csv=tmp_path / "a.csv",
            occultation_target_csv=tmp_path / "o.csv",
            output_dir=tmp_path / "output",
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
        )

        with pytest.raises(FileNotFoundError, match="Nope b"):
            _initialize_tracker(
                inputs, config,
                datetime(2026, 2, 5), datetime(2027, 2, 5),
                datetime(2026, 2, 5), datetime(2027, 2, 5),
            )

    def test_tracker_handles_archive(self, tmp_path):
        """Archive file reduces Transits Needed."""
        data_dir = tmp_path / "data"
        targets_dir = data_dir / "targets"

        target_list = pd.DataFrame(
            {
                "Planet Name": ["Archived b"],
                "Star Name": ["Archived"],
                "Primary Target": [True],
                "RA": [60.0],
                "DEC": [20.0],
                "Number of Transits to Capture": [5],
            }
        )

        planet_dir = targets_dir / "Archived" / "Archived b"
        planet_dir.mkdir(parents=True)
        self._make_planet_parquet(
            planet_dir / "Visibility for Archived b.parquet",
            n_transits=5,
            coverage=0.5,
        )

        # Write archive
        archive = pd.DataFrame({"Target": ["Archived b"], "Count": [1]})
        data_dir.mkdir(parents=True, exist_ok=True)
        archive.to_csv(data_dir / "Pandora_archive.csv", index=False)

        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=targets_dir,
            aux_targets_dir=data_dir / "aux_targets",
            baseline_dir=data_dir / "baseline",
        )
        inputs = SchedulerInputs(
            pandora_start=datetime(2026, 2, 5),
            pandora_stop=datetime(2027, 2, 5),
            sched_start=datetime(2026, 2, 5),
            sched_stop=datetime(2027, 2, 5),
            target_list=target_list,
            paths=paths,
            target_definition_files=["primary"],
            primary_target_csv=tmp_path / "p.csv",
            auxiliary_target_csv=tmp_path / "a.csv",
            occultation_target_csv=tmp_path / "o.csv",
            output_dir=tmp_path / "output",
        )
        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 2, 5),
            window_end=datetime(2027, 2, 5),
            transit_coverage_min=0.2,
        )

        tracker = _initialize_tracker(
            inputs, config,
            datetime(2026, 2, 5), datetime(2027, 2, 5),
            datetime(2026, 2, 5), datetime(2027, 2, 5),
        )

        # Archive subtracted 1 from Transits Needed (5 → 4)
        assert tracker["Transits Needed"].iloc[0] == 4
        assert tracker["Transits Acquired"].iloc[0] == 1
