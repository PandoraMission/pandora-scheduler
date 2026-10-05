"""Tests for Phase 2: Schedule Output Quality.

Tests that:
- _primary_transit_comment labels primary/secondary correctly
- Auxiliary scheduling produces Comments column
- Below-min-visibility targets are scheduled with a comment
"""

from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.scheduler import (
    AuxiliaryObservationStats,
    SchedulerInputs,
    SchedulerPaths,
    SchedulerState,
    _primary_transit_comment,
    _schedule_auxiliary_target,
)


# ---------------------------------------------------------------------------
# _primary_transit_comment
# ---------------------------------------------------------------------------


class TestPrimaryTransitComment:
    def test_primary_target_column_true(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Primary Target": [True],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "primary exoplanet transit"

    def test_primary_target_column_false(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Primary Target": [False],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "secondary exoplanet transit"

    def test_primary_target_column_numeric_1(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Primary Target": [1],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "primary exoplanet transit"

    def test_primary_target_column_numeric_0(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Primary Target": [0],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "secondary exoplanet transit"

    def test_primary_target_column_nan_treated_as_secondary(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Primary Target": [np.nan],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "secondary exoplanet transit"

    def test_fallback_to_transits_to_capture_10(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Number of Transits to Capture": [10],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "primary exoplanet transit"

    def test_fallback_to_transits_to_capture_3(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Number of Transits to Capture": [3],
        })
        assert _primary_transit_comment(tl, "PlanetA") == "secondary exoplanet transit"

    def test_planet_not_found(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "Primary Target": [True],
        })
        assert _primary_transit_comment(tl, "PlanetB") == ""

    def test_empty_target_list(self):
        tl = pd.DataFrame()
        assert _primary_transit_comment(tl, "PlanetA") == ""

    def test_no_relevant_columns(self):
        tl = pd.DataFrame({
            "Planet Name": ["PlanetA"],
            "RA": [90.0],
        })
        assert _primary_transit_comment(tl, "PlanetA") == ""


# ---------------------------------------------------------------------------
# Auxiliary scheduling — Comments column present
# ---------------------------------------------------------------------------


def _write_vis_parquet(vis_dir: Path, name: str, start: datetime, stop: datetime, visible: bool = True):
    """Create a visibility parquet file for a star with all-visible or all-invisible samples."""
    target_dir = vis_dir / name
    target_dir.mkdir(parents=True, exist_ok=True)
    n = 60  # one-minute samples
    times = [start + timedelta(minutes=i) for i in range(n)]
    mjd_times = Time(times, scale="utc").to_value("mjd")
    flags = [1 if visible else 0] * n
    utc_times = pd.to_datetime([t.isoformat() for t in times])
    pd.DataFrame({
        "Time(MJD_UTC)": mjd_times,
        "Time_UTC": utc_times,
        "Visible": flags,
    }).to_parquet(target_dir / f"Visibility for {name}.parquet", index=False)


class TestAuxiliaryCommentsColumn:
    """Verify the auxiliary scheduler produces a Comments column."""

    def _make_fixtures(self, tmp_path: Path):
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "targets").mkdir()
        (data_dir / "baseline").mkdir()

        # auxiliary target list
        aux_csv = data_dir / "auxiliary-standard_targets.csv"
        pd.DataFrame({
            "Star Name": ["AuxStar1"],
            "RA": [45.0],
            "DEC": [30.0],
            "Priority": [5.0],
            "Number of Hours Requested": [100.0],
        }).to_csv(aux_csv, index=False)

        # visibility data (fully visible)
        aux_vis_dir = data_dir / "aux_targets"
        _write_vis_parquet(aux_vis_dir, "AuxStar1",
                           datetime(2026, 1, 1, 0, 0), datetime(2026, 1, 1, 1, 0))

        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 2, 1),
            output_dir=str(tmp_path / "output"),
            min_sequence_minutes=5,
            std_obs_frequency_days=999,  # skip STD scheduling
        )
        state = SchedulerState(
            tracker=pd.DataFrame(),
            all_target_obs_time={},
            non_primary_obs_time={},
            last_std_obs=datetime(2026, 1, 1),
        )
        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=data_dir / "targets",
            aux_targets_dir=aux_vis_dir,
            baseline_dir=data_dir / "baseline",
        )
        inputs = SchedulerInputs(
            pandora_start=datetime(2026, 1, 1),
            pandora_stop=datetime(2026, 2, 1),
            sched_start=datetime(2026, 1, 1),
            sched_stop=datetime(2026, 2, 1),
            target_list=pd.DataFrame({
                "Planet Name": ["PlanetA"],
                "Star Name": ["StarA"],
                "RA": [90.0],
                "DEC": [45.0],
            }),
            target_definition_files=["exoplanet", "auxiliary-standard"],
            paths=paths,
            primary_target_csv=data_dir / "exoplanet_targets.csv",
            auxiliary_target_csv=aux_csv,
            occultation_target_csv=data_dir / "occultation-standard_targets.csv",
            output_dir=tmp_path / "output",
        )
        return config, state, inputs

    def test_comments_column_exists(self, tmp_path):
        config, state, inputs = self._make_fixtures(tmp_path)
        start = datetime(2026, 1, 1, 0, 0)
        stop = datetime(2026, 1, 1, 0, 30)
        result_df, _ = _schedule_auxiliary_target(start, stop, config, state, inputs)
        assert "Comments" in result_df.columns

    def test_full_visibility_empty_comment(self, tmp_path):
        config, state, inputs = self._make_fixtures(tmp_path)
        start = datetime(2026, 1, 1, 0, 0)
        stop = datetime(2026, 1, 1, 0, 30)
        result_df, _ = _schedule_auxiliary_target(start, stop, config, state, inputs)
        # Target should be scheduled with empty comment (full visibility)
        non_free = result_df[result_df["Target"] != "Free Time"]
        if not non_free.empty:
            assert non_free.iloc[0]["Comments"] == ""


class TestBelowMinVisibilityScheduled:
    """Below-min-visibility targets should be scheduled with a Visibility comment."""

    def _make_fixtures(self, tmp_path: Path, vis_fraction: float = 0.2):
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "targets").mkdir()
        (data_dir / "baseline").mkdir()

        # auxiliary target list
        aux_csv = data_dir / "auxiliary-standard_targets.csv"
        pd.DataFrame({
            "Star Name": ["LowVisStar"],
            "RA": [45.0],
            "DEC": [30.0],
            "Priority": [5.0],
            "Number of Hours Requested": [100.0],
        }).to_csv(aux_csv, index=False)

        # visibility data with partial visibility
        aux_vis_dir = data_dir / "aux_targets"
        target_dir = aux_vis_dir / "LowVisStar"
        target_dir.mkdir(parents=True)
        n = 60
        start = datetime(2026, 1, 1, 0, 0)
        times = [start + timedelta(minutes=i) for i in range(n)]
        mjd_times = Time(times, scale="utc").to_value("mjd")
        # First vis_fraction of samples visible, rest not
        visible_count = int(n * vis_fraction)
        flags = [1] * visible_count + [0] * (n - visible_count)
        utc_times = pd.to_datetime([t.isoformat() for t in times])
        pd.DataFrame({
            "Time(MJD_UTC)": mjd_times,
            "Time_UTC": utc_times,
            "Visible": flags,
        }).to_parquet(target_dir / "Visibility for LowVisStar.parquet", index=False)

        config = PandoraSchedulerConfig(
            window_start=datetime(2026, 1, 1),
            window_end=datetime(2026, 2, 1),
            output_dir=str(tmp_path / "output"),
            min_sequence_minutes=5,
            min_visibility=0.5,  # 50% threshold
            std_obs_frequency_days=999,
        )
        state = SchedulerState(
            tracker=pd.DataFrame(),
            all_target_obs_time={},
            non_primary_obs_time={},
            last_std_obs=datetime(2026, 1, 1),
        )
        paths = SchedulerPaths(
            package_dir=tmp_path,
            data_dir=data_dir,
            targets_dir=data_dir / "targets",
            aux_targets_dir=aux_vis_dir,
            baseline_dir=data_dir / "baseline",
        )
        inputs = SchedulerInputs(
            pandora_start=datetime(2026, 1, 1),
            pandora_stop=datetime(2026, 2, 1),
            sched_start=datetime(2026, 1, 1),
            sched_stop=datetime(2026, 2, 1),
            target_list=pd.DataFrame({
                "Planet Name": ["PlanetA"],
                "Star Name": ["StarA"],
                "RA": [90.0],
                "DEC": [45.0],
            }),
            target_definition_files=["exoplanet", "auxiliary-standard"],
            paths=paths,
            primary_target_csv=data_dir / "exoplanet_targets.csv",
            auxiliary_target_csv=aux_csv,
            occultation_target_csv=data_dir / "occultation-standard_targets.csv",
            output_dir=tmp_path / "output",
        )
        return config, state, inputs

    def test_below_min_vis_still_scheduled(self, tmp_path):
        """A target with 20% visibility should be scheduled (not Free Time)."""
        config, state, inputs = self._make_fixtures(tmp_path, vis_fraction=0.2)
        start = datetime(2026, 1, 1, 0, 0)
        stop = datetime(2026, 1, 1, 0, 30)
        result_df, log_info = _schedule_auxiliary_target(start, stop, config, state, inputs)
        non_free = result_df[result_df["Target"] != "Free Time"]
        assert not non_free.empty, "Below-min-vis target should be scheduled, not Free Time"

    def test_below_min_vis_has_comment(self, tmp_path):
        """The scheduled row should have a Visibility comment."""
        config, state, inputs = self._make_fixtures(tmp_path, vis_fraction=0.2)
        start = datetime(2026, 1, 1, 0, 0)
        stop = datetime(2026, 1, 1, 0, 30)
        result_df, _ = _schedule_auxiliary_target(start, stop, config, state, inputs)
        non_free = result_df[result_df["Target"] != "Free Time"]
        if not non_free.empty:
            comment = str(non_free.iloc[0]["Comments"])
            assert "Visibility" in comment
