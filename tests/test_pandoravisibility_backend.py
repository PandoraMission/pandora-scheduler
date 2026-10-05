"""Tests for the optional pandoravisibility adapter's input handling."""

from datetime import datetime

import numpy as np
import pytest

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.pandoravisibility_backend import (
    earth_center_to_limb_threshold_deg,
    resolve_tle,
)


def test_earth_center_keepout_conversion_is_altitude_dependent():
    limb_angles = np.deg2rad([24.0, 26.0])

    thresholds = earth_center_to_limb_threshold_deg(110.0, limb_angles)

    np.testing.assert_allclose(thresholds, [44.0, 46.0])


def test_resolve_tle_from_named_three_line_file(tmp_path):
    tle_file = tmp_path / "pandora.tle"
    tle_file.write_text(
        "PANDORA\n"
        "1 99999U 26001A   26231.00000000  .00000000  00000-0  00000-0 0  9999\n"
        "2 99999  97.0000 100.0000 0001000   0.0000   0.0000 15.00000000    01\n"
    )
    config = PandoraSchedulerConfig(
        window_start=datetime(2026, 8, 10),
        window_end=datetime(2026, 8, 11),
        visibility_backend="pandoravisibility",
        visibility_tle_file=tle_file,
    )

    line1, line2 = resolve_tle(config)

    assert line1.startswith("1 ")
    assert line2.startswith("2 ")


def test_resolve_tle_rejects_incomplete_file(tmp_path):
    tle_file = tmp_path / "pandora.tle"
    tle_file.write_text("PANDORA\n1 incomplete\n")
    config = PandoraSchedulerConfig(
        window_start=datetime(2026, 8, 10),
        window_end=datetime(2026, 8, 11),
        visibility_backend="pandoravisibility",
        visibility_tle_file=tle_file,
    )

    with pytest.raises(ValueError, match="must contain lines"):
        resolve_tle(config)
