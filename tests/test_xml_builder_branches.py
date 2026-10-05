"""Tests for xml/builder.py uncovered _to_datetime and _build_observational_parameters branches."""

from datetime import datetime

import pytest

from pandorascheduler_rework.xml.builder import (
    _build_observational_parameters,
    _duration_in_seconds,
    observation_sequence,
)


class TestBuildObservationalParametersDatetimeParsing:
    """Branch coverage for datetime parsing and RA/DEC fallback."""

    def test_string_iso_format(self):
        """String timestamps in ISO format are parsed correctly."""
        params = _build_observational_parameters(
            "Target", "High", "2026-03-01T12:00:00Z", "2026-03-01T18:00:00Z", 120.0, 30.0
        )
        assert params["Timing"][2] == "2026-03-01T12:00:00Z"
        assert params["Timing"][3] == "2026-03-01T18:00:00Z"

    def test_string_space_format(self):
        """String timestamps in 'YYYY-MM-DD HH:MM:SS' format."""
        params = _build_observational_parameters(
            "T", "Low", "2026-03-01 12:00:00", "2026-03-01 18:00:00", 0.0, 0.0
        )
        assert "2026-03-01T12:00:00Z" in params["Timing"][2]

    def test_string_date_only_format(self):
        """String timestamps in 'YYYY-MM-DD' format."""
        params = _build_observational_parameters(
            "T", "Low", "2026-03-01", "2026-03-02", 0.0, 0.0
        )
        assert "2026-03-01" in params["Timing"][2]

    def test_datetime_objects(self):
        """Native datetime objects pass through."""
        params = _build_observational_parameters(
            "T", "Low",
            datetime(2026, 3, 1, 12), datetime(2026, 3, 1, 18),
            120.0, 30.0,
        )
        assert "2026-03-01T12:00:00Z" in params["Timing"][2]

    def test_unparseable_timestamps_fallback(self):
        """Unparseable strings are used as-is with a warning."""
        params = _build_observational_parameters(
            "T", "Low", "NOPE", "ALSO_NOPE", 120.0, 30.0
        )
        assert params["Timing"][2] == "NOPE"
        assert params["Timing"][3] == "ALSO_NOPE"

    def test_invalid_ra_dec_sentinel(self):
        """Non-numeric RA/DEC falls back to sentinel -999.0."""
        params = _build_observational_parameters(
            "T", "Low",
            datetime(2026, 3, 1), datetime(2026, 3, 2),
            "bad_ra", "bad_dec",
        )
        assert params["Boresight"][2] == "-999.0"
        assert params["Boresight"][3] == "-999.0"

    def test_none_ra_dec_sentinel(self):
        """None RA/DEC falls back to sentinel."""
        params = _build_observational_parameters(
            "T", "Low",
            datetime(2026, 3, 1), datetime(2026, 3, 2),
            None, None,
        )
        assert params["Boresight"][2] == "-999.0"


class TestDurationInSeconds:
    """Branch coverage for _duration_in_seconds."""

    def test_datetime_objects(self):
        secs = _duration_in_seconds(
            datetime(2026, 3, 1, 0, 0), datetime(2026, 3, 1, 1, 0)
        )
        assert secs == pytest.approx(3600.0)

    def test_iso_strings(self):
        secs = _duration_in_seconds("2026-03-01T00:00:00Z", "2026-03-01T01:00:00Z")
        assert secs == pytest.approx(3600.0)

    def test_unparseable_returns_zero(self):
        secs = _duration_in_seconds("NOPE", "ALSO_NOPE")
        assert secs == 0.0

    def test_none_returns_zero(self):
        secs = _duration_in_seconds(None, None)
        assert secs == 0.0


class TestObservationSequenceWithStringTimestamps:
    """Covers observation_sequence with string timestamps."""

    def test_string_timestamps_produce_valid_xml(self):
        visit = None  # not used in typical calls; the function expects an ET.Element
        import xml.etree.ElementTree as ET

        visit = ET.Element("Visit")
        import pandas as pd

        targ_info = pd.DataFrame(
            {"Star Name": ["TestStar"], "Planet Name": ["TestStar b"], "RA": [120.0], "DEC": [30.0]}
        )
        result = observation_sequence(
            visit, 1, "TestStar b", "High",
            "2026-03-01T00:00:00Z", "2026-03-01T01:00:00Z",
            120.0, 30.0, targ_info,
        )
        assert result.tag == "Observation_Sequence"
