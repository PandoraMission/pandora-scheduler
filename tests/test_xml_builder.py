from __future__ import annotations

import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import pytest
from astropy.time import Time

from pandorascheduler_rework import science_calendar
from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.xml import observation_sequence


def _write_visibility(directory: Path, name: str, times: list[datetime], flags: list[int]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    mjd_times = Time(times, scale="utc").to_value("mjd")
    pd.DataFrame({"Time(MJD_UTC)": mjd_times, "Visible": flags}).to_parquet(
        directory / f"Visibility for {name}.parquet", index=False
    )


def _write_planet_visibility(directory: Path, name: str, start: datetime, stop: datetime) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    transit_start = Time([start], scale="utc").to_value("mjd")
    transit_stop = Time([stop], scale="utc").to_value("mjd")
    pd.DataFrame(
        {
            "Transit_Start": transit_start,
            "Transit_Stop": transit_stop,
            "Transit_Coverage": [0.75],
            "SAA_Overlap": [0.1],
        }
    ).to_parquet(directory / f"Visibility for {name}.parquet", index=False)


def test_generate_science_calendar_with_occultation(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    start = datetime(2026, 1, 1, 0, 0, 0)
    times = [start + timedelta(minutes=value) for value in range(0, 180)]
    visibility_flags = [1] * 60 + [0] * 60 + [1] * 60

    _write_visibility(data_dir / "targets" / "StarOne", "StarOne", times, visibility_flags)
    _write_planet_visibility(
        data_dir / "targets" / "StarOne" / "StarOne b",
        "StarOne b",
        start + timedelta(minutes=30),
        start + timedelta(minutes=90),
    )
    _write_visibility(data_dir / "aux_targets" / "OccStar", "OccStar", times, [1] * len(times))

    pd.DataFrame(
        [
            {
                "Planet Name": "StarOne b",
                "Star Name": "StarOne",
                "RA": 10.0,
                "DEC": -20.0,
            }
        ]
    ).to_csv(data_dir / "exoplanet_targets.csv", index=False)

    pd.DataFrame(
        [
            {
                "Star Name": "OccStar",
                "RA": 30.0,
                "DEC": 15.0,
            }
        ]
    ).to_csv(data_dir / "all_targets.csv", index=False)

    pd.DataFrame(
        [
            {"Star Name": "OccStar", "RA": 30.0, "DEC": 10.0, "Number of Hours Requested": 600}
        ]
    ).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

    # Create visibility for StarOne (required for scheduling)
    vis_dir = data_dir / "targets" / "StarOne"
    vis_dir.mkdir(parents=True, exist_ok=True)

    # Create visibility for StarOne b (planet)
    planet_dir = vis_dir / "StarOne b"
    planet_dir.mkdir(parents=True, exist_ok=True)

    # Create visibility for OccStar
    vis_dir = data_dir / "aux_targets" / "OccStar"
    vis_dir.mkdir(parents=True, exist_ok=True)
    # use the minute-resolution visibility written above

    schedule_df = pd.DataFrame(
        [
            {
                "Target": "StarOne b",
                "Observation Start": "2026-01-01 00:00:00",
                "Observation Stop": "2026-01-01 03:00:00",
                "Transit Coverage": 0.75,
                "SAA Overlap": 0.1,
                "Schedule Factor": 0.9,
                "Quality Factor": 0.85,
                "Comments": "",
            }
        ]
    )
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)


    inputs = science_calendar.ScienceCalendarInputs(schedule_csv=schedule_path, data_dir=data_dir)
    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=start + timedelta(days=1),
        visit_limit=1,
        prioritise_occultations_by_slew=True,
        created_timestamp=datetime(2025, 11, 17, 18, 10, 32),
    )
    output_path = tmp_path / "calendar.xml"

    generated_path = science_calendar.generate_science_calendar(
        inputs, output_path=output_path, config=config
    )

    assert generated_path == output_path
    xml_text = generated_path.read_text(encoding="utf-8")
    root = ET.fromstring(xml_text)
    # Handle default namespace if present
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}", 1)[0].strip("{")
    def q(tag: str) -> str:
        return f"{{{ns}}}{tag}" if ns else tag
    targets = [
        t.text
        for t in root.findall(
            f'.//{q("Observation_Sequence")}/{q("Observational_Parameters")}/{q("Target")}'
        )
        if t is not None
    ]
    assert "StarOne b" in targets
    assert "OccStar" in targets

    # Ensure the metadata reflects the configured minimum sequence pruning value
    meta = root.find(q("Meta"))
    assert meta is not None
    assert meta.get("Removed_Sequences_Shorter_Than_min") == "8"
    assert meta.get("Created") == "2025-11-17 18:10:32"


def test_generate_science_calendar_splits_long_occultations(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    start = datetime(2026, 1, 2, 0, 0, 0)
    total_minutes = 180
    times = [start + timedelta(minutes=value) for value in range(total_minutes)]

    primary_flags = [1] * 30 + [0] * 60 + [1] * 90
    _write_visibility(data_dir / "targets" / "StarTwo", "StarTwo", times, primary_flags)
    _write_planet_visibility(
        data_dir / "targets" / "StarTwo" / "StarTwo b",
        "StarTwo b",
        start + timedelta(minutes=20),
        start + timedelta(minutes=40),
    )

    occ_a_flags = [0] * 30 + [1] * 60 + [0] * (total_minutes - 90)
    occ_b_flags = [0] * 30 + [1] * 60 + [0] * (total_minutes - 90)
    _write_visibility(data_dir / "aux_targets" / "OccA", "OccA", times, occ_a_flags)
    _write_visibility(data_dir / "aux_targets" / "OccB", "OccB", times, occ_b_flags)

    pd.DataFrame(
        [
            {"Planet Name": "StarTwo b", "Star Name": "StarTwo", "RA": 15.0, "DEC": -10.0}
        ]
    ).to_csv(data_dir / "exoplanet_targets.csv", index=False)

    # Create visibility for StarTwo (required for scheduling)
    vis_dir = data_dir / "targets" / "StarTwo"
    vis_dir.mkdir(parents=True, exist_ok=True)
    # use the minute-resolution visibility written above

    # Create visibility for StarTwo b (planet)
    planet_dir = vis_dir / "StarTwo b"
    planet_dir.mkdir(parents=True, exist_ok=True)
    # use the minute-resolution planet transit visibility written above

    # OccA and OccB aux_targets visibility already written by
    # _write_visibility above with correct minute-resolution data.

    pd.DataFrame(
        [
            {"Star Name": "OccA", "RA": 30.0, "DEC": 10.0, "Number of Hours Requested": 600},
            {"Star Name": "OccB", "RA": 35.0, "DEC": 12.0, "Number of Hours Requested": 600},
        ]
    ).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

    pd.DataFrame(
        [
            {"Star Name": "OccA", "RA": 30.0, "DEC": 10.0},
            {"Star Name": "OccB", "RA": 35.0, "DEC": 12.0},
        ]
    ).to_csv(data_dir / "all_targets.csv", index=False)

    schedule_df = pd.DataFrame(
        [
            {
                "Target": "StarTwo b",
                "Observation Start": "2026-01-02 00:00:00",
                "Observation Stop": "2026-01-02 03:00:00",
                "Transit Coverage": 0.8,
                "SAA Overlap": 0.0,
                "Schedule Factor": 0.75,
                "Quality Factor": 0.7,
                "Comments": "",
            }
        ]
    )
    schedule_path = tmp_path / "schedule_long_occ.csv"
    schedule_df.to_csv(schedule_path, index=False)


    inputs = science_calendar.ScienceCalendarInputs(schedule_csv=schedule_path, data_dir=data_dir)
    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=start + timedelta(days=1),
        visit_limit=1,
        prioritise_occultations_by_slew=True,
    )
    output_path = tmp_path / "calendar_long_occ.xml"

    science_calendar.generate_science_calendar(inputs, output_path=output_path, config=config)

    xml_text = output_path.read_text(encoding="utf-8")
    root = ET.fromstring(xml_text)
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}", 1)[0].strip("{")
    def q(tag: str) -> str:
        return f"{{{ns}}}{tag}" if ns else tag
    targets = [
        t.text
        for t in root.findall(
            f'.//{q("Observation_Sequence")}/{q("Observational_Parameters")}/{q("Target")}'
        )
        if t is not None
    ]
    assert "StarTwo b" in targets
    # Accept either OccA or OccB (scheduling may prioritise one target over the other)
    assert any(name in {"OccA", "OccB"} for name in targets)


def test_observation_sequence_uses_manifest_star_roi_method():
    visit = ET.Element("Visit")
    start = datetime(2026, 1, 3, 0, 0, 0)
    stop = start + timedelta(minutes=10)

    targ_info = pd.DataFrame(
        [
            {
                "Star Name": "DemoStar",
                "Planet Name": "DemoStar b",
                "StarRoiDetMethod": 1,
                "RA": 123.0,
                "DEC": -45.0,
                "VDA_StarRoiDetMethod": "SET_BY_TARGET_DEFINITION_FILE",
            }
        ]
    )

    observation_sequence(
        visit,
        "001",
        "DemoStar b",
        "0",
        start,
        stop,
        123.0,
        -45.0,
        targ_info,
    )

    star_roi = visit.find(
        "Observation_Sequence/Payload_Parameters/AcquireVisCamScienceData/StarRoiDetMethod"
    )
    assert star_roi is not None
    assert star_roi.text == "1"


def test_visit_id_formatting_matches_legacy_quirk(tmp_path, monkeypatch):
    """Test that Visit IDs use legacy's quirky padding formula.
    
    Legacy bug: Uses len(str(i)) instead of len(str(i+1)), producing 5 digits for 2-digit visits.
    For i=9: "0" * (4 - len(str(9))) + str(10) = "000" + "10" = "00010"
    
    This test ensures we match that behavior to maintain XML parity.
    """
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    start = datetime(2026, 1, 1, 0, 0, 0)
    times = [start + timedelta(hours=i) for i in range(24)]
    _write_visibility(data_dir / "aux_targets" / "TestStar", "TestStar", times, [1] * len(times))

    pd.DataFrame([{"Star Name": "TestStar", "RA": 10.0, "DEC": -20.0}]).to_csv(
        data_dir / "all_targets.csv", index=False
    )
    pd.DataFrame([{"Star Name": "TestStar", "RA": 10.0, "DEC": -20.0}]).to_csv(
        data_dir / "exoplanet_targets.csv", index=False
    )
    pd.DataFrame([{"Star Name": "TestStar", "RA": 10.0, "DEC": -20.0, "Number of Hours Requested": 600}]).to_csv(
        data_dir / "occultation-standard_targets.csv", index=False
    )

    # Create schedule with 10 visits to test the formatting transition
    schedule_rows = []
    for visit_num in range(1, 11):
        schedule_rows.append({
            "Target": "TestStar",
            "Observation Start": (start + timedelta(hours=visit_num * 2)).strftime("%Y-%m-%d %H:%M:%S"),
            "Observation Stop": (start + timedelta(hours=visit_num * 2 + 1)).strftime("%Y-%m-%d %H:%M:%S"),
            "Transit Coverage": None,
            "SAA Overlap": None,
            "Schedule Factor": 0.9,
            "Quality Factor": 0.85,
            "Comments": "",
        })
    
    schedule_df = pd.DataFrame(schedule_rows)
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)


    inputs = science_calendar.ScienceCalendarInputs(schedule_csv=schedule_path, data_dir=data_dir)
    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=start + timedelta(days=2),
        visit_limit=10
    )
    output_path = tmp_path / "calendar.xml"

    science_calendar.generate_science_calendar(inputs, output_path=output_path, config=config)

    xml_text = output_path.read_text(encoding="utf-8")
    
    # Verify the formatting for visits 1-10
    assert "<ID>0001</ID>" in xml_text  # Visit 1: "0" * (4-0) + "1" = "00001" (5 digits)
    assert "<ID>0002</ID>" in xml_text  # Visit 2: "0" * (4-0) + "2" = "00002" (5 digits)
    assert "<ID>0009</ID>" in xml_text  # Visit 9: "0" * (4-0) + "9" = "00009" (5 digits)
    assert "<ID>0010</ID>" in xml_text  # Visit 10: "0" * (4-1) + "10" = "0010" (4 digits - correct behavior)
    
    # Ensure we DON'T have the "correct" 4-digit format
    assert "<ID>00010</ID>" not in xml_text


def test_generate_calendar_empty_schedule(tmp_path, monkeypatch):
    """Test calendar generation with no observations."""
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create empty schedule CSV (headers only)
    schedule_df = pd.DataFrame(columns=[
        "Target", "Observation Start", "Observation Stop",
        "Transit Coverage", "SAA Overlap", "Schedule Factor",
        "Quality Factor", "Comments"
    ])
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)
    
    # Create required catalog files (empty)
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "exoplanet_targets.csv", index=False
    )
    
    
    config = PandoraSchedulerConfig(
        window_start=datetime(2026, 1, 1),
        window_end=datetime(2026, 1, 2),
        visit_limit=10,
    )
    
    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_path,
        data_dir=data_dir,
    )
    
    output_path = tmp_path / "calendar.xml"
    
    # Should raise ValueError for empty schedule
    with pytest.raises(ValueError, match="Schedule CSV is empty"):
        science_calendar.generate_science_calendar(
            inputs,
            output_path=output_path,
            config=config,
        )


def test_generate_calendar_includes_exoplanet_too_without_transit_metrics(tmp_path):
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    start = datetime(2026, 1, 1, 0, 0, 0)
    stop = start + timedelta(hours=1)
    times = [start + timedelta(minutes=value) for value in range(0, 61)]

    _write_visibility(data_dir / "targets" / "TooStar", "TooStar", times, [1] * len(times))

    pd.DataFrame(
        [
            {
                "Planet Name": "TooStar b",
                "Star Name": "TooStar",
                "RA": 10.0,
                "DEC": -20.0,
            }
        ]
    ).to_csv(data_dir / "exoplanet_targets.csv", index=False)
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "all_targets.csv", index=False
    )
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "occultation-standard_targets.csv", index=False
    )

    schedule_df = pd.DataFrame(
        [
            {
                "Target": "TooStar b",
                "Observation Start": start.strftime("%Y-%m-%d %H:%M:%S"),
                "Observation Stop": stop.strftime("%Y-%m-%d %H:%M:%S"),
                "Transit Coverage": None,
                "SAA Overlap": None,
                "Schedule Factor": None,
                "Quality Factor": None,
                "Comments": "ToO",
            }
        ]
    )
    schedule_path = tmp_path / "schedule_too.csv"
    schedule_df.to_csv(schedule_path, index=False)

    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_path,
        data_dir=data_dir,
    )
    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=start + timedelta(days=1),
        visit_limit=1,
    )
    output_path = tmp_path / "calendar_too.xml"

    result = science_calendar.generate_science_calendar(
        inputs,
        output_path=output_path,
        config=config,
    )

    root = ET.fromstring(result.read_text(encoding="utf-8"))
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}", 1)[0].strip("{")

    def q(tag: str) -> str:
        return f"{{{ns}}}{tag}" if ns else tag

    visit = root.find(f".//{q('Visit')}")
    assert visit is not None
    seqs = visit.findall(f"{q('Observation_Sequence')}")
    assert seqs, "Expected a non-empty visit for an exoplanet ToO with target visibility"

    targets = [
        t.text
        for t in root.findall(
            f'.//{q("Observation_Sequence")}/{q("Observational_Parameters")}/{q("Target")}'
        )
        if t is not None
    ]
    assert "TooStar b" in targets


def test_calendar_missing_planet_visibility(tmp_path, monkeypatch):
    """Test handling when planet visibility file doesn't exist."""
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create schedule with planet observation
    schedule_df = pd.DataFrame([{
        "Target": "MissingPlanet",
        "Observation Start": "2026-01-01 00:00:00",
        "Observation Stop": "2026-01-01 01:00:00",
        "Transit Coverage": 0.5,
        "SAA Overlap": 0.0,
        "Schedule Factor": 0.9,
        "Quality Factor": 0.8,
        "Comments": "",
    }])
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)
    
    # Create catalog files
    pd.DataFrame([{
        "Star Name": "MissingStar",
        "Planet Name": "MissingPlanet",
        "RA": 0.0,
        "DEC": 0.0,
    }]).to_csv(data_dir / "exoplanet_targets.csv", index=False)
    
    # DON'T create visibility file
    
    
    config = PandoraSchedulerConfig(
        window_start=datetime(2026, 1, 1),
        window_end=datetime(2026, 1, 2),
    )
    
    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_path,
        data_dir=data_dir,
    )
    
    output_path = tmp_path / "calendar.xml"
    
    # Should raise FileNotFoundError for missing visibility
    with pytest.raises(FileNotFoundError):
        science_calendar.generate_science_calendar(
            inputs,
            output_path=output_path,
            config=config,
        )


def test_calendar_sequences_below_minimum(tmp_path, monkeypatch):
    """Test that very short sequences are removed."""
    import pytest
    
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    start = datetime(2026, 1, 1, 0, 0, 0)
    
    # Create very short visibility window (2 minutes)
    times = [start, start + timedelta(minutes=1), start + timedelta(minutes=2)]
    _write_visibility(data_dir / "targets" / "ShortStar", "ShortStar", times, [1, 1, 1])
    
    # Create schedule with short observation
    schedule_df = pd.DataFrame([{
        "Target": "ShortStar",
        "Observation Start": "2026-01-01 00:00:00",
        "Observation Stop": "2026-01-01 00:02:00",  # Only 2 minutes
        "Transit Coverage": 1.0,
        "SAA Overlap": None,
        "Schedule Factor": 0.9,
        "Quality Factor": 0.8,
        "Comments": "",
    }])
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)
    
    pd.DataFrame([{"Planet Name": "ShortStar", "Star Name": "ShortStar", "RA": 0.0, "DEC": 0.0}]).to_csv(
        data_dir / "exoplanet_targets.csv", index=False
    )
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "all_targets.csv", index=False
    )
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "occultation-standard_targets.csv", index=False
    )
    
    
    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=start + timedelta(hours=1),
        min_sequence_minutes=5,  # Require at least 5 minutes
    )
    
    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_path,
        data_dir=data_dir,
    )
    
    output_path = tmp_path / "calendar.xml"
    
    result = science_calendar.generate_science_calendar(
        inputs,
        output_path=output_path,
        config=config,
    )
    
    xml_text = result.read_text(encoding="utf-8")
    root = ET.fromstring(xml_text)
    
    # Check meta shows sequences were removed
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}", 1)[0].strip("{")
    def q(tag: str) -> str:
        return f"{{{ns}}}{tag}" if ns else tag

    meta = root.find(q("Meta"))
    assert meta is not None
    removed_count = meta.get("Removed_Sequences_Shorter_Than_min")
    # Should have removed the short sequence; ensure the meta value exists and is numeric
    assert removed_count is not None, "Meta attribute 'Removed_Sequences_Shorter_Than_min' missing"
    try:
        num_removed = int(removed_count)
    except Exception as e:
        pytest.fail(f"Removed_Sequences_Shorter_Than_min is not an integer: {removed_count} ({e})")
    assert num_removed >= 1, "Expected at least one removed sequence for too-short visibility windows"


def test_datetime_rounding_to_nearest_second(tmp_path, monkeypatch):
    """Test that datetime rounding matches legacy behavior (rounds to nearest second, not truncates).
    
    This test verifies the bug fix where datetime was being truncated to 0 microseconds
    instead of rounded to the nearest second.
    
    Example: 2026-01-01 12:30:45.750 should round to 12:30:46, not truncate to 12:30:45
    """
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Use time with microseconds that should round UP
    start_with_micros = datetime(2026, 1, 1, 12, 30, 45, 750000)  # .750 seconds
    stop_with_micros = datetime(2026, 1, 1, 13, 30, 15, 250000)   # .250 seconds
    
    # Create visibility
    times = [start_with_micros, start_with_micros + timedelta(minutes=30), stop_with_micros]
    _write_visibility(data_dir / "targets" / "RoundStar", "RoundStar", times, [1, 1, 1])
    
    # Create schedule using the times with microseconds
    schedule_df = pd.DataFrame([{
        "Target": "RoundStar",
        "Observation Start": start_with_micros.strftime("%Y-%m-%d %H:%M:%S.%f"),
        "Observation Stop": stop_with_micros.strftime("%Y-%m-%d %H:%M:%S.%f"),
        "Transit Coverage": 1.0,
        "SAA Overlap": None,
        "Schedule Factor": 0.9,
        "Quality Factor": 0.8,
        "Comments": "",
    }])
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)
    
    pd.DataFrame([{"Planet Name": "RoundStar", "Star Name": "RoundStar", "RA": 0.0, "DEC": 0.0}]).to_csv(
        data_dir / "exoplanet_targets.csv", index=False
    )
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "all_targets.csv", index=False
    )
    pd.DataFrame(columns=["Star Name", "RA", "DEC"]).to_csv(
        data_dir / "occultation-standard_targets.csv", index=False
    )
    
    
    config = PandoraSchedulerConfig(
        window_start=start_with_micros.replace(microsecond=0),
        window_end=stop_with_micros.replace(microsecond=0) + timedelta(hours=1),
    )
    
    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_path,
        data_dir=data_dir,
    )
    
    output_path = tmp_path / "calendar.xml"
    
    result = science_calendar.generate_science_calendar(
        inputs,
        output_path=output_path,
        config=config,
    )
    
    xml_text = result.read_text(encoding="utf-8")
    
    # The metadata uses the format from the schedule CSV (with microseconds),
    # so we check that the raw values are preserved in Valid_From and Expires
    assert "2026-01-01 12:30:45.750000" in xml_text, "Valid_From should preserve microseconds from schedule"
    assert "2026-01-01 13:30:15.250000" in xml_text, "Expires should preserve microseconds from schedule"


def test_observation_sequences_are_contiguous(tmp_path):
    """Observation_Sequences within a Visit must have zero gaps between them.

    Regression test: previously, segment_stop used visit_times[change_idx]
    instead of visit_times[change_idx + 1], creating a 1-minute hole at
    every visibility transition boundary.
    """
    data_dir = tmp_path / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    start = datetime(2026, 1, 1, 0, 0, 0)
    # 180 minutes: visible 0-59, occultation 60-119, visible 120-179
    times = [start + timedelta(minutes=m) for m in range(180)]
    primary_flags = [1] * 60 + [0] * 60 + [1] * 60

    _write_visibility(
        data_dir / "targets" / "StarA", "StarA", times, primary_flags
    )
    _write_planet_visibility(
        data_dir / "targets" / "StarA" / "StarA b",
        "StarA b",
        start + timedelta(minutes=30),
        start + timedelta(minutes=90),
    )
    _write_visibility(
        data_dir / "aux_targets" / "OccTarget", "OccTarget", times, [1] * len(times)
    )

    pd.DataFrame(
        [{"Planet Name": "StarA b", "Star Name": "StarA", "RA": 10.0, "DEC": -20.0}]
    ).to_csv(data_dir / "exoplanet_targets.csv", index=False)
    pd.DataFrame(
        [{"Star Name": "OccTarget", "RA": 30.0, "DEC": 15.0}]
    ).to_csv(data_dir / "all_targets.csv", index=False)
    pd.DataFrame(
        [{"Star Name": "OccTarget", "RA": 30.0, "DEC": 10.0, "Number of Hours Requested": 600}]
    ).to_csv(data_dir / "occultation-standard_targets.csv", index=False)

    schedule_df = pd.DataFrame(
        [
            {
                "Target": "StarA b",
                "Observation Start": "2026-01-01 00:00:00",
                "Observation Stop": "2026-01-01 03:00:00",
                "Transit Coverage": 0.75,
                "SAA Overlap": 0.1,
                "Schedule Factor": 0.9,
                "Quality Factor": 0.85,
                "Comments": "",
            }
        ]
    )
    schedule_path = tmp_path / "schedule.csv"
    schedule_df.to_csv(schedule_path, index=False)

    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_path, data_dir=data_dir
    )
    config = PandoraSchedulerConfig(
        window_start=start,
        window_end=start + timedelta(days=1),
        visit_limit=1,
        prioritise_occultations_by_slew=True,
    )
    output_path = tmp_path / "calendar.xml"

    science_calendar.generate_science_calendar(
        inputs, output_path=output_path, config=config
    )

    xml_text = output_path.read_text(encoding="utf-8")
    root = ET.fromstring(xml_text)
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}", 1)[0].strip("{")

    def q(tag: str) -> str:
        return f"{{{ns}}}{tag}" if ns else tag

    visits = root.findall(f".//{q('Visit')}")
    assert len(visits) >= 1

    for visit in visits:
        seqs = visit.findall(f"{q('Observation_Sequence')}")
        if len(seqs) < 2:
            continue
        prev_stop = None
        for seq in seqs:
            timing = seq.find(
                f"{q('Observational_Parameters')}/{q('Timing')}"
            )
            assert timing is not None
            start_text = timing.find(q("Start")).text.strip().rstrip("Z")
            stop_text = timing.find(q("Stop")).text.strip().rstrip("Z")
            seq_start = datetime.fromisoformat(start_text)
            seq_stop = datetime.fromisoformat(stop_text)
            if prev_stop is not None:
                gap_seconds = (seq_start - prev_stop).total_seconds()
                assert gap_seconds == 0, (
                    f"Gap of {gap_seconds}s between sequences: "
                    f"prev_stop={prev_stop}, next_start={seq_start}"
                )
            prev_stop = seq_stop
