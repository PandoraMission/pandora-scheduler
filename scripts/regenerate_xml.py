#!/usr/bin/env python3
"""Regenerate science calendar XML from existing schedule CSV."""

from pathlib import Path

from pandorascheduler_rework.science_calendar import (
    ScienceCalendarInputs,
    generate_science_calendar,
)

inputs = ScienceCalendarInputs(
    schedule_csv=Path("output_standalone/Pandora_Schedule_0.8_0.0_0.2_2026-02-05_to_2027-02-05.csv"),
    data_dir=Path("output_standalone/data"),
)

print("Regenerating science calendar XML from existing schedule CSV...")
xml_path = generate_science_calendar(inputs)
print(f"✅ XML generated: {xml_path}")
