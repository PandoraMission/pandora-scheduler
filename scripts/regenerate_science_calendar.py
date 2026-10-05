import logging
import sys
from datetime import datetime
from pathlib import Path

from pandorascheduler_rework import science_calendar
from pandorascheduler_rework.config import PandoraSchedulerConfig


def regenerate_calendar(output_dir_str: str, use_legacy_mode: bool = False):
    output_dir = Path(output_dir_str).resolve()
    # Find the schedule CSV
    schedule_csv = list(output_dir.glob("Pandora_Schedule_*.csv"))[0]
    data_dir = output_dir / "data"
    
    print("Regenerating Science Calendar...")
    print(f"  Schedule: {schedule_csv}")
    print(f"  Data Dir: {data_dir}")
    print(f"  Output:   {output_dir}")
    print(f"  Legacy mode: {use_legacy_mode}")

    if not data_dir.exists():
        print(f"ERROR: Data directory not found at {data_dir}")
        sys.exit(1)

    logging.basicConfig(level=logging.INFO)

    inputs = science_calendar.ScienceCalendarInputs(
        schedule_csv=schedule_csv,
        data_dir=data_dir,
    )

    # Dummy dates - XML generation extracts dates from the schedule CSV
    config = PandoraSchedulerConfig(
        window_start=datetime(2026, 2, 1),
        window_end=datetime(2027, 2, 2),
        use_legacy_mode=use_legacy_mode,
    )

    output_path = data_dir / "Pandora_science_calendar.xml"

    science_calendar.generate_science_calendar(
        inputs,
        config=config,
        output_path=output_path,
    )
    print(f"SUCCESS: Generated {output_path}")

if __name__ == "__main__":
    use_legacy = "--legacy-mode" in sys.argv or "-l" in sys.argv
    regenerate_calendar("output_standalone_test", use_legacy_mode=use_legacy)
