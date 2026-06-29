# Pandora Scheduler

This repository builds Pandora scheduling products from target manifests, target
definition files, and GMAT ephemerides.

Main outputs:
- `Pandora_Schedule_*.csv`
- `Pandora_science_calendar.xml`
- `run_config_manifest.json`
- `data_*/targets/.../Visibility for *.parquet`

For detailed JSON key descriptions, see `docs/EXAMPLE_SCHEDULER_CONFIG.md`.

## Prerequisites

From the repository root:

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler
poetry install
```

If you prefer not to use `poetry run`, the commands below also work with:

```bash
.venv/bin/python ...
```

The examples below use `.venv/bin/python` explicitly.

## Common Inputs

Main pipeline runs typically need:
- a JSON config, for example `scheduler_config_20260622_20260629_with_too.json`
- a target-definition base, usually provided in the JSON as:
  `extra_inputs.target_definition_base`
- a GMAT ephemeris, usually provided in the JSON as:
  `extra_inputs.visibility_gmat`

Common output layout:
- `output_*/Pandora_Schedule_*.csv`
- `output_*/Pandora_science_calendar.xml`
- `output_*/run_config_manifest.json`
- `output_*/data_<sun>_<moon>_<earth>/`

## Run Modes

### 1. Full Pipeline From A Config

This is the normal end-to-end run: manifests, visibility, schedule, XML, and
reports.

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_scheduler.py \
  --start "2026-06-22 00:00:00" \
  --end "2026-06-29 00:00:00" \
  --output output_20260622_20260629_EarthDay111 \
  --config scheduler_config_20260622_20260629_with_too.json
```
If the above doesn't work, try with: 
python3 run_scheduler.py \
  --start "2026-07-06 00:00:00" \
  --end "2026-07-13 00:00:00" \
  --output output_20260706_20260713_EarthDay111 \
  --config scheduler_config_20260706_20260713_with_too.json

Use this mode whenever the output directory does not already contain the
schedule CSV and visibility products you need.

### 2. Full Pipeline With Explicit Target Definitions And Visibility Generation

Use this when you want the command line to supply the target-definition base and
GMAT ephemeris directly instead of relying on the JSON.

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_scheduler.py \
  --start "2026-06-29 00:00:00" \
  --end "2026-07-06 00:00:00" \
  --output output_20260629_20260706_EarthDay111 \
  --config scheduler_config.json \
  --target-definitions /Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files \
  --generate-visibility \
  --gmat-ephemeris /Users/vkostov/Documents/GitHub/pandora-scheduler/gmat/PAN-GMAT-COM-20260610-VF20260610-EX20270710.txt
```

### 3. XML Only From An Existing Schedule CSV

Use this when the schedule already exists and you only want to regenerate the
science calendar XML.

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_scheduler.py \
  --schedule-csv output_20260622_20260629_EarthDay111/Pandora_Schedule_0.8_0.0_0.2_2026-06-22_to_2026-06-29.csv \
  --xml-data-dir output_20260622_20260629_EarthDay111/data_91_20_111 \
  --output output_20260622_20260629_EarthDay111 \
  --config scheduler_config_20260622_20260629_with_too.json
```

Notes:
- `--start` and `--end` are optional in this mode.
- `--output` defaults to the schedule CSV parent if omitted.
- `--xml-data-dir` should point at the matching `data_*` directory used to
  build the schedule.

### 4. XML Only For A Slice Of The Schedule CSV

Useful for debugging a subset of rows.

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_scheduler.py \
  --schedule-csv output_20260622_20260629_EarthDay111/Pandora_Schedule_0.8_0.0_0.2_2026-06-22_to_2026-06-29.csv \
  --schedule-row-start 1 \
  --schedule-row-end 20 \
  --xml-data-dir output_20260622_20260629_EarthDay111/data_91_20_111 \
  --output output_20260622_20260629_EarthDay111 \
  --config scheduler_config_20260622_20260629_with_too.json
```

### 5. Visualizer Only For An Existing XML

Use this when XML already exists and you want the plot without rerunning the
pipeline.

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python scripts/visualizer.py \
  output_20260622_20260629_EarthDay111/Pandora_science_calendar.xml \
  --mode visibility \
  --data-dir output_20260622_20260629_EarthDay111/data_91_20_111 \
  --out output_20260622_20260629_EarthDay111/visualizer_visibility.png
```

Supported `--mode` values:
- `priority`
- `timeline`
- `target-time`
- `simple`
- `visibility`

### 6. Auto-Run The Visualizer As Part Of The Pipeline

Set these keys in the JSON config:

```json
{
  "generate_xml": true,
  "run_visualizer_after_pipeline": true,
  "visualizer_mode": "visibility"
}
```

Then run the normal full pipeline command:

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_scheduler.py \
  --start "2026-06-22 00:00:00" \
  --end "2026-06-29 00:00:00" \
  --output output_20260622_20260629_EarthDay111 \
  --config scheduler_config_20260622_20260629_with_too.json
```

### 7. Validate XML Against Visibility Products

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python scripts/validate_xml_visibility.py \
  output_20260622_20260629_EarthDay111/Pandora_science_calendar.xml \
  --data-dir output_20260622_20260629_EarthDay111/data_91_20_111
```

This writes:
- `output_20260622_20260629_EarthDay111/Pandora_science_calendar_visibility_validation.csv`

### 8. Export Minute-By-Minute Visit Visibility Diagnostics

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python scripts/export_visit_visibility_diagnostics.py \
  --start "2026-06-04 15:30:00" \
  --stop "2026-06-07 15:30:00" \
  --target-name DS_Tuc_Ab \
  --schedule-csv output_20260601_20260608_EarthDay111/Pandora_Schedule_0.8_0.0_0.2_2026-06-01_to_2026-06-08.csv \
  --config scheduler_config_20260601_20260608_with_too.json \
  --output output_20260601_20260608_EarthDay111/DS_Tuc_Ab_2026-06-04_1530_2026-06-07_1530_visibility_diagnostics.csv
```

Use this for minute-by-minute Sun, Moon, Earth, and star-tracker keepout
diagnostics for a visit window.

## ToO Lists

To inject explicit ToO windows into the main scheduler pipeline, set:

```json
{
  "extra_inputs": {
    "too_list_csv": "/absolute/path/to/ToO_list.csv"
  }
}
```

Required CSV columns:
- `Target`
- `Obs Window Start`
- `Obs Window Stop`

Then run the normal full pipeline command:

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_scheduler.py \
  --start "2026-06-22 00:00:00" \
  --end "2026-06-29 00:00:00" \
  --output output_20260622_20260629_EarthDay111 \
  --config scheduler_config_20260622_20260629_with_too.json
```

Notes:
- ToO target names must match `Planet Name` in the active
  `exoplanet_targets.csv`.
- On rerun, the configured source file overwrites `output_*/ToO_list.csv`.

## Experimental 10 HJs + ToOs Wrapper

This is separate from the main `run_scheduler.py` pipeline.

```bash
cd /Users/vkostov/Documents/GitHub/pandora-scheduler

.venv/bin/python run_10_hjs_toos.py \
  --start "2026-04-27 00:00:00" \
  --end "2026-05-04 00:00:00" \
  --output output_1_week_10_HJs_5_ToOs \
  --config scheduler_config_10_HJs_5_ToOs.json
```

## Files To Read Next

If you need more detail:
- `QUICK_START.md`
- `docs/EXAMPLE_SCHEDULER_CONFIG.md`
- `scripts/export_visit_visibility_diagnostics.py`
- `scripts/validate_xml_visibility.py`
