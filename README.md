# Pandora Scheduler (rework)

Brief overview and quick links for developers working on the rework.

- See `QUICK_START.md` for runnable examples and common workflows.
- Example JSON configuration: `example_scheduler_config.json` (root) — use with `--config`.
- Detailed keys: `docs/EXAMPLE_SCHEDULER_CONFIG.md`.

## Quick Start

```bash
# 1) Install dependencies (poetry environment assumed)
poetry install

# 2) Run a quick test (assumes manifests/visibility already present)
poetry run python run_scheduler.py \
    --start "2026-02-05" \
    --end "2026-02-07" \
    --output ./output_test

# 3) Full pipeline with target definitions
poetry run python run_scheduler.py \
    --start "2026-02-05" \
    --end "2026-02-12" \
    --output ./output \
    --target-definitions /path/to/PandoraTargetList/target_definition_files \
    --generate-visibility \
    --gmat-ephemeris /path/to/ephemeris.txt \
    --show-progress
```

## Run The Pipeline

With the example JSON config:

```bash
poetry run python run_scheduler.py \
    --start "2026-04-01" \
    --end "2026-04-05" \
    --output output_standalone_test7_short \
    --config example_scheduler_config.json
```

This writes:

- schedule CSVs under `output_standalone_test7_short/`
- visibility parquet files under `output_standalone_test7_short/data_*`
- science calendar XML at `output_standalone_test7_short/Pandora_science_calendar.xml` unless `skip_xml` is `true`

If you want the visualizer to run automatically after the pipeline:

```json
"skip_xml": false,
"run_visualizer_after_pipeline": true,
"visualizer_mode": "priority"
```

## XML-Only From Existing Schedule
------------------------------
If you already have a schedule CSV and only want to regenerate the science calendar XML, you can skip the scheduling pipeline entirely:

```bash
poetry run python run_scheduler.py \
  --schedule-csv output_test1_short/Pandora_Schedule_0.8_0.0_0.2_2026-04-01_to_2026-04-15.csv \
  --config example_scheduler_config.json
```

Notes:
- `--start` / `--end` are optional in this mode; they are inferred from the schedule CSV if omitted.
- `--output` is optional in this mode; it defaults to the schedule CSV parent directory.
- The XML is still written to `output_*/Pandora_science_calendar.xml`.

Generate Visualizations

All visualization commands read the science calendar XML:

```bash
poetry run python scripts/visualizer.py \
    output_standalone_test7_short/Pandora_science_calendar.xml \
    --mode priority \
    --out output_standalone_test7_short/visualizer_priority.png
```

Available modes:

- `priority`: main Gantt-style plot colored by sequence priority

```bash
poetry run python scripts/visualizer.py \
    output_standalone_test7_short/Pandora_science_calendar.xml \
    --mode priority \
    --out output_standalone_test7_short/visualizer_priority.png
```

- `timeline`: simple chronological timeline

```bash
poetry run python scripts/visualizer.py \
    output_standalone_test7_short/Pandora_science_calendar.xml \
    --mode timeline \
    --out output_standalone_test7_short/visualizer_timeline.png
```

- `target-time`: total scheduled time per target

```bash
poetry run python scripts/visualizer.py \
    output_standalone_test7_short/Pandora_science_calendar.xml \
    --mode target-time \
    --out output_standalone_test7_short/visualizer_target_time.png
```

- `simple`: lighter-weight priority timeline

```bash
poetry run python scripts/visualizer.py \
    output_standalone_test7_short/Pandora_science_calendar.xml \
    --mode simple \
    --out output_standalone_test7_short/visualizer_simple.png
```

- `visibility`: priority plot with non-visible intervals overlaid from the run's parquet visibility files

```bash
poetry run python scripts/visualizer.py \
    output_standalone_test7_short/Pandora_science_calendar.xml \
    --mode visibility \
    --data-dir output_standalone_test7_short/data_91_20_115 \
    --out output_standalone_test7_short/visualizer_visibility.png
```

Optional flag for supported modes:

```bash
--show-sequence-labels
```

## Flowchart Diagrams

Generated reference diagrams live under `docs/`:

- [docs/pipeline_flowchart.png](/Users/vkostov/Documents/GitHub/pandora-scheduler/docs/pipeline_flowchart.png): end-to-end pipeline overview
- [docs/visibility_checks_flowchart.png](/Users/vkostov/Documents/GitHub/pandora-scheduler/docs/visibility_checks_flowchart.png): detailed visibility-generation and keepout logic
- [docs/earth_daynight_logic_flowchart.png](/Users/vkostov/Documents/GitHub/pandora-scheduler/docs/earth_daynight_logic_flowchart.png): focused Earth day/night threshold selection
- [docs/xml_segment_assignment_flowchart.png](/Users/vkostov/Documents/GitHub/pandora-scheduler/docs/xml_segment_assignment_flowchart.png): focused science-vs-occultation XML segment assignment

Renderer scripts live under `scripts/` if you want to regenerate them:

- [scripts/render_pipeline_flowchart.py](/Users/vkostov/Documents/GitHub/pandora-scheduler/scripts/render_pipeline_flowchart.py)
- [scripts/render_visibility_flowchart.py](/Users/vkostov/Documents/GitHub/pandora-scheduler/scripts/render_visibility_flowchart.py)
- [scripts/render_earth_daynight_flowchart.py](/Users/vkostov/Documents/GitHub/pandora-scheduler/scripts/render_earth_daynight_flowchart.py)
- [scripts/render_xml_segment_flowchart.py](/Users/vkostov/Documents/GitHub/pandora-scheduler/scripts/render_xml_segment_flowchart.py)

## Config Notes

- `daynight_mode` accepts two values:
  - `subsatellite`: classify day/night from whether the subsatellite point is sunlit
  - `limb`: classify day/night from whether the nearest Earth limb point in the target direction is sunlit
- `generate_occultation_xml` can be set in the JSON config to turn occultation filling in the science calendar XML on or off.
- `min_science_sequence_minutes` sets the minimum standalone science-visible fragment. Shorter science fragments are merged into a contiguous preceding science chunk when possible; otherwise they are handed to occultation filling.
- `min_occultation_sequence_minutes` sets the minimum standalone occultation tail. Short trailing occultation chunks are only absorbed into the preceding occultation chunk when the same occultation target remains visible there.
- `requested_occ_time_override` can be set in the JSON config to allow occultation scheduling to continue when requested-hours bookkeeping is incomplete or would otherwise block assignment.
- `run_visualizer_after_pipeline` can be set in the JSON config to generate a plot automatically after the XML is written.
- `visualizer_mode` accepts:
  - `priority`: main Gantt-style plot colored by sequence priority
  - `timeline`: simple chronological timeline
  - `target-time`: bar chart of total observation time per target
  - `simple`: lighter-weight priority timeline
  - `visibility`: priority plot with non-visible intervals overlaid from the run's visibility parquet files

If you need help, read `QUICK_START.md` for examples and troubleshooting tips.
