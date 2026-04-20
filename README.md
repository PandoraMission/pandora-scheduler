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
- science calendar XML at `output_standalone_test7_short/Pandora_science_calendar.xml` when `generate_xml` is `true`
- run configuration manifest at `output_standalone_test7_short/run_config_manifest.json`

Notes:
- If you pass `--start`, `--end`, and `--output` explicitly, those force a full pipeline run even if the JSON also contains `schedule_csv`.
- If you keep `schedule_csv_row_limit` in the JSON for XML-only testing, set it to `null` before a full pipeline run unless you intentionally want to truncate the XML build.

If you want the visualizer to run automatically after the pipeline:

```json
"generate_xml": true,
"run_visualizer_after_pipeline": true,
"visualizer_mode": "priority"
```

## Select ToO Targets

The scheduler can now select fixed-window exoplanet Targets of Opportunity
inside a normal pipeline run. The main science target list remains `10_HJs/`,
while `extra_inputs.too_candidate_target_definition_base` points at the full
Pandora exoplanet target list used only for visibility generation and ToO
selection.

```bash
poetry run python run_scheduler.py \
  --start "2026-04-27" \
  --end "2026-05-04" \
  --output output_1_week_10_HJs_5_ToOs \
  --config example_scheduler_config_10_HJs_5_ToOs.json
```

The run writes all products under the same output data directory, for example
`output_1_week_10_HJs_5_ToOs/data_95_30_115/`:

- `exoplanet_targets.csv`: 10 HJs plus the selected ToO targets appended with
  `Number of Transits to Capture = 0`
- `all_exoplanet_targets.csv`: full exoplanet candidate manifest
- `targets/`: visibility products for the full exoplanet candidate list
- `all_too_candidates_ranked_by_transit_depth.csv`: all non-10-HJ candidates
  with visible transits, ranked by transit depth
- `all_too_windows_ranked_by_transit_depth.csv`: every eligible candidate
  transit window
- `ToO_list.csv`: top selected fixed ToO windows used by the scheduler
- `transit_depth_cache.csv`: cached NASA Exoplanet Archive depth lookups

Configure the selection in `example_scheduler_config_10_HJs_5_ToOs.json`:

```json
"auto_select_toos": true,
"too_candidate_target_definition_base": "/Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files",
"too_top_n": 5,
"too_depth_min_percent": 1.0,
"too_transit_coverage_min": 0.4,
"too_query_depths": true
```

## Run 10 HJs Plus 5 Visible Science Targets

As an alternative to fixed-window ToOs, `10_HJs_plus_5_visible/` treats the five
selected full-list exoplanets as ordinary science targets. Its priority file
puts those five first, ranked by best visible transit coverage, followed by the
10 HJs in their existing priority order:

```bash
poetry run python run_scheduler.py \
  --start "2026-04-27" \
  --end "2026-05-04" \
  --output output_1_week_10_HJs_plus_5_visible \
  --config example_scheduler_config_10_HJs_plus_5_visible.json
```

Use a fresh output directory for this run so no previous `ToO_list.csv` is
picked up.

## What's New On `soft_ST_test`

This branch adds a few XML-generation and provenance improvements that are helpful for science-tail experiments:

- science soft-ST tail extension for XML sequence generation only
  - applies only to science-visible segments
  - keeps hard boresight limits (`sun_avoidance_deg`, `moon_avoidance_deg`, `earth_avoidance_*`) unchanged
  - can relax or disable star-tracker constraints only in the final `science_soft_startracker_tail_minutes`
- one-run XML A/B generation when `allow_science_soft_startracker_tail` is `true`
  - baseline outputs:
    - `Pandora_science_calendar.xml`
    - `Pandora_science_calendar_sequence_provenance.csv`
  - soft-ST outputs:
    - `Pandora_science_calendar_soft_ST.xml`
    - `Pandora_science_calendar_soft_ST_sequence_provenance.csv`
- root-level run manifest
  - every run writes `output_*/run_config_manifest.json`
- explicit soft-tail provenance on soft-ST runs
  - `science_soft_tail_used`
  - `science_soft_tail_minutes`
  - these columns are omitted from the baseline provenance CSV
- dedicated plotting helper for soft-ST provenance comparisons:
  - `scripts/plot_soft_st_provenance_comparison.py`
- non-exoplanet target manifests no longer inject fake `Planet Name` values into monitoring/occultation standard catalogs

## XML-Only From Existing Schedule

If you already have a schedule CSV and only want to regenerate the science calendar XML, you can skip the scheduling pipeline entirely:

```bash
poetry run python run_scheduler.py \
  --schedule-csv output_test1_short/Pandora_Schedule_0.8_0.0_0.2_2026-04-01_to_2026-04-15.csv \
  --config example_scheduler_config.json
```

Notes:
- `--start` / `--end` are optional in this mode; they are inferred from the schedule CSV if omitted.
- `--output` is optional in this mode; it defaults to the schedule CSV parent directory.
- The XML is written to `output_*/Pandora_science_calendar.xml`.
- If `allow_science_soft_startracker_tail` is `true`, the same run also writes:
  - `output_*/Pandora_science_calendar_soft_ST.xml`
  - `output_*/Pandora_science_calendar_soft_ST_sequence_provenance.csv`
- If you only want to test the XML builder on the first few schedule rows, add `--schedule-row-limit N`.

Example:

```bash
poetry run python run_scheduler.py \
  --schedule-csv output_test1_short/Pandora_Schedule_0.8_0.0_0.2_2026-04-01_to_2026-04-15.csv \
  --schedule-row-limit 10 \
  --config example_scheduler_config.json
```

## Validate XML Visibility

To validate the final XML against the run's visibility parquet files:

```bash
python3 scripts/validate_xml_visibility.py \
  output_test1_short/Pandora_science_calendar.xml \
  --data-dir output_test1_short/data_91_20_115
```

This writes:

- `output_test1_short/Pandora_science_calendar_visibility_validation.csv`

If present, the validator also merges sequence-level provenance from:

- `output_test1_short/Pandora_science_calendar_sequence_provenance.csv`

The validation CSV includes:

- sequence identity (`visit_id`, `sequence_id`, `target`)
- schedule timing (`start_utc`, `stop_utc`, `duration_minutes`)
- visibility check results (`status`, `visible_minutes`, `non_visible_minutes`)
- provenance fields (`sequence_type`, `occultation_pass`, `sequence_visibility_fraction`)

The provenance CSV is written automatically during XML generation and records one row per emitted XML sequence.

When `allow_science_soft_startracker_tail` is `true`, the soft-ST provenance CSV also includes:

- `science_soft_tail_used`
- `science_soft_tail_minutes`

These columns are intentionally omitted from the baseline provenance CSV.

## Plot Soft-ST Provenance

To visualise only the science rows affected by soft-ST tail extension:

```bash
MPLCONFIGDIR=/tmp/mplcache python3 scripts/plot_soft_st_provenance_comparison.py \
  --soft-provenance output_test1/Pandora_science_calendar_soft_ST_sequence_provenance.csv \
  --out output_test1/soft_ST_provenance_comparison.png
```

To generate the older two-panel baseline-vs-soft comparison:

```bash
MPLCONFIGDIR=/tmp/mplcache python3 scripts/plot_soft_st_provenance_comparison.py \
  --base-provenance output_test1/Pandora_science_calendar_sequence_provenance.csv \
  --soft-provenance output_test1/Pandora_science_calendar_soft_ST_sequence_provenance.csv \
  --layout two-panel \
  --out output_test1/soft_ST_provenance_comparison_two_panel.png
```

To compare `science_soft_st_required = 0, 1, 2` side by side:

```bash
MPLCONFIGDIR=/tmp/mplcache python3 scripts/plot_soft_st_required_comparison.py \
  --provenance-0 output_test1/Pandora_science_calendar_soft_ST_sequence_provenance_0.csv \
  --provenance-1 output_test1/Pandora_science_calendar_soft_ST_sequence_provenance_1.csv \
  --provenance-2 output_test1/Pandora_science_calendar_soft_ST_sequence_provenance_2.csv \
  --out output_test1/soft_ST_required_comparison.png
```

## Debug a Visit

If you want to inspect how occultation scheduling behaved for one specific XML visit, use the visit debug helper:

```bash
poetry run python scripts/debug_occultation_visit.py \
  --config example_scheduler_config.json \
  --schedule-csv output_test_6months/Pandora_Schedule_0.8_0.0_0.2_2026-04-01_to_2026-07-01.csv \
  --data-dir output_test1_short/data_91_20_115 \
  --visit-id 0145 \
  --validation-csv output_test_6months/Pandora_science_calendar_visibility_validation.csv
```

This prints:
- the visit window
- the science/occultation segments seen by the XML builder
- the occultation intervals passed into `schedule_occultation_targets(...)`
- the raw `occ_df` returned for that visit
- any validation failures for that visit

This is useful when a visit mixes:
- `scheduled_occultation` rows from Pass 1/2/3/4
- `catalog_fallback` rows emitted later by the XML builder

## Generate Visualizations

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
- `generate_xml` turns XML generation on or off.
- `include_occultation_sequences_in_xml` controls whether occultation filling is included in the science calendar XML.
- `allow_science_soft_startracker_tail` enables the science-tail XML A/B mode and writes both baseline and `_soft_ST` XML/provenance outputs in one run.
- `science_soft_startracker_tail_minutes` sets the maximum science-tail extension window.
- `science_soft_st_required` controls how many star trackers must pass inside the soft tail:
  - `2`: both trackers
  - `1`: at least one tracker
  - `0`: disable ST constraints for the tail and use hard boresight limits only
- `min_science_sequence_minutes` sets the minimum standalone science-visible fragment. Shorter science fragments are merged into a contiguous preceding science chunk when possible; otherwise they are handed to occultation filling.
- `min_occultation_sequence_minutes` sets the minimum standalone occultation chunk. Scheduled occultation rows are now chunked and merged upstream in the visit-level occultation planner; XML-time fallback chunks still use the local short-chunk merge rule when needed.
- `occultation_nonvisible_tolerance_minutes` sets how many non-visible minutes are tolerated inside an occultation interval before a later occultation pass or validation failure is triggered. This applies to occultation only.
- `requested_occ_time_override` defaults to `true` and allows occultation scheduling to continue when requested-hours bookkeeping is incomplete or would otherwise block assignment.
- Scheduled occultation XML rows now come directly from the visit-level occultation planner: Step A chunks and merges the occultation intervals first, then Step B emits those exact rows and only falls back to catalog search for genuinely uncovered time.
- `run_visualizer_after_pipeline` can be set in the JSON config to generate a plot automatically after the XML is written.
- `visualizer_mode` accepts:
  - `priority`: main Gantt-style plot colored by sequence priority
  - `timeline`: simple chronological timeline
  - `target-time`: bar chart of total observation time per target
  - `simple`: lighter-weight priority timeline
  - `visibility`: priority plot with non-visible intervals overlaid from the run's visibility parquet files

If you need help, read `QUICK_START.md` for examples and troubleshooting tips.
