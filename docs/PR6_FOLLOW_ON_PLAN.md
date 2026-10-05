# PR #6 Follow-On Implementation Plan

## Scope

This document tracks what we should carrying forward from the original PR #6 (`vbkostov:tb_main` → `mrtommyb:main`, 48 commits, 18 files changed, +1705/−230) after reviewing the current `pr6-minimal-fixes` branch state.

It excludes parity-only work and excludes large tooling additions such as the Streamlit visibility app, notebooks, devcontainer, and `requirements.txt`.

Status in this document is updated as of 2026-03-13 after Phase 4 implementation is complete. All phases are done.

## Current Status

Completed and already present on `pr6-minimal-fixes`:

1. Run-scoped data directory routing (`data_subdir` / `data_<sun>_<moon>_<earth>`) across the runner, pipeline, scheduler paths, visibility output, and XML input lookup.
2. XML predefined ROI slot 0 alignment so the first ROI uses the target RA/DEC (`xml/parameters.py`).
3. Regression coverage for both of those fixes.
4. Example config and config documentation updates.
5. `primary_only_mode` end to end, including config field, CLI `--primary-only` flag, scheduler gating, tests, and docs.
6. Name-based target manifest mapping (`"exoplanet"`, `"auxiliary-standard"`, etc.) instead of positional index.
7. `create_aux_list` backward-compatible path handling (accepts data dir directly or package root containing `data/`).
8. `use_legacy_mode` fast-path interval arithmetic in `check_if_transits_in_obs_window` and `_handle_targets_of_opportunity` to avoid building minute ranges.
9. Degenerate-window filtering (`interval_mask.sum() == 0` checks) in `observation_utils.py`.
10. **Phase 1: Occultation config toggles and CLI exposure** — three new config fields (`enable_occultation_xml`, `enable_occultation_pass1`, `strict_occultation_time_limits`), seven CLI flags, `_as_bool` helper, all three toggles wired through `science_calendar.py` and `observation_utils.py`, per-pass logging, `uncovered_minutes` tracking in Pass 4, startup logging, duplicate `base_path` line removed, docs and config updated, 7 new tests in `tests/test_occultation_toggles.py`.
11. **Phase 2: Schedule output quality** — `Comments` column added to `row_columns` in `_schedule_auxiliary_target` (all append calls updated), `_primary_transit_comment` function labels primary/secondary exoplanet transits, below-min-visibility targets now scheduled with `"Visibility = X.XX%"` comment instead of becoming Free Time, 13 new tests in `tests/test_schedule_output_quality.py`.
12. **Phase 3: Plumbing quality** — `_write_visibility_parquet` helper embeds keepout-angle metadata (`pandora.visibility_{sun,moon,earth}_deg`) in parquet schema, all 4 inline parquet writes replaced, `config` parameter added to `_apply_transit_overlaps`, three-way `generate_visibility` logic (explicit true/false/unset) in `_maybe_generate_visibility`, manifest.py already well-factored (no changes needed), 11 new tests in `tests/test_plumbing_quality.py`.
13. **Phase 4: Debug occultation visit script** — `scripts/debug_occultation_visit.py` provides per-visit occultation scheduling diagnostics: visit metadata, 4-pass analysis (coverage ranking, minute-resolution run breakdown), occultation list summary, and outcome prediction. Reuses production helpers (`resolve_star_visibility_file`, `load_visibility_arrays_cached`). Supports `--visit`, `--target`, and `--list-visits` modes. 4 new tests in `tests/test_debug_occultation_visit.py`.

Primary files already updated across that work:

- `run_scheduler.py`
- `src/pandorascheduler_rework/config.py`
- `src/pandorascheduler_rework/pipeline.py`
- `src/pandorascheduler_rework/scheduler.py`
- `src/pandorascheduler_rework/observation_utils.py`
- `src/pandorascheduler_rework/science_calendar.py`
- `src/pandorascheduler_rework/visibility/catalog.py`
- `src/pandorascheduler_rework/xml/parameters.py`
- `tests/test_pipeline_visibility.py`
- `tests/test_xml.py`
- `tests/test_config_propagation.py`
- `tests/test_scheduler_short_gaps.py`
- `tests/test_primary_only_mode.py`
- `tests/test_occultation_toggles.py`
- `tests/test_schedule_output_quality.py`
- `tests/test_plumbing_quality.py`
- `tests/test_debug_occultation_visit.py`
- `example_scheduler_config.json`
- `docs/EXAMPLE_SCHEDULER_CONFIG.md`

## Review Of The Original PR #6

After a full re-review of all 48 commits and the complete file-level diff of PR #6, the useful remaining work splits into four buckets.

### Already absorbed or superseded

These do not need to remain as planned follow-on items:

1. The original occultation multi-pass assignment work (Passes 1–4) is already reflected in the current `observation_utils.py` and `science_calendar.py` implementations.
2. Occultation time tracking and deprioritization have already landed, including the later improvement to use per-target `Number of Hours Requested` values rather than a single global limit.
3. The visibility-handling pieces from the original PR have already been folded into the current XML builder behavior and test coverage.
4. Broad XML comparison tooling is already present in `scripts/detailed_xml_diff.py` and `scripts/comprehensive_xml_analysis.py`.
5. Data sub-directory routing (`data_<sun>_<moon>_<earth>`) is fully implemented.
6. The PredefinedStarRoiRa[0] alignment fix is implemented in `xml/parameters.py`.
7. Name-based target manifest mapping and `create_aux_list` backward-compat path handling are present.
8. `use_legacy_mode` fast-path interval arithmetic avoids minute-grid intersection.
9. `primary_only_mode` is fully wired (config, CLI, scheduler, tests, docs).
10. Degenerate-window filtering in Passes 1 and 2 of `schedule_occultation_targets`.

### ~~Still valuable — occultation robustness improvements~~ ✅ DONE (Phase 1)

All three config toggles and their supporting infrastructure are now implemented:

1. **`enable_occultation_xml`** — gates occultation-target calculations during XML generation in `science_calendar.py`.
2. **`enable_occultation_pass1`** — makes Pass 1 optional via `use_pass1` parameter threaded through `_find_occultation_target` → `_build_occultation_schedule` → `schedule_occultation_targets`.
3. **`strict_occultation_time_limits`** — relaxed mode logs a warning and returns a large fallback limit instead of raising.

Supporting changes also landed:
- `_as_bool` helper in `run_scheduler.py` for robust boolean parsing.
- `uncovered_minutes` tracking in Pass 4 with warning log.
- Per-pass logging (assigned/remaining counts) across all four passes.
- Seven new CLI flags for all occultation knobs.
- Startup logging for toggle values.
- Duplicate `base_path` line removed from `observation_utils.py`.

### Still valuable — schedule output quality

~~The PR improves the quality and diagnostic value of schedule output. These are moderate-value improvements:~~ ✅ DONE (Phase 2)

All three items are now implemented:

1. **Comments column** in `_schedule_auxiliary_target` — `row_columns` extended to include `"Comments"`, all `scheduled_rows.append()` calls updated.
2. **`_primary_transit_comment`** function — labels primary vs secondary exoplanet transits using `Primary Target` column (with `Number of Transits to Capture >= 10` fallback).
3. **Below-min-visibility handling** — the best available target is now scheduled with a `"Visibility = X.XX%"` comment instead of being skipped.

### Still valuable — plumbing quality

~~Lower-priority improvements that strengthen the codebase:~~ ✅ DONE (Phase 3)

All items implemented:

1. **Visibility parquet metadata** — `_write_visibility_parquet` embeds keepout angles in parquet schema metadata. All 4 inline writes replaced.
2. **`_apply_transit_overlaps` config parameter** — accepts optional `config` to embed metadata when rewriting planet parquets.
3. **Explicit `generate_visibility` override logic** — three-way precedence: explicit `true` forces, explicit `false` disables even with GMAT, unset defaults to GMAT presence.
4. **`manifest.py` refactoring** — not needed; already well-factored with 13+ focused helper functions.

### Not worth carrying forward as separate work

These either add little operational value now or are better handled by the existing codebase shape:

1. Streamlit visibility app, notebooks, devcontainer, `requirements.txt`.
2. Another general-purpose XML diff script (we already have two).
3. Any broad reimplementation of occultation logic that already exists locally in a more evolved form.
4. The PR's `compare_earth_avoidance.py` and `compare_target_visibility.py` scripts (useful for one-off analysis but not core).
5. Logging format simplification (removing timestamps) — a style choice better left alone.

## Remaining Recommended Work

All four phases are now complete:\n\n1. ~~Occultation config toggles and CLI exposure.~~ ✅ COMPLETE\n2. ~~Schedule output quality (Comments column, primary/secondary labels, below-min-vis).~~ ✅ COMPLETE\n3. ~~Plumbing quality (visibility metadata, generate_visibility override, manifest.py refactor).~~ ✅ COMPLETE\n4. ~~Single-visit occultation debug tool.~~ ✅ COMPLETE

## ~~Phase 1: Occultation Config Toggles And CLI Exposure~~ ✅ COMPLETE

### Goal

Add the three new occultation behavior toggles from PR #6, wire them end-to-end through config → CLI → scheduling logic, and also expose the existing occultation knobs as CLI flags.

### Why

- These toggles are the single highest-value gap between PR #6 and our code.
- They provide fine-grained operational control over occultation scheduling.
- `enable_occultation_xml` lets users skip expensive occultation calculations entirely.
- `enable_occultation_pass1` lets users bypass the strict "one target for all intervals" requirement.
- `strict_occultation_time_limits` lets exploratory runs proceed with incomplete catalog data.
- The existing occultation knobs (`use_target_list_for_occultations`, `prioritise_occultations_by_slew`, `break_occultation_sequences`) still lack CLI flags.

### Current State ✅ IMPLEMENTED

All items below are now present:

- CLI flags for `use_target_list_for_occultations`, `prioritise_occultations_by_slew`, `break_occultation_sequences`
- Config fields, CLI parsing, and scheduling logic for `enable_occultation_xml`, `enable_occultation_pass1`, `strict_occultation_time_limits`
- `_as_bool` helper for robust boolean parsing
- `use_pass1` parameter on `schedule_occultation_targets`
- Pass 4 `uncovered_minutes` tracking with warning log
- Per-pass logging across all four passes

### Implementation

1. ✅ Add `enable_occultation_xml: bool = True`, `enable_occultation_pass1: bool = True`, and `strict_occultation_time_limits: bool = True` to `PandoraSchedulerConfig`. (Default `True` preserves existing strict behavior.)
2. ✅ Add `_as_bool` helper to `run_scheduler.py`.
3. ✅ Add CLI flags in `run_scheduler.py`:
   - `--use-target-list-for-occultations`
   - `--prioritise-occultations-by-slew`
   - `--no-break-occultation-sequences`
   - `--no-occultation-xml`
   - `--skip-occultation-pass1`
   - `--relaxed-occultation-time-limits`
4. ✅ Wire `enable_occultation_xml` into `science_calendar.py` to gate occultation calculations.
5. ✅ Wire `enable_occultation_pass1` through to `schedule_occultation_targets` via a `use_pass1` parameter.
6. ✅ Wire `strict_occultation_time_limits` into `science_calendar.py` so relaxed mode logs a warning and returns a very large effective limit instead of raising.
7. ✅ Add `uncovered_minutes` tracking in Pass 4 of `schedule_occultation_targets` with warning log.
8. ✅ Add per-pass logging (pass counts, remaining windows, escalation notes) to all passes.
9. ✅ Update `docs/EXAMPLE_SCHEDULER_CONFIG.md` and `example_scheduler_config.json`.
10. ✅ Add startup logging for the new toggle values (matching the PR's `logger.info("GENERATE_OCCULTATION_XML=...")` pattern).

### Tests ✅

All tests implemented in `tests/test_occultation_toggles.py` (7 tests):

1. ✅ Config defaults: all three fields default to `True`.
2. ✅ Config explicit `False`: can explicitly disable all three.
3. ✅ Relaxed time limits: missing target returns large fallback instead of raising.
4. ✅ Strict time limits: missing target raises `ValueError`.
5. ✅ Relaxed time limits: empty catalog returns large fallback.
6. ✅ `use_pass1=True` assigns targets via Pass 1.
7. ✅ `use_pass1=False` skips Pass 1, later passes still assign.

### Acceptance Criteria ✅

- ✅ All three new toggles configurable via JSON and CLI.
- ✅ Existing three occultation knobs configurable via CLI.
- ✅ Default behavior is unchanged (strict limits, Pass 1 enabled, occultation XML generated).
- ✅ A user can reproduce occultation strategy from CLI alone.

## ~~Phase 2: Schedule Output Quality~~ ✅ COMPLETE

### Goal

Improve the diagnostic value of schedule CSV output by adding a Comments column for auxiliary observations, labeling primary vs secondary transits, and scheduling below-min-visibility targets with a warning comment instead of skipping them.

### Why

- Makes schedule output self-documenting.
- Below-min-visibility targets being silently skipped is a debugging pain point.
- Primary/secondary labeling is needed for downstream analysis.
- All three changes are low-risk and operationally useful.

### Current Behavior ✅ IMPLEMENTED

- `_schedule_auxiliary_target` now has `row_columns = ["Target", "Observation Start", "Observation Stop", "RA", "DEC", "Comments"]`.
- `_primary_transit_comment` labels primary/secondary transits.
- Below-min-visibility targets scheduled with `"Visibility = X.XX%"` comment.

### Implementation

1. ✅ Add `"Comments"` to `row_columns` in `_schedule_auxiliary_target`.
2. ✅ Update all `scheduled_rows.append([...])` calls in `_schedule_auxiliary_target` to include a comments value.
3. ✅ Add `_primary_transit_comment(target_list, planet_name)` function using the PR's logic: check `Primary Target` column first, fall back to `Number of Transits to Capture == 10`.
4. ✅ Use `_primary_transit_comment` when building the primary-target schedule row.
5. ✅ Change below-min-visibility handling: schedule the best available target and add `"Visibility = X.XX%"` comment instead of skipping.

### Tests ✅

All tests implemented in `tests/test_schedule_output_quality.py` (13 tests):

1. ✅ `_primary_transit_comment` with Primary Target=True/False, numeric 1/0, transits fallback, missing planet, empty list, no columns.
2. ✅ Auxiliary schedule output includes Comments column.
3. ✅ Full-visibility auxiliary target has empty comment.
4. ✅ Below-min-visibility target is scheduled (not Free Time).
5. ✅ Below-min-visibility target has a Visibility comment.

### Acceptance Criteria ✅

- ✅ Schedule CSV includes Comments for auxiliary and primary observations.
- ✅ Primary targets labeled as `"primary exoplanet transit"` or `"secondary exoplanet transit"`.
- ✅ Below-min-visibility targets scheduled with diagnostic comment.

## ~~Phase 3: Plumbing Quality Improvements~~ ✅ COMPLETE

### Goal

Strengthen codebase quality with visibility parquet metadata, improved `generate_visibility` override semantics, and optional manifest refactoring.

### Why

- Parquet metadata makes visibility outputs self-describing and reproducible.
- Improved `generate_visibility` semantics match the PR's "explicit false disables even when GMAT is present" pattern.
- Manifest refactoring improves review-ability but is optional.

### Implementation

1. ✅ Add `_write_visibility_parquet(dataframe, output_path, config)` to `visibility/catalog.py` that embeds `pandora.visibility_{sun,moon,earth}_deg` in parquet schema metadata.
2. ✅ Replace all inline parquet writes with the new function.
3. ✅ Add `config` parameter to `_apply_transit_overlaps`.
4. ✅ Update `_maybe_generate_visibility` in `pipeline.py` to use explicit three-way logic: `true` → generate, `false` → skip even with GMAT, `unset` → default to GMAT presence.
5. ✅ Skipped — manifest.py already well-factored with 13+ focused helper functions.

### Tests ✅

All tests implemented in `tests/test_plumbing_quality.py` (11 tests):

1. ✅ Parquet file contains keepout metadata after `_write_visibility_parquet`.
2. ✅ BytesIO buffer contains keepout metadata.
3. ✅ Data round-trips correctly through new write path.
4. ✅ Existing pandas metadata not clobbered.
5. ✅ `_apply_transit_overlaps` with config embeds metadata.
6. ✅ `_apply_transit_overlaps` without config still works (backward compat).
7. ✅ `generate_visibility='true'` without GMAT forces generation.
8. ✅ `generate_visibility='false'` with GMAT skips generation.
9. ✅ Unset with GMAT generates.
10. ✅ Unset without GMAT skips.
11. ✅ `generate_visibility='no'` treated as explicit false.

### Acceptance Criteria ✅

- ✅ Visibility parquet files contain keepout angle metadata in schema.
- ✅ `generate_visibility=false` reliably disables visibility generation.

## ~~Phase 4: Single-Visit Occultation Debug Tool~~ ✅ COMPLETE

### Goal

Make it easy to explain why one visit produced a particular occultation schedule.

### Why

- This is still the hardest remaining debugging surface.
- It is more useful now than adding another broad comparison script, because broad XML diff tools already exist.
- The PR includes both `scripts/debug_occultation_chunks.py` and `scripts/run_one_visit_xml.py` which serve this purpose.
- Our version should reuse production logic rather than cloning it.

### Implementation ✅

1. ✅ Added `scripts/debug_occultation_visit.py` with three modes:
   - `--list-visits`: tabular listing of all schedule entries
   - `--visit N`: debug by zero-based row index
   - `--target NAME`: debug by target name (exact or partial match)
2. ✅ Output includes:
   - Visit metadata (target, start/stop, RA/DEC, duration, comments)
   - Candidate count from aux_targets/ visibility files
   - Pass 1 analysis: fully-visible candidates
   - Pass 2 analysis: single-interval equivalence note
   - Pass 3 analysis: partial-coverage ranking (top 15)
   - Pass 4 analysis: minute-resolution coverage with contiguous-run breakdown (top 5)
   - Occultation list summary: targets in/out of catalog
   - Final outcome prediction
3. ✅ Reuses production helpers: `resolve_star_visibility_file`, `load_visibility_arrays_cached`
4. ✅ Validated against real schedule output (`output_standalone/`).

### Tests ✅

All tests implemented in `tests/test_debug_occultation_visit.py` (4 tests):

1. ✅ `_load_schedule` correctly loads CSV and parses datetime columns.
2. ✅ `_find_occultation_candidates` finds stars with parquet files.
3. ✅ `_find_occultation_candidates` returns empty list when no aux dir.
4. ✅ `debug_visit` prints structured output with visit info.

### Acceptance Criteria ✅

- ✅ A single problematic visit can be explained without instrumenting production code.
- ✅ Output is detailed enough to compare segmentation and candidate selection decisions.

## Files Likely To Change

Primary implementation files:

- `run_scheduler.py`
- `src/pandorascheduler_rework/config.py`
- `src/pandorascheduler_rework/science_calendar.py`
- `src/pandorascheduler_rework/observation_utils.py`
- `src/pandorascheduler_rework/scheduler.py`
- `src/pandorascheduler_rework/pipeline.py`
- `src/pandorascheduler_rework/visibility/catalog.py`
- `src/pandorascheduler_rework/targets/manifest.py` (optional refactor)

Likely test files:

- `tests/test_config_behavior.py`
- `tests/test_config_propagation.py`
- `tests/test_occultation_deprioritization.py`
- `tests/test_xml_builder.py`
- optionally a new targeted script/helper test file

Likely docs/config files:

- `docs/EXAMPLE_SCHEDULER_CONFIG.md`
- `example_scheduler_config.json`
- `README.md`
- `QUICK_START.md`

Likely new script:

- `scripts/debug_occultation_visit.py`

## Suggested Execution Order

1. ~~Phase 1: Occultation config toggles + CLI exposure (highest value, most changes).~~ ✅ COMPLETE
2. ~~Phase 2: Schedule output quality (Comments, labels, below-min-vis).~~ ✅ COMPLETE
3. ~~Phase 3: Plumbing quality (visibility metadata, generate_visibility, manifest).~~ ✅ COMPLETE
4. ~~Phase 4: Debug tool.~~ ✅ COMPLETE

All phases complete.

Run targeted tests after each phase.

## Validation Plan

For each remaining phase:

1. Run directly affected unit tests.
2. Add one regression test for the new behavior.
3. Keep one default-path regression to ensure existing behavior did not drift.

End-of-tranche verification:

1. Run targeted config and propagation tests.
2. Run occultation-related tests.
3. Run XML-builder tests.
4. Optionally run one short real pipeline case using existing output data.

## Explicitly Out Of Scope

The following original PR #6 items are not recommended for additional implementation in this repo at this stage:

- Streamlit visibility app (`vis_app/`) and associated deployment docs
- Jupyter notebooks (`compare_target_visibility.ipynb`)
- devcontainer configuration (`.devcontainer/devcontainer.json`)
- `requirements.txt`
- Another broad XML diff tool
- `scripts/compare_earth_avoidance.py` and `scripts/compare_target_visibility.py` (useful for one-off analysis but not core)
- Re-porting older occultation helper structures that are already superseded by the current implementation
- Logging format changes (removing timestamps) — a style preference
- Default value changes (`occ_sequence_limit_min` 50→20, `visibility_earth_deg` 86→96) — these are site-specific operational choices, not code changes