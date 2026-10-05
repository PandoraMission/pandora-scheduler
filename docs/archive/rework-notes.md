# Pandora Scheduler Rework Notes

## Current Status (2025-11-16)

- Rework scheduler (`pandorascheduler_rework.scheduler`) produces identical outputs to the legacy implementation for the February 2026 – February 2027 window.
- `scripts/run_schedule_comparison.py` runs cleanly and reports MATCH for schedule, tracker CSV, observation report, and tracker pickle.
- Legacy helper modules are still imported dynamically; no behavioural divergences detected in the comparison run.

## Guardrails

- Preserve the golden rule: every meaningful change must be validated against the legacy scheduler via `poetry run python scripts/run_schedule_comparison.py`.
- Do not detach the rework from the legacy helpers until equivalent, well-tested replacements are ready; document assumptions as each function is ported.
- Summarize insights from the critical legacy modules (`scheduler_deprioritize_102925.py`, `transits.py`, `sched2xml_WIP.py`, etc.) here as they are understood so future work starts with shared context.
- Day-to-day terminal is fish; keep inline commands single-line (`python -c "..."`) and avoid heredocs or Bash-specific syntax.
- Target-definition fixtures currently live under `comparison_outputs/target_definition_files_limited/`; this is a temporary subset copied from the production `PandoraTargetList` repository. They can be used to understand the format of the files. Plan for the rework code to consume that external repo in-place rather than assuming the files ship with this project. The full files currently live in /Users/tsbarcl2/gitcode/PandoraTargetList/target_definition_files/.

## Outstanding Work

1. **Visibility Generation Port**
   - Legacy entrypoint: `scripts/vis_calc.py`, backed by `pandorascheduler.transits` and `helper_codes_aux`.
   - Outputs required by scheduler: per-target visibility CSVs under `data/targets`, auxiliary visibility under `data/aux_targets`, SAA overlap tables, partner overlap checks.
   - Need a rework module to orchestrate geometry calculations and emit the same directory structure.

2. **Science Calendar / XML Builder**
   - Legacy logic in `sched2xml_WIP.py` still unported.
   - `pandorascheduler_rework/xml_builder.py` currently raises `NotImplementedError`.
   - Requires porting observation sequencing helpers and XML serialization.

3. **Helper Porting and Cleanup**
   - `pandorascheduler_rework.helper_codes` and `_aux` are façade wrappers over the legacy modules.
   - Plan to incrementally reimplement the required helper functions with tests, then retire dynamic imports.

4. **Post-processing Tools**
   - Legacy scripts (`check_lost_days.py`, `confirm_XML_visibility.py`, `xml_to_visibility_figure.py`) remain untouched.
   - Decide which should be modernized alongside the rework.
5. **Run-time UX Improvements**
   - ✅ Scheduler now exposes a progress bar when `show_progress` is enabled (wired through `scripts/run_schedule_comparison.py`).
   - ✅ XML builder emits a `tqdm` progress bar during legacy comparisons to surface visit processing.

## Near-Term Plan

1. Draft a design for the visibility pipeline reimplementation (inputs, outputs, math checks, dependency needs).
2. Scaffold a `pandorascheduler_rework/visibility` module with a CLI to reproduce `vis_calc.py` behaviour, backed by unit tests on a small fixture dataset.
3. Once visibility generation is in place, revisit the XML builder and associated helper functions.

> Note: Re-run `poetry run python scripts/run_schedule_comparison.py` after substantial changes to ensure parity with the legacy scheduler remains intact.

## XML Parity Debug Funnel (2025-11-17)

1. **Confirm comparison baselines** – run `poetry run python scripts/run_schedule_comparison.py --window-days 7` and capture the diff targets reported by `scripts/debug_science_calendar.py`.
2. **Map mismatched visits** – note the schedule row index and visit ID for each target discrepancy (currently TOI-776, L_98-59, GJ_3470, and G6869…0304).
3. **Replicate legacy visibility windows** – for each affected visit:
   - load the relevant visibility CSVs;
   - run the legacy rounding + filtering (via `_extract_visibility_segment` mirror) to make sure start/stop times align.
4. **Compare occultation window slices** – dump `_occultation_windows` output for both legacy and rework paths to confirm the same number of segments and identical boundaries.
5. **Diff helper results** – use `scripts/debug_occultation.py --row-index <idx>` to compare the rework helper against `sched2xml_WIP.sch_occ_new` for the same intervals.
6. **Inspect visit assembly** – temporarily instrument `_ScienceCalendarBuilder._add_visit` to log sequence additions for the target under test, and compare to the equivalent legacy visit block.
7. **Rinse & repeat** – after each fix rerun the comparison, ensure the diff shrinks, and remove instrumentation once the XML diff is empty.

## Legacy Follow-up Notes

- `Transit_Coverage` and `SAA_Overlap` currently mirror the legacy minute-resolution intersection. This introduces up to ~0.09 fractional error; once parity is locked, consider a higher-fidelity integration.
- SIMBAD network lookups remain the default in legacy; evaluate whether to favor manifest/catalog coordinates or add caching so the rework is resilient to external failures.
- Legacy visibility skips rewriting files unless `force` is set. Explore hashing or freshness checks so the rework can detect stale artifacts automatically.
