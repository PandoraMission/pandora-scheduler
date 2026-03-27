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

## Config Notes

- `daynight_mode` accepts two values:
  - `subsatellite`: classify day/night from whether the subsatellite point is sunlit
  - `limb`: classify day/night from whether the nearest Earth limb point in the target direction is sunlit
- `generate_occultation_xml` can be set in the JSON config to turn occultation filling in the science calendar XML on or off.
- `requested_occ_time_override` can be set in the JSON config to allow occultation scheduling to continue when requested-hours bookkeeping is incomplete or would otherwise block assignment.

If you need help, read `QUICK_START.md` for examples and troubleshooting tips.
