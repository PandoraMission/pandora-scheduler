# Visibility Pipeline Reimplementation

## Goals

- Port the legacy visibility generation (`scripts/vis_calc.py` + `pandorascheduler/transits.py`) into the rework package while preserving output parity.
- Provide a typed, testable API that can regenerate the star/planet visibility artifacts consumed by the scheduler.
- Maintain the golden-rule regression guardrail by diffing generated files against legacy outputs before replacing them.

## Legacy Behaviour Summary

1. **`transits.star_vis`**
   - Loads Pandora orbital state vectors from GMAT output (CSV/TXT).
   - Interpolates spacecraft, Earth, Sun, Moon positions to 1-minute cadence across the observing window.
   - Applies sun/moon/earth avoidance angles to determine visibility for each target star; also tags South Atlantic Anomaly (SAA) crossings.
   - Writes per-star visibility CSV under `data/targets/<Star>/Visibility for <Star>.csv`.
2. **`transits.transit_timing`**
   - For each planet, loads star visibility CSV and target list metadata; derives mid-transit times and coverage using the visibility mask.
   - Writes per-planet visibility CSV under `data/targets/<Star>/<Planet>/Visibility for <Planet>.csv`.
3. **`transits.Transit_overlap`**
   - Compares transit windows among planets in the same system (target + partner lists) to compute overlap fractions.
   - Updates each planet’s visibility CSV with a `Transit_Overlap` column.
4. **`transits.SAA_overlap`**
   - Cross-references planet transit windows with SAA crossings from the star visibility file.
   - Updates the planet visibility CSV with `SAA_Overlap`.
5. **`scripts/vis_calc.py`** orchestration
   - Reads configuration constants (avoidance angles, observing window, GMAT file, etc.).
   - Runs `star_vis`, iterates through targets to call `transit_timing`, `SAA_overlap`, and partner-overlap checks.

## Proposed Rework Structure

```text
pandorascheduler_rework/
├── visibility/
│   ├── __init__.py
│   ├── geometry.py        # Vector math & interpolation utilities
│   ├── catalog.py         # Public API for generating visibility artifacts
│   ├── serializers.py     # File-writing helpers and layout management
│   └── config.py          # Dataclasses describing inputs (GMAT file, target lists, constraints)
```

### Key Components

- **Config dataclasses** capturing:
  - Observing window (`window_start`, `window_end`).
  - Avoidance angles (`sun_deg`, `moon_deg`, `earth_deg`).
  - Paths for GMAT ephemeris, target/partner manifests, output root.
- **Geometry utilities**
  - Load and interpolate GMAT ephemeris to minute cadence (`numpy.interp` parity check).
  - Compute angular separations using `astropy.coordinates.SkyCoord`.
  - Evaluate SAA footprint and return boolean masks.
- **Star visibility generator**
  - Iterate targets, normalise Simbad names (Gaia DR3 mapping), skip existing outputs unless forced.
  - Produce DataFrame matching legacy columns/rounding.
- **Planet transit generator**
  - Derive transit windows with `astropy.time.Time`; ensure truncation to minute resolution matches legacy behaviour.
  - Attach transit coverage, SAA overlap, and partner overlap metrics.
- **I/O helpers**
  - Maintain legacy folder/file structure.
  - Offer optional in-memory returns for tests (e.g., return DataFrames alongside paths).

## Testing Strategy

- Build fixtures under `tests/fixtures/visibility` with small target lists and trimmed GMAT samples.
- Write unit tests for:
  - Geometry interpolation (minute cadence, SAA mask).
  - Star visibility output parity (compare DataFrame columns and rounding rules).
  - Planet transit coverage and overlap calculations.
- Add integration test executing the catalog to reproduce a known set of outputs and diff against legacy files.

## Implementation Steps

1. [x] Scaffold `pandorascheduler_rework/visibility` package with config + loader stubs.
2. [x] Implement GMAT loader/interpolator and unit tests.
3. [x] Port star visibility logic, confirm parity on a single target fixture.
4. [x] Port planet transit timing; validate coverage and SAA overlap on fixtures.
5. [ ] Implement partner overlap calculations.
6. [ ] Port manifest generation (JSON → `*_targets.csv`) from `helper_codes.process_target_files` into the rework pipeline.
7. [ ] Build CLI entrypoint (e.g., `scripts/run_visibility_build.py`) to mirror `vis_calc.py`.
8. [ ] Update `docs/rework-notes.md` as components are ported; run comparison script after each milestone to verify scheduling outputs remain unchanged.
9. [x] Expose a pipeline flag that can regenerate visibility artifacts before scheduling.

## Open Questions

- Should we support partial recomputation (e.g., target subsets) natively? Legacy script already skips existing files—retain behaviour via `force` flag.
- Confirm the authoritative GMAT files and whether multiple versions (450km vs 600km) need configuration support at runtime.
- Determine how to expose visibility data to downstream XML builder—likely via shared config paths in the rework pipeline.
