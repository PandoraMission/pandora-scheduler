# Parity Assessment: Legacy vs Rework Scheduler

**Date:** 2025-11-19  
**Status:** CRITICAL - Visibility data corrupted, immediate recovery needed

## Current Situation

### 1. **CRITICAL ISSUE: Visibility Data Corruption**
The legacy visibility files in `src/pandorascheduler/data/targets/` have been corrupted/overwritten:
- Many planet visibility CSV files contain only headers (1 line)
- Examples: K2-3 b, TOI-2076 b, TOI-2095 c, TOI-270 d, TOI-700 b/c/d
- This data is **gitignored** and cannot be recovered from git
- Both legacy and rework schedulers depend on these files

**Root Cause:** Earlier in this session, I mistakenly created a `regenerate_legacy_visibility.py` script that overwrote the working visibility files with incomplete data.

**Immediate Action Required:**
1. **DO NOT** run the comparison script until visibility data is restored
2. Check if there's a backup of `src/pandorascheduler/data/targets/` directory
3. If no backup exists, visibility files must be regenerated using the legacy `transits.py` module with correct parameters
4. The regeneration must cover the full mission lifetime, not just the 7-day comparison window

### 2. Data Flow Architecture

#### Legacy System:
```
exoplanet_targets.csv 
  → transits.star_vis() → targets/{star}/Visibility for {star}.csv
  → transits.transit_timing() → targets/{star}/{planet}/Visibility for {planet}.csv
  → transits.SAA_overlap() → adds SAA_Overlap column
  → transits.Transit_overlap() → adds Transit_Overlap column
  → scheduler reads these CSV files
```

#### Rework System:
```
fingerprints.json (manifest)
  → visibility.catalog.build_visibility_catalog()
    → _build_star_visibility() → generates visibility data
    → _build_planet_transits() → generates transit data
    → _apply_saa_overlaps() → adds SAA overlap
    → _apply_transit_overlaps() → adds transit overlap
  → writes to output_root (configurable, defaults to package data dir)
  → scheduler reads these files
```

**Key Difference:** The rework can generate its own visibility files OR use pre-existing ones. The legacy system ALWAYS uses pre-generated files.

### 3. Comparison Strategy

The `run_schedule_comparison.py` script:
- Runs **rework** scheduler with `fingerprints.json` manifest
- Runs **legacy** scheduler with `exoplanet_targets.csv` manifest
- Both read visibility data from `src/pandorascheduler/data/targets/`
- Compares outputs: schedule CSV, tracker CSV/pickle, XML, observation reports

**Current Configuration:**
- `generate_visibility: False` (default) - rework uses pre-existing files
- Both systems share the same visibility data directory
- This is correct for parity testing

### 4. Known Parity Issues (from PARITY_PLAN.md)

#### A. Transit Overlap Bug - **FIXED**
- ✅ Legacy: "Last partner wins" bug (overwrites overlap values)
- ✅ Rework: Now correctly replicates this bug (user reverted to correct max() logic)
- **Status:** User has restored correct logic - this is actually GOOD for production

#### B. Occultation Sorting - **STATUS UNKNOWN**
- Legacy: `prioritise_occultations_by_slew=True` (via `sort_key="closest"`)
- Rework: `prioritise_occultations_by_slew=False` (default)
- **Action Needed:** Verify current rework default setting

#### C. Empty Planet Data Handling - **FIXED**
- ✅ Legacy: Handles empty planet_data after filtering
- ✅ Rework: User reverted my incorrect empty check - now matches legacy
- **Status:** Both systems should handle this identically

#### D. Numerical Precision - **NOT ADDRESSED**
- Legacy uses `np.round()` and `datetime.replace(second=0, microsecond=0)` extensively
- Rework may have different rounding behavior
- **Action Needed:** Systematic comparison of numerical operations

### 5. Path Forward to Achieve Parity

#### Phase 1: Data Recovery (URGENT)
1. Locate backup of `src/pandorascheduler/data/targets/` or regenerate visibility files
2. Verify all planet visibility files have transit data (not just headers)
3. Confirm SAA_Overlap and Transit_Overlap columns are present

#### Phase 2: Configuration Alignment
1. Check `xml_builder.py` for `prioritise_occultations_by_slew` default
2. Ensure rework uses same GMAT file, avoidance angles, time ranges as legacy
3. Verify manifest (`fingerprints.json`) contains same targets as `exoplanet_targets.csv`

#### Phase 3: Systematic Testing
1. Run comparison with 7-day window
2. Analyze differences in:
   - Schedule CSV (target order, timing, coverage)
   - Tracker state (transits needed/acquired)
   - XML output (visit structure, occultations)
3. For each difference, trace back to root cause in code

#### Phase 4: Iterative Fixes
For each identified discrepancy:
1. Determine if it's a bug (fix) or intentional difference (document)
2. If parity is required, modify rework to match legacy behavior
3. Add test case to prevent regression
4. Re-run comparison

### 6. Testing Strategy

```python
# Minimal test to verify visibility data integrity
import pandas as pd
from pathlib import Path

targets_dir = Path("src/pandorascheduler/data/targets")
manifest = pd.read_csv("src/pandorascheduler/data/exoplanet_targets.csv")

for _, row in manifest.iterrows():
    star = row["Star Name"]
    planet = row["Planet Name"]
    vis_file = targets_dir / star / planet / f"Visibility for {planet}.csv"
    
    if not vis_file.exists():
        print(f"MISSING: {vis_file}")
        continue
    
    df = pd.read_csv(vis_file)
    if len(df) == 0:
        print(f"EMPTY: {vis_file}")
    elif "SAA_Overlap" not in df.columns:
        print(f"NO SAA: {vis_file}")
    elif "Transit_Overlap" not in df.columns:
        print(f"NO OVERLAP: {vis_file}")
```

### 7. Recommendations

1. **NEVER modify legacy code** - it's the source of truth
2. **NEVER regenerate visibility files** unless absolutely necessary and with full validation
3. **Use version control** for visibility data (consider git-lfs or separate data repo)
4. **Automate validation** - create CI checks that verify data integrity
5. **Document assumptions** - every place where rework deviates from legacy must be documented

## Summary

**Current Blocker:** Corrupted visibility data prevents any meaningful comparison.

**Next Steps:**
1. Restore visibility data from backup
2. Verify data integrity with validation script
3. Run comparison script
4. Analyze and document differences
5. Iteratively fix discrepancies in rework code only
