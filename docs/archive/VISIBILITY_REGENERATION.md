# Visibility Data Regeneration - In Progress

**Date:** 2025-11-19  
**Status:** REGENERATING

## What Happened

The legacy visibility files in `src/pandorascheduler/data/targets/` were corrupted during this session when I mistakenly created and ran a flawed regeneration script. This overwrote working visibility data with incomplete files.

## Recovery Process

Running the legacy scheduler directly to regenerate all visibility files:

```bash
python src/pandorascheduler/scheduler_deprioritize_102925.py
```

This executes `Schedule_all_scratch()` which:
1. Generates star visibility for all targets using `transits.star_vis()`
2. Generates planet transit data using `transits.transit_timing()`
3. Adds SAA overlap data using `transits.SAA_overlap()`
4. Adds transit overlap data using `transits.Transit_overlap()`
5. Runs the full scheduler

## Parameters Used

- **Dates:** 2026-02-05 to 2027-02-05 (1 year)
- **Avoidance Angles:** Sun=91°, Moon=25°, Earth=86°
- **GMAT File:** Pandora-600km-withoutdrag-20251018.txt
- **Mode:** `vis_and_schedule` (generates visibility AND runs scheduler)

## Progress

The script is currently running and generating:
- Star visibility files (✅ Complete)
- Planet transit files (✅ Complete)  
- SAA overlap data (✅ Complete)
- Transit overlap data (🔄 In Progress - ~13% complete, ~22 minutes remaining)

## Next Steps

Once regeneration completes:

1. **Validate Data Integrity**
   ```bash
   python scripts/validate_visibility_data.py
   ```
   Should show all 30 planets with complete data (Transit_Start, Transit_Stop, SAA_Overlap, Transit_Overlap columns)

2. **Run Comparison**
   ```bash
   python scripts/run_schedule_comparison.py --window-days 7
   ```
   This will compare rework vs legacy outputs

3. **Analyze Differences**
   - Review `comparison_outputs/science_calendar.diff`
   - Identify discrepancies in schedule, tracker, XML
   - Trace each difference to root cause in code

4. **Iterative Fixes**
   - Modify **rework code only** to match legacy behavior
   - Never modify legacy code - it's the source of truth
   - Re-run comparison after each fix

## Lessons Learned

1. **Never regenerate visibility files** unless absolutely necessary
2. **Never modify legacy code** - it's the reference implementation
3. **Always validate** before and after any data operations
4. **Use the legacy scheduler's built-in functions** for data generation
5. **Keep backups** of critical data files (consider git-lfs or separate data repo)

## Key Files

- **Legacy Scheduler:** `src/pandorascheduler/scheduler_deprioritize_102925.py`
- **Legacy Transits:** `src/pandorascheduler/transits.py`
- **Validation Script:** `scripts/validate_visibility_data.py`
- **Comparison Script:** `scripts/run_schedule_comparison.py`
- **Visibility Data:** `src/pandorascheduler/data/targets/{star}/{planet}/`
