# Bug Fixes for Legacy Parity

This document describes bugs discovered during parity assessment and the fixes applied to match legacy behavior exactly.

## Bug 1: Zero-Duration Occultation Window Rejection

### Symptom
Visit 10 (K2-321, 2026-02-10 19:25 to 2026-02-11 17:49) was completely empty in rework output but populated in legacy. The `schedule_occultation_targets` function returned `False` despite successfully filling 27 of 28 occultation windows.

### Root Cause
The last occultation window had zero duration (start_mjd == stop_mjd at 60328.74167), creating an empty `interval_mask` when filtering the visibility time series. The rework code had an additional `np.any(interval_mask)` check:

```python
# BUGGY CODE (rework before fix)
if np.any(interval_mask) and np.all(visibility[interval_mask] == 1):
    # schedule target
```

Legacy code only checked:
```python
# LEGACY CODE
if np.all(visibility[interval_mask] == 1):
    # schedule target
```

Because `np.any([])` returns `False` but `np.all([])` returns `True` (vacuous truth), the rework rejected zero-duration windows while legacy accepted them.

### Fix
Removed the `np.any(interval_mask)` check from line 381 of `src/pandorascheduler_rework/observation_utils.py` (formerly `helper_codes.py`):

```python
# FIXED CODE (matches legacy)
if np.all(visibility[interval_mask] == 1):
    # schedule target
```

### Test Coverage
Added two regression tests in `tests/test_helper_codes.py`:
- `test_schedule_occultation_targets_handles_zero_duration_window`: Validates zero-duration windows accept first candidate
- `test_schedule_occultation_targets_fills_multiple_windows_with_one_unfilled`: Validates 28-window scenario with mixed durations (simulates K2-321)

### Impact
This bug affected any visit with occultation windows that land exactly on time series boundaries due to rounding. The fix ensures all windows are filled when a visible candidate exists, matching legacy behavior.

---

## Bug 2: Visit ID Padding Quirk

### Symptom
Visit 10 displayed as `<ID>0010</ID>` (4 digits) in rework but `<ID>00010</ID>` (5 digits) in legacy, causing 11-line XML diff.

### Root Cause
Legacy uses a quirky padding formula with an off-by-one error:

```python
# LEGACY CODE (in sched2xml_WIP.py)
visit_id = f'{("0"*(4-len(str(i))))+str(i+1)}'
```

Where `i` is a zero-indexed loop counter. For visit 10 (i=9):
- `len(str(9))` = 1
- `4 - 1` = 3 zeros
- Result: "000" + "10" = "00010" (5 digits!)

The rework used correct zero-padding:
```python
# BUGGY CODE (rework before fix)
visit_id = f"{visit_counter:04d}"  # Always produces 4 digits
```

### Fix
Replicated legacy's quirky formula in `src/pandorascheduler_rework/science_calendar.py` (formerly `xml_builder.py`) line ~152:

```python
# FIXED CODE (matches legacy bug)
visit_id = f"{'0' * max(0, 4 - len(str(visit_counter - 1)))}{visit_counter}"
```

### Test Coverage
Added regression test in `tests/test_xml_builder.py`:
- `test_visit_id_formatting_matches_legacy_quirk`: Creates 10 visits and validates formatting including "00010" for visit 10

### Impact
This is technically a bug in legacy code (should produce consistent 4-digit IDs), but we replicate it for byte-exact XML output matching. The quirk affects all visits with 2-digit numbers (10-99).

---

## Parity Validation

After applying both fixes, achieved full parity across all comparison categories:
- ✅ `schedule` - CSV schedule matches
- ✅ `tracker_csv` - Tracker CSV matches  
- ✅ `observation_report` - Observation time report matches
- ✅ `tracker_pickle` - Tracker pickle matches
- ✅ `science_calendar_xml` - Science calendar XML matches (0-line diff)

Validated on window sizes: 2-day, 7-day, and 21-day.

---

## Development Notes

### NumPy Vacuous Truth Behavior
- `np.all([])` returns `True` (universal quantification over empty set)
- `np.any([])` returns `False` (existential quantification over empty set)
- This behavior is critical for edge case handling in conditional logic

### Zero-Duration Windows
Zero-duration occultation windows arise when:
1. Occultation periods are broken into 31-minute segments
2. Segment boundaries coincide with time series sample points
3. Floating-point rounding produces start_mjd == stop_mjd

The legacy behavior is to accept the first candidate target for such windows, treating them as valid scheduling opportunities.

### Legacy Bug Replication
When achieving parity with legacy code, bugs in legacy must be replicated if they affect output format. The Visit ID padding quirk is an example where "correctness" conflicts with "byte-exact matching".
