# Parity Plan: Ensuring Identical Outputs

**Date:** 2025-11-19
**Objective:** Achieve bit-for-bit (or near bit-for-bit) parity between `pandorascheduler` (legacy) and `pandorascheduler_rework` outputs.

## 1. Identified Discrepancies

Analysis of the source code has revealed the following critical divergences:

### A. Transit Overlap "Last Partner Wins" Bug
*   **Legacy Behavior:** In `transits.py`, the `Transit_overlap` function iterates through partner planets. For each partner, it calculates the overlap and **overwrites** the `Transit_Overlap` column. This means if a planet has multiple partners, only the overlap with the *last* partner in the list is preserved.
*   **Rework Behavior:** `catalog.py` correctly calculates the **maximum** overlap across all partners.
*   **Action:** To achieve parity, we must temporarily downgrade the rework logic to replicate this bug.

### B. Occultation Target Sorting
*   **Legacy Behavior:** `sched2xml_WIP.py` calls `sch_occ_new` with `sort_key="closest"`, which prioritizes occultation targets by slew distance.
*   **Rework Behavior:** `science_calendar.py` (formerly `xml_builder.py`) defaults `prioritise_occultations_by_slew` to `False`.
*   **Action:** Enable `prioritise_occultations_by_slew=True` in the default configuration or runtime arguments.

### C. Occultation Window Slicing
*   **Legacy Behavior:** `sched2xml_WIP.py` has specific logic for handling the first and last visibility flags (`v_flag[0]`, `v_flag[-1]`) and appending to `oc_starts`/`oc_stops`.
*   **Rework Behavior:** `science_calendar.py` uses a simplified `_occultation_windows` function.
*   **Action:** Verify if the simplified logic produces identical boundaries. If not, port the legacy logic verbatim.

### D. Coordinate Resolution
*   **Legacy Behavior:** Tries to read `RA`/`DEC` from the input DataFrame row. If it fails, it falls back to `SkyCoord.from_name` (SIMBAD).
*   **Rework Behavior:** Tries to read from DataFrame row. If it fails, it falls back to `SkyCoord.from_name`.
*   **Action:** Ensure the input DataFrames (`target_info`) in rework contain the exact same columns and data types as legacy to ensure the "try" block succeeds/fails in the same way.

## 2. Implementation Steps

### Step 1: Configure for Parity
Update `ScienceCalendarConfig` defaults in `src/pandorascheduler_rework/science_calendar.py`:
```python
@dataclass(frozen=True)
class ScienceCalendarConfig:
    ...
    prioritise_occultations_by_slew: bool = True  # Changed from False
    ...
```

### Step 2: Replicate Transit Overlap Bug
Modify `_apply_transit_overlaps` in `src/pandorascheduler_rework/visibility/catalog.py`:
*   Instead of `best_overlap = max(...)`, implement a loop that overwrites `best_overlap` with the current partner's overlap, matching the legacy iteration order.

### Step 3: Verify Occultation Slicing
*   Create a unit test with a specific visibility pattern (e.g., `[0, 1, 0, 0, 1]`) that exercises the edge cases.
*   Run both legacy `sch_occ_new` logic (extracted) and rework `_occultation_windows` against it.
*   If they differ, replace rework logic with legacy logic.

### Step 4: Run Comparison
*   Execute `scripts/run_schedule_comparison.py`.
*   Analyze diffs.

## 3. Long-Term Strategy

Once parity is confirmed and the rework is trusted:
1.  **Fix the Overlap Bug:** Remove the "last partner wins" logic and restore the `max()` logic.
2.  **Standardize Sorting:** Decide if slew sorting is desired and make it explicit.
3.  **Clean Up Slicing:** Use the cleaner rework slicing logic if it is functionally equivalent or better.

## 4. Immediate Action Items

1.  [ ] Modify `science_calendar.py` to default `prioritise_occultations_by_slew=True`.
2.  [ ] Modify `catalog.py` to replicate "Last Partner Wins" bug.
3.  [ ] Run comparison script.
