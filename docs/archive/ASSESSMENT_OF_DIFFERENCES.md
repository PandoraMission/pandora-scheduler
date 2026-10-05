# Assessment of Output Differences: Legacy vs. Rework

**Date:** 2025-11-19
**Assessment Type:** Static Analysis (No code execution)

## Executive Summary

This document outlines the expected differences between the legacy `pandorascheduler` and the new `pandorascheduler_rework` pipelines. While the core scheduling logic has achieved high parity, several architectural and implementation improvements in the rework will lead to divergent outputs in specific areas.

## 1. Science Calendar XML (High Divergence)

The `ScienceCalendar.xml` output is the area with the most significant expected differences.

### A. Visit Mismatches
*   **Issue:** Known mismatches in visit generation for specific targets (e.g., `TOI-776`, `L_98-59`, `GJ_3470`).
*   **Cause:** Differences in floating-point rounding and filtering logic within visibility window extraction.
*   **Impact:** Some visits may be shifted by one minute or split differently compared to the legacy output.

### B. Occultation Window Slicing
*   **Issue:** The rework uses a cleaner, but different, logic for segmenting occultation windows (`_occultation_windows`).
*   **Cause:** The legacy `sched2xml_WIP.sch_occ_new` function has complex edge-case handling that was simplified in the rework.
*   **Impact:** Start and stop times for occultation observations may differ, and the number of segments in a block may vary.

### C. Coordinate Resolution
*   **Issue:** Differences in Right Ascension (RA) and Declination (DEC) values.
*   **Cause:**
    *   **Legacy:** Defaults to **SIMBAD network lookups** for coordinates.
    *   **Rework:** Prioritizes **catalog/manifest values** first, falling back to SIMBAD only if necessary.
*   **Impact:** Small numerical discrepancies in coordinate values if the input catalog differs from SIMBAD's current data.

## 2. Visibility CSVs (Medium Divergence)

### A. Overlap Calculations
*   **Issue:** Numerical differences in `Transit_Coverage` and `SAA_Overlap` columns.
*   **Cause:** The rework uses a minute-resolution intersection method which introduces a documented ~0.09 fractional error compared to the legacy method.
*   **Impact:** Planet visibility files will have slightly different overlap metrics, which could marginally affect scheduling priority for edge-case targets.

### B. File Freshness & Regeneration
*   **Issue:** The rework may regenerate files that the legacy code skips.
*   **Cause:** Legacy code skips generation if the file exists (unless forced). The rework implements smarter freshness checks.
*   **Impact:** If underlying GMAT or config data has changed, the rework will produce up-to-date files, while the legacy code might use stale cached files, leading to divergence.

## 3. Target Manifests (Medium Divergence)

*   **Issue:** Potential differences in `*_targets.csv` inputs.
*   **Cause:** The logic for converting JSON target definitions to CSV manifests is currently being ported.
*   **Impact:** If the rework generates these CSVs differently (e.g., different column order, formatting, or filtering), the scheduler inputs will differ.

## 4. Operational Differences

*   **Console Output:** The rework uses `tqdm` progress bars and structured logging, replacing the legacy print statements.
*   **Error Handling:** The rework is more robust against external failures (e.g., SIMBAD lookups), meaning it may successfully schedule a target that the legacy code would crash on or skip.

## Conclusion

The **Schedule CSV** and **Tracker CSV** are expected to match closely (verified by current tests), but the **Science Calendar XML** and **Visibility CSVs** will show differences due to improved logic and architectural changes. These differences are largely intentional improvements but must be accounted for during validation.
