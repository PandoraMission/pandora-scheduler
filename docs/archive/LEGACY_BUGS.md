# Legacy Code Bugs & Issues

**Date:** 2025-11-19
**Module:** `pandorascheduler.transits`

## 1. Critical Logic Bug: "Last Partner Wins"

**Location:** `transits.py`, lines 451-498 (specifically line 488)
**Function:** `Transit_overlap`

### Description
When calculating transit overlaps for a target planet against multiple partner planets in the same system, the code iterates through each partner sequentially. For each partner, it calculates the overlap fraction and **assigns** it to the `overlap` DataFrame, overwriting any value set by a previous partner.

### Code Snippet
```python
for m in range(len(All_start_transits.columns)): # Loop over partners
    ...
    for n in range(len(All_start_transits[planet_name].dropna())): # Loop over target transits
        ...
        # Calculate overlap with CURRENT partner
        transit_overlap = len(overlap_times)/len(transit_rng)
        
        # BUG: Overwrites previous value!
        overlap.loc[n, 'Transit_Overlap'] = transit_overlap
```

### Impact
If a planet has multiple partners (e.g., Partner A and Partner B):
1.  Overlap with Partner A is calculated (e.g., 0.5). `overlap` is set to 0.5.
2.  Overlap with Partner B is calculated (e.g., 0.0). `overlap` is overwritten to 0.0.
3.  The final result ignores Partner A entirely.

**Consequence:** The scheduler may incorrectly schedule a transit that is actually blocked by an earlier-processed partner, leading to a conflict or wasted observation time.

---

## 2. Logging Syntax Errors

**Location:** `transits.py`, lines 419, 421, 464
**Function:** `Transit_overlap`

### Description
The code uses `print`-style arguments for `logging.info`, which is incorrect for the standard Python `logging` module.

### Code Snippet
```python
logging.info('Only one planet in list around ', star_name)
```

### Impact
This will likely raise a `TypeError` (if the logger expects a format string) or simply log the first string and ignore the rest, meaning the `star_name` variable is never logged. This makes debugging difficult.

---

## 3. Inefficient Transit Matching

**Location:** `transits.py`, lines 466-473
**Function:** `Transit_overlap`

### Description
The code uses a nested loop to compare every transit of the target planet against every transit of the partner planet.

### Impact
*   **Complexity:** O(N * M) where N and M are the number of transits.
*   **Performance:** For short-period planets (e.g., 1-day period) over a 1-year mission, this results in ~365 * 365 = ~133,000 comparisons per pair. While not catastrophic for small N, it scales poorly.

---

## 4. Path Handling Assumptions

**Location:** `transits.py`, line 402

### Description
The code assumes that `target_list` and `partner_list` arguments are filenames residing specifically in `PACKAGEDIR/data/`.

### Code Snippet
```python
target_data = pd.read_csv(f'{PACKAGEDIR}/data/' + target_list, sep=',')
```

### Impact
Passing an absolute path or a path relative to the current working directory will cause a `FileNotFoundError` because the code prepends the package data directory path blindly.

---

## 5. Hardcoded GMAT Columns

**Location:** `transits.py`, lines 70-98

### Description
The code has a large `if/elif` block to handle specific GMAT file names, hardcoding the time column name (`Earth.A1ModJulian` vs `Earth.UTCModJulian`) based on the filename.

### Impact
If a user renames a GMAT file or uses a new one that isn't in this hardcoded list, the code may fail to detect the correct time column or default to the wrong one, leading to incorrect ephemeris interpolation.
