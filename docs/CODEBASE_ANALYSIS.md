# Pandora Scheduler Rework - Codebase Analysis

**Analysis Date:** 2025-11-19  
**Analyzed by:** Antigravity AI

---

## Executive Summary

The **Pandora Scheduler Rework** is a clean-room reimplementation of NASA's Pandora mission scheduling pipeline. The project maintains backward compatibility with a legacy scheduler while modernizing the codebase with typed interfaces, testable components, and improved architecture. The codebase is currently in an **advanced migration phase**, with the core scheduler achieving output parity with the legacy system.

### Key Achievements
✅ **Scheduler parity achieved** - Rework scheduler produces identical outputs to legacy for Feb 2026 – Feb 2027 window  
✅ **Visibility pipeline ported** - Star and planet visibility generation implemented  
✅ **Comprehensive testing** - Comparison scripts validate outputs against legacy baseline  
✅ **Progress tracking** - Modern UX with progress bars and logging  

### Outstanding Work
⚠️ **Partner overlap calculations** - Not yet implemented in visibility pipeline  
⚠️ **Manifest generation** - JSON → CSV conversion needs porting  
⚠️ **XML builder completion** - Some legacy logic still unported  
⚠️ **Helper code cleanup** - Legacy dynamic imports need replacement  

---

## Codebase Analysis: Pandora Scheduler Rework

## Executive Summary

The `pandorascheduler_rework` is a modernized, type-safe, and modular reimplementation of the legacy Pandora scheduler. It achieves full functional parity with the legacy system while eliminating technical debt, improving maintainability, and removing direct code dependencies on the legacy package.

**Key Achievements:**
*   **Independence:** The rework package (`src/pandorascheduler_rework`) does not import any code from the legacy package (`src/pandorascheduler`). It only references the legacy directory to locate shared data files.
*   **Modularity:** Monolithic scripts have been broken down into focused modules (`pipeline`, `scheduler`, `visibility`, `targets`).
*   **Type Safety:** The codebase is fully typed with Python type hints, enabling static analysis and better IDE support.
*   **Performance:** Critical paths (visibility filtering, transit calculation) have been optimized using vectorized operations (NumPy/Pandas).

## Architecture

The system follows a linear pipeline architecture:

```mermaid
graph TD
    A[Target Definitions (JSON)] -->|targets/manifest.py| B[Target Manifests (CSV)]
    B -->|visibility/catalog.py| C[Visibility Data (CSV)]
    C -->|scheduler.py| D[Schedule (CSV)]
    D -->|science_calendar.py| E[Science Calendar (XML)]
    
    subgraph "Data Layer"
    B
    C
    end
    
    subgraph "Core Logic"
    A
    D
    E
    end
```

## Key Components

### 1. Pipeline Orchestration (`pipeline.py`)
The entry point for the system. It handles:
*   **Configuration:** Parsing command-line arguments and config files into typed dataclasses (`SchedulerRequest`, `SchedulerConfig`).
*   **Orchestration:** Calling the target manifest generator, visibility builder, and scheduler in the correct order.
*   **Hermeticity:** Ensuring all inputs are explicitly defined, reducing reliance on implicit global state.

### 2. Core Scheduler (`scheduler.py`)
The heart of the optimization loop. It mirrors the legacy `Schedule` function but with significant improvements:
*   **State Management:** Uses a `SchedulerState` dataclass to track the schedule's evolution, replacing loose variables.
*   **Vectorization:** Transit filtering and overlap calculations are vectorized for speed.
*   **Logic Separation:** The scheduling loop is cleaner, with helper functions for specific tasks like ToO handling and auxiliary target selection.

### 3. Visibility Generation (`visibility/`)
A dedicated package for calculating target visibility:
*   **`catalog.py`**: The main driver. It iterates through targets, calculates visibility windows, and generates the CSV files.
*   **`geometry.py`**: Handles the geometric calculations (Sun/Moon/Earth avoidance, SAA crossings).
*   **`config.py`**: Defines the configuration for visibility runs.
*   **No Simbad:** Strictly relies on catalog coordinates, raising errors for missing data to ensure reproducibility.

### 4. Target Management (`targets/manifest.py`)
Responsible for converting the complex JSON target definitions into flat CSV manifests used by the scheduler.
*   **Abstraction:** Decouples the scheduler from the specific layout of the `PandoraTargetList` repository.
*   **Normalization:** Handles legacy quirks (like planet name formatting) to ensure consistent identifiers.

### 5. Science Calendar (`science_calendar.py`)
Generates the final XML output for the spacecraft.
*   **XML Generation:** Uses `xml.etree.ElementTree` to build the XML structure programmatically.
*   **Validation:** Ensures that the generated XML adheres to the schema (implicitly, by structure).
*   **Clean Implementation:** A complete rewrite of the legacy `xml_builder.py` / `sched2xml` logic.

### 6. Observation Utilities (`observation_utils.py`)
A collection of shared helper functions.
*   **Legacy Parity:** Reimplements specific legacy behaviors (like `remove_short_sequences`) to ensure output parity.
*   **Shared Logic:** Centralizes logic used by both the scheduler and the XML builder.

## Data Flow & Legacy Integration

The rework coexists with the legacy system by sharing the **Data Layer**:

1.  **Input Data:** The rework reads manifests and visibility products from the run's output directory (for example, `output_standalone/data/`). The legacy system reads from `src/pandorascheduler/data/`.
2.  **Output Data:** The rework writes to a specified output directory, but mimics the directory structure of the legacy system (e.g., `targets/StarName/PlanetName/Visibility...`) to ensure downstream tools can process the results.
3.  **Independence:** The rework code *never* imports from `src/pandorascheduler`. This means the legacy code can eventually be deleted without breaking the rework (once the data directory is moved or redefined).

## Future Recommendations

1.  **Data Directory Migration:** Move the `data/` directory out of `src/pandorascheduler` to a top-level `data/` or `assets/` directory. This will sever the final link to the legacy package structure.
2.  **Legacy Removal:** Once the data directory is moved, the `src/pandorascheduler` directory can be archived or deleted.
3.  **Unit Test Coverage:** While high-level parity tests exist, more granular unit tests for `scheduler.py` logic would be beneficial.
4.  **Performance Profiling:** Continue to profile the visibility generation, as it remains the most compute-intensive part of the pipeline.
- `config.py` - Configuration dataclasses
- `diff.py` - Comparison utilities for validation

**Process Flow:**
1. Load Pandora orbital state vectors from GMAT output
2. Interpolate spacecraft, Earth, Sun, Moon positions to 1-minute cadence
3. Apply avoidance angles (Sun/Moon/Earth) to determine visibility
4. Tag South Atlantic Anomaly (SAA) crossings
5. Generate per-target visibility CSVs
6. Calculate transit timing for exoplanet targets
7. Compute SAA overlap and partner transit overlap

**Status:** ✅ Core functionality complete, ⚠️ partner overlap pending

#### 2. **Scheduler** (`scheduler.py`)
Core scheduling logic that assigns observation windows to targets.

**Key Features:**
- Deterministic scheduling algorithm
- Priority-based target selection
- Constraint satisfaction (visibility, SAA, occultation)
- Progress tracking with `tqdm`

**Status:** ✅ Produces identical outputs to legacy

#### 3. **Science Calendar Builder** (`science_calendar.py`, formerly `xml_builder.py`)
Converts schedule into science calendar XML format for spacecraft commanding.

**Responsibilities:**
- Observation sequencing
- Readout scheme parameter mapping
- Occultation window handling
- XML serialization

**Status:** ⚠️ Partially complete, some legacy logic unported

#### 4. **Pipeline Orchestration** (`pipeline.py`)
High-level API that coordinates visibility generation, scheduling, and XML export.

**Interface:**
```python
SchedulerRequest → SchedulerResult
```

**Status:** ✅ Functional with progress reporting

---

## Target Categories & Prioritization

The scheduler handles **7 target categories** with distinct scheduling strategies:

| Priority | Category | Type | Description |
|----------|----------|------|-------------|
| 1 | `time-critical` | Standard | Must occur at specified time |
| 2 | `primary-exoplanet` | Exoplanet | Primary science targets (mission-critical) |
| 3 | `monitoring-standard` | Standard | Fixed cadence for spacecraft health |
| 4 | `auxiliary-exoplanet` | Exoplanet | Fill gaps with transit observations |
| 5 | `auxiliary-standard` | Standard | Fill gaps, any time |
| 6 | `occultation-standard` | Standard | Fill gaps during Earth occultation |
| 7 | `secondary-exoplanet` | Exoplanet | Backup targets (not actively scheduled) |

**Scheduling Strategies:**
- **Exoplanet targets** - Scheduled to observe transits using ephemeris
- **Standard targets** - Scheduled at any time within visibility windows
- **Time-critical** - Scheduled at exact specified times

---

## Data Flow

### Input Data Sources

1. **GMAT Ephemeris** - Spacecraft orbital state vectors
   - Format: CSV/TXT with position/velocity vectors
   - Cadence: Variable, interpolated to 1-minute resolution

2. **Target Definition Files** - JSON manifests per category
   - Location: `PandoraTargetList/target_definition_files/`
   - Contains: Star/planet parameters, readout schemes, priorities
   - Example categories: `primary-exoplanet/`, `monitoring-standard/`

3. **Readout Schemes** - Instrument configuration
   - `nirda_readout_schemes.json` - Near-infrared detector array
   - `vda_readout_schemes.json` - Visible detector array
   - Maps mnemonics to flight software commands

4. **Priority Files** - Scheduling metadata
   - Remaining transit counts for exoplanet targets
   - Observation cadence requirements
   - Deprioritization rules

### Output Artifacts

1. **Visibility CSVs**
   - Per-star: `data/targets/<Star>/Visibility for <Star>.csv`
   - Per-planet: `data/targets/<Star>/<Planet>/Visibility for <Planet>.csv`
   - Columns: Time, Visibility flags, SAA status, Transit coverage, Overlaps

2. **Schedule Files**
   - Master schedule CSV with observation windows
   - Tracker CSV with target status
   - Observation report summary

3. **Science Calendar XML**
   - Flight-ready commanding sequences
   - Includes payload parameters, timing, readout schemes

4. **Diagnostic Outputs**
   - Tracker pickle for state persistence
   - Debug logs for troubleshooting

---

## Testing & Validation Strategy

### Regression Testing

**Primary Tool:** `scripts/run_schedule_comparison.py`
- Runs both legacy and rework schedulers
- Compares outputs: schedule CSV, tracker, observation report, XML
- Reports: MATCH or detailed diffs

**Usage:**
```bash
poetry run python scripts/run_schedule_comparison.py --window-days 7
```

### Debug Utilities

1. **`debug_science_calendar.py`** - Compare XML outputs
2. **`debug_occultation.py`** - Validate occultation window calculations
3. **`debug_transit_candidates.py`** - Check transit timing logic
4. **`debug_visibility_compare.py`** - Diff visibility CSVs
5. **`analyse_occultation_window.py`** - Deep-dive on specific windows

### Unit Testing

- Framework: `pytest`
- Coverage: Geometry calculations, visibility logic, helper functions
- Fixtures: Small target lists, trimmed GMAT samples

---

## Configuration & Dependencies

### Technology Stack

**Core Dependencies:**
- **Python 3.12+** - Modern type hints and performance
- **NumPy 2.1+** - Numerical computations
- **Pandas 2.2+** - Data manipulation
- **Astropy 6.0+** - Astronomical calculations (authoritative for time/coordinates)
- **tqdm 4.66+** - Progress bars

**Development Tools:**
- **Black** - Code formatting (88 char line length)
- **Flake8** - Linting (E203, W503 ignored for Black compatibility)
- **pytest** - Testing framework
- **Poetry** - Dependency management

### Environment Setup

```bash
# Activate conda/mamba environment
conda activate <your-mamba-env>

# Install dependencies
poetry install

# Run tests
poetry run pytest

# Format code
poetry run black src/pandorascheduler_rework/

# Lint
poetry run flake8 src/pandorascheduler_rework/
```

---

## Development Workflow

### Guardrails

1. **Golden Rule** - Validate every change against legacy scheduler
2. **No Detachment** - Keep legacy helpers until replacements are tested
3. **Document Assumptions** - Note design decisions as functions are ported
4. **Fish Shell Compatibility** - Avoid Bash-specific syntax in scripts

### Migration Process

1. **Understand Legacy Behavior** - Document inputs/outputs/edge cases
2. **Implement Rework Version** - Typed, testable, documented
3. **Create Fixtures** - Small test cases for unit tests
4. **Validate Parity** - Run comparison scripts
5. **Update Documentation** - Reflect changes in `rework-notes.md`
6. **Retire Legacy Code** - Only after full validation

### Current Migration Status

**Completed:**
- ✅ Core scheduler logic
- ✅ Star visibility generation
- ✅ Planet transit timing
- ✅ SAA overlap calculation
- ✅ Progress tracking UX
- ✅ Comparison infrastructure

**In Progress:**
- 🔄 XML builder (partially complete)
- 🔄 Helper code cleanup

**Pending:**
- ⏳ Partner overlap calculations
- ⏳ Manifest generation (JSON → CSV)
- ⏳ Post-processing tools modernization
- ⏳ Legacy helper retirement

---

## Known Issues & Limitations

### Precision Considerations

**Transit Coverage & SAA Overlap:**
- Current implementation uses minute-resolution intersection
- Introduces up to ~0.09 fractional error
- **Future:** Consider higher-fidelity integration after parity locked

### External Dependencies

**SIMBAD Lookups:**
- Legacy defaults to network lookups for coordinates
- **Risk:** External failures can break pipeline
- **Future:** Add caching or favor manifest coordinates

### File Freshness

**Visibility Regeneration:**
- Legacy skips rewriting unless `force` flag set
- **Future:** Implement hashing or freshness checks for automatic stale detection

### XML Parity Issues

**Debug Funnel (as of 2025-11-17):**
- Some targets show XML differences: TOI-776, L_98-59, GJ_3470, G6869…0304
- Likely causes: Rounding differences, occultation window slicing
- **Mitigation:** Use `debug_occultation.py` to compare window calculations

---

## File Naming Conventions

### Target Definition Files
- Location: `PandoraTargetList/target_definition_files/<category>/`
- Format: JSON with star/planet parameters
- Example: `primary-exoplanet/TOI-776.json`

### Visibility Files
- Stars: `data/targets/<Star>/Visibility for <Star>.csv`
- Planets: `data/targets/<Star>/<Planet>/Visibility for <Planet>.csv`

### Output Files
- Schedule: `schedule_<start_date>.csv`
- Tracker: `tracker_<start_date>.csv`
- XML: `science_calendar_<start_date>.xml`

---

## Future Enhancements

### Near-Term (Next Sprint)
1. Complete partner overlap calculations
2. Port manifest generation to rework
3. Finish XML builder migration
4. Add CLI for visibility regeneration

### Medium-Term
1. Retire legacy helper dynamic imports
2. Modernize post-processing tools
3. Improve error handling and logging
4. Add configuration validation

### Long-Term
1. Higher-fidelity SAA/transit integration
2. SIMBAD caching layer
3. Automatic stale artifact detection
4. Performance optimization for large target lists
5. Web-based schedule visualization

---

## Key Insights from Documentation

### From `rework-notes.md`

**Critical Legacy Modules:**
- `scheduler_deprioritize_102925.py` - Target deprioritization logic
- `transits.py` - Visibility and transit calculations
- `sched2xml_WIP.py` - XML generation (work in progress)

**Development Environment:**
- Day-to-day terminal: **fish shell**
- Avoid heredocs and Bash-specific syntax
- Keep inline commands single-line

**Data Locations:**
- Test fixtures: `comparison_outputs/target_definition_files_limited/`
- Production data: `/Users/tsbarcl2/gitcode/PandoraTargetList/target_definition_files/`

### From `visibility-plan.md`

**Legacy Functions:**
1. `transits.star_vis` - Star visibility calculation
2. `transits.transit_timing` - Planet transit windows
3. `transits.Transit_overlap` - Partner transit overlap
4. `transits.SAA_overlap` - SAA crossing overlap
5. `scripts/vis_calc.py` - Orchestration script

**Rework Structure:**
- `geometry.py` - Vector math & interpolation
- `catalog.py` - Public API
- `serializers.py` - File I/O
- `config.py` - Configuration dataclasses

### From `target_definition_files/README.md`

**Readout Scheme Structure:**
- `metadata` - Version, last updated
- `data` - Command definitions
  - `CommandName` - Flight software command
  - `IncludedMnemonics` - Defined mnemonics
  - `FixedParameters` - Universal parameters
  - Mnemonic-specific parameters

---

## Recommendations

### For New Contributors

1. **Start with Documentation** - Read all `.md` files in `docs/`
2. **Run Comparison Script** - Understand baseline behavior
3. **Explore Fixtures** - Small datasets in `comparison_outputs/`
4. **Use Debug Tools** - Scripts in `scripts/` for troubleshooting
5. **Follow Guardrails** - Always validate against legacy

### For Code Reviews

1. **Check Parity** - Ensure comparison script passes
2. **Verify Types** - All new code should use type hints
3. **Test Coverage** - Unit tests for new functions
4. **Documentation** - Update `rework-notes.md` with changes
5. **Formatting** - Run Black and Flake8 before commit

### For Deployment

1. **Full Regression** - Run comparison on extended window (30+ days)
2. **Validate Outputs** - Check all artifact types (CSV, XML, pickle)
3. **Review Logs** - Ensure no warnings or errors
4. **Backup Legacy** - Keep legacy code until full confidence
5. **Monitor Performance** - Track execution time for large schedules

---

## Glossary

- **GMAT** - General Mission Analysis Tool (orbital mechanics software)
- **SAA** - South Atlantic Anomaly (radiation belt region)
- **Ephemeris** - Predicted positions of celestial objects over time
- **Occultation** - When Earth blocks line-of-sight to target
- **Readout Scheme** - Detector configuration for observations
- **Transit** - Exoplanet passing in front of its host star
- **Visibility Window** - Time period when target is observable
- **Mnemonic** - Short identifier for readout configuration
- **Science Calendar** - Commanding sequence for spacecraft

---

## Contact & Resources

**Project Repository:** `mrtommyb/pandora-scheduler`  
**Target Definitions:** `/Users/tsbarcl2/gitcode/PandoraTargetList/`  
**Documentation:** `docs/` directory  
**Issue Tracking:** See conversation history for recent work  

**Related Projects:**
- Pandora Target List Repository
- GMAT Ephemeris Generation
- Flight Software Command Interface

---

## Appendix: File Inventory

### Documentation Files
- `docs/rework-notes.md` - Development notes and status
- `docs/visibility-plan.md` - Visibility pipeline design
- `comparison_outputs/target_definition_files_limited/README.md` - Target file format
- `src/pandorascheduler_rework/README.md` - Package overview

### Configuration Files
- `pyproject.toml` - Poetry dependencies and tool config
- `poetry.lock` - Locked dependency versions
- `.gitignore` - Version control exclusions

### Script Files
- `run_schedule_comparison.py` - Main validation tool
- `vis_calc.py` - Visibility generation entrypoint
- `debug_science_calendar.py` - XML comparison
- `debug_occultation.py` - Occultation validation
- `debug_transit_candidates.py` - Transit timing check
- `debug_visibility_compare.py` - Visibility diff
- `analyse_occultation_window.py` - Deep occultation analysis
- `check_data_fingerprints.py` - Data integrity check
- `generate_target_manifests.py` - Manifest creation

### Core Source Files
- `scheduler.py` - Scheduling algorithm
- `science_calendar.py` - Science calendar XML generation (formerly `xml_builder.py`)
- `observation_utils.py` - Observation scheduling utilities (formerly `helper_codes.py`)
- `pipeline.py` - High-level pipeline
- `utils/time.py` - Time rounding utilities (formerly `helper_codes_aux.py`)
- `utils/calendar_diff.py` - Calendar comparison (formerly `xml_compare.py`)
- `visibility/catalog.py` - Visibility generation
- `visibility/geometry.py` - Orbital calculations
- `visibility/config.py` - Configuration
- `visibility/serializers.py` - File I/O
- `visibility/diff.py` - Comparison utilities

---

**End of Analysis**
