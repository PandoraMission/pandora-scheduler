"""Unified configuration system for Pandora Scheduler.

This module consolidates the scattered configuration classes into a single,
hierarchical system that's easier to understand and maintain.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class PandoraSchedulerConfig:
    """Master configuration for the Pandora Scheduler pipeline.

    This consolidates SchedulerConfig, ScienceCalendarConfig, and VisibilityConfig
    into a single, coherent configuration object.
    """

    # ============================================================================
    # TIMING & WINDOWS
    # ============================================================================

    window_start: datetime
    """Start of the scheduling window."""

    window_end: datetime
    """End of the scheduling window."""

    schedule_step: timedelta = timedelta(hours=24)
    """Rolling scheduling step size (default: 24 hours).

    This controls how far the scheduler advances its rolling window each
    iteration. It is not a per-target visit duration.
    """

    commissioning_days: int = 0
    """Number of commissioning days at start of mission."""

    # ============================================================================
    # PATHS & DATA SOURCES
    # ============================================================================

    targets_manifest: Optional[Path] = None
    """Path to target definition manifest/directory."""

    gmat_ephemeris: Optional[Path] = None
    """Path to GMAT ephemeris file (for visibility generation)."""

    output_dir: Optional[Path] = None
    """Output directory for generated files."""

    # ============================================================================
    # SCHEDULING THRESHOLDS
    # ============================================================================

    transit_coverage_min: float = 0.2
    """Minimum transit coverage to schedule (0-1). Lower = more transits scheduled."""

    min_visibility: float = 0.5
    """Minimum visibility fraction to consider observable (0–1).

    Default matches the ``--min-visibility`` CLI default.  0.0 schedules
    any target regardless of how little of the visit window is visible.
    """

    # ============================================================================
    # TRANSIT EDGE BUFFER PARAMETERS
    # ============================================================================

    short_visit_threshold_hours: float = 12.0
    """Visits shorter than this use short_visit_edge_buffer_hours."""

    short_visit_edge_buffer_hours: float = 1.5
    """Edge buffer (pre/post transit) for visits < short_visit_threshold_hours."""

    long_visit_edge_buffer_hours: float = 4.0
    """Edge buffer (pre/post transit) for visits >= short_visit_threshold_hours."""

    # ============================================================================
    # WEIGHTING FACTORS (must sum to 1.0)
    # ============================================================================

    transit_scheduling_weights: Tuple[float, float, float] = (0.8, 0.0, 0.2)
    """Unified transit scheduling weights: (coverage, saa, schedule).

    - weights[0]: transit_coverage (fraction of transit visible)
    - weights[1]: SAA overlap penalty — **pre-computed from GMAT ephemeris**.
      Setting this non-zero penalises transits that overlap the South Atlantic
      Anomaly, but the overlap cannot be avoided dynamically at scheduling time
      because it is fixed in the visibility catalog.  A ``UserWarning`` is
      emitted when this weight is non-zero.
    - weights[2]: schedule_factor = 1 − (gap_to_window_start / visit_duration)

    This triple is also recorded in the science calendar XML metadata.
    """

    # ============================================================================
    # KEEPOUT ANGLES (degrees)
    # ============================================================================

    sun_avoidance_deg: float = 91.0
    """Minimum angle from Sun (degrees)."""

    moon_avoidance_deg: float = 25.0
    """Minimum angle from Moon (degrees)."""

    earth_avoidance_deg: float = 110.0
    """Default Earth-center avoidance angle for the **boresight** (degrees).

    Measured as the angular separation between the boresight direction and the
    Earth's *centre* (centre-separation geometry).  This is intentionally
    different from the star-tracker Earth-limb keepout angles
    (``st_earthlimb_min_deg`` / ``st1_earthlimb_min_deg`` / ``st2_earthlimb_min_deg``)
    which use limb-angle geometry — see note on dual Earth models below.

    Used uniformly when both day/night overrides are None.  When either
    ``earth_avoidance_day_deg`` or ``earth_avoidance_night_deg`` is set,
    those values take precedence for the corresponding orbital phase.

    Dual Earth-avoidance model note
    --------------------------------
    Boresight keepout uses **centre-separation** (angle to Earth's centre).
    Star-tracker keepout uses **limb-angle** (angular elevation above the limb,
    computed via ``fast_limb_deg()`` in ``visibility/constraints.py``).
    These two quantities have different physical interpretations and cannot be
    compared directly.  At 600 km LEO the Earth's limb subtends ~20° from
    nadir, so a 110° centre-sep threshold corresponds roughly to a target
    being ~20° above the limb on the nadir side.
    """

    earth_avoidance_day_deg: Optional[float] = None
    """Boresight Earth-center avoidance when the nearest limb is sunlit (degrees).

    Uses centre-separation geometry (same as ``earth_avoidance_deg``).
    Set to ``None`` to use ``earth_avoidance_deg`` uniformly.
    Recommended value when enabled: 110.0.
    """

    earth_avoidance_night_deg: Optional[float] = None
    """Boresight Earth-center avoidance when the nearest limb is in shadow (degrees).

    Uses centre-separation geometry (same as ``earth_avoidance_deg``).
    Set to ``None`` to use ``earth_avoidance_deg`` uniformly.
    Recommended value when enabled: 80.0.
    """

    twilight_margin_deg: float = 0.0
    """Degrees past the geometric terminator to classify as sunlit (degrees).

    Expands the "sunlit" classification for day/night Earth-limb keepout.
    0 (default) gives a sharp terminator matching prior behaviour.  18 is
    analogous to astronomical twilight.  Only has an effect when
    ``earth_avoidance_day_deg`` or ``earth_avoidance_night_deg`` is set.
    """

    daynight_mode: str = "subsatellite"
    """Algorithm for determining day vs night for Earth-avoidance switching.

    * ``"subsatellite"`` (default): day/night is based on whether the subsatellite
      point (ground directly below the spacecraft) is sunlit.
    * ``"limb"`` : day/night is based on whether the nearest
      Earth limb point *in the target direction* is sunlit.

    Only has an effect when ``earth_avoidance_day_deg`` or
    ``earth_avoidance_night_deg`` is set.
    """

    # ============================================================================
    # STAR TRACKER KEEPOUT ANGLES (degrees)
    # ============================================================================

    st_sun_min_deg: float = 0.0
    """Minimum star-tracker–Sun separation (degrees). 0 = disabled."""

    st_moon_min_deg: float = 0.0
    """Minimum star-tracker–Moon separation (degrees). 0 = disabled."""

    st_earthlimb_min_deg: float = 0.0
    """Minimum star-tracker–Earth-limb separation (degrees). 0 = disabled.

    Uses **limb-angle geometry** (angular elevation above the Earth's limb),
    not centre-separation.  See the dual Earth-avoidance model note on
    ``earth_avoidance_deg`` for the distinction.
    """

    st1_earthlimb_min_deg: Optional[float] = None
    """Per-tracker override for ST1 Earth-limb keepout. None = use shared."""

    st2_earthlimb_min_deg: Optional[float] = None
    """Per-tracker override for ST2 Earth-limb keepout. None = use shared."""

    st_required: int = 1
    """Number of star trackers required to pass: 0 (skip), 1 (OR), or 2 (AND)."""

    # ============================================================================
    # ROLL OPTIMISATION
    # ============================================================================

    roll_step_deg: float = 2.0
    """Roll sweep step size (degrees). Smaller = more accurate but slower."""

    min_power_frac: float = 0.7
    """Minimum solar power fraction to accept a roll angle."""

    # ============================================================================
    # XML GENERATION PARAMETERS
    # ============================================================================

    obs_sequence_duration_min: int = 90
    """Observation sequence duration in minutes."""

    occ_sequence_limit_min: int = 50
    """Maximum occultation sequence duration in minutes."""

    min_sequence_minutes: int = 8
    """Minimum sequence length to include in XML (shorter sequences merged)."""

    min_science_sequence_minutes: Optional[int] = None
    """Minimum science-visible fragment length in minutes for XML handling.

    When ``None``, falls back to ``min_sequence_minutes``.
    """

    min_occultation_sequence_minutes: Optional[int] = None
    """Minimum occultation fragment length in minutes for XML tail handling.

    When ``None``, falls back to ``min_sequence_minutes``.
    """

    occultation_nonvisible_tolerance_minutes: int = 3
    """Allowed non-visible minutes inside an occultation interval.

    This tolerance applies only to occultation scheduling/validation and is
    used to absorb small boundary mismatches without forcing later occultation
    passes.
    """

    allow_science_soft_startracker_tail: bool = False
    """Allow science-visible XML segments to extend into a soft-ST tail.

    When enabled, the XML builder may extend the end of a science-visible
    sequence by up to ``science_soft_startracker_tail_minutes`` into the
    immediately following non-visible gap if the extension passes the hard
    boresight constraints and a softened star-tracker-only visibility check.
    This does not affect target selection or transit admission.
    """

    science_soft_startracker_tail_minutes: int = 10
    """Maximum tail length, in minutes, for soft-ST science extension."""

    science_soft_st_sun_min_deg: Optional[float] = None
    """Soft-tail override for star-tracker Sun keepout. None = use nominal."""

    science_soft_st_moon_min_deg: Optional[float] = None
    """Soft-tail override for star-tracker Moon keepout. None = use nominal."""

    science_soft_st_earthlimb_min_deg: Optional[float] = None
    """Soft-tail shared override for tracker Earth-limb keepout. None = nominal."""

    science_soft_st1_earthlimb_min_deg: Optional[float] = None
    """Soft-tail override for ST1 Earth-limb keepout. None = shared/nominal."""

    science_soft_st2_earthlimb_min_deg: Optional[float] = None
    """Soft-tail override for ST2 Earth-limb keepout. None = shared/nominal."""

    science_soft_st_required: Optional[int] = None
    """Soft-tail override for number of passing trackers required. None = nominal."""

    @property
    def effective_min_science_sequence_minutes(self) -> int:
        """Resolved science-sequence minimum in minutes."""
        if self.min_science_sequence_minutes is None:
            return self.min_sequence_minutes
        return self.min_science_sequence_minutes

    @property
    def effective_min_occultation_sequence_minutes(self) -> int:
        """Resolved occultation-sequence minimum in minutes."""
        if self.min_occultation_sequence_minutes is None:
            return self.min_sequence_minutes
        return self.min_occultation_sequence_minutes

    break_occultation_sequences: bool = True
    """Break long occultation sequences into chunks."""

    # ============================================================================
    # STANDARD OBSERVATIONS
    # ============================================================================

    std_obs_duration_hours: float = 0.5
    """Duration of standard star observations in hours."""

    std_obs_frequency_days: float = 3.0
    """Frequency of standard star observations in days."""

    # ============================================================================
    # BEHAVIOR FLAGS
    # ============================================================================

    show_progress: bool = False
    """Show progress bars during processing."""

    force_regenerate: bool = False
    """Force regeneration of files even if they exist."""

    primary_only_mode: bool = False
    """Disable non-primary gap-filling observations."""

    use_target_list_for_occultations: bool = False
    """Use target list for occultation scheduling (vs. separate list)."""

    prioritise_occultations_by_slew: bool = False
    """Prioritize occultation targets by slew angle."""

    enable_occultation_xml: bool = True
    """Enable occultation-target calculations during XML generation."""

    enable_occultation_pass1: bool = True
    """Enable Pass 1 in occultation assignment (single target covers all intervals)."""

    requested_occ_time_override: bool = False
    """When true, allow scheduling occultation targets beyond requested-hour limits."""

    allow_occ_startracker_violation: bool = False
    """Allow occultation targets that violate only star-tracker keepout.

    When True, occultation targets whose *boresight* constraints (Sun, Moon,
    Earth) all pass but whose star-tracker keepout fails are still accepted.
    Among such targets, candidates with shorter star-tracker violation are
    preferred.  Targets that violate any boresight constraint are still
    rejected.
    """

    # ============================================================================
    # PARALLELISM
    # ============================================================================

    parallel_workers: int = 0
    """Number of parallel workers for visibility generation.

    0  = auto (use all available CPUs).
    1  = serial (no multiprocessing overhead, useful for debugging).
    N  = use exactly N worker processes.
    """

    # ============================================================================
    # LEGACY COMPATIBILITY
    # ============================================================================

    use_legacy_mode: bool = False
    """Enable legacy scheduling behavior for validation against old outputs.
    
    When True, uses legacy algorithms that match the original scheduler exactly.
    When False (default), uses improved algorithms that may produce slightly
    different but equally valid (or better) results.
    
    Legacy behaviors controlled by this flag:
    - Visibility filtering: Uses MJD-based filtering (legacy) vs datetime-based
      filtering (modern). MJD filtering can exclude boundary points due to
      floating-point precision, while datetime filtering is more precise.
    
    Set to True when validating against historical baseline outputs.
    Set to False for production use with improved algorithms.
    """

    # ============================================================================
    # AUXILIARY SORTING
    # ============================================================================

    aux_sort_key: str = "sort_by_tdf_priority"
    """Key for sorting auxiliary targets."""

    # ============================================================================
    # METADATA
    # ============================================================================

    author: Optional[str] = None
    """Author name for XML metadata."""

    created_timestamp: Optional[datetime | str] = None
    """Creation timestamp for XML metadata."""

    visit_limit: Optional[int] = None
    """Limit number of visits (for testing). None = no limit."""

    target_filters: Sequence[str] = field(default_factory=tuple)
    """Target name filters for visibility generation."""

    extra_inputs: Dict[str, object] = field(default_factory=dict)
    """Additional input files (auxiliary lists, etc.)."""

    # ============================================================================
    # VALIDATION
    # ============================================================================

    def __post_init__(self) -> None:
        """Validate configuration consistency."""
        # Validate transit_scheduling_weights sum to 1.0
        if not np.isclose(sum(self.transit_scheduling_weights), 1.0):
            raise ValueError(
                "transit_scheduling_weights must sum to 1.0, got %s"
                % (sum(self.transit_scheduling_weights),)
            )

        # Warn when SAA weight is non-zero: SAA_Overlap is pre-computed and
        # cannot be avoided dynamically at scheduling time.
        if self.transit_scheduling_weights[1] != 0.0:
            warnings.warn(
                f"transit_scheduling_weights[1] (SAA weight) is "
                f"{self.transit_scheduling_weights[1]:.2f}. SAA_Overlap is "
                "pre-computed from the GMAT ephemeris and fixed in the "
                "visibility catalog — it cannot be avoided dynamically at "
                "scheduling time. The weight penalises transits that overlap "
                "the SAA but does not route around them.",
                UserWarning,
                stacklevel=2,
            )

        # Validate transit_coverage_min in range
        if not 0.0 <= self.transit_coverage_min <= 1.0:
            raise ValueError(
                "transit_coverage_min must be in [0, 1], got %s"
                % (self.transit_coverage_min,)
            )

        # Validate star tracker required count
        if self.st_required not in (0, 1, 2):
            raise ValueError(
                "st_required must be 0, 1, or 2, got %s" % (self.st_required,)
            )

        # Validate roll step
        if self.roll_step_deg <= 0:
            raise ValueError(
                "roll_step_deg must be > 0, got %s" % (self.roll_step_deg,)
            )

        # Validate min_power_frac
        if not 0.0 <= self.min_power_frac <= 1.0:
            raise ValueError(
                "min_power_frac must be in [0, 1], got %s"
                % (self.min_power_frac,)
            )

        for field_name, value in (
            ("min_sequence_minutes", self.min_sequence_minutes),
            (
                "min_science_sequence_minutes",
                self.min_science_sequence_minutes,
            ),
            (
                "min_occultation_sequence_minutes",
                self.min_occultation_sequence_minutes,
            ),
            (
                "occultation_nonvisible_tolerance_minutes",
                self.occultation_nonvisible_tolerance_minutes,
            ),
            (
                "science_soft_startracker_tail_minutes",
                self.science_soft_startracker_tail_minutes,
            ),
        ):
            if value is not None and value < 0:
                raise ValueError(
                    f"{field_name} must be >= 0, got {value}"
                )

        if (
            self.science_soft_st_required is not None
            and self.science_soft_st_required not in (0, 1, 2)
        ):
            raise ValueError(
                "science_soft_st_required must be 0, 1, or 2 when set, got %s"
                % (self.science_soft_st_required,)
            )

        # Validate daynight_mode
        if self.daynight_mode not in ("limb", "subsatellite"):
            raise ValueError(
                "daynight_mode must be 'limb' or 'subsatellite', got %r"
                % (self.daynight_mode,)
            )

        # Validate parallel worker count
        if self.parallel_workers < 0:
            raise ValueError(
                "parallel_workers must be >= 0, got %s"
                % (self.parallel_workers,)
            )

def build_default_data_subdir(
    sun_avoidance_deg: float,
    moon_avoidance_deg: float,
    earth_avoidance_deg: float,
    earth_avoidance_day_deg: Optional[float] = None,
) -> str:
    """Build the default run data directory name from keepout angles."""

    earth_label = (
        earth_avoidance_day_deg
        if earth_avoidance_day_deg is not None
        else earth_avoidance_deg
    )
    return (
        f"data_{int(float(sun_avoidance_deg))}_"
        f"{int(float(moon_avoidance_deg))}_"
        f"{int(float(earth_label))}"
    )


def resolve_data_subdir(
    extra_inputs: Mapping[str, object] | None,
    *,
    sun_avoidance_deg: float,
    moon_avoidance_deg: float,
    earth_avoidance_deg: float,
    earth_avoidance_day_deg: Optional[float] = None,
) -> str:
    """Resolve the run data directory name.

    When ``extra_inputs.data_subdir`` is not provided, derive the directory name
    from the keepout angles so multiple runs under one output root can coexist.
    """

    raw_value = None if extra_inputs is None else extra_inputs.get("data_subdir")
    if raw_value is None or str(raw_value).strip() == "":
        return build_default_data_subdir(
            sun_avoidance_deg,
            moon_avoidance_deg,
            earth_avoidance_deg,
            earth_avoidance_day_deg,
        )

    candidate = str(raw_value).strip()
    path_candidate = Path(candidate)
    if path_candidate.is_absolute():
        raise ValueError("extra_inputs.data_subdir must be a relative directory name")
    if path_candidate.name != candidate:
        raise ValueError(
            "extra_inputs.data_subdir must not include path separators"
        )
    if candidate in {"", ".", ".."}:
        raise ValueError("extra_inputs.data_subdir is invalid")
    return candidate
