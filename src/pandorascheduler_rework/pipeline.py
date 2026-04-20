"""High-level entry points for the reworked scheduler pipeline."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import pandas as pd

from pandorascheduler_rework import observation_utils as rework_helper
from pandorascheduler_rework.config import PandoraSchedulerConfig, resolve_data_subdir
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    run_scheduler,
)
from pandorascheduler_rework.too_selection import select_depth_ranked_toos
from pandorascheduler_rework.utils.io import read_csv_cached
from pandorascheduler_rework.visibility.catalog import build_visibility_catalog

LOGGER = logging.getLogger(__name__)


def _compute_visibility_window_end(
    target_list,
    config: PandoraSchedulerConfig,
) -> datetime:
    """Return the visibility horizon needed to finish late-starting visits."""
    max_extension = max(
        timedelta(hours=float(config.std_obs_duration_hours)),
        timedelta(minutes=float(config.occ_sequence_limit_min)),
    )

    if {"Planet Name", "Obs Window (hrs)"}.issubset(
        set(getattr(target_list, "columns", []))
    ):
        for _, row in target_list.iterrows():
            planet = str(row.get("Planet Name", "") or "").strip()
            if not planet:
                continue
            try:
                visit_duration = rework_helper.get_target_visit_duration(
                    planet, target_list
                )
            except Exception:
                continue
            if visit_duration > max_extension:
                max_extension = visit_duration

    return config.window_end + max_extension


@dataclass
class SchedulerResult:
    """Outputs produced by a scheduling run.

    Maintaining a structured return value simplifies comparisons with the
    legacy scheduler while leaving room to attach additional diagnostics later
    (plots, metrics, logs, ...).
    """

    schedule_csv: Optional[Path] = None
    reports: Dict[str, Path] = field(default_factory=dict)
    diagnostics: Dict[str, object] = field(default_factory=dict)

    def iter_output_files(self) -> Iterable[Path]:
        """Yield every concrete file generated during the run."""

        if self.schedule_csv:
            yield self.schedule_csv
        for path in self.reports.values():
            yield path


def build_schedule(config: PandoraSchedulerConfig) -> SchedulerResult:
    """Run the modern scheduler and persist outputs alongside diagnostics.

    Args:
        config: Unified configuration object

    Returns:
        SchedulerResult with paths to generated files
    """
    if config.output_dir is None:
        raise ValueError(
            "config.output_dir must be set (no default output directory fallback)"
        )

    output_dir = config.output_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    data_subdir = resolve_data_subdir(
        config.extra_inputs,
        sun_avoidance_deg=config.sun_avoidance_deg,
        moon_avoidance_deg=config.moon_avoidance_deg,
        earth_avoidance_deg=config.earth_avoidance_deg,
        earth_avoidance_day_deg=config.earth_avoidance_day_deg,
    )

    # Filesystem layout for this run. Use the run's `output_dir` as the
    # package root so generated data lives under `<output_dir>/<data_subdir>`.
    paths = SchedulerPaths.from_package_root(output_dir, data_dir_name=data_subdir)

    # Prepare extra_inputs
    extra_inputs = config.extra_inputs

    # When generating target manifests from a provided target definition base,
    # write the generated CSV manifests into the run's output data directory so
    # subsequent steps (visibility & calendar generation) can find them there.
    out_data = output_dir / data_subdir

    # We need to resolve these paths now to pass to target manifest generation
    # IMPORTANT: Use absolute paths to avoid legacy folder lookups
    primary_target_csv = _coerce_path(
        extra_inputs.get("primary_target_csv"),
        out_data / "exoplanet_targets.csv",
    ).resolve()
    auxiliary_target_csv = _coerce_path(
        extra_inputs.get("auxiliary_target_csv"),
        out_data / "auxiliary-standard_targets.csv",
    ).resolve()
    monitoring_target_csv = _coerce_path(
        extra_inputs.get("monitoring_target_csv"),
        out_data / "monitoring-standard_targets.csv",
    ).resolve()
    occultation_target_csv = _coerce_path(
        extra_inputs.get("occultation_target_csv"),
        out_data / "occultation-standard_targets.csv",
    ).resolve()

    # Handle target definition files if provided
    target_definition_files_raw = extra_inputs.get("target_definition_files")
    if target_definition_files_raw:
        target_definition_files = _resolve_target_definition_files(
            target_definition_files_raw,
            [
                _target_definition_from_csv(primary_target_csv),
                _target_definition_from_csv(auxiliary_target_csv),
                _target_definition_from_csv(monitoring_target_csv),
                _target_definition_from_csv(occultation_target_csv),
            ],
        )

        target_definition_base = _coerce_optional_path(
            extra_inputs.get("target_definition_base")
        )

        if target_definition_base is not None:
            # Allow callers to opt-out of regenerating manifests (use existing CSVs)
            skip_manifests = _as_bool(extra_inputs.get("skip_manifests"), False)
            if skip_manifests:
                LOGGER.info(
                    "Skipping generation of target manifests (skip_manifests=True)"
                )
            else:
                _generate_target_manifests(
                    target_definition_files,
                    target_definition_base,
                    primary_target_csv,
                    auxiliary_target_csv,
                    monitoring_target_csv,
                    occultation_target_csv,
                )
    else:
        # Fallback if not explicitly provided (legacy behavior relied on implicit paths)
        # But for V2 we prefer explicit. If missing, we assume the CSVs exist.
        target_definition_files = [
            _target_definition_from_csv(primary_target_csv),
            _target_definition_from_csv(auxiliary_target_csv),
            _target_definition_from_csv(monitoring_target_csv),
            _target_definition_from_csv(occultation_target_csv),
        ]

    too_candidate_target_csv = _maybe_generate_too_candidate_manifest(
        config,
        out_data,
    )

    if config.targets_manifest and not config.targets_manifest.exists():
        raise FileNotFoundError(
            f"Provided targets_manifest not found: {config.targets_manifest}"
        )

    # Load primary manifest early so we can validate before generating visibility
    # (visibility generation can take a long time).
    target_list = read_csv_cached(str(primary_target_csv))
    if target_list is None:
        raise FileNotFoundError(
            f"Primary target manifest not found: {primary_target_csv}"
        )

    _validate_primary_visit_windows(target_list, config)

    visibility_window_end = _compute_visibility_window_end(target_list, config)
    visibility_config = (
        replace(config, window_end=visibility_window_end)
        if visibility_window_end > config.window_end
        else config
    )
    if visibility_window_end > config.window_end:
        LOGGER.info(
            "Extending visibility generation horizon from %s to %s so boundary visits can complete",
            config.window_end,
            visibility_window_end,
        )

    _maybe_generate_visibility(
        visibility_config,
        paths,
        visibility_config.window_start,
        visibility_config.window_end,
        primary_target_csv,
        auxiliary_target_csv,
        monitoring_target_csv,
        occultation_target_csv,
        primary_visibility_csv=too_candidate_target_csv,
    )

    if too_candidate_target_csv is not None:
        too_result = _maybe_select_and_apply_auto_toos(
            config,
            out_data,
            primary_target_csv,
            too_candidate_target_csv,
        )
        if too_result is not None:
            LOGGER.info(
                "Selected %d automatic ToO target(s): %s",
                len(too_result.selected_targets),
                ", ".join(too_result.selected_targets) or "none",
            )
            target_list = read_csv_cached(str(primary_target_csv))
            if target_list is None:
                raise FileNotFoundError(
                    f"Primary target manifest not found after ToO selection: {primary_target_csv}"
                )

    scheduler_inputs = SchedulerInputs(
        pandora_start=config.window_start,
        pandora_stop=config.window_end,
        sched_start=config.window_start,
        sched_stop=config.window_end,
        target_list=target_list,
        paths=paths,
        target_definition_files=target_definition_files,
        primary_target_csv=primary_target_csv,
        auxiliary_target_csv=auxiliary_target_csv,
        occultation_target_csv=occultation_target_csv,
        output_dir=output_dir,
        tracker_pickle_path=_coerce_optional_path(extra_inputs.get("tracker_pickle")),
    )

    # scheduler_config = config.to_scheduler_config() (Removed)

    outputs = run_scheduler(scheduler_inputs, config)

    reports: Dict[str, Path] = {}
    if outputs.observation_report_path is not None:
        reports["observation_time"] = outputs.observation_report_path
    if outputs.tracker_csv_path is not None:
        reports["tracker_csv"] = outputs.tracker_csv_path
    if outputs.tracker_pickle_path is not None:
        reports["tracker_pickle"] = outputs.tracker_pickle_path

    diagnostics: Dict[str, Any] = {
        "schedule_dataframe": outputs.schedule,
        "tracker_dataframe": outputs.tracker,
    }

    return SchedulerResult(
        schedule_csv=outputs.schedule_path,
        reports=reports,
        diagnostics=diagnostics,
    )


def _validate_primary_visit_windows(
    target_list, config: PandoraSchedulerConfig
) -> None:
    """Validate that required visit windows can accommodate transit scheduling.

    We require per-target 'Obs Window (hrs)' and fail fast if any target's visit
    duration is too short to fit its transit duration plus required edge buffers.
    """
    required_columns = {
        "Planet Name",
        "Number of Transits to Capture",
        "Transit Duration (hrs)",
    }
    missing = required_columns.difference(set(getattr(target_list, "columns", [])))
    if missing:
        raise ValueError(
            "Primary target manifest is missing required columns for validation: "
            + ", ".join(sorted(missing))
        )

    problems: list[str] = []
    for _, row in target_list.iterrows():
        planet = str(row.get("Planet Name", "") or "")
        if not planet:
            continue

        transits_needed = row.get("Number of Transits to Capture")
        try:
            if transits_needed is None or float(transits_needed) <= 0:
                continue
        except (TypeError, ValueError):
            # If it can't be parsed, let downstream strictness catch it elsewhere.
            continue

        visit = rework_helper.get_target_visit_duration(planet, target_list)

        td_hours = row.get("Transit Duration (hrs)")
        if td_hours is None:
            raise ValueError(
                f"Target '{planet}' is missing required 'Transit Duration (hrs)'"
            )
        try:
            transit_duration = timedelta(hours=float(td_hours))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Target '{planet}' has invalid 'Transit Duration (hrs)': {td_hours!r}"
            ) from exc

        edge_buffer = rework_helper.compute_edge_buffer(
            visit,
            config.short_visit_threshold_hours,
            config.short_visit_edge_buffer_hours,
            config.long_visit_edge_buffer_hours,
        )
        required = transit_duration + 2 * edge_buffer
        if required > visit:
            problems.append(
                f"{planet}: visit={visit} transit={transit_duration} "
                f"edge_buffer={edge_buffer} required={required}"
            )

    if problems:
        sample = "\n".join(problems[:20])
        extra = "" if len(problems) <= 20 else f"\n... and {len(problems) - 20} more"
        raise ValueError(
            "One or more targets have Obs Window (hrs) too short to schedule the transit "
            "with required edge buffers. Fix the target definition file(s) by increasing "
            "Obs Window (hrs), or adjust edge buffer parameters if that is intended.\n"
            + sample
            + extra
        )


def _maybe_generate_too_candidate_manifest(
    config: PandoraSchedulerConfig,
    out_data: Path,
) -> Path | None:
    extra_inputs = config.extra_inputs
    if not _as_bool(extra_inputs.get("auto_select_toos"), False):
        return None

    raw_base = extra_inputs.get("too_candidate_target_definition_base")
    if raw_base is None:
        raise ValueError(
            "extra_inputs.auto_select_toos requires "
            "extra_inputs.too_candidate_target_definition_base"
        )

    candidate_base = Path(str(raw_base)).expanduser().resolve()
    candidate_csv = _coerce_path(
        extra_inputs.get("too_candidate_exoplanet_csv"),
        out_data / "all_exoplanet_targets.csv",
    ).resolve()
    LOGGER.info(
        "Building automatic ToO candidate manifest from %s",
        candidate_base,
    )
    manifest = rework_helper.process_target_files(
        "exoplanet",
        base_path=candidate_base,
    )
    candidate_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(candidate_csv, index=False)
    return candidate_csv


def _maybe_select_and_apply_auto_toos(
    config: PandoraSchedulerConfig,
    out_data: Path,
    primary_target_csv: Path,
    too_candidate_target_csv: Path,
):
    extra_inputs = config.extra_inputs
    if not _as_bool(extra_inputs.get("auto_select_toos"), False):
        return None

    primary_targets = pd.read_csv(primary_target_csv)
    main_target_names = set(primary_targets["Planet Name"].dropna().astype(str))

    top_n = int(extra_inputs.get("too_top_n", 5))
    if top_n < 0:
        raise ValueError("extra_inputs.too_top_n must be >= 0")
    depth_min_percent = float(extra_inputs.get("too_depth_min_percent", 1.0))
    coverage_min = float(
        extra_inputs.get("too_transit_coverage_min", config.transit_coverage_min)
    )
    query_depths = _as_bool(extra_inputs.get("too_query_depths"), True)

    result = select_depth_ranked_toos(
        data_dir=out_data,
        candidate_manifest=too_candidate_target_csv,
        main_targets=main_target_names,
        start=config.window_start,
        end=config.window_end,
        coverage_min=coverage_min,
        depth_min_percent=depth_min_percent,
        top_n=top_n,
        query_depths=query_depths,
    )
    _append_selected_toos_to_primary_manifest(
        primary_target_csv,
        too_candidate_target_csv,
        result.selected_targets,
    )
    LOGGER.info("Wrote ranked automatic ToO candidates: %s", result.ranked_targets_path)
    LOGGER.info("Wrote automatic ToO list: %s", result.too_list_path)
    return result


def _append_selected_toos_to_primary_manifest(
    primary_target_csv: Path,
    candidate_target_csv: Path,
    selected_targets: Sequence[str],
) -> None:
    if not selected_targets:
        return

    primary = pd.read_csv(primary_target_csv)
    candidate = pd.read_csv(candidate_target_csv)
    existing = set(primary["Planet Name"].dropna().astype(str))

    rows = []
    for target in selected_targets:
        if target in existing:
            continue
        match = candidate.loc[candidate["Planet Name"].astype(str) == str(target)]
        if match.empty:
            LOGGER.warning("Selected ToO %s missing from candidate manifest", target)
            continue
        row = match.iloc[0].copy()
        if "Number of Transits to Capture" in row.index:
            row["Number of Transits to Capture"] = 0
        if "Primary Target" in row.index:
            row["Primary Target"] = 0
        rows.append(row)

    if not rows:
        return

    appended = pd.concat([primary, pd.DataFrame(rows)], ignore_index=True)
    appended.to_csv(primary_target_csv, index=False)
    rework_helper.create_aux_list(
        [
            "exoplanet",
            "auxiliary-standard",
            "monitoring-standard",
            "occultation-standard",
        ],
        primary_target_csv.parent,
    )
    LOGGER.info(
        "Appended %d selected ToO target(s) to scheduler manifest %s",
        len(rows),
        primary_target_csv,
    )


def _maybe_generate_visibility(
    config: PandoraSchedulerConfig,
    paths: SchedulerPaths,
    pandora_start: datetime,
    pandora_stop: datetime,
    primary_target_csv: Path,
    auxiliary_target_csv: Path,
    monitoring_target_csv: Path,
    occultation_target_csv: Path,
    *,
    primary_visibility_csv: Path | None = None,
) -> None:
    # Three-way generate_visibility logic:
    #   explicit "true"/"yes"/"1"  -> always generate
    #   explicit "false"/"no"/"0"  -> never generate (even with GMAT)
    #   unset / empty              -> default to GMAT ephemeris presence
    raw = str(config.extra_inputs.get("generate_visibility", "")).strip().lower()
    if raw in {"1", "true", "yes", "y"}:
        generate_visibility = True
    elif raw in {"0", "false", "no", "n"}:
        generate_visibility = False
    else:
        generate_visibility = bool(config.gmat_ephemeris)

    if not generate_visibility:
        return

    # 1. Primary Targets -> <data_subdir>/targets
    LOGGER.info(
        "Generating visibility for Primary Targets in %s",
        paths.targets_dir if config.output_dir else Path("output") / "data" / "targets",
    )
    build_visibility_catalog(
        config,
        target_list=primary_visibility_csv or primary_target_csv,
        partner_list=auxiliary_target_csv,
        output_subpath="targets",
    )

    # 2. Auxiliary Targets -> <data_subdir>/aux_targets
    LOGGER.info(
        "Generating visibility for Auxiliary Targets in %s",
        paths.aux_targets_dir if config.output_dir else Path("output") / "data" / "aux_targets",
    )
    build_visibility_catalog(
        config,
        target_list=auxiliary_target_csv,
        partner_list=None,
        output_subpath="aux_targets",
    )

    # 3. Monitoring Targets -> <data_subdir>/aux_targets
    LOGGER.info(
        "Generating visibility for Monitoring Targets in %s",
        paths.aux_targets_dir if config.output_dir else Path("output") / "data" / "aux_targets",
    )
    build_visibility_catalog(
        config,
        target_list=monitoring_target_csv,
        partner_list=None,
        output_subpath="aux_targets",
    )

    # 4. Occultation Targets -> <data_subdir>/aux_targets
    LOGGER.info(
        "Generating visibility for Occultation Targets in %s",
        paths.aux_targets_dir if config.output_dir else Path("output") / "data" / "aux_targets",
    )
    build_visibility_catalog(
        config,
        target_list=occultation_target_csv,
        partner_list=None,
        output_subpath="aux_targets",
    )


def _generate_target_manifests(
    target_definition_files: Sequence[str],
    base_dir: Path,
    primary_target_csv: Path,
    auxiliary_target_csv: Path,
    monitoring_target_csv: Path,
    occultation_target_csv: Path,
) -> None:
    mapping = {
        "exoplanet": primary_target_csv,
        "auxiliary-standard": auxiliary_target_csv,
        "monitoring-standard": monitoring_target_csv,
        "occultation-standard": occultation_target_csv,
    }

    for category in target_definition_files:
        destination = mapping.get(str(category))
        if destination is None:
            continue

        manifest = rework_helper.process_target_files(
            category,
            base_path=base_dir,
        )
        destination.parent.mkdir(parents=True, exist_ok=True)
        manifest.to_csv(destination, index=False)

    rework_helper.create_aux_list(
        target_definition_files,
        primary_target_csv.parent,
    )




def _coerce_path(value: object, default: Path) -> Path:
    if value is None:
        return default
    return Path(str(value)).expanduser().resolve()


def _coerce_optional_path(value: object) -> Optional[Path]:
    if value is None:
        return None
    return Path(str(value)).expanduser().resolve()


def _resolve_target_definition_files(
    value: object, fallback: Sequence[str]
) -> List[str]:
    if value is None:
        return list(fallback)
    if isinstance(value, (list, tuple)):
        return [str(item) for item in value]
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    raise TypeError(
        "target_definition_files must be a sequence or comma-separated string"
    )


def _target_definition_from_csv(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_targets"):
        stem = stem[: -len("_targets")]
    return stem


def _as_bool(value: object, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized == "":
            return default
        return normalized in {"1", "true", "yes", "y", "on"}
    if isinstance(value, (int, float)):
        return bool(value)
    return default
