"""Standalone 10 HJs plus automatic exoplanet-ToO experiment.

This module deliberately lives outside ``pandorascheduler_rework.pipeline`` so
the production calendar pipeline does not carry experiment-specific ToO logic.
It reuses the core manifest, visibility, and scheduler primitives, but owns the
experimental steps:

* build the main 10 HJ manifests;
* build a full-exoplanet visibility/ToO candidate manifest;
* select depth-ranked visible ToOs;
* append selected ToOs to the scheduler manifest;
* run the normal scheduler on the experiment manifest.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from dataclasses import fields, replace
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any, Optional, Sequence

import pandas as pd

from pandorascheduler_rework import observation_utils as rework_helper
from pandorascheduler_rework.config import PandoraSchedulerConfig, resolve_data_subdir
from pandorascheduler_rework.pipeline import (
    SchedulerResult,
    _as_bool,
    _coerce_optional_path,
    _coerce_path,
    _compute_visibility_window_end,
    _generate_target_manifests,
    _maybe_generate_visibility,
    _resolve_target_definition_files,
    _stage_too_list,
    _target_definition_from_csv,
    _validate_primary_visit_windows,
)
from pandorascheduler_rework.scheduler import (
    SchedulerInputs,
    SchedulerPaths,
    run_scheduler,
)
from pandorascheduler_rework.science_calendar import (
    ScienceCalendarInputs,
    generate_science_calendar,
)
from pandorascheduler_experiments.too_selection import (
    TooSelectionResult,
    select_depth_ranked_toos,
)
from pandorascheduler_rework.utils.io import read_csv_cached

LOGGER = logging.getLogger(__name__)


def build_schedule(config: PandoraSchedulerConfig) -> SchedulerResult:
    """Run the standalone 10 HJs + ToOs experiment."""

    if config.output_dir is None:
        raise ValueError("config.output_dir must be set")

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    data_subdir = resolve_data_subdir(
        config.extra_inputs,
        sun_avoidance_deg=config.sun_avoidance_deg,
        moon_avoidance_deg=config.moon_avoidance_deg,
        earth_avoidance_deg=config.earth_avoidance_deg,
        earth_avoidance_day_deg=config.earth_avoidance_day_deg,
        earth_keepouts=config.earth_keepouts,
        earth_avoidance_night_deg=config.earth_avoidance_night_deg,
        earth_avoidance_day_deg_occultation=config.earth_avoidance_day_deg_occultation,
        earth_avoidance_night_deg_occultation=config.earth_avoidance_night_deg_occultation,
    )
    paths = SchedulerPaths.from_package_root(output_dir, data_dir_name=data_subdir)
    out_data = output_dir / data_subdir
    extra_inputs = config.extra_inputs
    write_all_targets_csv = _as_bool(
        extra_inputs.get("write_all_targets_csv"),
        False,
    )

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

    target_definition_files = _resolve_target_definition_files(
        extra_inputs.get("target_definition_files"),
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
    if target_definition_base is None:
        raise ValueError(
            "The 10 HJs + ToOs experiment requires "
            "extra_inputs.target_definition_base"
        )

    if _as_bool(extra_inputs.get("skip_manifests"), False):
        LOGGER.info("Skipping generation of 10 HJ manifests")
    else:
        _generate_target_manifests(
            target_definition_files,
            target_definition_base,
            primary_target_csv,
            auxiliary_target_csv,
            monitoring_target_csv,
            occultation_target_csv,
            write_all_targets_csv=write_all_targets_csv,
        )

    too_candidate_target_csv = _generate_too_candidate_manifest(
        config,
        out_data,
        primary_target_csv,
    )

    target_list = read_csv_cached(str(primary_target_csv))
    if target_list is None:
        raise FileNotFoundError(f"Primary target manifest not found: {primary_target_csv}")

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

    too_result = _select_and_apply_toos(
        config,
        out_data,
        primary_target_csv,
        too_candidate_target_csv,
        write_all_targets_csv=write_all_targets_csv,
    )
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

    outputs = run_scheduler(scheduler_inputs, config)

    reports: dict[str, Path] = {}
    if outputs.observation_report_path is not None:
        reports["observation_time"] = outputs.observation_report_path
    if outputs.tracker_csv_path is not None:
        reports["tracker_csv"] = outputs.tracker_csv_path
    if outputs.tracker_pickle_path is not None:
        reports["tracker_pickle"] = outputs.tracker_pickle_path

    return SchedulerResult(
        schedule_csv=outputs.schedule_path,
        reports=reports,
        diagnostics={
            "schedule_dataframe": outputs.schedule,
            "tracker_dataframe": outputs.tracker,
            "too_selection": too_result,
        },
    )


def _generate_too_candidate_manifest(
    config: PandoraSchedulerConfig,
    out_data: Path,
    primary_target_csv: Path,
) -> Path:
    extra_inputs = config.extra_inputs
    raw_base = extra_inputs.get("too_candidate_target_definition_base")
    if raw_base is None:
        raise ValueError(
            "The 10 HJs + ToOs experiment requires "
            "extra_inputs.too_candidate_target_definition_base"
        )

    candidate_base = Path(str(raw_base)).expanduser().resolve()
    candidate_csv = _coerce_path(
        extra_inputs.get("too_candidate_exoplanet_csv"),
        out_data / "all_exoplanet_targets.csv",
    ).resolve()
    LOGGER.info("Building ToO candidate manifest from %s", candidate_base)
    manifest = rework_helper.process_target_files(
        "exoplanet",
        base_path=candidate_base,
    )
    manifest = _merge_primary_targets_for_visibility(manifest, primary_target_csv)
    candidate_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(candidate_csv, index=False)
    return candidate_csv


def _merge_primary_targets_for_visibility(
    manifest: pd.DataFrame,
    primary_target_csv: Path,
) -> pd.DataFrame:
    """Return full exoplanets with primary science rows taking precedence."""

    primary = pd.read_csv(primary_target_csv)
    if "Planet Name" not in primary.columns or "Planet Name" not in manifest.columns:
        raise ValueError(
            "Primary and ToO candidate manifests must contain a 'Planet Name' column"
        )

    primary_names = set(primary["Planet Name"].dropna().astype(str))
    non_primary_candidates = manifest.loc[
        ~manifest["Planet Name"].fillna("").astype(str).isin(primary_names)
    ]
    replaced_count = len(manifest) - len(non_primary_candidates)
    LOGGER.info(
        "Using %d primary science row(s) from the main target list in the "
        "ToO visibility manifest: %s",
        len(primary),
        ", ".join(primary["Planet Name"].dropna().astype(str)),
    )
    if replaced_count:
        LOGGER.info(
            "Replaced %d matching row(s) from the ToO candidate manifest with "
            "main science definitions",
            replaced_count,
        )
    return pd.concat([non_primary_candidates, primary], ignore_index=True)


def _select_and_apply_toos(
    config: PandoraSchedulerConfig,
    out_data: Path,
    primary_target_csv: Path,
    too_candidate_target_csv: Path,
    *,
    write_all_targets_csv: bool = False,
) -> TooSelectionResult:
    extra_inputs = config.extra_inputs
    primary_targets = pd.read_csv(primary_target_csv)
    main_target_names = set(primary_targets["Planet Name"].dropna().astype(str))
    existing_too_list = _stage_too_list(extra_inputs, out_data)
    if existing_too_list is not None:
        candidate_targets = pd.read_csv(too_candidate_target_csv)
        selected_targets = _selected_targets_from_existing_too_list(
            existing_too_list,
            primary_targets,
            candidate_targets,
        )
        _append_selected_toos_to_primary_manifest(
            primary_target_csv,
            too_candidate_target_csv,
            selected_targets,
            write_all_targets_csv=write_all_targets_csv,
        )
        LOGGER.info(
            "Using existing ToO list instead of automatic ToO selection: %s",
            existing_too_list,
        )
        return TooSelectionResult(
            ranked_targets_path=out_data / "all_too_candidates_ranked_by_transit_depth.csv",
            ranked_windows_path=out_data / "all_too_windows_ranked_by_transit_depth.csv",
            too_list_path=existing_too_list,
            selected_targets=selected_targets,
        )

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
        write_all_targets_csv=write_all_targets_csv,
    )
    LOGGER.info("Wrote ranked automatic ToO candidates: %s", result.ranked_targets_path)
    LOGGER.info("Wrote automatic ToO list: %s", result.too_list_path)
    return result


def _selected_targets_from_existing_too_list(
    too_list_csv: Path,
    primary_targets: pd.DataFrame,
    candidate_targets: pd.DataFrame,
) -> list[str]:
    """Validate and read selected targets from an existing ToO list."""

    too_table = pd.read_csv(too_list_csv)
    required = {"Target", "Obs Window Start", "Obs Window Stop"}
    missing = required.difference(too_table.columns)
    if missing:
        raise ValueError(
            f"ToO list {too_list_csv} is missing required columns: "
            + ", ".join(sorted(missing))
        )

    starts = pd.to_datetime(too_table["Obs Window Start"], errors="coerce")
    stops = pd.to_datetime(too_table["Obs Window Stop"], errors="coerce")
    bad_times = too_table.loc[starts.isna() | stops.isna()]
    if not bad_times.empty:
        raise ValueError(
            f"ToO list {too_list_csv} contains unparseable time rows: "
            + str(
                bad_times[
                    ["Target", "Obs Window Start", "Obs Window Stop"]
                ]
                .head(10)
                .to_dict(orient="records")
            )
        )
    bad_windows = too_table.loc[stops <= starts]
    if not bad_windows.empty:
        raise ValueError(
            f"ToO list {too_list_csv} contains non-positive windows: "
            + str(
                bad_windows[
                    ["Target", "Obs Window Start", "Obs Window Stop"]
                ]
                .head(10)
                .to_dict(orient="records")
            )
        )

    selected_targets = [
        str(target)
        for target in too_table["Target"].dropna().astype(str).tolist()
        if str(target).strip()
    ]
    known_targets = set(primary_targets["Planet Name"].dropna().astype(str))
    known_targets.update(candidate_targets["Planet Name"].dropna().astype(str))
    unknown = sorted(set(selected_targets).difference(known_targets))
    if unknown:
        raise ValueError(
            "ToO list contains target(s) that are not present in the 10 HJ "
            "manifest or full exoplanet candidate manifest: "
            + ", ".join(unknown)
        )

    return list(dict.fromkeys(selected_targets))


def _append_selected_toos_to_primary_manifest(
    primary_target_csv: Path,
    candidate_target_csv: Path,
    selected_targets: Sequence[str],
    *,
    write_all_targets_csv: bool = False,
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
    if write_all_targets_csv:
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


def config_from_json(
    config_path: Path,
    *,
    start: str,
    end: str,
    output: Path,
) -> tuple[PandoraSchedulerConfig, dict[str, Any]]:
    """Build a scheduler config for this experiment from a JSON file."""

    raw = json.loads(config_path.read_text(encoding="utf-8"))
    extra_inputs = dict(raw.get("extra_inputs") or {})
    gmat_ephemeris = raw.get("gmat_ephemeris") or extra_inputs.get("visibility_gmat")

    kwargs: dict[str, Any] = {
        "window_start": _parse_datetime(start),
        "window_end": _parse_datetime(end),
        "output_dir": output.expanduser().resolve(),
        "extra_inputs": extra_inputs,
    }
    if gmat_ephemeris is not None:
        kwargs["gmat_ephemeris"] = Path(str(gmat_ephemeris)).expanduser().resolve()
    if "schedule_step_hours" in raw:
        kwargs["schedule_step"] = timedelta(hours=float(raw["schedule_step_hours"]))
    if "transit_scheduling_weights" in raw:
        kwargs["transit_scheduling_weights"] = tuple(raw["transit_scheduling_weights"])
    if "include_occultation_sequences_in_xml" in raw:
        kwargs["enable_occultation_xml"] = _as_bool(
            raw["include_occultation_sequences_in_xml"],
            True,
        )
    if "generate_occultation_xml" in raw:
        kwargs["enable_occultation_xml"] = _as_bool(
            raw["generate_occultation_xml"],
            True,
        )

    field_names = {field.name for field in fields(PandoraSchedulerConfig)}
    special = {
        "window_start",
        "window_end",
        "output_dir",
        "extra_inputs",
        "gmat_ephemeris",
        "schedule_step_hours",
        "transit_scheduling_weights",
        "include_occultation_sequences_in_xml",
        "generate_occultation_xml",
    }
    path_fields = {"targets_manifest", "gmat_ephemeris", "output_dir"}
    tuple_fields = {"transit_scheduling_weights", "target_filters"}

    for key, value in raw.items():
        if key in special or key not in field_names or value is None:
            continue
        if key in path_fields:
            kwargs[key] = Path(str(value)).expanduser().resolve()
        elif key == "created_timestamp":
            kwargs[key] = value
        elif key in tuple_fields:
            kwargs[key] = tuple(value)
        else:
            kwargs[key] = value

    return PandoraSchedulerConfig(**kwargs), raw


def _parse_datetime(value: str) -> datetime:
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(value, fmt)
        except ValueError:
            pass
    return datetime.fromisoformat(value)


def _write_run_config_manifest(
    output_dir: Path,
    config_path: Path,
    raw_config: dict[str, Any],
) -> Path:
    manifest_path = output_dir / "run_config_manifest.json"
    payload = {
        "generated_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "experiment": "10_hjs_toos",
        "source_config_path": str(config_path.resolve()),
        "json_config": raw_config,
    }
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def _generate_xml_if_requested(
    config: PandoraSchedulerConfig,
    raw_config: dict[str, Any],
    result: SchedulerResult,
) -> Optional[Path]:
    if not _as_bool(raw_config.get("generate_xml"), False):
        return None
    if result.schedule_csv is None or config.output_dir is None:
        return None

    data_dir = config.output_dir / resolve_data_subdir(
        config.extra_inputs,
        sun_avoidance_deg=config.sun_avoidance_deg,
        moon_avoidance_deg=config.moon_avoidance_deg,
        earth_avoidance_deg=config.earth_avoidance_deg,
        earth_avoidance_day_deg=config.earth_avoidance_day_deg,
        earth_keepouts=config.earth_keepouts,
        earth_avoidance_night_deg=config.earth_avoidance_night_deg,
        earth_avoidance_day_deg_occultation=config.earth_avoidance_day_deg_occultation,
        earth_avoidance_night_deg_occultation=config.earth_avoidance_night_deg_occultation,
    )
    return generate_science_calendar(
        inputs=ScienceCalendarInputs(
            schedule_csv=result.schedule_csv,
            data_dir=data_dir,
        ),
        config=config,
        output_path=config.output_dir / "Pandora_science_calendar.xml",
        progress_label="Building science calendar",
    )


def _run_visualizer_if_requested(
    config: PandoraSchedulerConfig,
    raw_config: dict[str, Any],
    xml_path: Optional[Path],
) -> Optional[Path]:
    if not _as_bool(raw_config.get("run_visualizer_after_pipeline"), False):
        return None
    if xml_path is None:
        LOGGER.warning(
            "run_visualizer_after_pipeline is enabled, but no XML was generated. "
            "Set generate_xml=true to enable automatic visualization."
        )
        return None
    if config.output_dir is None:
        return None

    visualizer_mode = str(raw_config.get("visualizer_mode", "priority")).strip().lower()
    visualizer_script = Path(__file__).resolve().parents[2] / "scripts" / "visualizer.py"
    visualizer_output = config.output_dir / f"visualizer_{visualizer_mode}.png"
    visualizer_cmd = [
        sys.executable,
        str(visualizer_script),
        str(xml_path),
        "--mode",
        visualizer_mode,
        "--out",
        str(visualizer_output),
    ]
    if visualizer_mode == "visibility":
        data_dir = config.output_dir / resolve_data_subdir(
            config.extra_inputs,
            sun_avoidance_deg=config.sun_avoidance_deg,
            moon_avoidance_deg=config.moon_avoidance_deg,
            earth_avoidance_deg=config.earth_avoidance_deg,
            earth_avoidance_day_deg=config.earth_avoidance_day_deg,
            earth_keepouts=config.earth_keepouts,
            earth_avoidance_night_deg=config.earth_avoidance_night_deg,
            earth_avoidance_day_deg_occultation=config.earth_avoidance_day_deg_occultation,
            earth_avoidance_night_deg_occultation=config.earth_avoidance_night_deg_occultation,
        )
        visualizer_cmd.extend(["--data-dir", str(data_dir)])

    try:
        completed = subprocess.run(
            visualizer_cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else str(exc)
        LOGGER.warning("Visualizer generation failed: %s", stderr)
        return None

    if completed.stdout.strip():
        LOGGER.info(completed.stdout.strip())
    LOGGER.info("Visualizer written to: %s", visualizer_output)
    return visualizer_output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the standalone 10 HJs + automatic ToOs experiment."
    )
    parser.add_argument("--start", required=True, help="Start date/time")
    parser.add_argument("--end", required=True, help="End date/time")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--config", required=True, type=Path, help="JSON config")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    config, raw_config = config_from_json(
        args.config,
        start=args.start,
        end=args.end,
        output=args.output,
    )
    if config.output_dir is not None:
        config.output_dir.mkdir(parents=True, exist_ok=True)
        manifest = _write_run_config_manifest(config.output_dir, args.config, raw_config)
        LOGGER.info("Wrote run config manifest: %s", manifest)

    LOGGER.info("Starting standalone 10 HJs + ToOs experiment")
    result = build_schedule(config)
    xml_path = _generate_xml_if_requested(config, raw_config, result)
    visualizer_path = _run_visualizer_if_requested(config, raw_config, xml_path)
    if result.schedule_csv is not None:
        print(f"Schedule generated: {result.schedule_csv}")
    if xml_path is not None:
        print(f"Science calendar generated: {xml_path}")
    if visualizer_path is not None:
        print(f"Visualizer generated: {visualizer_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
