#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

import pandas as pd

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.science_calendar import (
    ScienceCalendarInputs,
    _ScienceCalendarBuilder,
    _extract_visibility_segment,
    _normalise_target_name,
    _parse_datetime,
    _read_visibility,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Dump occultation scheduling details for a single XML visit."
    )
    parser.add_argument("--config", type=Path, required=True, help="JSON config path")
    parser.add_argument(
        "--schedule-csv", type=Path, required=True, help="Schedule CSV used to build the XML"
    )
    parser.add_argument(
        "--data-dir", type=Path, required=True, help="Resolved data_* directory used by the XML builder"
    )
    parser.add_argument(
        "--visit-id",
        type=str,
        required=True,
        help="Visit ID as shown in XML/provenance, e.g. 0145 or 145",
    )
    parser.add_argument(
        "--validation-csv",
        type=Path,
        default=None,
        help="Optional visibility validation CSV for highlighting problematic chunks",
    )
    return parser.parse_args()


def _json_value(cfg: dict, key: str, default):
    value = cfg.get(key, default)
    return default if value is None else value


def build_config(config_path: Path, data_dir: Path, output_dir: Path) -> PandoraSchedulerConfig:
    cfg = json.loads(config_path.read_text())
    extra_inputs = dict(cfg.get("extra_inputs") or {})

    start_raw = cfg.get("window_start") or cfg.get("start") or "2026-04-01"
    end_raw = cfg.get("window_end") or cfg.get("end") or "2026-07-01"
    window_start = datetime.fromisoformat(str(start_raw).replace("Z", ""))
    window_end = datetime.fromisoformat(str(end_raw).replace("Z", ""))

    gmat = extra_inputs.get("visibility_gmat")
    gmat_path = Path(gmat) if gmat else None

    return PandoraSchedulerConfig(
        window_start=window_start,
        window_end=window_end,
        schedule_step=timedelta(hours=float(_json_value(cfg, "schedule_step_hours", 24.0))),
        targets_manifest=data_dir,
        gmat_ephemeris=gmat_path,
        output_dir=output_dir,
        transit_coverage_min=float(_json_value(cfg, "transit_coverage_min", 0.4)),
        min_visibility=float(_json_value(cfg, "min_visibility", 0.5)),
        transit_scheduling_weights=tuple(_json_value(cfg, "transit_scheduling_weights", [0.8, 0.0, 0.2])),
        sun_avoidance_deg=float(_json_value(cfg, "sun_avoidance_deg", 91.0)),
        moon_avoidance_deg=float(_json_value(cfg, "moon_avoidance_deg", 25.0)),
        earth_avoidance_deg=float(_json_value(cfg, "earth_avoidance_deg", 110.0)),
        earth_avoidance_day_deg=cfg.get("earth_avoidance_day_deg"),
        earth_avoidance_night_deg=cfg.get("earth_avoidance_night_deg"),
        twilight_margin_deg=float(_json_value(cfg, "twilight_margin_deg", 0.0)),
        daynight_mode=str(_json_value(cfg, "daynight_mode", "subsatellite")),
        st_sun_min_deg=float(_json_value(cfg, "st_sun_min_deg", 0.0)),
        st_moon_min_deg=float(_json_value(cfg, "st_moon_min_deg", 0.0)),
        st_earthlimb_min_deg=float(_json_value(cfg, "st_earthlimb_min_deg", 0.0)),
        st_required=int(_json_value(cfg, "st_required", 1)),
        roll_step_deg=float(_json_value(cfg, "roll_step_deg", 2.0)),
        min_power_frac=float(_json_value(cfg, "min_power_frac", 0.7)),
        obs_sequence_duration_min=int(_json_value(cfg, "obs_sequence_duration_min", 90)),
        occ_sequence_limit_min=int(_json_value(cfg, "occ_sequence_limit_min", 50)),
        min_sequence_minutes=int(_json_value(cfg, "min_sequence_minutes", 8)),
        min_science_sequence_minutes=cfg.get("min_science_sequence_minutes"),
        min_occultation_sequence_minutes=cfg.get("min_occultation_sequence_minutes"),
        break_occultation_sequences=bool(_json_value(cfg, "break_occultation_sequences", True)),
        show_progress=False,
        primary_only_mode=bool(_json_value(cfg, "primary_only_mode", False)),
        use_target_list_for_occultations=bool(_json_value(cfg, "use_target_list_for_occultations", False)),
        prioritise_occultations_by_slew=bool(_json_value(cfg, "prioritise_occultations_by_slew", False)),
        enable_occultation_xml=bool(_json_value(cfg, "generate_occultation_xml", True)),
        enable_occultation_pass1=bool(_json_value(cfg, "enable_occultation_pass1", True)),
        requested_occ_time_override=bool(_json_value(cfg, "requested_occ_time_override", False)),
        allow_occ_startracker_violation=bool(_json_value(cfg, "allow_occ_startracker_violation", False)),
        parallel_workers=int(_json_value(cfg, "parallel_workers", 0)),
        extra_inputs=extra_inputs,
    )


def shorten(name: str) -> str:
    if name.startswith("G") and len(name) > 4 and name[1:4].isdigit():
        return f"{name[:4]}..."
    return name


def resolve_visit_row(builder: _ScienceCalendarBuilder, visit_id: str) -> tuple[int, pd.Series]:
    target_visit = int(str(visit_id).lstrip("0") or "0")
    visit_counter = 0
    for _, row in builder.schedule.iterrows():
        target_label = str(row.get("Target", ""))
        if not target_label or target_label == "Free Time" or target_label.startswith("WARNING"):
            continue
        visit_counter += 1
        if visit_counter == target_visit:
            return visit_counter, row
    raise ValueError(f"Visit {visit_id} not found in schedule")


def print_validation_excerpt(validation_csv: Optional[Path], visit_id: str) -> None:
    if validation_csv is None or not validation_csv.exists():
        return
    df = pd.read_csv(validation_csv, dtype={"visit_id": str, "sequence_id": str})
    vid = f"{int(str(visit_id).lstrip('0') or '0'):04d}"
    sub = df[(df["visit_id"] == vid) & (df["status"] != "ok")].copy()
    if sub.empty:
        print("\nValidation: no failing chunks for this visit")
        return
    sub["target"] = sub["target"].astype(str).map(shorten)
    cols = [
        c
        for c in [
            "sequence_id",
            "target",
            "sequence_type",
            "occultation_pass",
            "duration_minutes",
            "visible_minutes",
            "nonvisible_minutes",
            "sequence_visibility_fraction",
            "status",
        ]
        if c in sub.columns
    ]
    print("\nValidation failures for this visit:")
    print(sub[cols].to_string(index=False))


def load_validation_failures(validation_csv: Optional[Path], visit_id: str) -> pd.DataFrame:
    rows = load_visit_validation_rows(validation_csv, visit_id)
    if rows.empty:
        return rows
    return rows[rows["status"] != "ok"].copy()


def load_visit_validation_rows(validation_csv: Optional[Path], visit_id: str) -> pd.DataFrame:
    if validation_csv is None or not validation_csv.exists():
        return pd.DataFrame()
    df = pd.read_csv(validation_csv, dtype={"visit_id": str, "sequence_id": str})
    vid = f"{int(str(visit_id).lstrip('0') or '0'):04d}"
    sub = df[df["visit_id"] == vid].copy()
    if sub.empty:
        return sub
    sub["_start_dt"] = pd.to_datetime(sub["start_utc"], utc=True, errors="coerce").dt.tz_localize(None)
    sub["_stop_dt"] = pd.to_datetime(sub["stop_utc"], utc=True, errors="coerce").dt.tz_localize(None)
    return sub


def prepare_occ_df(occ_df: pd.DataFrame) -> pd.DataFrame:
    occ_df = occ_df.copy()
    if "Target" in occ_df.columns:
        occ_df["Target"] = occ_df["Target"].astype(str)
    if "start" in occ_df.columns:
        occ_df["_start_dt"] = pd.to_datetime(occ_df["start"], utc=True, errors="coerce").dt.tz_localize(None)
    else:
        occ_df["_start_dt"] = pd.NaT
    if "stop" in occ_df.columns:
        occ_df["_stop_dt"] = pd.to_datetime(occ_df["stop"], utc=True, errors="coerce").dt.tz_localize(None)
    else:
        occ_df["_stop_dt"] = pd.NaT
    return occ_df


def _format_occ_rows(rows: pd.DataFrame) -> str:
    cols = [
        c
        for c in ["start", "stop", "Target", "Occultation Pass", "RA", "DEC"]
        if c in rows.columns
    ]
    display = rows.loc[:, cols].copy()
    if "Target" in display.columns:
        display["Target"] = display["Target"].astype(str).map(shorten)
    return display.to_string(index=False)


def _describe_scheduled_row_decision(
    builder: _ScienceCalendarBuilder,
    occ_row: pd.Series,
    seg_start: datetime,
    seg_stop: datetime,
) -> list[str]:
    occ_target = str(occ_row.get("Target", "")).strip()
    if not occ_target or occ_target.lower() == "no target":
        return [
            f"scheduled target: {occ_target or '<blank>'}",
            "time-limit reject: not applicable",
            "occ_visibility_score reject: not applicable",
            "would enter catalog_fallback immediately: yes (blank / No target row)",
        ]

    current_occ_time = builder.occultation_obs_time.get(occ_target, timedelta())
    target_time_limit = builder._get_occultation_time_limit(occ_target)
    time_limit_reject = current_occ_time >= target_time_limit
    acceptable, st_frac = builder._occ_visibility_score(occ_target, seg_start, seg_stop)
    visibility_reject = not acceptable

    return [
        f"scheduled target: {shorten(occ_target)}",
        (
            "time-limit reject: "
            f"{'yes' if time_limit_reject else 'no'} "
            f"({current_occ_time.total_seconds() / 3600:.2f} / "
            f"{target_time_limit.total_seconds() / 3600:.2f} hr used)"
        ),
        (
            "occ_visibility_score reject: "
            f"{'yes' if visibility_reject else 'no'} "
            f"(st_violation_frac={st_frac:.2f})"
        ),
        (
            "would enter catalog_fallback immediately: no"
            if not time_limit_reject and not visibility_reject
            else "would enter catalog_fallback immediately: no "
            "(current code skips/retries first; it does not jump straight to fallback)"
        ),
    ]


def print_occ_df_for_problematic_segments(
    builder: _ScienceCalendarBuilder,
    occ_df: pd.DataFrame,
    validation_csv: Optional[Path],
    visit_id: str,
) -> None:
    failures = load_validation_failures(validation_csv, visit_id)
    if failures.empty:
        return

    occ_index = prepare_occ_df(occ_df)
    print("\nocc_df rows for problematic segments:")

    for _, row in failures.iterrows():
        sequence_id = str(row["sequence_id"])
        target = shorten(str(row["target"]))
        seq_type = str(row.get("sequence_type", ""))
        start_dt = row["_start_dt"]
        stop_dt = row["_stop_dt"]
        visible = row.get("visible_minutes")
        nonvisible = row.get("nonvisible_minutes")
        frac = row.get("sequence_visibility_fraction")

        print(
            f"\nSequence {sequence_id}  {target}  {start_dt} -> {stop_dt}  "
            f"({row['duration_minutes']} min, visible={visible}, non-visible={nonvisible}, frac={frac})"
        )

        if seq_type != "occultation":
            print("  Not an occultation sequence; no occ_df row is expected.")
            continue

        exact_mask = (
            occ_index["_start_dt"].notna()
            & occ_index["_stop_dt"].notna()
            & (occ_index["_start_dt"] <= start_dt)
            & (occ_index["_stop_dt"] >= stop_dt)
        )
        overlap_mask = (
            occ_index["_start_dt"].notna()
            & occ_index["_stop_dt"].notna()
            & (occ_index["_start_dt"] < stop_dt)
            & (occ_index["_stop_dt"] > start_dt)
        )

        exact_rows = occ_index.loc[exact_mask]
        overlap_rows = occ_index.loc[overlap_mask]

        if not exact_rows.empty:
            print("  Exact/covering occ_df row(s):")
            print(_format_occ_rows(exact_rows))
            for line in _describe_scheduled_row_decision(
                builder, exact_rows.iloc[0], start_dt, stop_dt
            ):
                print(f"  {line}")
            continue

        if not overlap_rows.empty:
            print("  Overlapping occ_df row(s):")
            print(_format_occ_rows(overlap_rows))
            continue

        prev_rows = occ_index.loc[
            occ_index["_stop_dt"].notna() & (occ_index["_stop_dt"] <= start_dt)
        ].sort_values("_stop_dt", ascending=False).head(2)
        next_rows = occ_index.loc[
            occ_index["_start_dt"].notna() & (occ_index["_start_dt"] >= stop_dt)
        ].sort_values("_start_dt", ascending=True).head(2)

        print("  No matching occ_df row for this problematic segment.")
        if not prev_rows.empty:
            print("  Nearest previous occ_df row(s):")
            print(_format_occ_rows(prev_rows))
        if not next_rows.empty:
            print("  Nearest next occ_df row(s):")
            print(_format_occ_rows(next_rows))


def print_all_occultation_intervals(
    occ_starts: list[datetime],
    occ_stops: list[datetime],
    occ_df: pd.DataFrame,
    validation_csv: Optional[Path],
    visit_id: str,
) -> None:
    occ_index = prepare_occ_df(occ_df)
    visit_rows = load_visit_validation_rows(validation_csv, visit_id)
    if not visit_rows.empty:
        visit_rows = visit_rows[visit_rows.get("sequence_type", "").astype(str) == "occultation"].copy()

    print("\nAll occultation intervals in this visit:")
    for idx, (seg_start, seg_stop) in enumerate(zip(occ_starts, occ_stops), start=1):
        minutes = int((seg_stop - seg_start).total_seconds() / 60)
        print(f"\nInterval {idx:02d}  {seg_start} -> {seg_stop}  ({minutes} min)")

        exact_mask = (
            occ_index["_start_dt"].notna()
            & occ_index["_stop_dt"].notna()
            & (occ_index["_start_dt"] <= seg_start)
            & (occ_index["_stop_dt"] >= seg_stop)
        )
        overlap_mask = (
            occ_index["_start_dt"].notna()
            & occ_index["_stop_dt"].notna()
            & (occ_index["_start_dt"] < seg_stop)
            & (occ_index["_stop_dt"] > seg_start)
        )
        exact_rows = occ_index.loc[exact_mask]
        overlap_rows = occ_index.loc[overlap_mask]

        if not exact_rows.empty:
            print("  Matching occ_df row(s):")
            print(_format_occ_rows(exact_rows))
        elif not overlap_rows.empty:
            print("  Overlapping occ_df row(s):")
            print(_format_occ_rows(overlap_rows))
        else:
            print("  No matching occ_df row.")

        if visit_rows.empty:
            continue

        final_mask = (
            visit_rows["_start_dt"].notna()
            & visit_rows["_stop_dt"].notna()
            & (visit_rows["_start_dt"] < seg_stop)
            & (visit_rows["_stop_dt"] > seg_start)
        )
        final_rows = visit_rows.loc[final_mask].copy()
        if final_rows.empty:
            print("  No final XML occultation sequence recorded for this interval.")
            continue

        final_display_cols = [
            c
            for c in [
                "sequence_id",
                "target",
                "start_utc",
                "stop_utc",
                "duration_minutes",
                "occultation_pass",
                "status",
                "visible_minutes",
                "nonvisible_minutes",
                "sequence_visibility_fraction",
            ]
            if c in final_rows.columns
        ]
        final_rows.loc[:, "target"] = final_rows["target"].astype(str).map(shorten)
        print("  Final XML occultation sequence(s):")
        print(final_rows.loc[:, final_display_cols].to_string(index=False))


def main() -> int:
    args = parse_args()
    config = build_config(args.config, args.data_dir, args.schedule_csv.parent)
    inputs = ScienceCalendarInputs(schedule_csv=args.schedule_csv, data_dir=args.data_dir)
    builder = _ScienceCalendarBuilder(inputs, config)
    visit_counter, row = resolve_visit_row(builder, args.visit_id)

    target_label = str(row.get("Target", ""))
    target_name, star_name = _normalise_target_name(target_label)
    start = _parse_datetime(row.get("Observation Start"))
    stop = _parse_datetime(row.get("Observation Stop"))
    if start is None or stop is None:
        raise ValueError("Could not parse visit start/stop")

    visibility_df = _read_visibility(builder.data_dir / "targets" / star_name, star_name)
    if visibility_df is None or visibility_df.empty:
        raise ValueError(f"No science visibility data for {target_name}/{star_name}")

    visit_times, visibility_flags = _extract_visibility_segment(
        visibility_df,
        start,
        stop,
        config.min_sequence_minutes,
    )
    final_time = stop
    raw_segments = builder._raw_visit_segments(visit_times, visibility_flags, start, final_time)
    visit_segments = builder._apply_science_fragment_policy(raw_segments)
    oc_starts, oc_stops = builder._occultation_windows_from_segments(visit_segments)

    print(f"Visit {visit_counter:04d}: {target_name}")
    print(f"Window: {start} -> {stop}")
    print("\nVisit segments:")
    for seg_start, seg_stop, is_visible in visit_segments:
        kind = "science" if is_visible else "occultation"
        minutes = int((seg_stop - seg_start).total_seconds() / 60)
        print(f"- {kind:7s} {seg_start} -> {seg_stop}  ({minutes} min)")

    print("\nOccultation intervals passed into schedule_occultation_targets:")
    for seg_start, seg_stop in zip(oc_starts, oc_stops):
        minutes = int((seg_stop - seg_start).total_seconds() / 60)
        print(f"- {seg_start} -> {seg_stop}  ({minutes} min)")

    occultation_info = builder._find_occultation_target(
        oc_starts,
        oc_stops,
        start,
        final_time,
        float(row["RA"]),
        float(row["DEC"]),
    )
    print(f"\n_find_occultation_target returned: {occultation_info is not None}")
    if occultation_info is None:
        print("No occ_df returned")
        print_validation_excerpt(args.validation_csv, args.visit_id)
        return 0

    occ_df, scheduled = occultation_info
    print(f"scheduled = {scheduled}")
    print(f"occ_df rows = {len(occ_df)}")
    print("\nocc_df:")
    print(_format_occ_rows(occ_df) if not occ_df.empty else "<empty>")

    print_all_occultation_intervals(oc_starts, oc_stops, occ_df, args.validation_csv, args.visit_id)
    print_validation_excerpt(args.validation_csv, args.visit_id)
    print_occ_df_for_problematic_segments(builder, occ_df, args.validation_csv, args.visit_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
