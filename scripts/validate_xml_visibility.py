#!/usr/bin/env python3
"""Validate XML observation sequences against visibility parquet files."""

from __future__ import annotations

import argparse
import csv
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
from astropy.time import Time


@dataclass(frozen=True)
class XmlSequence:
    visit_id: str
    sequence_id: str
    target: str
    priority: str
    start: pd.Timestamp
    stop: pd.Timestamp

    @property
    def duration_minutes(self) -> int:
        return int(round((self.stop - self.start).total_seconds() / 60.0))


def _child_text(parent: ET.Element, path: str, default: str = "") -> str:
    node = parent.find(path)
    if node is None or node.text is None:
        return default
    return node.text.strip()


def parse_science_calendar(xml_path: Path) -> list[XmlSequence]:
    tree = ET.parse(xml_path)
    root = tree.getroot()

    sequences: list[XmlSequence] = []
    for visit_index, visit_el in enumerate(root.findall(".//{*}Visit"), start=1):
        visit_id = _child_text(visit_el, "./{*}ID", str(visit_index))
        for seq_index, seq_el in enumerate(
            visit_el.findall("./{*}Observation_Sequence"),
            start=1,
        ):
            seq_id = _child_text(seq_el, "./{*}ID", f"{visit_id}-{seq_index}")
            target = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Target",
                "",
            )
            priority = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Priority",
                "0",
            )
            start_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Timing/{*}Start",
            )
            stop_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Timing/{*}Stop",
            )
            if not target or not start_text or not stop_text:
                continue
            sequences.append(
                XmlSequence(
                    visit_id=visit_id,
                    sequence_id=seq_id,
                    target=target,
                    priority=priority,
                    start=pd.Timestamp(Time(start_text, format="isot", scale="utc").datetime).tz_localize(None),
                    stop=pd.Timestamp(Time(stop_text, format="isot", scale="utc").datetime).tz_localize(None),
                )
            )
    return sequences


def _load_planet_to_star(data_dir: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    exoplanet_csv = data_dir / "exoplanet_targets.csv"
    if not exoplanet_csv.exists():
        return mapping

    try:
        df = pd.read_csv(exoplanet_csv, usecols=["Planet Name", "Star Name"])
    except Exception:
        return mapping

    for _, row in df.iterrows():
        planet = str(row.get("Planet Name", "")).strip()
        star = str(row.get("Star Name", "")).strip()
        if planet and star:
            mapping[planet] = star
    return mapping


def _candidate_visibility_paths(
    data_dir: Path,
    target_name: str,
    planet_to_star: dict[str, str],
) -> list[Path]:
    candidates = [
        data_dir / "targets" / target_name / f"Visibility for {target_name}.parquet",
        data_dir / "aux_targets" / target_name / f"Visibility for {target_name}.parquet",
    ]
    star_name = planet_to_star.get(target_name)
    if star_name:
        candidates.extend(
            [
                data_dir / "targets" / star_name / f"Visibility for {star_name}.parquet",
                data_dir / "aux_targets" / star_name / f"Visibility for {star_name}.parquet",
            ]
        )
    return candidates


def _prepare_visibility_df(path: Path) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_parquet(path, columns=["Time_UTC", "Visible"])
        times = pd.to_datetime(df["Time_UTC"], errors="coerce")
    except Exception:
        try:
            df = pd.read_parquet(path, columns=["Time(MJD_UTC)", "Visible"])
        except Exception:
            return None
        times = pd.to_datetime(
            Time(df["Time(MJD_UTC)"].to_numpy(dtype=float), format="mjd", scale="utc").to_datetime()
        )

    index = pd.DatetimeIndex(times)
    if getattr(index, "tz", None) is not None:
        index = index.tz_localize(None)
    index = index.round("min")

    prepared = pd.DataFrame(
        {"Visible": (df["Visible"].to_numpy(dtype=float) > 0.5)},
        index=index,
    )
    prepared = prepared.groupby(level=0)["Visible"].max().to_frame()
    return prepared


def _resolve_visibility_df(
    data_dir: Path,
    target_name: str,
    planet_to_star: dict[str, str],
    cache: dict[str, Optional[pd.DataFrame]],
) -> tuple[Optional[pd.DataFrame], Optional[Path]]:
    if target_name in cache:
        return cache[target_name], None

    for candidate in _candidate_visibility_paths(data_dir, target_name, planet_to_star):
        if not candidate.exists():
            continue
        prepared = _prepare_visibility_df(candidate)
        if prepared is not None:
            cache[target_name] = prepared
            return prepared, candidate

    cache[target_name] = None
    return None, None


def _load_sequence_provenance(xml_path: Path) -> pd.DataFrame:
    provenance_path = xml_path.with_name(f"{xml_path.stem}_sequence_provenance.csv")
    if not provenance_path.exists():
        return pd.DataFrame()

    try:
        df = pd.read_csv(
            provenance_path,
            dtype={"visit_id": str, "sequence_id": str},
        )
    except Exception:
        return pd.DataFrame()

    needed = {"visit_id", "sequence_id"}
    if not needed.issubset(df.columns):
        return pd.DataFrame()

    keep_cols = [
        "visit_id",
        "sequence_id",
        "sequence_type",
        "occultation_pass",
        "visibility_fraction",
    ]
    available = [col for col in keep_cols if col in df.columns]
    return df[available].copy()


def validate_sequences(
    sequences: list[XmlSequence],
    data_dir: Path,
    provenance_df: Optional[pd.DataFrame] = None,
    occultation_nonvisible_tolerance_minutes: int = 3,
) -> list[dict[str, object]]:
    planet_to_star = _load_planet_to_star(data_dir)
    cache: dict[str, Optional[pd.DataFrame]] = {}
    results: list[dict[str, object]] = []
    provenance_lookup: dict[tuple[str, str], dict[str, object]] = {}
    if provenance_df is not None and not provenance_df.empty:
        for _, row in provenance_df.iterrows():
            provenance_lookup[(str(row.get("visit_id", "")), str(row.get("sequence_id", "")))] = row.to_dict()

    for seq in sequences:
        provenance = provenance_lookup.get((seq.visit_id, seq.sequence_id), {})
        visibility_df, source_path = _resolve_visibility_df(
            data_dir,
            seq.target,
            planet_to_star,
            cache,
        )
        n_minutes = seq.duration_minutes
        minute_index = (
            pd.date_range(seq.start, periods=n_minutes, freq="min")
            if n_minutes > 0
            else pd.DatetimeIndex([])
        )

        if visibility_df is None:
            results.append(
                {
                    "visit_id": seq.visit_id,
                    "sequence_id": seq.sequence_id,
                    "target": seq.target,
                    "priority": seq.priority,
                    "start_utc": seq.start.isoformat(),
                    "stop_utc": seq.stop.isoformat(),
                    "duration_minutes": n_minutes,
                    "status": "missing_visibility_file",
                    "visible_minutes": 0,
                    "nonvisible_minutes": n_minutes,
                    "sequence_type": provenance.get("sequence_type", ""),
                    "occultation_pass": provenance.get("occultation_pass", ""),
                    "sequence_visibility_fraction": provenance.get("visibility_fraction", ""),
                }
            )
            continue

        if n_minutes <= 0:
            results.append(
                {
                    "visit_id": seq.visit_id,
                    "sequence_id": seq.sequence_id,
                    "target": seq.target,
                    "priority": seq.priority,
                    "start_utc": seq.start.isoformat(),
                    "stop_utc": seq.stop.isoformat(),
                    "duration_minutes": n_minutes,
                    "status": "zero_duration",
                    "visible_minutes": 0,
                    "nonvisible_minutes": 0,
                    "sequence_type": provenance.get("sequence_type", ""),
                    "occultation_pass": provenance.get("occultation_pass", ""),
                    "sequence_visibility_fraction": provenance.get("visibility_fraction", ""),
                }
            )
            continue

        aligned = visibility_df["Visible"].reindex(minute_index, fill_value=False)
        visible_minutes = int(aligned.sum())
        nonvisible_minutes = int((~aligned).sum())
        first_problem = ""
        if nonvisible_minutes:
            first_problem = aligned.index[(~aligned).to_numpy()][0].isoformat()

        sequence_type = str(provenance.get("sequence_type", "") or "")
        tolerated_occultation_miss = (
            sequence_type == "occultation"
            and nonvisible_minutes <= occultation_nonvisible_tolerance_minutes
        )

        results.append(
            {
                "visit_id": seq.visit_id,
                "sequence_id": seq.sequence_id,
                "target": seq.target,
                "priority": seq.priority,
                "start_utc": seq.start.isoformat(),
                "stop_utc": seq.stop.isoformat(),
                "duration_minutes": n_minutes,
                "status": (
                    "ok"
                    if nonvisible_minutes == 0 or tolerated_occultation_miss
                    else "nonvisible_minutes_found"
                ),
                "visible_minutes": visible_minutes,
                "nonvisible_minutes": nonvisible_minutes,
                "sequence_type": sequence_type,
                "occultation_pass": provenance.get("occultation_pass", ""),
                "sequence_visibility_fraction": provenance.get("visibility_fraction", ""),
            }
        )

    return results


def write_csv(rows: list[dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "visit_id",
        "sequence_id",
        "target",
        "priority",
        "start_utc",
        "stop_utc",
        "duration_minutes",
        "status",
        "visible_minutes",
        "nonvisible_minutes",
        "sequence_type",
        "occultation_pass",
        "sequence_visibility_fraction",
    ]
    rounded_rows: list[dict[str, object]] = []
    for row in rows:
        row_copy = dict(row)
        value = row_copy.get("sequence_visibility_fraction")
        try:
            if value not in {"", None}:
                row_copy["sequence_visibility_fraction"] = round(float(value), 2)
        except (TypeError, ValueError):
            pass
        rounded_rows.append(row_copy)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rounded_rows)




def print_summary(rows: list[dict[str, object]], csv_path: Path) -> int:
    total = len(rows)
    ok = sum(row["status"] == "ok" for row in rows)
    missing = sum(row["status"] == "missing_visibility_file" for row in rows)
    failed = sum(row["status"] == "nonvisible_minutes_found" for row in rows)
    zero = sum(row["status"] == "zero_duration" for row in rows)
    nonvisible_minutes = sum(int(row["nonvisible_minutes"]) for row in rows)

    print("\nXML VISIBILITY VALIDATION")
    print("=" * 60)
    print(f"Total sequences:           {total}")
    print(f"Fully visible sequences:   {ok}")
    print(f"Missing visibility files:  {missing}")
    print(f"Sequences with violations: {failed}")
    print(f"Zero-duration sequences:   {zero}")
    print(f"Total non-visible minutes: {nonvisible_minutes}")
    print(f"CSV report:                {csv_path}")

    problem_rows = [
        row for row in rows
        if row["status"] in {"missing_visibility_file", "nonvisible_minutes_found"}
    ]
    if problem_rows:
        print("\nFirst issues:")
        for row in problem_rows[:10]:
            pass_info = str(row.get("occultation_pass", "")).strip()
            pass_suffix = f", source={pass_info}" if pass_info else ""
            seq_frac = row.get("sequence_visibility_fraction", "")
            frac_suffix = f", seq_vis_frac={seq_frac}" if seq_frac not in {"", None} else ""
            print(
                f"- {row['visit_id']} / {row['sequence_id']} / {row['target']}: "
                f"{row['status']} "
                f"(non-visible min={row['nonvisible_minutes']}{pass_suffix}{frac_suffix})"
            )
        return 1

    print("\nValidation passed: all XML sequences are visible.")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate final XML observation sequences against visibility parquet files.",
    )
    parser.add_argument("xml", type=Path, help="Path to Pandora_science_calendar.xml")
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Run data directory containing visibility parquet files.",
    )
    parser.add_argument(
        "--csv-out",
        type=Path,
        help="Optional CSV output path. Defaults next to the XML file.",
    )
    parser.add_argument(
        "--occultation-nonvisible-tolerance-minutes",
        type=int,
        default=3,
        help=(
            "Treat occultation XML sequences with up to this many non-visible "
            "minutes as acceptable. Default: 3."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sequences = parse_science_calendar(args.xml)
    provenance_df = _load_sequence_provenance(args.xml)
    csv_path = args.csv_out or args.xml.with_name(
        f"{args.xml.stem}_visibility_validation.csv"
    )
    results = validate_sequences(
        sequences,
        args.data_dir,
        provenance_df=provenance_df,
        occultation_nonvisible_tolerance_minutes=args.occultation_nonvisible_tolerance_minutes,
    )
    write_csv(results, csv_path)
    return print_summary(results, csv_path)


if __name__ == "__main__":
    raise SystemExit(main())
