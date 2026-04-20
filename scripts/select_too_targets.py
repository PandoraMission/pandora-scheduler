#!/usr/bin/env python3
"""Select candidate exoplanet Targets of Opportunity from visibility products.

This script scans a scheduler ``data_*`` directory, ranks non-main-science
exoplanets by their best visible transit coverage in a requested time window,
adds transit depths from NASA Exoplanet Archive when available, and writes:

* a ranked per-target CSV
* a ranked per-window CSV
* a scheduler-readable ``ToO_list.csv`` for the top N distinct targets
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd


NASA_TAP_SYNC = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"


@dataclass(frozen=True)
class TransitDepth:
    target: str
    depth_percent: float | None
    err_plus_percent: float | None
    err_minus_percent: float | None
    source_table: str
    archive_name: str

    @property
    def depth_ppm(self) -> float | None:
        if self.depth_percent is None:
            return None
        return self.depth_percent * 10_000.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank eligible exoplanet ToOs and write ToO_list.csv."
    )
    parser.add_argument(
        "--data-dir",
        required=True,
        type=Path,
        help="Scheduler data_* directory containing exoplanet_targets.csv and targets/.",
    )
    parser.add_argument(
        "--start",
        required=True,
        help="Inclusive UTC window start, e.g. 2026-04-27 or '2026-04-27 00:00:00'.",
    )
    parser.add_argument(
        "--end",
        required=True,
        help="Exclusive UTC window end, e.g. 2026-05-04 or '2026-05-04 00:00:00'.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.4,
        help="Minimum Transit_Coverage to include (default: 0.4).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of distinct targets to write to ToO_list.csv (default: 5).",
    )
    parser.add_argument(
        "--exclude-priorities",
        type=Path,
        help=(
            "Optional exoplanet_priorities.csv whose target column should be "
            "excluded, e.g. the current 10_HJs priorities file."
        ),
    )
    parser.add_argument(
        "--ranked-out",
        type=Path,
        help="Output ranked per-target CSV (default: <data-dir>/eligible_too_targets_ranked_by_coverage.csv).",
    )
    parser.add_argument(
        "--windows-out",
        type=Path,
        help="Output ranked per-window CSV (default: <data-dir>/eligible_too_windows_ranked_by_coverage.csv).",
    )
    parser.add_argument(
        "--too-out",
        type=Path,
        help="Output scheduler ToO_list.csv (default: <data-dir>/ToO_list.csv).",
    )
    parser.add_argument(
        "--depth-cache",
        type=Path,
        help=(
            "Transit-depth cache CSV path. Existing values are reused; fetched "
            "values are written back. Default: <data-dir>/transit_depth_cache.csv."
        ),
    )
    parser.add_argument(
        "--no-depth-query",
        action="store_true",
        help="Do not query NASA Exoplanet Archive; use only --depth-cache if present.",
    )
    parser.add_argument(
        "--no-too-list",
        action="store_true",
        help="Do not write ToO_list.csv; only write ranked CSV outputs.",
    )
    parser.add_argument(
        "--main-target-definition-dir",
        type=Path,
        help=(
            "Optional target-definition base for the main science targets, e.g. "
            "10_HJs. Used with --candidate-target-definition-dir and "
            "--target-definition-out to build a combined main+ToO target list."
        ),
    )
    parser.add_argument(
        "--candidate-target-definition-dir",
        type=Path,
        help=(
            "Optional full target-definition base to copy selected ToO target "
            "definitions from, e.g. PandoraTargetList/target_definition_files."
        ),
    )
    parser.add_argument(
        "--target-definition-out",
        type=Path,
        help=(
            "Optional output target-definition base containing main targets plus "
            "selected ToOs. Selected ToOs are written with transits_req=0 so "
            "they do not compete as ordinary science targets."
        ),
    )
    args = parser.parse_args()
    if args.threshold < 0 or args.threshold > 1:
        parser.error("--threshold must be between 0 and 1")
    if args.top_n < 0:
        parser.error("--top-n must be >= 0")
    if args.target_definition_out and (
        not args.main_target_definition_dir or not args.candidate_target_definition_dir
    ):
        parser.error(
            "--target-definition-out requires --main-target-definition-dir and "
            "--candidate-target-definition-dir"
        )
    return args


def read_excluded_targets(path: Path | None) -> set[str]:
    if path is None:
        return set()
    if not path.exists():
        raise SystemExit(f"--exclude-priorities not found: {path}")
    table = pd.read_csv(path, comment="#")
    if "target" not in table.columns:
        raise SystemExit(f"--exclude-priorities is missing target column: {path}")
    return set(table["target"].dropna().astype(str))


def find_eligible_windows(
    data_dir: Path,
    start: pd.Timestamp,
    end: pd.Timestamp,
    threshold: float,
    excluded_targets: set[str],
) -> pd.DataFrame:
    manifest_path = data_dir / "exoplanet_targets.csv"
    if not manifest_path.exists():
        raise SystemExit(f"Missing manifest: {manifest_path}")
    manifest = pd.read_csv(manifest_path)

    if end <= start:
        raise SystemExit(f"--end must be after --start: {start} >= {end}")

    rows: list[dict[str, object]] = []
    for _, target in manifest.iterrows():
        planet = str(target["Planet Name"])
        star = str(target["Star Name"])
        if planet in excluded_targets:
            continue

        visibility_path = (
            data_dir
            / "targets"
            / star
            / planet
            / f"Visibility for {planet}.parquet"
        )
        if not visibility_path.exists():
            continue

        visibility = pd.read_parquet(visibility_path)
        required = {"Transit_Start_UTC", "Transit_Stop_UTC", "Transit_Coverage"}
        if not required.issubset(set(visibility.columns)):
            missing = sorted(required.difference(visibility.columns))
            raise SystemExit(f"{visibility_path} missing columns: {missing}")

        transit_start = pd.to_datetime(visibility["Transit_Start_UTC"])
        transit_stop = pd.to_datetime(visibility["Transit_Stop_UTC"])
        mask = (
            (transit_start < end)
            & (transit_stop > start)
            & (visibility["Transit_Coverage"] >= threshold)
        )
        for _, transit in visibility.loc[mask].iterrows():
            rows.append(
                {
                    "target": planet,
                    "star": star,
                    "priority": float(target.get("Priority", 0.0)),
                    "period_days": _float_or_none(target.get("Period (days)")),
                    "obs_window_start": transit["Transit_Start_UTC"],
                    "obs_window_stop": transit["Transit_Stop_UTC"],
                    "transit_coverage": float(transit["Transit_Coverage"]),
                    "saa_overlap": float(transit.get("SAA_Overlap", 0.0)),
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                "target",
                "star",
                "priority",
                "period_days",
                "obs_window_start",
                "obs_window_stop",
                "transit_coverage",
                "saa_overlap",
            ]
        )

    return pd.DataFrame(rows).sort_values(
        ["transit_coverage", "priority", "obs_window_start"],
        ascending=[False, False, True],
        ignore_index=True,
    )


def build_ranked_targets(windows: pd.DataFrame) -> pd.DataFrame:
    if windows.empty:
        ranked = windows.copy()
        ranked.insert(0, "rank_by_transit_coverage", [])
        ranked["eligible_window_count"] = []
        return ranked
    ranked = windows.drop_duplicates("target", keep="first").copy()
    ranked.insert(0, "rank_by_transit_coverage", range(1, len(ranked) + 1))
    counts = windows.groupby("target").size().rename("eligible_window_count")
    ranked = ranked.merge(counts, on="target", how="left")
    return ranked


def load_depth_cache(path: Path) -> dict[str, TransitDepth]:
    if not path.exists():
        return {}
    table = pd.read_csv(path)
    depths: dict[str, TransitDepth] = {}
    for _, row in table.iterrows():
        target = str(row.get("target", "") or "")
        if not target:
            continue
        depths[target] = TransitDepth(
            target=target,
            depth_percent=_float_or_none(row.get("transit_depth_percent")),
            err_plus_percent=_float_or_none(row.get("transit_depth_err_plus_percent")),
            err_minus_percent=_float_or_none(row.get("transit_depth_err_minus_percent")),
            source_table=str(row.get("transit_depth_source_table", "") or ""),
            archive_name=str(row.get("transit_depth_archive_name", "") or ""),
        )
    return depths


def write_depth_cache(path: Path, depths: dict[str, TransitDepth]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "target",
                "transit_depth_percent",
                "transit_depth_ppm",
                "transit_depth_err_plus_percent",
                "transit_depth_err_minus_percent",
                "transit_depth_source_table",
                "transit_depth_archive_name",
            ],
        )
        writer.writeheader()
        for target in sorted(depths):
            depth = depths[target]
            writer.writerow(depth_to_row(depth))


def depth_to_row(depth: TransitDepth) -> dict[str, object]:
    return {
        "target": depth.target,
        "transit_depth_percent": depth.depth_percent,
        "transit_depth_ppm": depth.depth_ppm,
        "transit_depth_err_plus_percent": depth.err_plus_percent,
        "transit_depth_err_minus_percent": depth.err_minus_percent,
        "transit_depth_source_table": depth.source_table,
        "transit_depth_archive_name": depth.archive_name,
    }


def enrich_with_depths(
    ranked: pd.DataFrame,
    depth_cache: Path,
    query_archive: bool,
) -> pd.DataFrame:
    if ranked.empty:
        return ranked

    depths = load_depth_cache(depth_cache)
    missing_targets = [
        target for target in ranked["target"].astype(str).tolist()
        if target not in depths
        or (
            query_archive
            and depths[target].depth_percent is None
            and depths[target].source_table in {"", "not found"}
        )
    ]
    if query_archive and missing_targets:
        fetched = fetch_transit_depths(ranked, missing_targets)
        depths.update(fetched)
        write_depth_cache(depth_cache, depths)

    enriched = ranked.copy()
    for column in [
        "transit_depth_percent",
        "transit_depth_ppm",
        "transit_depth_err_plus_percent",
        "transit_depth_err_minus_percent",
        "transit_depth_source_table",
        "transit_depth_archive_name",
    ]:
        enriched[column] = pd.NA

    for index, row in enriched.iterrows():
        depth = depths.get(str(row["target"]))
        if depth is None:
            continue
        for key, value in depth_to_row(depth).items():
            if key != "target":
                enriched.at[index, key] = value
    return enriched


def fetch_transit_depths(
    ranked: pd.DataFrame,
    targets: Iterable[str],
) -> dict[str, TransitDepth]:
    requested = set(targets)
    archive_names = {target: pandora_to_archive_planet_name(target) for target in requested}
    depths = fetch_pscomppars_depths(archive_names)

    still_missing = sorted(requested.difference(depths))
    if still_missing:
        depths.update(fetch_toi_depths_by_period(ranked, still_missing))

    for target in still_missing:
        if target not in depths:
            depths[target] = TransitDepth(
                target=target,
                depth_percent=None,
                err_plus_percent=None,
                err_minus_percent=None,
                source_table="not found",
                archive_name=archive_names[target],
            )
    return depths


def fetch_pscomppars_depths(
    archive_names: dict[str, str],
) -> dict[str, TransitDepth]:
    if not archive_names:
        return {}
    quoted = ",".join(sql_quote(name) for name in archive_names.values())
    query = (
        "select pl_name,hostname,pl_trandep,pl_trandeperr1,pl_trandeperr2 "
        f"from pscomppars where pl_name in ({quoted})"
    )
    try:
        table = tap_query(query)
    except Exception as exc:  # noqa: BLE001 - show a useful warning and continue.
        print(f"WARNING: pscomppars depth query failed: {exc}", file=sys.stderr)
        return {}

    reverse = {name: target for target, name in archive_names.items()}
    depths: dict[str, TransitDepth] = {}
    for _, row in table.iterrows():
        archive_name = str(row.get("pl_name", "") or "")
        target = reverse.get(archive_name)
        if target is None:
            continue
        depths[target] = TransitDepth(
            target=target,
            depth_percent=_float_or_none(row.get("pl_trandep")),
            err_plus_percent=_float_or_none(row.get("pl_trandeperr1")),
            err_minus_percent=_float_or_none(row.get("pl_trandeperr2")),
            source_table="NASA Exoplanet Archive pscomppars",
            archive_name=archive_name,
        )
    return depths


def fetch_toi_depths_by_period(
    ranked: pd.DataFrame,
    targets: Iterable[str],
) -> dict[str, TransitDepth]:
    candidates = ranked[ranked["target"].isin(list(targets))].copy()
    toi_prefixes = sorted(
        {
            match.group(1)
            for target in candidates["target"].astype(str)
            for match in [re.match(r"^TOI-(\d+)", target)]
            if match
        }
    )
    if not toi_prefixes:
        return {}

    query = (
        "select toi,toidisplay,toipfx,pl_orbper,pl_trandep,pl_trandeperr1,pl_trandeperr2 "
        "from toi where toipfx in ("
        + ",".join(sql_quote(prefix) for prefix in toi_prefixes)
        + ")"
    )
    try:
        toi_table = tap_query(query)
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: TOI depth query failed: {exc}", file=sys.stderr)
        return {}

    depths: dict[str, TransitDepth] = {}
    for _, target_row in candidates.iterrows():
        target = str(target_row["target"])
        prefix_match = re.match(r"^TOI-(\d+)", target)
        period = _float_or_none(target_row.get("period_days"))
        if prefix_match is None or period is None:
            continue
        prefix = prefix_match.group(1)
        subset = toi_table.loc[toi_table["toipfx"].astype(str) == prefix].copy()
        if subset.empty or "pl_orbper" not in subset:
            continue
        subset["period_delta"] = (pd.to_numeric(subset["pl_orbper"]) - period).abs()
        best = subset.sort_values("period_delta").head(1)
        if best.empty:
            continue
        row = best.iloc[0]
        depth_ppm = _float_or_none(row.get("pl_trandep"))
        err_plus_ppm = _float_or_none(row.get("pl_trandeperr1"))
        err_minus_ppm = _float_or_none(row.get("pl_trandeperr2"))
        depths[target] = TransitDepth(
            target=target,
            depth_percent=None if depth_ppm is None else depth_ppm / 10_000.0,
            err_plus_percent=None if err_plus_ppm is None else err_plus_ppm / 10_000.0,
            err_minus_percent=None if err_minus_ppm is None else err_minus_ppm / 10_000.0,
            source_table=(
                "NASA Exoplanet Archive toi; matched by TOI prefix and period"
            ),
            archive_name=str(row.get("toidisplay", row.get("toi", "")) or ""),
        )
    return depths


def tap_query(query: str) -> pd.DataFrame:
    url = NASA_TAP_SYNC + "?" + urlencode({"query": query, "format": "csv"})
    with urlopen(url, timeout=60) as response:
        payload = response.read().decode("utf-8")
    return pd.read_csv(io.StringIO(payload))


def pandora_to_archive_planet_name(target: str) -> str:
    name = target.replace("_", " ")
    # Split the final planet letter to match Exoplanet Archive names:
    # LTT_3780b -> LTT 3780 b, HD_93963_Ab -> HD 93963 A b,
    # 55_Cnce -> 55 Cnc e.
    name = re.sub(r"(?<! )([bcdefgh])$", r" \1", name)
    return name


def sql_quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _float_or_none(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def write_too_list(path: Path, top_targets: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Target", "Obs Window Start", "Obs Window Stop"],
        )
        writer.writeheader()
        for _, row in top_targets.iterrows():
            writer.writerow(
                {
                    "Target": row["target"],
                    "Obs Window Start": str(row["obs_window_start"]),
                    "Obs Window Stop": str(row["obs_window_stop"]),
                }
            )


def write_combined_target_definitions(
    main_dir: Path,
    candidate_dir: Path,
    out_dir: Path,
    top_targets: pd.DataFrame,
) -> None:
    main_dir = main_dir.expanduser().resolve()
    candidate_dir = candidate_dir.expanduser().resolve()
    out_dir = out_dir.expanduser().resolve()

    main_exoplanet_dir = main_dir / "exoplanet"
    candidate_exoplanet_dir = candidate_dir / "exoplanet"
    out_exoplanet_dir = out_dir / "exoplanet"

    for path in (main_exoplanet_dir, candidate_exoplanet_dir):
        if not path.is_dir():
            raise SystemExit(f"Target definition directory not found: {path}")

    out_exoplanet_dir.mkdir(parents=True, exist_ok=True)

    selected_targets = top_targets["target"].astype(str).tolist()
    copied_targets: set[str] = set()

    for json_path in sorted(main_exoplanet_dir.glob("*_target_definition.json")):
        shutil.copy2(json_path, out_exoplanet_dir / json_path.name)
        copied_targets.add(json_path.name.removesuffix("_target_definition.json"))

    for target in selected_targets:
        json_path = candidate_exoplanet_dir / f"{target}_target_definition.json"
        if not json_path.is_file():
            raise SystemExit(f"Selected ToO target definition not found: {json_path}")
        shutil.copy2(json_path, out_exoplanet_dir / json_path.name)
        copied_targets.add(target)

    main_priorities = _read_priority_rows(main_exoplanet_dir / "exoplanet_priorities.csv")
    candidate_priorities = _read_priority_rows(
        candidate_exoplanet_dir / "exoplanet_priorities.csv"
    )
    candidate_by_target = {
        str(row["target"]): row for _, row in candidate_priorities.iterrows()
    }

    rows: list[dict[str, object]] = []
    for _, row in main_priorities.iterrows():
        target = str(row["target"])
        if target not in copied_targets:
            continue
        rows.append(
            {
                "target": target,
                "priority": _float_or_none(row.get("priority")) or 0.0,
                "transits_req": _float_or_none(row.get("transits_req")) or 0.0,
                "transits_obs": _float_or_none(row.get("transits_obs")) or 0.0,
                "transits_rem": _float_or_none(row.get("transits_rem")) or 0.0,
            }
        )

    main_targets = {str(row["target"]) for row in rows}
    for target in selected_targets:
        if target in main_targets:
            continue
        source = candidate_by_target.get(target, {})
        rows.append(
            {
                "target": target,
                "priority": _float_or_none(source.get("priority")) or 0.0,
                "transits_req": 0.0,
                "transits_obs": 0.0,
                "transits_rem": 0.0,
            }
        )

    priority_path = out_exoplanet_dir / "exoplanet_priorities.csv"
    with priority_path.open("w", newline="") as handle:
        handle.write("# Prioritization file for main science targets plus selected ToOs\n")
        handle.write("# Version: 1.0.0\n")
        handle.write("# Updated by: select_too_targets.py\n")
        handle.write("# Selected ToOs use transits_req=0 and are scheduled only via ToO_list.csv\n")
        handle.write("#\n")
        handle.write("# Column 1: Rank\n")
        handle.write("# Column 2: Target Name\n")
        handle.write("# Column 3: Priority\n")
        handle.write("# Column 4: Number of Transits Requested\n")
        handle.write("# Column 5: Number of Transits Observed\n")
        handle.write("# Column 6: Number of Transits Remaining\n")
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "target",
                "priority",
                "transits_req",
                "transits_obs",
                "transits_rem",
            ],
        )
        writer.writeheader()
        for rank, row in enumerate(rows):
            writer.writerow({"rank": rank, **row})

    for name in [
        "auxiliary-standard",
        "monitoring-standard",
        "occultation-standard",
        "nirda_readout_schemes.json",
        "vda_readout_schemes.json",
    ]:
        destination = out_dir / name
        if destination.exists() or destination.is_symlink():
            continue
        source = main_dir / name
        if not source.exists():
            source = candidate_dir / name
        if source.exists():
            destination.symlink_to(source)


def _read_priority_rows(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise SystemExit(f"Priority table not found: {path}")
    return pd.read_csv(path, comment="#")


def main() -> int:
    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    if not data_dir.is_dir():
        raise SystemExit(f"--data-dir not found: {data_dir}")

    ranked_out = (
        args.ranked_out.expanduser().resolve()
        if args.ranked_out
        else data_dir / "eligible_too_targets_ranked_by_coverage.csv"
    )
    windows_out = (
        args.windows_out.expanduser().resolve()
        if args.windows_out
        else data_dir / "eligible_too_windows_ranked_by_coverage.csv"
    )
    too_out = (
        args.too_out.expanduser().resolve()
        if args.too_out
        else data_dir / "ToO_list.csv"
    )
    depth_cache = (
        args.depth_cache.expanduser().resolve()
        if args.depth_cache
        else data_dir / "transit_depth_cache.csv"
    )

    for output_path in (ranked_out, windows_out, too_out, depth_cache):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    excluded = read_excluded_targets(args.exclude_priorities)
    windows = find_eligible_windows(
        data_dir=data_dir,
        start=pd.Timestamp(args.start),
        end=pd.Timestamp(args.end),
        threshold=args.threshold,
        excluded_targets=excluded,
    )
    windows.to_csv(windows_out, index=False)

    ranked = build_ranked_targets(windows)
    ranked = enrich_with_depths(
        ranked,
        depth_cache=depth_cache,
        query_archive=not args.no_depth_query,
    )
    ranked.to_csv(ranked_out, index=False)

    if not args.no_too_list:
        write_too_list(too_out, ranked.head(args.top_n))

    if args.target_definition_out:
        write_combined_target_definitions(
            args.main_target_definition_dir,
            args.candidate_target_definition_dir,
            args.target_definition_out,
            ranked.head(args.top_n),
        )

    print(f"Eligible windows: {len(windows)}")
    print(f"Eligible targets: {len(ranked)}")
    print(f"Wrote ranked targets: {ranked_out}")
    print(f"Wrote ranked windows: {windows_out}")
    if not args.no_too_list:
        print(f"Wrote ToO list: {too_out}")
    if args.target_definition_out:
        print(f"Wrote target definitions: {args.target_definition_out}")
    if len(ranked) < args.top_n:
        print(
            f"WARNING: only {len(ranked)} eligible targets for top-n={args.top_n}",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
