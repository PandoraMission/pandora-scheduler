"""Automatic depth-ranked Target-of-Opportunity selection."""

from __future__ import annotations

import csv
import io
import logging
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd

from pandorascheduler_rework.utils.io import build_visibility_path

LOGGER = logging.getLogger(__name__)
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


@dataclass(frozen=True)
class TooSelectionResult:
    ranked_targets_path: Path
    ranked_windows_path: Path
    too_list_path: Path
    selected_targets: list[str]


def select_depth_ranked_toos(
    *,
    data_dir: Path,
    candidate_manifest: Path,
    main_targets: Iterable[str],
    start: datetime,
    end: datetime,
    coverage_min: float,
    depth_min_percent: float,
    top_n: int,
    query_depths: bool = True,
) -> TooSelectionResult:
    """Rank visible non-main exoplanets by transit depth and write ToO files."""

    data_dir = data_dir.resolve()
    ranked_targets_path = data_dir / "all_too_candidates_ranked_by_transit_depth.csv"
    ranked_windows_path = data_dir / "all_too_windows_ranked_by_transit_depth.csv"
    too_list_path = data_dir / "ToO_list.csv"
    depth_cache_path = data_dir / "transit_depth_cache.csv"

    excluded_targets = {str(target) for target in main_targets}
    windows = find_eligible_windows(
        data_dir=data_dir,
        manifest_path=candidate_manifest,
        start=pd.Timestamp(start),
        end=pd.Timestamp(end),
        coverage_min=coverage_min,
        excluded_targets=excluded_targets,
    )
    ranked_windows = windows.sort_values(
        ["target", "transit_coverage", "obs_window_start"],
        ascending=[True, False, True],
        ignore_index=True,
    )
    ranked_windows_path.parent.mkdir(parents=True, exist_ok=True)
    ranked_windows.to_csv(ranked_windows_path, index=False)

    ranked_targets = build_ranked_targets(windows)
    ranked_targets = enrich_with_depths(
        ranked_targets,
        depth_cache=depth_cache_path,
        query_archive=query_depths,
    )
    ranked_targets = rank_by_depth(ranked_targets)
    ranked_targets.to_csv(ranked_targets_path, index=False)

    if ranked_targets.empty:
        selected = ranked_targets
    else:
        selected = ranked_targets.loc[
            pd.to_numeric(ranked_targets["transit_depth_percent"], errors="coerce")
            >= depth_min_percent
        ].head(top_n)
    write_too_list(too_list_path, selected)

    return TooSelectionResult(
        ranked_targets_path=ranked_targets_path,
        ranked_windows_path=ranked_windows_path,
        too_list_path=too_list_path,
        selected_targets=selected["target"].astype(str).tolist(),
    )


def find_eligible_windows(
    *,
    data_dir: Path,
    manifest_path: Path,
    start: pd.Timestamp,
    end: pd.Timestamp,
    coverage_min: float,
    excluded_targets: set[str],
) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path)
    rows: list[dict[str, object]] = []

    if end <= start:
        raise ValueError(f"ToO end must be after start: {start} >= {end}")

    for _, target in manifest.iterrows():
        planet = str(target["Planet Name"])
        star = str(target["Star Name"])
        if planet in excluded_targets:
            continue

        visibility_path = build_visibility_path(data_dir / "targets", star, planet)
        if not visibility_path.exists():
            continue

        visibility = pd.read_parquet(visibility_path)
        required = {"Transit_Start_UTC", "Transit_Stop_UTC", "Transit_Coverage"}
        missing = required.difference(visibility.columns)
        if missing:
            raise ValueError(
                f"{visibility_path} missing required columns: {sorted(missing)}"
            )

        transit_start = pd.to_datetime(visibility["Transit_Start_UTC"])
        transit_stop = pd.to_datetime(visibility["Transit_Stop_UTC"])
        mask = (
            (transit_start < end)
            & (transit_stop > start)
            & (visibility["Transit_Coverage"] >= coverage_min)
        )
        for _, transit in visibility.loc[mask].iterrows():
            rows.append(
                {
                    "target": planet,
                    "star": star,
                    "priority": _float_or_none(target.get("Priority")) or 0.0,
                    "period_days": _float_or_none(target.get("Period (days)")),
                    "obs_window_start": transit["Transit_Start_UTC"],
                    "obs_window_stop": transit["Transit_Stop_UTC"],
                    "transit_coverage": float(transit["Transit_Coverage"]),
                    "saa_overlap": float(transit.get("SAA_Overlap", 0.0)),
                }
            )

    columns = [
        "target",
        "star",
        "priority",
        "period_days",
        "obs_window_start",
        "obs_window_stop",
        "transit_coverage",
        "saa_overlap",
    ]
    if not rows:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(rows, columns=columns).sort_values(
        ["transit_coverage", "priority", "obs_window_start"],
        ascending=[False, False, True],
        ignore_index=True,
    )


def build_ranked_targets(windows: pd.DataFrame) -> pd.DataFrame:
    if windows.empty:
        ranked = windows.copy()
        ranked.insert(0, "rank_by_transit_depth", [])
        ranked["eligible_window_count"] = []
        return ranked

    ranked = windows.drop_duplicates("target", keep="first").copy()
    counts = windows.groupby("target").size().rename("eligible_window_count")
    return ranked.merge(counts, on="target", how="left")


def rank_by_depth(ranked: pd.DataFrame) -> pd.DataFrame:
    if ranked.empty:
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
        return enriched

    sortable = ranked.copy()
    sortable["_depth_sort"] = pd.to_numeric(
        sortable["transit_depth_percent"],
        errors="coerce",
    ).fillna(-1.0)
    sortable = sortable.sort_values(
        ["_depth_sort", "transit_coverage", "priority", "obs_window_start"],
        ascending=[False, False, False, True],
        ignore_index=True,
    ).drop(columns=["_depth_sort"])
    sortable.insert(0, "rank_by_transit_depth", range(1, len(sortable) + 1))
    return sortable


def enrich_with_depths(
    ranked: pd.DataFrame,
    *,
    depth_cache: Path,
    query_archive: bool,
) -> pd.DataFrame:
    if ranked.empty:
        return ranked

    depths = load_depth_cache(depth_cache)
    targets = ranked["target"].astype(str).tolist()
    missing_targets = [
        target
        for target in targets
        if target not in depths
        or (
            query_archive
            and depths[target].depth_percent is None
            and depths[target].source_table in {"", "not found"}
        )
    ]
    if query_archive and missing_targets:
        depths.update(fetch_transit_depths(ranked, missing_targets))
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
            writer.writerow(depth_to_row(depths[target]))


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
        "select pl_name,pl_trandep,pl_trandeperr1,pl_trandeperr2 "
        f"from pscomppars where pl_name in ({quoted})"
    )
    try:
        table = tap_query(query)
    except Exception as exc:
        LOGGER.warning("pscomppars depth query failed: %s", exc)
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
    except Exception as exc:
        LOGGER.warning("TOI depth query failed: %s", exc)
        return {}

    depths: dict[str, TransitDepth] = {}
    for _, target_row in candidates.iterrows():
        target = str(target_row["target"])
        prefix_match = re.match(r"^TOI-(\d+)", target)
        period = _float_or_none(target_row.get("period_days"))
        if prefix_match is None or period is None:
            continue
        subset = toi_table.loc[toi_table["toipfx"].astype(str) == prefix_match.group(1)].copy()
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
            source_table="NASA Exoplanet Archive toi; matched by TOI prefix and period",
            archive_name=str(row.get("toidisplay", row.get("toi", "")) or ""),
        )
    return depths


def tap_query(query: str) -> pd.DataFrame:
    url = NASA_TAP_SYNC + "?" + urlencode({"query": query, "format": "csv"})
    with urlopen(url, timeout=60) as response:
        payload = response.read().decode("utf-8")
    return pd.read_csv(io.StringIO(payload))


def write_too_list(path: Path, selected: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Target", "Obs Window Start", "Obs Window Stop"],
        )
        writer.writeheader()
        for _, row in selected.iterrows():
            writer.writerow(
                {
                    "Target": row["target"],
                    "Obs Window Start": str(row["obs_window_start"]),
                    "Obs Window Stop": str(row["obs_window_stop"]),
                }
            )


def pandora_to_archive_planet_name(target: str) -> str:
    name = target.replace("_", " ")
    return re.sub(r"(?<! )([bcdefgh])$", r" \1", name)


def sql_quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _float_or_none(value: object) -> float | None:
    if value is None or pd.isna(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
