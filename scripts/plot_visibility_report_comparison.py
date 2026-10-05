#!/usr/bin/env python3
"""Plot a high-signal comparison between two visibility report text files."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROW_RE = re.compile(r"^\s*(\d{4})\s+(\d{3})\s+(.+?)\s{2,}([DN])\s+")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare two minute-by-minute visibility report text files and "
            "render a summary plot."
        )
    )
    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Short-term sched visibility report, e.g. PAN-SCICAL-...-violations.txt",
    )
    parser.add_argument(
        "--candidate",
        type=Path,
        required=True,
        help=(
            "Long-term sched visibility report, e.g. "
            "Pandora_science_calendar_sequence_provenance_visibility_report.txt"
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output PNG path",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional path to write per-sequence comparison metrics",
    )
    return parser.parse_args()


def _parse_report(path: Path, label: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            match = ROW_RE.match(line)
            if not match:
                continue
            visit_id, sequence_id, target, dn = match.groups()
            rows.append(
                {
                    "source": label,
                    "visit_id": int(visit_id),
                    "sequence_id": int(sequence_id),
                    "target": target.strip(),
                    "dn": dn,
                    "violation": "*** VIOLATION" in line,
                    "line_no": line_no,
                }
            )
    if not rows:
        raise ValueError(f"No sequence rows found in {path}")
    return pd.DataFrame(rows)


def _aggregate(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    grouped = (
        df.groupby(["visit_id", "sequence_id", "target"], as_index=False)
        .agg(
            minute_count=("violation", "size"),
            violation_count=("violation", "sum"),
            day_minutes=("dn", lambda s: int((s == "D").sum())),
            night_minutes=("dn", lambda s: int((s == "N").sum())),
        )
        .rename(
            columns={
                "minute_count": f"{prefix}_minutes",
                "violation_count": f"{prefix}_violations",
                "day_minutes": f"{prefix}_day_minutes",
                "night_minutes": f"{prefix}_night_minutes",
            }
        )
    )
    return grouped


def _label_sequence(row: pd.Series) -> str:
    return f"V{int(row['visit_id']):02d}-S{int(row['sequence_id']):03d} {row['target']}"


def _build_summary(reference: pd.DataFrame, candidate: pd.DataFrame) -> pd.DataFrame:
    ref = _aggregate(reference, "ref")
    cand = _aggregate(candidate, "cand")
    merged = ref.merge(
        cand,
        on=["visit_id", "sequence_id", "target"],
        how="outer",
    ).fillna(0)

    int_cols = [
        "ref_minutes",
        "ref_violations",
        "ref_day_minutes",
        "ref_night_minutes",
        "cand_minutes",
        "cand_violations",
        "cand_day_minutes",
        "cand_night_minutes",
    ]
    for col in int_cols:
        merged[col] = merged[col].astype(int)

    merged["minute_delta"] = merged["cand_minutes"] - merged["ref_minutes"]
    merged["violation_delta"] = merged["cand_violations"] - merged["ref_violations"]
    merged["sequence_label"] = merged.apply(_label_sequence, axis=1)
    return merged.sort_values(["visit_id", "sequence_id", "target"]).reset_index(
        drop=True
    )


def _plot_totals(ax: plt.Axes, reference: pd.DataFrame, candidate: pd.DataFrame) -> None:
    totals = pd.DataFrame(
        [
            {
                "report": "Short-term sched",
                "minutes": len(reference),
                "violations": int(reference["violation"].sum()),
            },
            {
                "report": "Long-term sched",
                "minutes": len(candidate),
                "violations": int(candidate["violation"].sum()),
            },
        ]
    )
    x = range(len(totals))
    width = 0.36
    ax.bar(
        [i - width / 2 for i in x],
        totals["minutes"],
        width=width,
        color="#9ecae1",
        label="Minutes",
    )
    ax.bar(
        [i + width / 2 for i in x],
        totals["violations"],
        width=width,
        color="#ef8a62",
        label="Violations",
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels(totals["report"])
    ax.set_title("Overall Counts")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    for i, row in totals.iterrows():
        ax.text(i - width / 2, row["minutes"] + 20, str(row["minutes"]), ha="center", fontsize=8)
        ax.text(i + width / 2, row["violations"] + 8, str(row["violations"]), ha="center", fontsize=8)


def _plot_visit_violations(ax: plt.Axes, reference: pd.DataFrame, candidate: pd.DataFrame) -> None:
    ref_visit = reference.groupby("visit_id")["violation"].sum().rename("Short-term sched")
    cand_visit = candidate.groupby("visit_id")["violation"].sum().rename("Long-term sched")
    visit_df = pd.concat([ref_visit, cand_visit], axis=1).fillna(0).astype(int).reset_index()
    x = range(len(visit_df))
    width = 0.42
    ax.bar(
        [i - width / 2 for i in x],
        visit_df["Short-term sched"],
        width=width,
        color="#6baed6",
        label="Short-term sched",
    )
    ax.bar(
        [i + width / 2 for i in x],
        visit_df["Long-term sched"],
        width=width,
        color="#fb6a4a",
        label="Long-term sched",
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{int(v):02d}" for v in visit_df["visit_id"]], rotation=90, fontsize=7)
    ax.set_title("Violations Per Visit")
    ax.set_xlabel("Visit ID")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)


def _plot_sequence_scatter(ax: plt.Axes, summary: pd.DataFrame) -> None:
    exact = (summary["minute_delta"] == 0) & (summary["violation_delta"] == 0)
    ax.scatter(
        summary.loc[exact, "ref_minutes"],
        summary.loc[exact, "cand_minutes"],
        s=24,
        alpha=0.7,
        color="#74c476",
        label="Exact match",
    )
    ax.scatter(
        summary.loc[~exact, "ref_minutes"],
        summary.loc[~exact, "cand_minutes"],
        s=28,
        alpha=0.85,
        color="#e34a33",
        label="Different",
    )
    min_val = 0
    max_val = max(summary["ref_minutes"].max(), summary["cand_minutes"].max()) + 2
    ax.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="#666666", linewidth=1)
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.set_xlabel("Short-term sched minutes per sequence")
    ax.set_ylabel("Long-term sched minutes per sequence")
    ax.set_title("Per-Sequence Minute Count")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="lower right")


def _plot_top_deltas(ax: plt.Axes, summary: pd.DataFrame) -> None:
    ranked = summary.copy()
    ranked["abs_delta"] = ranked["minute_delta"].abs() + ranked["violation_delta"].abs()
    ranked = ranked[ranked["abs_delta"] > 0].sort_values(
        ["abs_delta", "violation_delta", "minute_delta"], ascending=False
    )
    top = ranked.head(12).iloc[::-1]
    if top.empty:
        ax.text(0.5, 0.5, "No per-sequence differences", ha="center", va="center", fontsize=11)
        ax.axis("off")
        return
    y = range(len(top))
    ax.barh(y, top["minute_delta"], color="#9ecae1", label="Minute delta")
    ax.barh(y, top["violation_delta"], color="#ef8a62", alpha=0.85, label="Violation delta")
    ax.set_yticks(list(y))
    ax.set_yticklabels(top["sequence_label"], fontsize=7)
    ax.axvline(0, color="#666666", linewidth=1)
    ax.set_title("Largest Sequence-Level Deltas")
    ax.set_xlabel("Long-term sched - short-term sched")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(fontsize=8)


def _render_plot(
    reference: pd.DataFrame,
    candidate: pd.DataFrame,
    summary: pd.DataFrame,
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    _plot_totals(axes[0, 0], reference, candidate)
    _plot_visit_violations(axes[0, 1], reference, candidate)
    _plot_sequence_scatter(axes[1, 0], summary)
    _plot_top_deltas(axes[1, 1], summary)

    exact_matches = int(((summary["minute_delta"] == 0) & (summary["violation_delta"] == 0)).sum())
    total_sequences = len(summary)
    fig.suptitle(
        "Visibility Report Comparison\n"
        f"Exact sequence matches: {exact_matches}/{total_sequences}",
        fontsize=15,
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = _parse_args()
    reference = _parse_report(args.reference, "reference")
    candidate = _parse_report(args.candidate, "candidate")
    summary = _build_summary(reference, candidate)
    if args.summary_csv is not None:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary_csv, index=False)
    _render_plot(reference, candidate, summary, args.out)


if __name__ == "__main__":
    main()
