#!/usr/bin/env python3
"""Plot science-sequence provenance rows affected by soft-ST tail extension."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot science rows affected by soft-ST tail extension from XML "
            "sequence provenance CSVs."
        )
    )
    parser.add_argument(
        "--soft-provenance",
        type=Path,
        required=True,
        help="Path to Pandora_science_calendar_soft_ST_sequence_provenance.csv",
    )
    parser.add_argument(
        "--base-provenance",
        type=Path,
        help=(
            "Optional baseline Pandora_science_calendar_sequence_provenance.csv. "
            "Required only for --layout two-panel."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output PNG path",
    )
    parser.add_argument(
        "--layout",
        choices=("right-only", "two-panel"),
        default="right-only",
        help="Figure layout to render (default: right-only)",
    )
    return parser.parse_args()


def _load_affected_soft_rows(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df["start_dt"] = pd.to_datetime(df["start_utc"], utc=True)
    df["science_soft_tail_minutes"] = pd.to_numeric(
        df.get("science_soft_tail_minutes"), errors="coerce"
    ).fillna(0)
    df["visible_minutes"] = pd.to_numeric(
        df.get("visible_minutes"), errors="coerce"
    ).fillna(0)
    df["non_visible_minutes"] = pd.to_numeric(
        df.get("non_visible_minutes"), errors="coerce"
    ).fillna(0)
    affected = df[
        (df["sequence_type"] == "science") & (df["science_soft_tail_minutes"] > 0)
    ].copy()
    affected = affected.sort_values(["visit_id", "start_dt"]).reset_index(drop=True)
    affected["label"] = affected.apply(
        lambda row: (
            f"V{int(row['visit_id']):02d} · {row['target']} · "
            f"{pd.Timestamp(row['start_dt']).strftime('%m-%d %H:%M')}"
        ),
        axis=1,
    )
    return affected


def _build_two_panel_match(
    base_path: Path,
    soft_affected: pd.DataFrame,
) -> pd.DataFrame:
    base = pd.read_csv(base_path, low_memory=False)
    base["start_dt"] = pd.to_datetime(base["start_utc"], utc=True)

    base_sci = base[base["sequence_type"] == "science"].copy()
    base_sci = base_sci.sort_values(["visit_id", "target", "start_dt"]).copy()
    base_sci["match_idx"] = base_sci.groupby(["visit_id", "target"]).cumcount()

    soft = soft_affected.copy()
    soft = soft.sort_values(["visit_id", "target", "start_dt"]).copy()
    soft["match_idx"] = soft.groupby(["visit_id", "target"]).cumcount()

    merged = soft.merge(
        base_sci[
            [
                "visit_id",
                "target",
                "match_idx",
                "duration_minutes",
                "visible_minutes",
                "non_visible_minutes",
            ]
        ].rename(
            columns={
                "duration_minutes": "base_duration_minutes",
                "visible_minutes": "base_visible_minutes",
                "non_visible_minutes": "base_non_visible_minutes",
            }
        ),
        on=["visit_id", "target", "match_idx"],
        how="left",
    )
    merged["duration_minutes"] = pd.to_numeric(
        merged["duration_minutes"], errors="coerce"
    )
    merged["base_duration_minutes"] = pd.to_numeric(
        merged["base_duration_minutes"], errors="coerce"
    )
    merged["delta_duration"] = (
        merged["duration_minutes"] - merged["base_duration_minutes"]
    )
    return merged.sort_values(["visit_id", "start_dt"]).reset_index(drop=True)


def _plot_right_only(affected: pd.DataFrame, out_path: Path) -> None:
    n_rows = len(affected)
    fig_h = max(8, 0.34 * n_rows + 1.8)
    fig, ax = plt.subplots(1, 1, figsize=(10.5, fig_h))
    y = list(range(n_rows))

    ax.barh(
        y,
        affected["visible_minutes"],
        height=0.56,
        color="#4caf50",
        label="Visible minutes",
    )
    ax.barh(
        y,
        affected["non_visible_minutes"],
        left=affected["visible_minutes"],
        height=0.56,
        color="#ffb74d",
        label="Nominal non-visible minutes",
    )

    for i, row in affected.iterrows():
        total = float(row["visible_minutes"]) + float(row["non_visible_minutes"])
        ax.text(
            total + 0.35,
            i,
            f"tail {int(round(float(row['science_soft_tail_minutes'])))}m",
            va="center",
            fontsize=8,
            color="#7a1f1f",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(affected["label"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Soft-ST Sequence Minutes")
    ax.set_title("Soft-ST Science Sequence Composition For Affected Rows")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(loc="lower right", fontsize=9)

    visit_summary = ", ".join(str(v) for v in affected["visit_id"].drop_duplicates())
    fig.suptitle(
        "Science Sequences Affected by Soft-ST Tail Extension\n"
        f"Affected visit_ids: {visit_summary}",
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_two_panel(merged: pd.DataFrame, out_path: Path) -> None:
    n_rows = len(merged)
    fig_h = max(8, 0.34 * n_rows + 2.2)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(17, fig_h),
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )
    ax0, ax1 = axes
    y = list(range(n_rows))

    ax0.barh(
        y,
        merged["base_duration_minutes"],
        height=0.34,
        color="#c7ced6",
        label="Baseline duration",
    )
    ax0.barh(
        [i + 0.18 for i in y],
        merged["duration_minutes"],
        height=0.34,
        color="#2b7a78",
        label="Soft-ST duration",
    )
    for i, row in merged.iterrows():
        base_val = row["base_duration_minutes"]
        soft_val = row["duration_minutes"]
        if pd.notna(base_val) and pd.notna(soft_val):
            ax0.plot([base_val, soft_val], [i, i + 0.18], color="#5b646d", linewidth=1)
            ax0.text(
                max(base_val, soft_val) + 0.35,
                i + 0.09,
                f"+{int(round(float(row['delta_duration'])))}m",
                va="center",
                fontsize=8,
                color="#1f2d3d",
            )
    ax0.set_yticks([i + 0.09 for i in y])
    ax0.set_yticklabels(merged["label"], fontsize=8)
    ax0.invert_yaxis()
    ax0.set_xlabel("Sequence Duration (min)")
    ax0.set_title("Baseline vs Soft-ST Science Sequence Duration")
    ax0.grid(axis="x", alpha=0.25)
    ax0.legend(loc="lower right", fontsize=9)

    ax1.barh(
        y,
        merged["visible_minutes"],
        height=0.52,
        color="#4caf50",
        label="Visible minutes",
    )
    ax1.barh(
        y,
        merged["non_visible_minutes"],
        left=merged["visible_minutes"],
        height=0.52,
        color="#ffb74d",
        label="Nominal non-visible minutes",
    )
    for i, row in merged.iterrows():
        total = float(row["visible_minutes"]) + float(row["non_visible_minutes"])
        ax1.text(
            total + 0.35,
            i,
            f"tail {int(round(float(row['science_soft_tail_minutes'])))}m",
            va="center",
            fontsize=8,
            color="#7a1f1f",
        )
    ax1.set_yticks(y)
    ax1.set_yticklabels([])
    ax1.invert_yaxis()
    ax1.set_xlabel("Soft-ST Sequence Minutes")
    ax1.set_title("Soft-ST Sequence Composition")
    ax1.grid(axis="x", alpha=0.25)
    ax1.legend(loc="lower right", fontsize=9)

    visit_summary = ", ".join(str(v) for v in merged["visit_id"].drop_duplicates())
    fig.suptitle(
        "Science Sequences Affected by Soft-ST Tail Extension\n"
        f"Affected visit_ids: {visit_summary}",
        fontsize=14,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = _parse_args()
    affected = _load_affected_soft_rows(args.soft_provenance)
    if affected.empty:
        raise ValueError(
            "No affected science rows found in the soft provenance file "
            "(science_soft_tail_minutes > 0)."
        )

    if args.layout == "two-panel":
        if args.base_provenance is None:
            raise ValueError("--base-provenance is required for --layout two-panel")
        merged = _build_two_panel_match(args.base_provenance, affected)
        _plot_two_panel(merged, args.out)
    else:
        _plot_right_only(affected, args.out)

    print(args.out)


if __name__ == "__main__":
    main()
