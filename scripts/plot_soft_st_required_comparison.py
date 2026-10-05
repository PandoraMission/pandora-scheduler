#!/usr/bin/env python3
"""Compare soft-ST provenance outputs across science_soft_st_required settings."""

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
            "Plot a comparison of soft-ST-affected science rows across multiple "
            "Pandora_science_calendar_soft_ST_sequence_provenance_*.csv files."
        )
    )
    parser.add_argument(
        "--provenance-0",
        type=Path,
        required=True,
        help="Soft-ST provenance CSV for science_soft_st_required = 0",
    )
    parser.add_argument(
        "--provenance-1",
        type=Path,
        required=True,
        help="Soft-ST provenance CSV for science_soft_st_required = 1",
    )
    parser.add_argument(
        "--provenance-2",
        type=Path,
        required=True,
        help="Soft-ST provenance CSV for science_soft_st_required = 2",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output PNG path",
    )
    return parser.parse_args()


def _load_affected_rows(path: Path, label: str) -> tuple[pd.DataFrame, dict[str, object]]:
    df = pd.read_csv(path, low_memory=False)
    df = df[df["sequence_type"] == "science"].copy()
    df["science_soft_tail_minutes"] = pd.to_numeric(
        df.get("science_soft_tail_minutes"), errors="coerce"
    ).fillna(0)
    df["start_dt"] = pd.to_datetime(df["start_utc"], utc=True)

    affected = df[df["science_soft_tail_minutes"] > 0].copy()
    affected["row_key"] = affected.apply(
        lambda row: (
            f"V{int(row['visit_id']):02d} · {row['target']} · "
            f"{pd.Timestamp(row['start_dt']).strftime('%m-%d %H:%M')}"
        ),
        axis=1,
    )

    summary = {
        "setting": label,
        "affected_rows": int(len(affected)),
        "total_tail_minutes": float(affected["science_soft_tail_minutes"].sum()),
        "targets": int(affected["target"].nunique()),
        "visits": int(affected["visit_id"].nunique()),
    }
    return affected, summary


def main() -> None:
    args = _parse_args()

    files = {
        "0": args.provenance_0,
        "1": args.provenance_1,
        "2": args.provenance_2,
    }

    rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []

    for label, path in files.items():
        affected, summary = _load_affected_rows(path, label)
        summaries.append(summary)
        for _, row in affected.iterrows():
            rows.append(
                {
                    "setting": label,
                    "row_key": row["row_key"],
                    "tail_minutes": float(row["science_soft_tail_minutes"]),
                }
            )

    summary_df = pd.DataFrame(summaries).sort_values("setting")
    rows_df = pd.DataFrame(rows)

    if rows_df.empty:
        raise ValueError("No affected science rows found across the provided provenance CSVs.")

    heat = (
        rows_df.pivot_table(
            index="row_key",
            columns="setting",
            values="tail_minutes",
            aggfunc="first",
        )
        .fillna(0.0)
    )
    heat["max_tail"] = heat.max(axis=1)
    heat = heat.sort_values(["max_tail"], ascending=False).drop(columns=["max_tail"])

    fig_h = max(8, 0.24 * len(heat) + 1.5)
    fig = plt.figure(figsize=(16, fig_h))
    gs = fig.add_gridspec(
        nrows=1, ncols=2, width_ratios=[1.0, 1.45]
    )

    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])

    x = range(len(summary_df))
    width = 0.35
    ax0.bar(
        [i - width / 2 for i in x],
        summary_df["affected_rows"],
        width=width,
        color="#2b7a78",
        label="Affected rows",
    )
    ax0.bar(
        [i + width / 2 for i in x],
        summary_df["total_tail_minutes"],
        width=width,
        color="#ffb74d",
        label="Total tail minutes",
    )
    for i, row in summary_df.reset_index(drop=True).iterrows():
        ax0.text(
            i - width / 2,
            row["affected_rows"] + 0.8,
            str(int(row["affected_rows"])),
            ha="center",
            va="bottom",
            fontsize=9,
        )
        ax0.text(
            i + width / 2,
            row["total_tail_minutes"] + 0.8,
            str(int(row["total_tail_minutes"])),
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax0.set_xticks(list(x))
    ax0.set_xticklabels([f"st_required={s}" for s in summary_df["setting"]])
    ax0.set_ylabel("Count / Minutes")
    ax0.set_title("Soft-ST Summary By science_soft_st_required")
    ax0.legend(loc="lower right")
    ax0.grid(axis="y", alpha=0.25)

    im = ax1.imshow(heat[["0", "1", "2"]].to_numpy(), aspect="auto", cmap="YlOrRd", vmin=0)
    ax1.set_xticks([0, 1, 2])
    ax1.set_xticklabels(["0", "1", "2"])
    ax1.set_xlabel("science_soft_st_required")
    ax1.set_yticks(range(len(heat)))
    ax1.set_yticklabels(heat.index.tolist(), fontsize=7)
    ax1.set_title("Tail Minutes Per Affected Science Row")

    for i in range(len(heat)):
        for j, col in enumerate(["0", "1", "2"]):
            value = float(heat.iloc[i][col])
            if value > 0:
                ax1.text(
                    j,
                    i,
                    str(int(value)),
                    ha="center",
                    va="center",
                    fontsize=7,
                    color="black",
                )

    cbar = fig.colorbar(im, ax=ax1, fraction=0.03, pad=0.02)
    cbar.set_label("Tail minutes")

    fig.suptitle(
        "Comparison Of Soft-ST Tail Behavior For science_soft_st_required = 0, 1, 2",
        fontsize=14,
        y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.985])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(args.out)


if __name__ == "__main__":
    main()
