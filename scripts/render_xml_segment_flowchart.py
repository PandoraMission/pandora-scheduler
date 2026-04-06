#!/usr/bin/env python3
"""Render a focused PNG flowchart for XML segment assignment logic."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


@dataclass(frozen=True)
class Node:
    key: str
    x: float
    y: float
    w: float
    h: float
    title: str
    body: str
    color: str


COLORS = {
    "input": "#DCEAF7",
    "science": "#DCEFD8",
    "decision": "#F9F2C8",
    "occ": "#F5D9D4",
    "output": "#E6E0F5",
}


def add_node(ax, node: Node) -> None:
    patch = FancyBboxPatch(
        (node.x, node.y),
        node.w,
        node.h,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.4,
        edgecolor="#243447",
        facecolor=node.color,
    )
    ax.add_patch(patch)
    ax.text(
        node.x + 0.16,
        node.y + node.h - 0.16,
        node.title,
        ha="left",
        va="top",
        fontsize=10.8,
        fontweight="bold",
        color="#102030",
    )
    ax.text(
        node.x + 0.16,
        node.y + node.h - 0.44,
        node.body,
        ha="left",
        va="top",
        fontsize=8.7,
        linespacing=1.2,
        color="#102030",
    )


def add_arrow(ax, start, end, label=None, rad=0.0):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.4,
        color="#334E68",
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arrow)
    if label:
        ax.text(
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2 + 0.1,
            label,
            fontsize=8.0,
            ha="center",
            va="bottom",
            color="#334E68",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.2, "alpha": 0.92},
        )


def render(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(22, 15), dpi=180)
    fig.patch.set_facecolor("#FBFCFE")
    ax.set_facecolor("#FBFCFE")
    ax.set_xlim(0, 22.5)
    ax.set_ylim(0, 15.5)
    ax.axis("off")

    nodes = [
        Node(
            "visit",
            0.5,
            12.8,
            5.0,
            2.0,
            "1. Start Visit",
            "_add_visit()\n- read schedule row\n- load target visibility parquet\n- extract visit_times + raw Visible flags\n- set final_time = CSV stop boundary",
            COLORS["input"],
        ),
        Node(
            "raw",
            0.5,
            9.9,
            5.0,
            2.1,
            "2. Raw Segments",
            "_raw_visit_segments()\n- compute visibility-change indices\n- build alternating (start, stop, is_visible) segments\n- no short-fragment cleanup here",
            COLORS["input"],
        ),
        Node(
            "xml_toggle",
            0.5,
            7.1,
            5.0,
            1.8,
            "3. Occultation XML Enabled?",
            "If include_occultation_sequences_in_xml / generate_occultation_xml / enable_occultation_xml is false:\n- only emit science segments\n- raw occultation intervals are ignored for occultation filling",
            COLORS["decision"],
        ),
        Node(
            "science_policy",
            6.4,
            11.8,
            5.2,
            2.8,
            "4. Science Fragment Policy",
            "_apply_science_fragment_policy()\nFor each raw segment:\n- if science-visible and duration >= min_science_sequence_minutes: keep as science\n- if shorter and contiguous preceding science exists: merge into preceding science chunk\n- else convert that interval into an occultation-fill candidate",
            COLORS["science"],
        ),
        Node(
            "coalesce",
            6.4,
            8.1,
            5.2,
            2.2,
            "5. Coalesce Adjacent Segments",
            "_coalesce_segments()\n- merge adjacent same-kind segments after policy rewrite\n- final visit_segments now define science vs occultation intervals for XML building",
            COLORS["science"],
        ),
        Node(
            "oc_windows",
            6.4,
            5.1,
            5.2,
            2.0,
            "6. Occultation Interval Extraction",
            "_occultation_windows_from_segments()\n- collect occultation segment starts/stops\n- these intervals are the occultation-fill request for this visit",
            COLORS["occ"],
        ),
        Node(
            "occ_source",
            12.4,
            11.5,
            5.0,
            3.0,
            "7. Resolve Occultation Source",
            "_find_occultation_target()\n- break occultation intervals by occ_sequence_limit_min if enabled\n- build scheduled occ_df from occultation catalog(s)\n- obey requested_occ_time_override\n- use pass1 / pass2 search logic\nIf no occ_df is available, fall back to per-segment catalog selection later",
            COLORS["occ"],
        ),
        Node(
            "science_emit",
            12.4,
            7.7,
            5.0,
            2.6,
            "8. Emit Science Segments",
            "_emit_science_sequences()\n- visible segments only\n- enforce min_science_sequence_minutes\n- chunk by obs_sequence_duration_min\n- absorb short science tail via _next_chunk_end()",
            COLORS["science"],
        ),
        Node(
            "occ_emit",
            12.4,
            3.6,
            5.0,
            3.2,
            "9. Emit Occultation Segments",
            "For each occultation segment:\n- use scheduled occ_df row if available, else fallback catalog target\n- visibility gate via _occ_visibility_score()\n- short occultation tail:\n  absorb only if same target stays visible across combined interval\n  else keep separate so another target can be chosen\n- chunk by occ_sequence_limit_min",
            COLORS["occ"],
        ),
        Node(
            "final",
            18.3,
            7.0,
            3.6,
            3.5,
            "10. XML Output",
            "Observation_Sequence entries written in visit order\n- science chunks keep target priority\n- occultation chunks get priority 0\n- if no visible occultation target exists for an occultation interval, a warning is logged and the interval remains uncovered",
            COLORS["output"],
        ),
    ]

    node_map = {n.key: n for n in nodes}
    for node in nodes:
        add_node(ax, node)

    def center_bottom(k):
        n = node_map[k]
        return (n.x + n.w / 2, n.y)

    def center_top(k):
        n = node_map[k]
        return (n.x + n.w / 2, n.y + n.h)

    def left_mid(k):
        n = node_map[k]
        return (n.x, n.y + n.h / 2)

    def right_mid(k):
        n = node_map[k]
        return (n.x + n.w, n.y + n.h / 2)

    add_arrow(ax, center_bottom("visit"), center_top("raw"))
    add_arrow(ax, center_bottom("raw"), center_top("xml_toggle"))
    add_arrow(ax, right_mid("xml_toggle"), left_mid("science_policy"), "XML occultation fill enabled", rad=0.06)
    add_arrow(ax, center_bottom("science_policy"), center_top("coalesce"))
    add_arrow(ax, center_bottom("coalesce"), center_top("oc_windows"))
    add_arrow(ax, right_mid("oc_windows"), left_mid("occ_source"), "occultation intervals", rad=0.03)
    add_arrow(ax, center_bottom("occ_source"), center_top("science_emit"), "visit_segments drive visible chunks", rad=-0.10)
    add_arrow(ax, center_bottom("science_emit"), center_top("occ_emit"))
    add_arrow(ax, right_mid("science_emit"), left_mid("final"), "science sequences", rad=0.03)
    add_arrow(ax, right_mid("occ_emit"), left_mid("final"), "occultation sequences", rad=-0.03)
    add_arrow(ax, center_bottom("xml_toggle"), (20.1, 10.5), "science-only path", rad=-0.40)

    ax.text(
        0.5,
        15.1,
        "XML Segment Assignment Logic",
        fontsize=20,
        fontweight="bold",
        ha="left",
        va="top",
        color="#12263A",
    )
    ax.text(
        0.5,
        14.68,
        "Focused view of how a scheduled visit is broken into science and occultation XML sequences.",
        fontsize=10.8,
        ha="left",
        va="top",
        color="#486581",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    output_path = repo_root / "docs" / "xml_segment_assignment_flowchart.png"
    render(output_path)
    print(output_path)


if __name__ == "__main__":
    main()
