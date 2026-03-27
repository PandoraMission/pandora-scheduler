#!/usr/bin/env python3
"""Render a focused PNG flowchart for Earth day/night keepout logic."""

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
    "decision": "#F9F2C8",
    "logic": "#DCEFD8",
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
        fontsize=11,
        fontweight="bold",
        color="#102030",
    )
    ax.text(
        node.x + 0.16,
        node.y + node.h - 0.44,
        node.body,
        ha="left",
        va="top",
        fontsize=8.8,
        linespacing=1.22,
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
            fontsize=8.1,
            ha="center",
            va="bottom",
            color="#334E68",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.2, "alpha": 0.92},
        )


def render(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(18, 12), dpi=180)
    fig.patch.set_facecolor("#FBFCFE")
    ax.set_facecolor("#FBFCFE")
    ax.set_xlim(0, 18.5)
    ax.set_ylim(0, 12.5)
    ax.axis("off")

    nodes = [
        Node(
            "inputs",
            0.5,
            9.7,
            4.5,
            1.8,
            "1. Inputs",
            "effective_earth_threshold()\n- target_unit\n- nadir_unit\n- sun_unit\n- limb_angle_rad\n- earth_avoidance_deg\n- earth_avoidance_day_deg\n- earth_avoidance_night_deg\n- twilight_margin_deg\n- daynight_mode",
            COLORS["input"],
        ),
        Node(
            "override",
            0.5,
            7.0,
            4.5,
            1.8,
            "2. Day/Night Override Active?",
            "If both day_deg and night_deg are None:\n- return scalar earth_avoidance_deg\n- twilight_margin has no effect\nElse continue to day/night classification.",
            COLORS["decision"],
        ),
        Node(
            "classify",
            6.0,
            9.0,
            5.0,
            2.5,
            "3. Choose Day/Night Classifier",
            "daynight_mode:\n- 'subsatellite' => subsatellite_is_sunlit()\n- 'limb' => earthlimb_is_sunlit()\nElse raise ValueError.\nBoth classifiers apply twilight_margin via threshold = -sin(twilight_margin_deg).",
            COLORS["decision"],
        ),
        Node(
            "subsat",
            6.0,
            5.8,
            5.0,
            2.2,
            "4a. Subsatellite Path",
            "subsatellite_is_sunlit()\n- zenith = -nadir\n- compute dot(zenith, sun)\n- sunlit when dot > -sin(twilight_margin)\nThis asks: is the ground point directly below Pandora sunlit?",
            COLORS["logic"],
        ),
        Node(
            "limb",
            11.9,
            5.8,
            5.0,
            2.8,
            "4b. Limb Path",
            "earthlimb_is_sunlit()\n- zenith = -nadir\n- project target into plane perpendicular to zenith\n- get nearest limb direction\n- with limb_angle_rad, build surface normal:\n  n = cos(la)*zenith + sin(la)*limb_dir\n- sunlit when dot(n, sun) > -sin(twilight_margin)",
            COLORS["logic"],
        ),
        Node(
            "select",
            6.0,
            2.4,
            5.0,
            2.3,
            "5. Select Threshold",
            "Resolve actual day/night thresholds:\n- day = earth_avoidance_day_deg or default\n- night = earth_avoidance_night_deg or default\n- return np.where(sunlit, day, night)\nThis is a per-minute boresight Earth-center threshold.",
            COLORS["logic"],
        ),
        Node(
            "final",
            11.9,
            2.1,
            5.0,
            2.6,
            "6. Final Earth Check",
            "compute_visibility_with_constraints()\n- Earth_Sep = angle(target, Earth center)\n- earth_ok = Earth_Sep > earth_threshold\n- boresight_visible = sun_ok & moon_ok & earth_ok\nImportant: this is Earth-center geometry, not star-tracker Earth-limb geometry.",
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

    def right_mid(k):
        n = node_map[k]
        return (n.x + n.w, n.y + n.h / 2)

    def left_mid(k):
        n = node_map[k]
        return (n.x, n.y + n.h / 2)

    add_arrow(ax, center_bottom("inputs"), center_top("override"))
    add_arrow(ax, right_mid("override"), left_mid("classify"), "overrides set", rad=0.05)
    add_arrow(ax, center_bottom("classify"), center_top("subsat"), "subsatellite", rad=0.05)
    add_arrow(ax, right_mid("classify"), center_top("limb"), "limb", rad=-0.18)
    add_arrow(ax, center_bottom("subsat"), center_top("select"))
    add_arrow(ax, center_bottom("limb"), right_mid("select"), rad=0.10)
    add_arrow(ax, center_bottom("select"), center_top("final"))
    add_arrow(ax, center_bottom("override"), (14.4, 4.7), "no overrides:\nuse default directly", rad=-0.36)

    ax.text(
        0.5,
        12.1,
        "Earth Day/Night Keepout Logic",
        fontsize=20,
        fontweight="bold",
        ha="left",
        va="top",
        color="#12263A",
    )
    ax.text(
        0.5,
        11.67,
        "Focused view of how the visibility engine chooses the Earth boresight threshold minute-by-minute.",
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
    output_path = repo_root / "docs" / "earth_daynight_logic_flowchart.png"
    render(output_path)
    print(output_path)


if __name__ == "__main__":
    main()
