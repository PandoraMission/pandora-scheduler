#!/usr/bin/env python3
"""Render a detailed PNG flowchart for the Pandora scheduler pipeline."""

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


PALETTE = {
    "setup": "#DCEAF7",
    "visibility": "#DCEFD8",
    "scheduler": "#F7E2C6",
    "xml": "#F5D9D4",
    "output": "#E6E0F5",
    "decision": "#F9F2C8",
}


def add_node(ax, node: Node) -> None:
    patch = FancyBboxPatch(
        (node.x, node.y),
        node.w,
        node.h,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.6,
        edgecolor="#243447",
        facecolor=node.color,
    )
    ax.add_patch(patch)
    ax.text(
        node.x + 0.2,
        node.y + node.h - 0.22,
        node.title,
        ha="left",
        va="top",
        fontsize=11.5,
        fontweight="bold",
        color="#102030",
    )
    ax.text(
        node.x + 0.2,
        node.y + node.h - 0.55,
        node.body,
        ha="left",
        va="top",
        fontsize=9.3,
        color="#102030",
        linespacing=1.25,
    )


def add_arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    label: str | None = None,
    rad: float = 0.0,
) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.5,
        color="#334E68",
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arrow)
    if label:
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(
            mid_x,
            mid_y + 0.12,
            label,
            fontsize=8.4,
            ha="center",
            va="bottom",
            color="#334E68",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.2, "alpha": 0.9},
        )


def build_nodes() -> list[Node]:
    return [
        Node(
            "cli",
            0.5,
            14.5,
            4.2,
            1.8,
            "1. Runner / Config",
            "run_scheduler.py\n- parse CLI + JSON config\n- resolve output/data_subdir\n- build PandoraSchedulerConfig\n- optional skip_xml / visualizer",
            PALETTE["setup"],
        ),
        Node(
            "pipeline",
            0.5,
            12.1,
            4.2,
            1.9,
            "2. Pipeline Setup",
            "pipeline.build_schedule()\n- create SchedulerPaths\n- derive manifest CSV paths\n- generate target manifests when requested\n- validate primary visit windows",
            PALETTE["setup"],
        ),
        Node(
            "vis_decision",
            0.5,
            9.8,
            4.2,
            1.55,
            "3. Generate Visibility?",
            "Decision gate\n- requires GMAT ephemeris\n- honors force_regenerate\n- writes under output/data_*",
            PALETTE["decision"],
        ),
        Node(
            "visibility",
            5.5,
            12.7,
            4.9,
            4.4,
            "4. Visibility Catalog",
            "visibility.catalog.build_visibility_catalog()\n- load target manifest(s)\n- build minute cadence\n- interpolate GMAT ephemeris\n- compute Sun/Moon/Earth/SAA geometry\n- generate star visibility parquet\n- generate planet transit parquet\n- apply transit overlap fields",
            PALETTE["visibility"],
        ),
        Node(
            "scheduler_loop",
            11.2,
            12.4,
            5.3,
            4.8,
            "5. Scheduler Loop",
            "scheduler.run_scheduler()\n- initialize tracker dataframe\n- preload transit windows\n- iterate rolling windows (schedule_step)\n- score transits by coverage/SAA/schedule weights\n- apply ToO / urgency logic\n- choose primary target or non-primary fill\n- update tracker + observation stats",
            PALETTE["scheduler"],
        ),
        Node(
            "scheduler_checks",
            11.2,
            7.4,
            5.3,
            4.3,
            "6. In-Window Checks",
            "Per-candidate checks\n- filter visibility by time window\n- min_visibility / transit_coverage_min\n- star-tracker keepout policy\n- requested-hours accounting\n- standard-star cadence\n- occultation fallback availability",
            PALETTE["scheduler"],
        ),
        Node(
            "scheduler_outputs",
            11.2,
            3.2,
            5.3,
            3.3,
            "7. Scheduler Outputs",
            "Write artifacts\n- Pandora_Schedule_*.csv\n- observation-time report\n- tracker CSV / pickle\n- diagnostics returned in SchedulerResult",
            PALETTE["output"],
        ),
        Node(
            "xml_decision",
            17.4,
            12.1,
            4.6,
            1.8,
            "8. Build XML?",
            "Decision gate\n- skip_xml=false required\n- generate_occultation_xml controls whether occultation intervals are filled, not whether XML exists",
            PALETTE["decision"],
        ),
        Node(
            "xml_builder",
            17.2,
            8.0,
            5.0,
            3.5,
            "9. Science Calendar Builder",
            "science_calendar.generate_science_calendar()\n- read schedule CSV + manifests\n- read visibility parquet for each visit\n- segment visit into science/occultation intervals\n- chunk science sequences\n- chunk occultation sequences\n- emit Pandora_science_calendar.xml",
            PALETTE["xml"],
        ),
        Node(
            "segment_logic",
            17.2,
            3.3,
            5.0,
            4.0,
            "10. Segment Assignment Logic",
            "Current flow\n- build raw visibility-change segments\n- short science fragment:\n  merge into contiguous preceding science chunk if possible\n  else convert to occultation-fill interval\n- occultation interval:\n  assign scheduled occultation target or fallback catalog target\n- short occultation tail:\n  absorb only if same target stays visible\n  else leave separate for alternate target selection",
            PALETTE["xml"],
        ),
        Node(
            "final_outputs",
            23.0,
            8.2,
            4.4,
            3.1,
            "11. Final Outputs",
            "Run output directory\n- XML at output_*/Pandora_science_calendar.xml\n- schedule_gaps / analyzer scripts\n- optional visualizer PNGs\n- reports + parquet remain under output/data_*",
            PALETTE["output"],
        ),
        Node(
            "visualizer",
            23.0,
            3.7,
            4.4,
            3.1,
            "12. Visualizations",
            "scripts/visualizer.py\n- priority / timeline / simple / target-time\n- visibility mode overlays XML with parquet visibility\n- optional auto-run after pipeline via JSON toggle",
            PALETTE["output"],
        ),
    ]


def render(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(24, 16), dpi=180)
    fig.patch.set_facecolor("#FBFCFE")
    ax.set_facecolor("#FBFCFE")
    ax.set_xlim(0, 28)
    ax.set_ylim(0, 17.5)
    ax.axis("off")

    nodes = build_nodes()
    for node in nodes:
        add_node(ax, node)

    node_map = {node.key: node for node in nodes}

    def center_bottom(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x + n.w / 2, n.y)

    def center_top(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x + n.w / 2, n.y + n.h)

    def right_mid(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x + n.w, n.y + n.h / 2)

    def left_mid(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x, n.y + n.h / 2)

    add_arrow(ax, center_bottom("cli"), center_top("pipeline"))
    add_arrow(ax, center_bottom("pipeline"), center_top("vis_decision"))
    add_arrow(ax, right_mid("vis_decision"), left_mid("visibility"), "yes", rad=0.05)
    add_arrow(ax, center_bottom("visibility"), center_top("scheduler_loop"), "catalog ready", rad=-0.05)
    add_arrow(ax, center_bottom("vis_decision"), (13.85, 12.4), "skip / reuse existing", rad=-0.25)
    add_arrow(ax, center_bottom("scheduler_loop"), center_top("scheduler_checks"))
    add_arrow(ax, center_bottom("scheduler_checks"), center_top("scheduler_outputs"))
    add_arrow(ax, right_mid("scheduler_outputs"), left_mid("xml_decision"), "schedule CSV", rad=0.03)
    add_arrow(ax, center_bottom("xml_decision"), center_top("xml_builder"), "skip_xml = false")
    add_arrow(ax, center_bottom("xml_builder"), center_top("segment_logic"))
    add_arrow(ax, right_mid("xml_builder"), left_mid("final_outputs"), "write XML", rad=0.02)
    add_arrow(ax, right_mid("segment_logic"), left_mid("visualizer"), "XML + data_dir", rad=-0.02)
    add_arrow(ax, center_bottom("final_outputs"), center_top("visualizer"), "optional post-run visualizer", rad=-0.25)

    ax.text(
        0.5,
        17.05,
        "Pandora Scheduler Pipeline Flowchart",
        fontsize=22,
        fontweight="bold",
        ha="left",
        va="top",
        color="#12263A",
    )
    ax.text(
        0.5,
        16.55,
        "Covers runner setup, visibility generation, scheduler loop, XML segment assignment, and visualization outputs.",
        fontsize=11.2,
        ha="left",
        va="top",
        color="#486581",
    )

    legend_items = [
        ("Setup / config", PALETTE["setup"]),
        ("Visibility generation", PALETTE["visibility"]),
        ("Scheduling / selection", PALETTE["scheduler"]),
        ("XML / segment logic", PALETTE["xml"]),
        ("Decision / output", PALETTE["decision"]),
    ]
    legend_x = 17.4
    legend_y = 16.7
    for idx, (label, color) in enumerate(legend_items):
        y = legend_y - idx * 0.34
        ax.add_patch(
            FancyBboxPatch(
                (legend_x, y - 0.14),
                0.35,
                0.2,
                boxstyle="round,pad=0.02,rounding_size=0.03",
                facecolor=color,
                edgecolor="#243447",
                linewidth=1.0,
            )
        )
        ax.text(
            legend_x + 0.48,
            y - 0.04,
            label,
            fontsize=9.0,
            ha="left",
            va="center",
            color="#243447",
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    output_path = repo_root / "docs" / "pipeline_flowchart.png"
    render(output_path)
    print(output_path)


if __name__ == "__main__":
    main()
