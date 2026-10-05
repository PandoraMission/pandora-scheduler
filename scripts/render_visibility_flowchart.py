#!/usr/bin/env python3
"""Render a detailed PNG flowchart for visibility-generation logic."""

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
    "geometry": "#DCEFD8",
    "check": "#F9F2C8",
    "roll": "#F7E2C6",
    "output": "#E6E0F5",
    "planet": "#F5D9D4",
}


def add_node(ax, node: Node) -> None:
    patch = FancyBboxPatch(
        (node.x, node.y),
        node.w,
        node.h,
        boxstyle="round,pad=0.02,rounding_size=0.06",
        linewidth=1.5,
        edgecolor="#243447",
        facecolor=node.color,
    )
    ax.add_patch(patch)
    ax.text(
        node.x + 0.18,
        node.y + node.h - 0.18,
        node.title,
        ha="left",
        va="top",
        fontsize=11.0,
        fontweight="bold",
        color="#102030",
    )
    ax.text(
        node.x + 0.18,
        node.y + node.h - 0.48,
        node.body,
        ha="left",
        va="top",
        fontsize=8.75,
        linespacing=1.22,
        color="#102030",
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
        linewidth=1.45,
        color="#334E68",
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arrow)
    if label:
        ax.text(
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2 + 0.12,
            label,
            fontsize=8.2,
            ha="center",
            va="bottom",
            color="#334E68",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 0.2},
        )


def nodes() -> list[Node]:
    return [
        Node(
            "entry",
            0.4,
            15.0,
            4.8,
            1.9,
            "1. Visibility Entry",
            "build_visibility_catalog()\n- load target manifest(s)\n- filter by target_filters\n- derive output_root under output/data_*\n- decide which star parquet files need generation",
            PALETTE["setup"],
        ),
        Node(
            "cadence",
            0.4,
            12.4,
            4.8,
            2.0,
            "2. Shared Minute Cadence",
            "geometry.build_minute_cadence()\n- generate 1-minute timestamps for window_start..window_end\n- store Time(MJD_UTC) and Time_UTC\n- one shared cadence reused for all stars",
            PALETTE["geometry"],
        ),
        Node(
            "ephem",
            0.4,
            9.6,
            4.8,
            2.3,
            "3. GMAT Interpolation",
            "geometry.interpolate_gmat_ephemeris()\n- detect time column + spacecraft name\n- interpolate Earth / spacecraft / Sun / Moon vectors\n- interpolate spacecraft latitude / longitude\n- build spacecraft-centric vectors",
            PALETTE["geometry"],
        ),
        Node(
            "payload",
            0.4,
            6.3,
            4.8,
            2.8,
            "4. Shared Geometry Payload",
            "catalog._build_base_payload()\n- compute SAA_Crossing from lat/lon box\n- normalize nadir / zenith / Sun / Moon unit vectors\n- compute observer distance and limb_angle_rad\n- detect orbit boundaries / orbit slices\n- cache Time_UTC for parquet output",
            PALETTE["geometry"],
        ),
        Node(
            "target",
            6.1,
            15.0,
            5.2,
            2.0,
            "5. Per-Star Target Vector",
            "catalog._resolve_star_coord()\n- resolve RA/Dec from manifest / metadata fallback\n- build fixed ICRS target unit vector\n- broadcast to every minute in cadence",
            PALETTE["setup"],
        ),
        Node(
            "boresight",
            6.1,
            12.1,
            5.2,
            2.5,
            "6. Boresight Checks",
            "constraints.compute_visibility_with_constraints()\nPhase A:\n- Sun_Sep = sep(target, Sun)\n- Moon_Sep = sep(target, Moon)\n- Earth_Sep = sep(target, Earth center)\n- apply sun_ok, moon_ok, earth_ok",
            PALETTE["check"],
        ),
        Node(
            "earth_mode",
            6.1,
            8.8,
            5.2,
            2.7,
            "7. Day/Night Earth Threshold",
            "constraints.effective_earth_threshold()\n- if no day/night overrides: use earth_avoidance_deg\n- else classify day/night by:\n  subsatellite_is_sunlit() or\n  earthlimb_is_sunlit()\n- choose day vs night Earth-center threshold\n- compare Earth_Sep > threshold",
            PALETTE["check"],
        ),
        Node(
            "st_gate",
            12.1,
            14.7,
            5.4,
            2.2,
            "8. Any ST Checks Active?",
            "constraints._st_thresholds_active()\n- st_required > 0\n- any ST Sun/Moon/Earth-limb threshold enabled?\n- no => Visible = boresight_visible directly",
            PALETTE["check"],
        ),
        Node(
            "roll",
            12.1,
            11.8,
            5.4,
            4.0,
            "9. Per-Orbit Roll Sweep",
            "constraints.find_best_roll_per_orbit()\nFor each orbit slice:\n- sweep roll angles 0..360 by roll_step_deg\n- construct fixed_roll_attitude()\n- reject rolls with mean solar power < min_power_frac\n- rotate ST1 / ST2 boresights to ECI\n- evaluate tracker keepouts at each minute",
            PALETTE["roll"],
        ),
        Node(
            "st_eval",
            12.1,
            7.1,
            5.4,
            4.0,
            "10. Star-Tracker Evaluation",
            "evaluate_star_tracker()\nPer tracker and per minute:\n- Sun separation >= st_sun_min_deg\n- Moon separation >= st_moon_min_deg\n- Earth-limb angle >= ST Earth-limb threshold\nCombine:\n- st_required == 1 => ST1 OR ST2\n- st_required == 2 => ST1 AND ST2\nPick roll with best visible-minute count, then best power",
            PALETTE["roll"],
        ),
        Node(
            "star_output",
            18.3,
            11.0,
            5.6,
            3.6,
            "11. Star Visibility Parquet",
            "catalog._build_star_visibility()\nWrite columns:\n- Time(MJD_UTC), Time_UTC\n- Visible\n- Earth_Sep, Moon_Sep, Sun_Sep\n- Roll_Deg\n- N_ST_Pass\n- SAA_Crossing\nThen _write_visibility_parquet()",
            PALETTE["output"],
        ),
        Node(
            "planet_entry",
            18.3,
            6.4,
            5.6,
            2.8,
            "12. Planet Transit Extraction",
            "catalog._build_planet_transits()\n- requires planet columns in manifest\n- reuse star visibility parquet\n- compute transit windows from period / epoch / duration\n- floor start/stop to minute precision",
            PALETTE["planet"],
        ),
        Node(
            "planet_checks",
            24.6,
            10.8,
            5.3,
            4.0,
            "13. Planet Coverage Checks",
            "catalog._compute_planet_transits()\nFor each transit:\n- build minute set across transit window\n- Transit_Coverage = visible minutes / total transit minutes\n- SAA_Overlap = SAA minutes / total transit minutes\n- write Transit_Start/Stop UTC + MJD-style values",
            PALETTE["planet"],
        ),
        Node(
            "overlap",
            24.6,
            6.0,
            5.3,
            3.6,
            "14. Transit Overlap Pass",
            "catalog._apply_transit_overlaps()\n- group planets by host star\n- if multiple planets share a star, add Transit_Overlap\n- skip recomputation when parquet schema already has overlap column",
            PALETTE["planet"],
        ),
        Node(
            "final",
            18.3,
            1.7,
            11.6,
            3.1,
            "15. Visibility Artifacts on Disk",
            "Run output/data_*/\n- targets/<star>/Visibility for <star>.parquet\n- targets/<star>/<planet>/Visibility for <planet>.parquet\n- aux_targets/... for standards\nThese files are then consumed by scheduler.py, science_calendar.py, and scripts/visualizer.py --mode visibility",
            PALETTE["output"],
        ),
    ]


def render(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(26, 18), dpi=180)
    fig.patch.set_facecolor("#FBFCFE")
    ax.set_facecolor("#FBFCFE")
    ax.set_xlim(0, 30.5)
    ax.set_ylim(0, 17.5)
    ax.axis("off")

    node_list = nodes()
    for node in node_list:
        add_node(ax, node)
    node_map = {node.key: node for node in node_list}

    def center_bottom(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x + n.w / 2, n.y)

    def center_top(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x + n.w / 2, n.y + n.h)

    def left_mid(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x, n.y + n.h / 2)

    def right_mid(key: str) -> tuple[float, float]:
        n = node_map[key]
        return (n.x + n.w, n.y + n.h / 2)

    add_arrow(ax, center_bottom("entry"), center_top("cadence"))
    add_arrow(ax, center_bottom("cadence"), center_top("ephem"))
    add_arrow(ax, center_bottom("ephem"), center_top("payload"))
    add_arrow(ax, right_mid("payload"), left_mid("target"), "shared payload", rad=-0.20)
    add_arrow(ax, center_bottom("target"), center_top("boresight"))
    add_arrow(ax, center_bottom("boresight"), center_top("earth_mode"))
    add_arrow(ax, right_mid("boresight"), left_mid("st_gate"), "boresight_visible", rad=0.05)
    add_arrow(ax, center_bottom("st_gate"), center_top("roll"), "ST active")
    add_arrow(ax, center_bottom("roll"), center_top("st_eval"))
    add_arrow(ax, right_mid("st_gate"), left_mid("star_output"), "no ST checks", rad=-0.08)
    add_arrow(ax, right_mid("st_eval"), left_mid("star_output"), "Visible + Roll_Deg + N_ST_Pass", rad=-0.02)
    add_arrow(ax, center_bottom("star_output"), center_top("planet_entry"), "host-star visibility", rad=0.10)
    add_arrow(ax, right_mid("planet_entry"), left_mid("planet_checks"), "per transit", rad=0.05)
    add_arrow(ax, center_bottom("planet_checks"), center_top("overlap"))
    add_arrow(ax, center_bottom("star_output"), (24.1, 4.8), "star parquet", rad=-0.25)
    add_arrow(ax, center_bottom("overlap"), (24.1, 4.8), "planet parquet + overlap", rad=0.15)

    ax.text(
        0.4,
        17.05,
        "Visibility Checks Flowchart",
        fontsize=22,
        fontweight="bold",
        ha="left",
        va="top",
        color="#12263A",
    )
    ax.text(
        0.4,
        16.58,
        "Detailed path from GMAT ephemeris to star and planet visibility parquet files, including all boresight and star-tracker checks.",
        fontsize=11.2,
        ha="left",
        va="top",
        color="#486581",
    )

    legend = [
        ("Setup / target prep", PALETTE["setup"]),
        ("Shared geometry", PALETTE["geometry"]),
        ("Boresight / threshold checks", PALETTE["check"]),
        ("Roll + star tracker logic", PALETTE["roll"]),
        ("Planet transit products", PALETTE["planet"]),
        ("Outputs", PALETTE["output"]),
    ]
    base_x = 20.2
    base_y = 16.75
    for idx, (label, color) in enumerate(legend):
        y = base_y - idx * 0.33
        ax.add_patch(
            FancyBboxPatch(
                (base_x, y - 0.13),
                0.35,
                0.19,
                boxstyle="round,pad=0.02,rounding_size=0.03",
                facecolor=color,
                edgecolor="#243447",
                linewidth=0.9,
            )
        )
        ax.text(
            base_x + 0.46,
            y - 0.035,
            label,
            fontsize=8.9,
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
    output_path = repo_root / "docs" / "visibility_checks_flowchart.png"
    render(output_path)
    print(output_path)


if __name__ == "__main__":
    main()
