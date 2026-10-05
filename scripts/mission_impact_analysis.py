#!/usr/bin/env python3
"""Mission-impact analysis: How do day/night Earth keepout and star tracker
constraints affect Pandora's ability to capture 10 transits per target?

This script generates visibility + transit parquets for all science targets
under 6 constraint configurations and compares the number of schedulable
transits against the 10-transit mission requirement.

Usage:
    poetry run python scripts/mission_impact_analysis.py

Runtime estimate: ~20-40 minutes total (the GMAT ephemeris interpolation is
done once per config, then reused for all targets).
"""

from __future__ import annotations

import gc
import logging
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from pandorascheduler_rework.config import PandoraSchedulerConfig
from pandorascheduler_rework.visibility.catalog import build_visibility_catalog

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    datefmt="%H:%M:%S",
)
LOG = logging.getLogger("mission_impact")

# ── File paths ───────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GMAT_PATH = (
    PROJECT_ROOT
    / "src"
    / "pandorascheduler"
    / "data"
    / "Pandora-600km-withoutdrag-20251018.txt"
)
MANIFEST_PATH = (
    PROJECT_ROOT / "src" / "pandorascheduler" / "data" / "exoplanet_targets.csv"
)
OUTPUT_DIR = PROJECT_ROOT / "output_standalone" / "mission_impact"

# ── Mission window ───────────────────────────────────────────────────────────

WINDOW_START = datetime(2026, 2, 5, 0, 0, 0)
WINDOW_END = datetime(2027, 2, 5, 0, 0, 0)
TRANSIT_COVERAGE_MIN = 0.2  # 20% minimum coverage to count a transit

# ── Scenario definitions ─────────────────────────────────────────────────────
# Each scenario: (name, description, config overrides)

SCENARIOS: list[tuple[str, str, dict]] = [
    (
        "A_baseline_86",
        "Old defaults: Earth=86° uniform, no ST",
        dict(
            earth_avoidance_deg=86.0,
            earth_avoidance_day_deg=None,
            earth_avoidance_night_deg=None,
            st_required=0,
        ),
    ),
    (
        "B_earth_110",
        "New Earth default: 110° uniform, no ST",
        dict(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=None,
            earth_avoidance_night_deg=None,
            st_required=0,
        ),
    ),
    (
        "C_daynight_110_80",
        "Day/night: 110°/80°, no ST",
        dict(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            st_required=0,
        ),
    ),
    (
        "D_daynight_ST_moderate",
        "Day/night 110°/80° + moderate ST (Sun=45°, Moon=20°, Limb=10°, 1 req)",
        dict(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            st_sun_min_deg=45.0,
            st_moon_min_deg=20.0,
            st_earthlimb_min_deg=10.0,
            st_required=1,
            roll_step_deg=5.0,
            min_power_frac=0.5,
        ),
    ),
    (
        "E_daynight_ST_strict",
        "Day/night 110°/80° + strict ST (Sun=60°, Moon=30°, Limb=20°, 1 req)",
        dict(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            st_sun_min_deg=60.0,
            st_moon_min_deg=30.0,
            st_earthlimb_min_deg=20.0,
            st_required=1,
            roll_step_deg=5.0,
            min_power_frac=0.5,
        ),
    ),
    (
        "F_daynight_ST_strict_both",
        "Day/night 110°/80° + strict ST (Sun=60°, Moon=30°, Limb=20°, 2 req)",
        dict(
            earth_avoidance_deg=110.0,
            earth_avoidance_day_deg=110.0,
            earth_avoidance_night_deg=80.0,
            st_sun_min_deg=60.0,
            st_moon_min_deg=30.0,
            st_earthlimb_min_deg=20.0,
            st_required=2,
            roll_step_deg=5.0,
            min_power_frac=0.5,
        ),
    ),
]


# ── Load target info ─────────────────────────────────────────────────────────


def load_target_info(manifest_path: Path) -> pd.DataFrame:
    """Load the exoplanet target manifest."""
    df = pd.read_csv(manifest_path)
    return df


# ── Collect transit results from generated parquets ──────────────────────────


def collect_transit_results(
    output_dir: Path,
    targets: pd.DataFrame,
    coverage_min: float,
) -> dict[str, tuple[int, int, float]]:
    """Read planet transit parquets and count schedulable transits.

    Returns dict: planet_name -> (total_transits, schedulable_transits, mean_coverage)
    """
    results: dict[str, tuple[int, int, float]] = {}
    targets_subdir = output_dir / "data" / "targets"

    for _, row in targets.iterrows():
        planet_name = str(row["Planet Name"])
        star_name = str(row["Star Name"])

        planet_parquet = (
            targets_subdir
            / star_name
            / planet_name
            / f"Visibility for {planet_name}.parquet"
        )

        if not planet_parquet.exists():
            LOG.warning("  Missing parquet for %s: %s", planet_name, planet_parquet)
            results[planet_name] = (0, 0, 0.0)
            continue

        try:
            df = pd.read_parquet(planet_parquet)
            if df.empty or "Transit_Coverage" not in df.columns:
                results[planet_name] = (0, 0, 0.0)
                continue

            total = len(df)
            schedulable = int((df["Transit_Coverage"] >= coverage_min).sum())
            mean_cov = float(df["Transit_Coverage"].mean())
            results[planet_name] = (total, schedulable, mean_cov)
        except Exception as e:
            LOG.warning("  Error reading parquet for %s: %s", planet_name, e)
            results[planet_name] = (0, 0, 0.0)

    return results


# ── Collect star visibility fractions from generated parquets ────────────────


def collect_star_visibility(
    output_dir: Path,
    star_names: list[str],
) -> dict[str, float]:
    """Read star visibility parquets and compute visible fraction."""
    vis_fracs: dict[str, float] = {}
    targets_subdir = output_dir / "data" / "targets"

    for star_name in star_names:
        star_parquet = (
            targets_subdir / star_name / f"Visibility for {star_name}.parquet"
        )
        if not star_parquet.exists():
            vis_fracs[star_name] = 0.0
            continue

        try:
            df = pd.read_parquet(star_parquet, columns=["Visible"])
            vis_fracs[star_name] = float((df["Visible"] > 0).mean())
        except Exception as e:
            LOG.warning("  Error reading star visibility for %s: %s", star_name, e)
            vis_fracs[star_name] = 0.0

    return vis_fracs


# ── Main analysis ────────────────────────────────────────────────────────────


def main() -> int:
    LOG.info("=" * 72)
    LOG.info("  MISSION IMPACT ANALYSIS")
    LOG.info("  Window: %s to %s", WINDOW_START.date(), WINDOW_END.date())
    LOG.info("  Transit coverage threshold: %.0f%%", TRANSIT_COVERAGE_MIN * 100)
    LOG.info("=" * 72)

    if not GMAT_PATH.exists():
        LOG.error("GMAT ephemeris not found: %s", GMAT_PATH)
        return 1
    if not MANIFEST_PATH.exists():
        LOG.error("Target manifest not found: %s", MANIFEST_PATH)
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load targets
    all_targets = load_target_info(MANIFEST_PATH)
    primary_targets = all_targets[all_targets["Primary Target"] == 1].copy()
    LOG.info("Loaded %d total targets (%d primary)", len(all_targets), len(primary_targets))

    # Unique stars
    unique_stars = sorted(all_targets["Star Name"].unique())
    primary_stars = sorted(primary_targets["Star Name"].unique())
    LOG.info("Unique stars: %d (primary: %d)", len(unique_stars), len(primary_stars))

    # Collect all results: scenario -> planet -> (total, schedulable, mean_cov)
    all_results: dict[str, dict[str, tuple[int, int, float]]] = {}
    all_vis_fracs: dict[str, dict[str, float]] = {}

    for scenario_idx, (scenario_name, scenario_desc, overrides) in enumerate(SCENARIOS):
        LOG.info("")
        LOG.info("-" * 72)
        LOG.info(
            "  Scenario %d/%d: %s",
            scenario_idx + 1,
            len(SCENARIOS),
            scenario_desc,
        )
        LOG.info("-" * 72)

        scenario_output = OUTPUT_DIR / scenario_name

        config = PandoraSchedulerConfig(
            window_start=WINDOW_START,
            window_end=WINDOW_END,
            gmat_ephemeris=GMAT_PATH,
            output_dir=scenario_output,
            sun_avoidance_deg=91.0,
            moon_avoidance_deg=25.0,
            force_regenerate=True,
            **overrides,
        )

        # Run the full visibility + transit pipeline
        t1 = time.time()
        build_visibility_catalog(
            config=config,
            target_list=MANIFEST_PATH,
        )
        elapsed = time.time() - t1
        LOG.info("  Visibility catalog built (%.1f s)", elapsed)

        # Collect results from generated parquets
        scenario_results = collect_transit_results(
            scenario_output, primary_targets, TRANSIT_COVERAGE_MIN
        )
        scenario_vis_fracs = collect_star_visibility(scenario_output, primary_stars)

        all_results[scenario_name] = scenario_results
        all_vis_fracs[scenario_name] = scenario_vis_fracs

        # Log per-target summary for this scenario
        for _, row in primary_targets.iterrows():
            pname = row["Planet Name"]
            needed = int(row["Number of Transits to Capture"])
            if pname in scenario_results:
                total, sched, mean_cov = scenario_results[pname]
                surplus = sched - needed
                status = "OK" if surplus >= 0 else f"SHORT by {-surplus}"
                LOG.info(
                    "    %-16s  total=%3d  sched=%3d  needed=%2d  surplus=%+3d  [%s]",
                    pname,
                    total,
                    sched,
                    needed,
                    surplus,
                    status,
                )

        gc.collect()

    # ── Build summary tables ─────────────────────────────────────────────
    LOG.info("")
    LOG.info("=" * 72)
    LOG.info("  BUILDING SUMMARY TABLES")
    LOG.info("=" * 72)

    # Per-target table
    planet_names = sorted(primary_targets["Planet Name"].unique())
    planet_info = {}
    for _, row in primary_targets.iterrows():
        planet_info[row["Planet Name"]] = {
            "Star": row["Star Name"],
            "Needed": int(row["Number of Transits to Capture"]),
            "Period_d": float(row["Period (days)"]),
            "Tdur_hr": float(row["Transit Duration (hrs)"]),
            "RA": float(row["RA"]),
            "DEC": float(row["DEC"]),
        }

    rows = []
    for pname in planet_names:
        info = planet_info[pname]
        row = {
            "Planet": pname,
            "Star": info["Star"],
            "Needed": info["Needed"],
            "Period_d": info["Period_d"],
            "Tdur_hr": info["Tdur_hr"],
        }
        for sname, _, _ in SCENARIOS:
            if pname in all_results.get(sname, {}):
                total, sched, mean_cov = all_results[sname][pname]
                row[f"{sname}_total"] = total
                row[f"{sname}_sched"] = sched
                row[f"{sname}_surplus"] = sched - info["Needed"]
                row[f"{sname}_mean_cov"] = round(mean_cov, 3)
            else:
                row[f"{sname}_total"] = 0
                row[f"{sname}_sched"] = 0
                row[f"{sname}_surplus"] = -info["Needed"]
                row[f"{sname}_mean_cov"] = 0.0
        rows.append(row)

    summary_df = pd.DataFrame(rows)
    summary_csv = OUTPUT_DIR / "mission_impact_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    LOG.info("Saved: %s", summary_csv)

    # Star visibility fraction table
    vis_rows = []
    for star_name in primary_stars:
        vrow = {"Star": star_name}
        for sname, _, _ in SCENARIOS:
            vrow[f"{sname}_vis%"] = round(
                all_vis_fracs.get(sname, {}).get(star_name, 0) * 100, 1
            )
        vis_rows.append(vrow)
    vis_df = pd.DataFrame(vis_rows)
    vis_csv = OUTPUT_DIR / "star_visibility_fractions.csv"
    vis_df.to_csv(vis_csv, index=False)
    LOG.info("Saved: %s", vis_csv)

    # ── Print console report ─────────────────────────────────────────────
    report_lines = []

    def rprint(s=""):
        report_lines.append(s)
        print(s)

    rprint("")
    rprint("=" * 100)
    rprint("  MISSION IMPACT REPORT")
    rprint(f"  Window: {WINDOW_START.date()} to {WINDOW_END.date()}")
    rprint(f"  Transit coverage threshold: {TRANSIT_COVERAGE_MIN * 100:.0f}%")
    rprint("=" * 100)

    # Scenario summary
    rprint("\n  SCENARIOS:")
    for sname, sdesc, _ in SCENARIOS:
        sched_counts = [
            all_results[sname][p][1]
            for p in planet_names
            if p in all_results[sname]
        ]
        needs = [
            planet_info[p]["Needed"]
            for p in planet_names
            if p in all_results[sname]
        ]
        shortfalls = [max(0, n - s) for s, n in zip(sched_counts, needs)]
        targets_met = sum(1 for s, n in zip(sched_counts, needs) if s >= n)
        total_deficit = sum(shortfalls)
        min_surplus = min(s - n for s, n in zip(sched_counts, needs)) if needs else 0
        rprint(f"\n    {sname}: {sdesc}")
        rprint(f"      Targets meeting 10-transit req: {targets_met}/{len(planet_names)}")
        rprint(f"      Total transit deficit: {total_deficit}")
        rprint(f"      Minimum surplus: {min_surplus:+d}")

    # Detailed per-target table (schedulable transits only)
    rprint(
        "\n\n  SCHEDULABLE TRANSITS PER TARGET (threshold >= {:.0f}% coverage):".format(
            TRANSIT_COVERAGE_MIN * 100
        )
    )
    header = f"  {'Planet':<18s} {'Nd':>3s}"
    for sname, _, _ in SCENARIOS:
        short = sname.split("_", 1)[0]  # A, B, C, D, E, F
        header += f"  {short:>5s}(+/-)"
    rprint(header)
    rprint("  " + "-" * (len(header) - 2))

    for pname in planet_names:
        info = planet_info[pname]
        line = f"  {pname:<18s} {info['Needed']:>3d}"
        for sname, _, _ in SCENARIOS:
            if pname in all_results[sname]:
                _, sched, _ = all_results[sname][pname]
                surplus = sched - info["Needed"]
                marker = " " if surplus >= 0 else "*"
                line += f"  {sched:>3d}{marker:1s}({surplus:+3d})"
            else:
                line += f"  {'N/A':>9s}"
        rprint(line)

    # Targets at risk
    rprint("\n\n  TARGETS AT RISK (negative surplus in any scenario):")
    at_risk = set()
    for sname, _, _ in SCENARIOS:
        for pname in planet_names:
            if pname in all_results[sname]:
                _, sched, _ = all_results[sname][pname]
                if sched < planet_info[pname]["Needed"]:
                    at_risk.add(pname)

    if at_risk:
        for pname in sorted(at_risk):
            info = planet_info[pname]
            rprint(
                f"\n    {pname} (Need {info['Needed']}, "
                f"Period={info['Period_d']:.2f}d, Star={info['Star']}):"
            )
            for sname, sdesc, _ in SCENARIOS:
                if pname in all_results[sname]:
                    total, sched, mean_cov = all_results[sname][pname]
                    surplus = sched - info["Needed"]
                    rprint(
                        f"      {sname}: {sched}/{total} schedulable "
                        f"(surplus={surplus:+d}, mean_cov={mean_cov:.2f})"
                    )
    else:
        rprint("    None — all targets meet requirements in all scenarios!")

    # Star visibility comparison
    rprint("\n\n  STAR VISIBILITY FRACTIONS (% of year visible):")
    header = f"  {'Star':<18s}"
    for sname, _, _ in SCENARIOS:
        short = sname.split("_", 1)[0]
        header += f"  {short:>7s}"
    rprint(header)
    rprint("  " + "-" * (len(header) - 2))
    for star_name in primary_stars:
        line = f"  {star_name:<18s}"
        for sname, _, _ in SCENARIOS:
            vfrac = all_vis_fracs.get(sname, {}).get(star_name, 0) * 100
            line += f"  {vfrac:>6.1f}%"
        rprint(line)

    # Key findings: compare each scenario to baseline
    rprint("\n\n  IMPACT vs BASELINE (A_baseline_86):")
    baseline = "A_baseline_86"
    for sname, sdesc, _ in SCENARIOS[1:]:
        improvements = 0
        degradations = 0
        same = 0
        total_sched_base = 0
        total_sched_new = 0
        for pname in planet_names:
            if pname in all_results[baseline] and pname in all_results[sname]:
                b_sched = all_results[baseline][pname][1]
                s_sched = all_results[sname][pname][1]
                total_sched_base += b_sched
                total_sched_new += s_sched
                if s_sched > b_sched:
                    improvements += 1
                elif s_sched < b_sched:
                    degradations += 1
                else:
                    same += 1
        delta = total_sched_new - total_sched_base
        rprint(f"\n    {sname}: {sdesc}")
        rprint(
            f"      Improved: {improvements}  Same: {same}  Degraded: {degradations}"
        )
        rprint(
            f"      Total schedulable transits: {total_sched_new} "
            f"(baseline: {total_sched_base}, delta: {delta:+d})"
        )

    # Day/night benefit: compare C vs B
    if "B_earth_110" in all_results and "C_daynight_110_80" in all_results:
        rprint("\n\n  DAY/NIGHT BENEFIT (C vs B — same day angle, relaxed night):")
        for pname in planet_names:
            if pname in all_results["B_earth_110"] and pname in all_results["C_daynight_110_80"]:
                b_sched = all_results["B_earth_110"][pname][1]
                c_sched = all_results["C_daynight_110_80"][pname][1]
                if c_sched != b_sched:
                    delta = c_sched - b_sched
                    rprint(f"    {pname:<18s}: B={b_sched}, C={c_sched} (delta={delta:+d})")

    rprint("\n\n" + "=" * 100)
    rprint("  ANALYSIS COMPLETE")
    rprint(f"  Results saved to: {OUTPUT_DIR}")
    rprint("=" * 100 + "\n")

    # Save report to file
    report_path = OUTPUT_DIR / "mission_impact_report.txt"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    LOG.info("Report saved to: %s", report_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
