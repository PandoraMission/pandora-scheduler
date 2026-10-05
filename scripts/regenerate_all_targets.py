#!/usr/bin/env python3
"""Regenerate all target lists with priority-based sorting."""

from pathlib import Path

from pandorascheduler_rework.targets.manifest import build_target_manifest

target_def_base = Path("../PandoraTargetList/target_definition_files")
output_dir = Path("output_standalone/data")
output_dir.mkdir(parents=True, exist_ok=True)

categories = [
    "exoplanet",
    "auxiliary-standard",
    "monitoring-standard",
    "occultation-standard",
]

for category in categories:
    print(f"Generating {category}_targets.csv...")
    df = build_target_manifest(category, target_def_base)
    output_path = output_dir / f"{category}_targets.csv"
    df.to_csv(output_path, index=False)
    print(f"  ✅ Generated {len(df)} targets")

print("\n✅ All target lists regenerated with priority-based sorting")
