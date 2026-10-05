#!/usr/bin/env python3
"""Build scheduler target CSV manifests from JSON definition sets."""

from __future__ import annotations

import argparse
from pathlib import Path

from pandorascheduler_rework.targets.manifest import (
    TargetDefinitionError,
    build_target_manifest,
)

DEFAULT_CATEGORIES = (
    "exoplanet",
    "auxiliary-standard",
    "monitoring-standard",
    "occultation-standard",
)


def _resolve_base_dir(explicit: str | None) -> Path:
    if explicit:
        resolved = Path(explicit).expanduser().resolve()
        if resolved.is_dir():
            return resolved
        raise SystemExit(f"Provided --base-dir does not exist: {explicit}")

    raise SystemExit(
        "Unable to locate target_definition_files. Provide --base-dir."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate *_targets.csv manifests")
    parser.add_argument(
        "--base-dir",
        help="Root directory containing target definition folders and readout scheme JSON",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parents[1] / "src" / "pandorascheduler" / "data"),
        help="Directory to write the generated CSV files (default: scheduler data dir)",
    )
    parser.add_argument(
        "--categories",
        nargs="+",
        help="Specific categories to export (default: standard scheduler set)",
    )
    args = parser.parse_args()

    base_dir = _resolve_base_dir(args.base_dir)
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    categories = tuple(args.categories) if args.categories else DEFAULT_CATEGORIES

    for category in categories:
        try:
            manifest = build_target_manifest(category, base_dir)
        except TargetDefinitionError as exc:
            raise SystemExit(str(exc)) from exc

        destination = output_dir / f"{category}_targets.csv"
        manifest.to_csv(destination, index=False)
        print(f"Wrote {destination}")


if __name__ == "__main__":
    main()
