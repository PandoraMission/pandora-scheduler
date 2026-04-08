"""Utilities for building scheduler target manifests from JSON definitions.

The rework pipeline loads target definition JSON, priority tables, and readout
scheme metadata from the external Pandora target list repository. The
``build_target_manifest`` helper assembles these artefacts into the flattened
CSV shape expected by the scheduler while remaining agnostic to where the files
are stored on disk.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, MutableMapping, Tuple

import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from erfa import ErfaWarning

warnings.filterwarnings("ignore", category=ErfaWarning)


_METADATA_FIELDS = {"Time Created", "Version", "Author", "Time Updated"}
_EXOPLANET_CATEGORIES = {
    "exoplanet",
    "auxiliary-exoplanet",
    "primary-exoplanet",
    "secondary-exoplanet",
}
_STANDARD_CATEGORIES = {
    "auxiliary-standard",
    "monitoring-standard",
    "occultation-standard",
}
_OCCULTATION_CATEGORY = "occultation-standard"
_DEFAULT_OBSERVATION_EPOCH = Time("2026-01-05")


class TargetDefinitionError(RuntimeError):
    """Raised when target definition assets are missing or inconsistent."""


_PLACEHOLDER_MARKERS = {"SET_BY_TARGET_DEFINITION_FILE", "SET_BY_SCHEDULER"}


@dataclass(frozen=True)
class _ReadoutSchemes:
    nirda_fixed: Mapping[str, object]
    nirda_schemes: Mapping[str, Mapping[str, object]]
    vda_fixed: Mapping[str, object]
    vda_schemes: Mapping[str, Mapping[str, object]]


def build_target_manifest(
    category: str,
    base_dir: Path,
    *,
    observation_epoch: Time | str = _DEFAULT_OBSERVATION_EPOCH,
) -> pd.DataFrame:
    """Return the flattened manifest for *category*.

    Parameters
    ----------
    category
        Target category directory name, e.g. ``"exoplanet"`` or
        ``"monitoring-standard"``.
    base_dir
        Path containing the category directory alongside the readout scheme
        JSON files.
    observation_epoch
        Epoch at which proper motion should be applied when propagating
        coordinates. Defaults to 2026-01-05 to match the legacy scheduler.
    """

    category = category.strip()
    base_dir = Path(base_dir)
    category_dir = base_dir / category
    if not category_dir.is_dir():
        raise TargetDefinitionError(
            f"Target definition directory not found: {category_dir}"
        )

    if isinstance(observation_epoch, str):
        observation_epoch = Time(observation_epoch)

    readouts = _load_readout_schemes(base_dir)
    priority_table = _load_priority_table(category_dir, category)

    rows: List[MutableMapping[str, object]] = []
    for json_path in sorted(category_dir.glob("*_target_definition.json")):
        row = _load_target_definition(json_path)
        row["Original Filename"] = json_path.name.replace("_target_definition.json", "")

        _validate_required_target_fields(row, category)

        _apply_priority(row, category, priority_table)
        _apply_identity_columns(row, category)
        _apply_readout_settings(row, readouts)
        _apply_proper_motion(row, observation_epoch)

        rows.append(row)

    if not rows:
        raise TargetDefinitionError(
            f"No target definition JSON files found in {category_dir}"
        )

    manifest = pd.DataFrame(rows)
    manifest = _normalise_manifest_columns(manifest, category)
    manifest = _standardise_dtypes(manifest)

    # Sort all target categories by priority (descending) to match Legacy behavior
    if "Priority" in manifest.columns:
        manifest = manifest.sort_values(by="Priority", ascending=False).reset_index(
            drop=True
        )

    return manifest


def _load_target_definition(path: Path) -> MutableMapping[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        payload: Dict[str, object] = json.load(handle)

    # Remove unwanted keys to match Legacy behavior
    for key in ["Time Created", "Version", "Author", "Time Updated"]:
        payload.pop(key, None)

    return _flatten_dict(payload)


def _flatten_dict(
    payload: Mapping[str, object],
    parent_key: str | None = None,
    *,
    sep: str = "_",
) -> MutableMapping[str, object]:
    flat: Dict[str, object] = {}
    for key, value in payload.items():
        new_key = f"{parent_key}{sep}{key}" if parent_key else key
        if isinstance(value, Mapping):
            flat.update(_flatten_dict(value, new_key, sep=sep))
        elif isinstance(value, list):
            for idx, item in enumerate(value):
                child_key = f"{new_key}{sep}{idx}"
                if isinstance(item, Mapping):
                    flat.update(_flatten_dict(item, child_key, sep=sep))
                else:
                    flat[child_key] = item
        else:
            flat[new_key] = value
    return flat


def _load_readout_schemes(base_dir: Path) -> _ReadoutSchemes:
    def load_file(
        filename: str,
    ) -> Tuple[Mapping[str, object], Mapping[str, Mapping[str, object]]]:
        path = base_dir / filename
        if not path.is_file():
            raise TargetDefinitionError(f"Readout scheme file missing: {path}")
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)["data"]
        fixed = data.get("FixedParameters", {})
        schemes = {
            key: value
            for key, value in data.items()
            if key not in {"FixedParameters", "CommandName", "IncludedMnemonics"}
        }
        return fixed, schemes

    nirda_fixed, nirda_schemes = load_file("nirda_readout_schemes.json")
    vda_fixed, vda_schemes = load_file("vda_readout_schemes.json")
    return _ReadoutSchemes(nirda_fixed, nirda_schemes, vda_fixed, vda_schemes)


def _load_priority_table(category_dir: Path, category: str) -> pd.DataFrame | None:
    priority_path = category_dir / f"{category}_priorities.csv"
    if not priority_path.is_file():
        # raise an error here. we need these files
        raise TargetDefinitionError(f"Priority table missing: {priority_path}")

    metadata, table = _read_priority_csv(priority_path)
    if table.empty:
        raise TargetDefinitionError(f"Priority table {priority_path} contained no rows")
    return table


def _read_priority_csv(path: Path) -> Tuple[Dict[str, str], pd.DataFrame]:
    metadata: Dict[str, str] = {}
    data_start = 0
    with path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or not stripped.startswith("#"):
            data_start = idx
            break
        if ":" in stripped:
            key, value = stripped.lstrip("# ").split(":", 1)
            metadata[key.strip()] = value.strip()

    table = pd.read_csv(path, skiprows=data_start)
    return metadata, table


def _apply_priority(
    row: MutableMapping[str, object],
    category: str,
    priority_table: pd.DataFrame | None,
) -> None:
    filename = row.get("Original Filename")

    if category in _EXOPLANET_CATEGORIES:
        if priority_table is None:
            raise TargetDefinitionError(
                f"Missing priority table for category '{category}'"
            )
        match = priority_table.loc[priority_table["target"] == filename]
        if match.empty:
            raise TargetDefinitionError(
                f"Target '{filename}' not present in {category} priority table"
            )
        row["Priority"] = float(match["priority"].iloc[0])
        row["Number of Transits to Capture"] = int(
            round(float(match["transits_req"].iloc[0]))
        )
    elif category in _STANDARD_CATEGORIES:
        if priority_table is None:
            raise TargetDefinitionError(
                f"Missing priority table for category '{category}'"
            )
        match = priority_table.loc[priority_table["target"] == filename]
        if match.empty:
            raise TargetDefinitionError(
                f"Target '{filename}' not present in {category} priority table"
            )
        row["Priority"] = float(match["priority"].iloc[0])
        # Add Number of Hours Requested from priority table (required for standard categories)
        if "hours_req" not in match.columns:
            raise TargetDefinitionError(
                f"Priority table for category '{category}' is missing required 'hours_req' column"
            )
        hours_req_value = match["hours_req"].iloc[0]
        if pd.isna(hours_req_value):
            raise TargetDefinitionError(
                f"Target '{filename}' in {category} has missing 'hours_req' value"
            )
        row["Number of Hours Requested"] = int(round(float(hours_req_value)))
    else:
        raise TargetDefinitionError(f"Unsupported target category '{category}'")


def _apply_identity_columns(row: MutableMapping[str, object], category: str) -> None:
    star_name = str(row.get("Star Name", "") or "")

    if category in _EXOPLANET_CATEGORIES:
        planet_name_raw = str(row.get("Planet Name", "") or "")
        if planet_name_raw:
            planet_name = _format_planet_name(planet_name_raw)
            row["Planet Name"] = planet_name
            row["Planet Simbad Name"] = planet_name
        else:
            row["Planet Simbad Name"] = planet_name_raw
        row["Star Simbad Name"] = star_name
        if "Transit Epoch (BJD_TDB)" in row:
            row["Transit Epoch (BJD_TDB-2400000.5)"] = (
                float(row["Transit Epoch (BJD_TDB)"]) - 2400000.5
            )
    elif category in _STANDARD_CATEGORIES:
        row["Planet Name"] = ""
        row["Planet Simbad Name"] = ""
        row["Star Simbad Name"] = star_name
    else:
        row["Star Simbad Name"] = star_name


def _format_planet_name(name: str) -> str:
    """Return planet name as-is (legacy doesn't modify planet names)."""
    return name


def _apply_readout_settings(
    row: MutableMapping[str, object], readouts: _ReadoutSchemes
) -> None:
    for key, value in readouts.nirda_fixed.items():
        row[f"NIRDA_{key}"] = value

    nirda_setting = row.get("NIRDA Setting")
    if not isinstance(nirda_setting, str) or not nirda_setting.strip():
        raise TargetDefinitionError(
            "Missing required 'NIRDA Setting' for target definition "
            f"'{row.get('Original Filename', '')}'"
        )
    scheme = readouts.nirda_schemes.get(nirda_setting)
    if scheme is None:
        raise TargetDefinitionError(
            "Unknown NIRDA Setting for target definition "
            f"'{row.get('Original Filename', '')}': {nirda_setting!r}"
        )
    for key, value in scheme.items():
        row[f"NIRDA_{key}"] = value

    for key, value in readouts.vda_fixed.items():
        row[f"VDA_{key}"] = value

    vda_setting = row.get("VDA Setting")
    if not isinstance(vda_setting, str) or not vda_setting.strip():
        raise TargetDefinitionError(
            "Missing required 'VDA Setting' for target definition "
            f"'{row.get('Original Filename', '')}'"
        )
    scheme = readouts.vda_schemes.get(vda_setting)
    if scheme is None:
        raise TargetDefinitionError(
            "Unknown VDA Setting for target definition "
            f"'{row.get('Original Filename', '')}': {vda_setting!r}"
        )
    for key, value in scheme.items():
        row[f"VDA_{key}"] = value


def _validate_required_target_fields(row: Mapping[str, object], category: str) -> None:
    """Validate required fields are present in the target definition JSON.

    Policy: we require StarRoiDetMethod and explicit readout scheme names.
    The readout scheme values are validated later (against scheme files) in
    _apply_readout_settings; here we enforce presence and basic shape.
    """
    filename = str(row.get("Original Filename", "") or "")

    for key in ("NIRDA Setting", "VDA Setting"):
        value = row.get(key)
        if not isinstance(value, str) or not value.strip():
            raise TargetDefinitionError(
                f"Target '{filename}' is missing required '{key}'"
            )

    star_roi = row.get("StarRoiDetMethod")
    if star_roi is None or (isinstance(star_roi, float) and pd.isna(star_roi)):
        raise TargetDefinitionError(
            f"Target '{filename}' is missing required 'StarRoiDetMethod'"
        )
    if isinstance(star_roi, str):
        if not star_roi.strip() or star_roi.strip() in _PLACEHOLDER_MARKERS:
            raise TargetDefinitionError(
                f"Target '{filename}' has unresolved 'StarRoiDetMethod': {star_roi!r}"
            )
        try:
            star_roi_int = int(star_roi)
        except ValueError as exc:
            raise TargetDefinitionError(
                f"Target '{filename}' has invalid 'StarRoiDetMethod': {star_roi!r}"
            ) from exc
    else:
        try:
            star_roi_int = int(star_roi)  # type: ignore[arg-type]
        except (TypeError, ValueError) as exc:
            raise TargetDefinitionError(
                f"Target '{filename}' has invalid 'StarRoiDetMethod': {star_roi!r}"
            ) from exc

    if star_roi_int not in {1, 2}:
        raise TargetDefinitionError(
            f"Target '{filename}' has unsupported 'StarRoiDetMethod': {star_roi_int}"
        )


def _apply_proper_motion(
    row: MutableMapping[str, object], observation_epoch: Time
) -> None:
    try:
        ra0 = float(row["RA"])
        dec0 = float(row["DEC"])
        pm_ra = float(row.get("pm_RA", 0.0))
        pm_dec = float(row.get("pm_DEC", 0.0))
    except (KeyError, TypeError, ValueError):
        return

    epoch_start = Time("J2016.0")

    coord = SkyCoord(
        ra=ra0 * u.deg,
        dec=dec0 * u.deg,
        pm_ra_cosdec=pm_ra * u.mas / u.yr,
        pm_dec=pm_dec * u.mas / u.yr,
        frame="icrs",
        obstime=epoch_start,
        distance=None,
    )

    propagated = coord.apply_space_motion(new_obstime=observation_epoch)
    row["RA"] = propagated.ra.deg
    row["DEC"] = propagated.dec.deg


def _normalise_manifest_columns(df: pd.DataFrame, category: str) -> pd.DataFrame:
    columns_order: List[str]

    if category in _STANDARD_CATEGORIES or category == _OCCULTATION_CATEGORY:
        columns_order = [col for col in ("Star Name", "Star Simbad Name") if col in df]
        remaining = [
            col for col in df.columns if col not in columns_order and col != "Priority"
        ]
        columns_order.extend(remaining)
        if "Priority" in df.columns:
            columns_order.append("Priority")
    else:
        base_cols = [
            "Planet Name",
            "Planet Simbad Name",
            "Star Name",
            "Star Simbad Name",
            "Number of Transits to Capture",
            "Priority",
            "Original Filename",
        ]
        columns_order = [col for col in base_cols if col in df.columns]
        extra_optional = "Transit Epoch (BJD_TDB-2400000.5)"
        if extra_optional in df.columns and extra_optional not in columns_order:
            columns_order.append(extra_optional)
        columns_order.extend([col for col in df.columns if col not in columns_order])

    nirda_vda_columns = [
        col
        for col in columns_order
        if col.startswith("NIRDA_") or col.startswith("VDA_")
    ]
    columns_order = [col for col in columns_order if col not in nirda_vda_columns]
    columns_order.extend(nirda_vda_columns)

    return df.reindex(columns=columns_order)


def _standardise_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    numeric_int_columns = [
        "Number of Transits to Capture",
        "Number of Hours Requested",
        "Primary Target",
        "numPredefinedStarRois",
    ]
    for column in numeric_int_columns:
        if column in result.columns:
            result[column] = pd.to_numeric(result[column], errors="coerce")
            if pd.api.types.is_float_dtype(result[column]):
                if result[column].isna().any():
                    result[column] = result[column].astype("Int64")
                else:
                    result[column] = result[column].round().astype(int)

    if "Priority" in result.columns:
        result["Priority"] = pd.to_numeric(result["Priority"], errors="coerce")

    return result


__all__ = ["build_target_manifest", "TargetDefinitionError"]
