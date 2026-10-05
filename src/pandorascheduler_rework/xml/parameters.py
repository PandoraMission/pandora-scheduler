"""XML parameter population logic for NIRDA and VDA."""

import ast
import logging
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from pandorascheduler_rework.utils.string_ops import target_identifier

_PLACEHOLDER_MARKERS = {"SET_BY_TARGET_DEFINITION_FILE", "SET_BY_SCHEDULER"}

LOGGER = logging.getLogger(__name__)


def populate_nirda_parameters(
    payload_parameters: ET.Element, targ_info: pd.DataFrame, diff_in_seconds: float
) -> None:
    """Populate NIRDA parameters in the XML payload."""
    if targ_info.empty:
        LOGGER.warning(
            "Missing target payload info (NIRDA): targ_info is empty; writing empty AcquireInfCamImages"
        )
        ET.SubElement(payload_parameters, "AcquireInfCamImages")
        return

    row = targ_info.iloc[0]

    nirda_columns = targ_info.columns[targ_info.columns.str.startswith("NIRDA_")]
    nirda_element = ET.SubElement(payload_parameters, "AcquireInfCamImages")

    if nirda_columns.empty:
        LOGGER.warning(
            "No NIRDA_* columns present for target %s; skipping NIRDA parameter population",
            target_identifier(row),
        )
        return

    columns_to_ignore = {
        "IncludeFieldSolnsInResp",
        "NIRDA_TargetID",
        "NIRDA_SC_Integrations",
        "NIRDA_FramesPerIntegration",
        "NIRDA_IntegrationTime_s",
    }

    for nirda_key, nirda_value in row[nirda_columns].items():
        column_name = str(nirda_key)
        if pd.isna(nirda_value):
            continue

        if column_name not in columns_to_ignore:
            ET.SubElement(nirda_element, column_name.replace("NIRDA_", "")).text = str(
                nirda_value
            )
            continue

        if column_name == "NIRDA_TargetID":
            ET.SubElement(nirda_element, "TargetID").text = target_identifier(row)
            continue

        if column_name == "NIRDA_SC_Integrations":
            integration_time = row.get("NIRDA_IntegrationTime_s")
            if pd.notna(integration_time) and integration_time:
                integrations = int(np.round(diff_in_seconds / integration_time))
                integrations = max(integrations, 0)
                ET.SubElement(nirda_element, "SC_Integrations").text = str(integrations)
            continue


def populate_vda_parameters(
    payload_parameters: ET.Element, targ_info: pd.DataFrame, diff_in_seconds: float
) -> None:
    """Populate VDA parameters in the XML payload."""
    if targ_info.empty:
        LOGGER.warning(
            "Missing target payload info (VDA): targ_info is empty; writing empty AcquireVisCamScienceData"
        )
        ET.SubElement(payload_parameters, "AcquireVisCamScienceData")
        return

    row = targ_info.iloc[0]

    vda_columns = targ_info.columns[targ_info.columns.str.startswith("VDA_")]
    vda_element = ET.SubElement(payload_parameters, "AcquireVisCamScienceData")

    if vda_columns.empty:
        LOGGER.warning(
            "No VDA_* columns present for target %s; skipping VDA parameter population",
            target_identifier(row),
        )
        return
    columns_to_ignore = {
        "VDA_NumExposuresMax",
        "VDA_NumTotalFramesRequested",
        "VDA_TargetID",
        "VDA_TargetRA",
        "VDA_TargetDEC",
        "VDA_StarRoiDetMethod",
        "VDA_numPredefinedStarRois",
        "VDA_PredefinedStarRoiRa",
        "VDA_PredefinedStarRoiDec",
        "VDA_IntegrationTime_s",
        "VDA_MaxNumStarRois",
    }

    warned: set[str] = set()

    for vda_key, vda_value in row[vda_columns].items():
        column_name = str(vda_key)
        if pd.isna(vda_value):
            continue

        shortened_key = column_name.replace("VDA_", "")

        if column_name not in columns_to_ignore:
            ET.SubElement(vda_element, shortened_key).text = str(vda_value)
            continue

        if column_name == "VDA_TargetID":
            ET.SubElement(vda_element, "TargetID").text = target_identifier(row)
            continue

        if column_name == "VDA_TargetRA":
            ra = row.get("RA")
            if pd.isna(ra) and column_name not in warned:
                LOGGER.warning(
                    "Missing RA for target %s; cannot populate VDA TargetRA reliably",
                    target_identifier(row),
                )
                warned.add(column_name)
            ET.SubElement(vda_element, "TargetRA").text = str(row.get("RA", vda_value))
            continue

        if column_name == "VDA_TargetDEC":
            dec = row.get("DEC")
            if pd.isna(dec) and column_name not in warned:
                LOGGER.warning(
                    "Missing DEC for target %s; cannot populate VDA TargetDEC reliably",
                    target_identifier(row),
                )
                warned.add(column_name)
            ET.SubElement(vda_element, "TargetDEC").text = str(
                row.get("DEC", vda_value)
            )
            continue

        if column_name == "VDA_StarRoiDetMethod":
            value = row.at[column_name]
            fallback = (
                row.get("StarRoiDetMethod")
                if isinstance(value, str) and value in _PLACEHOLDER_MARKERS
                else value
            )

            if fallback is None:
                raise ValueError(
                    "Missing required StarRoiDetMethod for target "
                    f"{target_identifier(row)}"
                )
            if isinstance(fallback, str) and fallback in _PLACEHOLDER_MARKERS:
                raise ValueError(
                    "Unresolved placeholder for StarRoiDetMethod for target "
                    f"{target_identifier(row)}: {fallback!r}"
                )
            if isinstance(fallback, float) and pd.isna(fallback):
                raise ValueError(
                    f"NaN StarRoiDetMethod for target {target_identifier(row)}"
                )

            try:
                fallback_value = int(fallback)
            except (TypeError, ValueError):
                fallback_value = fallback

            ET.SubElement(vda_element, "StarRoiDetMethod").text = str(fallback_value)
            continue

        if column_name == "VDA_MaxNumStarRois":
            method = row.get("StarRoiDetMethod")
            if method == 1:
                # Use numPredefinedStarRois value when method is 1 (matches Legacy)
                value = row.get("numPredefinedStarRois", 0)
            elif method == 2:
                value = 9
            else:
                value = vda_value
            ET.SubElement(vda_element, "MaxNumStarRois").text = str(int(value))
            continue

        if column_name == "VDA_numPredefinedStarRois":
            method = row.get("StarRoiDetMethod")
            if method == 2:
                continue
            field = row.get("numPredefinedStarRois")
            text_value = str(field) if pd.notna(field) else "-9999"
            ET.SubElement(vda_element, "numPredefinedStarRois").text = text_value
            continue

        if column_name in {"VDA_PredefinedStarRoiRa", "VDA_PredefinedStarRoiDec"}:
            method = row.get("StarRoiDetMethod")
            if method == 2:
                continue
            roi_coord_columns = [
                col
                for col in targ_info.columns
                if col.startswith("ROI_coord_") and col != "ROI_coord_epoch"
            ]
            roi_coord_values = targ_info[roi_coord_columns].dropna(axis=1)
            if roi_coord_values.empty:
                if column_name not in warned:
                    LOGGER.warning(
                        "Missing ROI_coord_* columns for target %s; skipping predefined ROI coordinates",
                        target_identifier(row),
                    )
                    warned.add(column_name)
                continue
            try:
                coordinates = np.asarray(
                    [ast.literal_eval(item) for item in roi_coord_values.iloc[0]]
                )
            except (ValueError, SyntaxError):
                if column_name not in warned:
                    LOGGER.warning(
                        "Failed to parse ROI_coord_* values for target %s; skipping predefined ROI coordinates",
                        target_identifier(row),
                    )
                    warned.add(column_name)
                continue

            try:
                target_ra = float(row.get("RA"))
                target_dec = float(row.get("DEC"))
            except (TypeError, ValueError):
                target_ra, target_dec = None, None

            coordinates = np.atleast_2d(coordinates).astype(float, copy=False)
            if (
                coordinates.ndim == 2
                and coordinates.shape[1] >= 2
                and target_ra is not None
                and target_dec is not None
            ):
                coordinates = coordinates.copy()
                coordinates[0, 0] = target_ra
                coordinates[0, 1] = target_dec

            element = ET.SubElement(vda_element, shortened_key)
            for index, coordinate in enumerate(coordinates):
                tag = "RA" if column_name == "VDA_PredefinedStarRoiRa" else "Dec"
                sub = ET.SubElement(element, f"{tag}{index + 1}")
                sub.text = f"{coordinate[0 if tag == 'RA' else 1]:.6f}"
            continue

        if column_name == "VDA_NumTotalFramesRequested":
            exposure_time_us = row.get("VDA_ExposureTime_us")
            frames_per_coadd = row.get("VDA_FramesPerCoadd")
            if (
                pd.notna(exposure_time_us)
                and pd.notna(frames_per_coadd)
                and frames_per_coadd
            ):
                exposure_seconds = 1e-6 * float(exposure_time_us)
                if exposure_seconds > 0:
                    coadd = int(frames_per_coadd)
                    frames = int(
                        np.floor(diff_in_seconds / exposure_seconds / coadd) * coadd
                    )
                    ET.SubElement(vda_element, "NumTotalFramesRequested").text = str(
                        max(frames, 0)
                    )
            continue
