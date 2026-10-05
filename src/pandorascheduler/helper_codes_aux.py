"""Utility helpers shared across Pandora scheduler scripts."""

import functools
import json
import os
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))


def general_parameters(obs_sequence_duration: int = 90, occ_sequence_limit: int = 30) -> Tuple[int, int]:
    """Return the default observation and occultation durations in minutes."""

    return obs_sequence_duration, occ_sequence_limit


def observation_sequence(visit, obs_seq_id: str, target: str, priority: str, start, stop, ra, dec):
    """Attach a fully populated observation sequence to the provided visit element."""

    sequence_element = ET.SubElement(visit, "Observation_Sequence")
    sequence_id = ET.SubElement(sequence_element, "ID")
    sequence_id.text = obs_seq_id

    observational_params, nirda_params, vda_params = params_obs_NIRDA_VDA(
        target, priority, start, stop, ra, dec
    )

    obs_parameters = ET.SubElement(sequence_element, "Observational_Parameters")
    for key, value in observational_params.items():
        obs_param_element = ET.SubElement(obs_parameters, key)
        if key in {"Timing", "Boresight"}:
            for idx in range(2):
                sub_element = ET.SubElement(obs_param_element, value[idx])
                sub_element.text = value[idx + 2]
            continue
        obs_param_element.text = str(value)

    payload_parameters = ET.SubElement(sequence_element, "Payload_Parameters")

    nirda = ET.SubElement(payload_parameters, "NIRDA")
    for key, value in nirda_params.items():
        nirda_element = ET.SubElement(nirda, key)
        nirda_element.text = str(value)

    vda = ET.SubElement(payload_parameters, "VDA")
    for key, value in vda_params.items():
        vda_element = ET.SubElement(vda, key)
        vda_element.text = str(value)

    return sequence_element


def params_obs_NIRDA_VDA(target: str, priority: str, start, stop, ra, dec):
    """Assemble observational and payload parameters for a single sequence."""

    try:
        start_format = datetime.strftime(start, "%Y-%m-%dT%H:%M:%SZ")
        stop_format = datetime.strftime(stop, "%Y-%m-%dT%H:%M:%SZ")
    except (TypeError, ValueError, AttributeError):
        start_format, stop_format = start, stop

    try:
        ra_tmp = f"{float(ra)}"
        dec_tmp = f"{float(dec)}"
    except (TypeError, ValueError):
        ra_tmp, dec_tmp = "-999.0", "-999.0"

    observational_parameters = {
        "Target": target,
        "Priority": f"{priority}",
        "Timing": ["Start", "Stop", start_format, stop_format],
        "Boresight": ["RA", "DEC", ra_tmp, dec_tmp],
    }

    params_nirda = {
        "AverageGroups": "1",
        "ROI_StartX": "0",
        "ROI_StartY": "824",
        "ROI_SizeX": "80",
        "ROI_SizeY": "400",
        "SC_Resets1": "1",
        "SC_Resets2": "1",
        "SC_DropFrames1": "0",
        "SC_DropFrames2": "16",
        "SC_DropFrames3": "0",
        "SC_ReadFrames": "4",
        "TargetID": target,
        "SC_Groups": "2",
        "SC_Integrations": "525",
    }

    params_vda = {
        "StartRoiDetMethod": 0,
        "FramesPerCoadd": 50,
        "NumTotalFramesRequested": 9000,
        "TargetRA": ra_tmp,
        "TargetDEC": dec_tmp,
        "IncludeFieldSolnsInResp": 1,
        "StarRoiDimension": 50,
        "MaxNumStarRois": 0,
        "numPredefinedStarRois": 5,
        "PredefinedStarRoiRa": [60.1, 60.2, 60.3, 60.4, 60.5],
        "PredefinedStarRoiDec": [-30.1, -30.2, -30.3, -30.4, -30.5],
        "TargetID": target,
        "NumExposuresMax": 1,
        "ExposureTime_us": 200000,
    }

    return observational_parameters, params_nirda, params_vda


def remove_short_sequences(array: Sequence[float], sequence_too_short: int):
    """Zero-out runs of ones shorter than ``sequence_too_short`` and report their spans."""

    cleaned = np.asarray(array).copy()
    start_index = None
    positions: List[Tuple[int, int]] = []

    for idx, value in enumerate(cleaned):
        if value == 1 and start_index is None:
            start_index = idx
        elif value == 0 and start_index is not None:
            if idx - start_index < sequence_too_short:
                positions.append((start_index, idx - 1))
            start_index = None

    if start_index is not None and len(cleaned) - start_index < sequence_too_short:
        positions.append((start_index, len(cleaned) - 1))

    for start_idx, stop_idx in positions:
        cleaned[start_idx : stop_idx + 1] = 0.0

    return cleaned, positions


def break_long_sequences(start, end, step) -> List[Tuple]:
    """Split a continuous interval into step-sized segments."""

    ranges: List[Tuple] = []
    current = start
    while current < end:
        next_val = min(current + step, end)
        ranges.append((current, next_val))
        current = next_val
    return ranges


@functools.lru_cache(maxsize=None)
def load_visibility_data(target_name: str, path: str):
    """Load cached visibility data for the provided target."""

    file_path = os.path.join(path, target_name, f"Visibility for {target_name}.csv")
    vis = pd.read_csv(file_path)
    return vis["Time(MJD_UTC)"], vis["Visible"]


def schedule_occultation_targets(
    v_names: Iterable[str],
    starts: Sequence[float],
    stops: Sequence[float],
    path: str,
    o_df: pd.DataFrame,
    o_list: pd.DataFrame,
    try_occ_targets: str,
):
    """Select occultation targets that remain visible for each interval."""

    schedule = pd.DataFrame(
        {"Stop": stops, "Target": np.nan, "Visibility": np.nan},
        index=pd.Index(starts, name="Start"),
    )

    if "Visibility" not in o_df.columns:
        o_df["Visibility"] = np.nan

    for v_name in tqdm(v_names, desc=f"Finding visible occultation target from {try_occ_targets}", leave=False):
        vis_times, visibility = load_visibility_data(v_name, path)

        mask = schedule["Target"].isna()
        for start in schedule.index[mask]:
            stop = schedule.loc[start, "Stop"]
            interval_mask = (vis_times >= start) & (vis_times <= stop)
            if np.all(visibility[interval_mask] == 1):
                schedule.loc[start, "Target"] = v_name
                schedule.loc[start, "Visibility"] = 1
                match = o_list.loc[o_list["Star Name"] == v_name]
                if match.empty:
                    continue
                match_row = match.iloc[0]
                o_df.loc[start, "Target"] = v_name
                o_df.loc[start, "RA"] = match_row["RA"]
                o_df.loc[start, "DEC"] = match_row["DEC"]
                o_df.loc[start, "Visibility"] = 1
            elif pd.isna(schedule.loc[start, "Visibility"]):
                schedule.loc[start, "Visibility"] = 0
                o_df.loc[start, "Visibility"] = 0

        if not schedule["Target"].isna().any():
            return o_df, True

    mask = schedule["Target"].isna()
    schedule.loc[mask, "Target"] = "No target"
    schedule.loc[mask, "Visibility"] = 0
    o_df.loc[o_df["Target"].isna(), "Target"] = "No target"
    o_df.loc[o_df["Visibility"].isna(), "Visibility"] = 0

    return o_df, False


def read_json_files(target_list: pd.DataFrame, filename: str) -> pd.DataFrame:
    """Merge JSON key/value data into the target list, expanding columns as needed."""

    target_list_copy = target_list.copy()
    with open(filename, "r", encoding="utf-8") as file:
        data = json.load(file)

    for key, value in data.items():
        if isinstance(value, (list, np.ndarray)):
            value = str(value)

        if key not in target_list_copy.columns:
            target_list_copy[key] = np.nan

        if value not in target_list_copy[key].values:
            nan_indices = target_list_copy[key].isna()
            if nan_indices.any():
                target_list_copy.loc[nan_indices.idxmax(), key] = value
            else:
                target_list_copy = pd.concat(
                    [target_list_copy, pd.DataFrame({key: [value]})],
                    ignore_index=True,
                )

    old_column_name = "Transit Epoch (BJD_TDB-ZZZZZ)"
    new_column_name = "Transit Epoch (BJD_TDB-2400000.5)"
    if old_column_name in target_list_copy.columns:
        target_list_copy[old_column_name] = (
            target_list_copy[old_column_name] - 2400000.5
        )
        target_list_copy = target_list_copy.rename(
            columns={old_column_name: new_column_name}
        )

    return target_list_copy


def round_to_nearest_second(timestamp: datetime) -> datetime:
    """Return ``timestamp`` rounded to the nearest second."""

    if timestamp.microsecond >= 500_000:
        return timestamp + timedelta(seconds=1) - timedelta(microseconds=timestamp.microsecond)
    return timestamp - timedelta(microseconds=timestamp.microsecond)


def print_element_from_xml(element: ET.Element, level: int = 0) -> None:
    """Recursively print an XML element tree for quick debugging."""

    indent = "  " * level
    text = element.text.strip() if element.text else ""
    print(f"{indent}{element.tag}: {text}")
    for child in element:
        print_element_from_xml(child, level + 1)