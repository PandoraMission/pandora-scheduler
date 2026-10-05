"""XML builder for observation sequences."""

import logging
import xml.etree.ElementTree as ET
from datetime import datetime

import pandas as pd

from pandorascheduler_rework.xml.parameters import (
    populate_nirda_parameters,
    populate_vda_parameters,
)

LOGGER = logging.getLogger(__name__)


def observation_sequence(
    visit,
    obs_seq_id: str,
    target_name: str,
    priority: str,
    start,
    stop,
    ra,
    dec,
    targ_info: pd.DataFrame,
):
    """Create an Observation_Sequence XML element.

    Args:
        visit: Parent XML element
        obs_seq_id: Sequence ID string
        target_name: Name of the target
        priority: Priority string
        start: Start time (datetime or str)
        stop: Stop time (datetime or str)
        ra: Right Ascension
        dec: Declination
        targ_info: DataFrame containing target parameters (NIRDA/VDA)

    Returns:
        The created sequence Element
    """
    sequence_element = ET.SubElement(visit, "Observation_Sequence")
    ET.SubElement(sequence_element, "ID").text = obs_seq_id

    observational_parameters = _build_observational_parameters(
        target_name, priority, start, stop, ra, dec
    )

    obs_parameters = ET.SubElement(sequence_element, "Observational_Parameters")
    for key, value in observational_parameters.items():
        obs_param_element = ET.SubElement(obs_parameters, key)
        if key in {"Timing", "Boresight"}:
            for index in range(2):
                sub_element = ET.SubElement(obs_param_element, value[index])
                sub_element.text = value[index + 2]
        else:
            obs_param_element.text = str(value)

    diff_in_seconds = _duration_in_seconds(start, stop)

    payload_parameters = ET.SubElement(sequence_element, "Payload_Parameters")
    populate_nirda_parameters(payload_parameters, targ_info, diff_in_seconds)
    populate_vda_parameters(payload_parameters, targ_info, diff_in_seconds)

    return sequence_element


def _build_observational_parameters(target_name, priority, start, stop, ra, dec):
    """Build dictionary of observational parameters."""

    def _to_datetime(val):
        if isinstance(val, datetime):
            return val
        if isinstance(val, str):
            # Try common ISO-like formats used in the project
            for fmt in ("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
                try:
                    return datetime.strptime(val, fmt)
                except ValueError:
                    continue
        return None

    start_dt = _to_datetime(start)
    stop_dt = _to_datetime(stop)

    if start_dt and stop_dt:
        start_format = start_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
        stop_format = stop_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    else:
        LOGGER.warning(
            "Unable to parse start/stop datetimes for target '%s' (start=%r stop=%r); "
            "writing times as-is",
            target_name,
            start,
            stop,
        )
        # Fallback: use provided values as-is
        start_format, stop_format = start, stop

    try:
        ra_value = f"{float(ra)}"
        dec_value = f"{float(dec)}"
    except (TypeError, ValueError):
        LOGGER.warning(
            "Unable to parse RA/DEC for target '%s' (ra=%r dec=%r); writing sentinel values",
            target_name,
            ra,
            dec,
        )
        ra_value, dec_value = "-999.0", "-999.0"

    return {
        "Target": target_name,
        "Priority": f"{priority}",
        "Timing": ["Start", "Stop", start_format, stop_format],
        "Boresight": ["RA", "DEC", ra_value, dec_value],
    }


def _duration_in_seconds(start, stop) -> float:
    """Calculate duration in seconds between start and stop times."""

    def _to_datetime(val):
        if isinstance(val, datetime):
            return val
        if isinstance(val, str):
            for fmt in ("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
                try:
                    return datetime.strptime(val, fmt)
                except ValueError:
                    continue
        return None

    start_dt = _to_datetime(start)
    stop_dt = _to_datetime(stop)
    if start_dt and stop_dt:
        return (stop_dt - start_dt).total_seconds()
    if start is not None or stop is not None:
        LOGGER.warning(
            "Unable to compute duration (seconds) from start/stop (start=%r stop=%r); using 0.0",
            start,
            stop,
        )
    return 0.0
