"""Data model for the shortschedule package.

This module defines the in-memory structures used to represent a
PAN-SCICAL science calendar: ObservationSequence, Visit and
ScienceCalendar. These classes include helper methods to read, update
and flatten payload XML elements and to compute useful summary
statistics used by the scheduler and visualizer.

Notes
-----
- Times are represented with `astropy.time.Time` objects and durations
    as `astropy.time.TimeDelta`/astropy units when appropriate.
- Payload parameters are stored as ElementTree elements to preserve XML
    structure; helpers are provided to convert to nested dicts or flattened
    dot-notation mappings.
"""

# Standard library
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

# Third-party
import astropy.units as u
import numpy as np
from astropy.time import Time, TimeDelta


class ObservationSequence:
    """Represents an observation sequence within a visit."""

    def __init__(
        self,
        id: str,
        target: str,
        priority: int,
        start_time: Time,
        stop_time: Time,
        ra: float,
        dec: float,
        payload_params: Dict[str, Any],
        roll: Optional[float] = None,
    ):
        self.id = id
        self.target = target
        self.priority = priority
        self.start_time = start_time
        self.stop_time = stop_time
        self.ra = ra
        self.dec = dec
        self.payload_params = payload_params
        self.roll = roll  # Spacecraft roll angle in degrees

    @property
    def duration(self) -> TimeDelta:
        """Return the sequence duration as an Astropy TimeDelta.

        Returns
        -------
        astropy.time.TimeDelta
            TimeDelta representing stop_time - start_time. The caller can
            convert to seconds/minutes/hours via `.sec` or astropy unit methods.
        """
        delta = self.stop_time - self.start_time  # TimeDelta
        return delta

    def copy(self) -> "ObservationSequence":
        """Create a deep copy of this observation sequence."""
        # Deep copy the payload parameters XML elements
        payload_copy = {}
        for key, element in self.payload_params.items():
            if isinstance(element, ET.Element):
                payload_copy[key] = self._copy_xml_element(element)
            else:
                payload_copy[key] = element

        return ObservationSequence(
            id=self.id,
            target=self.target,
            priority=self.priority,
            start_time=self.start_time,
            stop_time=self.stop_time,
            ra=self.ra,
            dec=self.dec,
            payload_params=payload_copy,
            roll=self.roll,
        )

    def _copy_xml_element(self, element):
        """Create a deep copy of an XML element."""
        # Convert to string and parse back to create a true copy
        xml_str = ET.tostring(element, encoding="unicode")
        return ET.fromstring(xml_str)

    @property
    def start_time_str(self) -> str:
        """Format start time as ISO string with Z suffix."""
        return self.start_time.isot

    @property
    def stop_time_str(self) -> str:
        """Format stop time as ISO string with Z suffix."""
        return self.stop_time.isot

    def get_payload_parameter(
        self, category: str, parameter_name: str, default: Any = None
    ) -> Any:
        """Retrieve a payload parameter from the stored XML payload element.

        Parameters
        ----------
        category : str
            Top-level payload category name (e.g. 'AcquireVisCamScienceData').
        parameter_name : str
            Child element name to retrieve under the category element.
        default : any, optional
            Value to return if the parameter or category is not present.

        Returns
        -------
        str or dict or None
            - If the element is a leaf with text, a stripped string is returned.
            - If the element has children, a dict mapping child tag -> value is returned.
            - If missing, `default` is returned.
        """
        if category not in self.payload_params:
            return default

        elem = self.payload_params[category].find(parameter_name)
        if elem is not None:
            # If element has text content and no children, return the text
            if elem.text and elem.text.strip() and len(elem) == 0:
                return elem.text.strip()
            # If element has children, return all child elements as a dictionary
            elif len(elem) > 0:
                children = {}
                for child in elem:
                    if child.text and child.text.strip():
                        children[child.tag] = child.text.strip()
                    else:
                        # Handle nested containers recursively if needed
                        children[child.tag] = None
                return children
            # If element exists but has no text and no children
            else:
                return None
        return default

    def get_nested_payload_parameter(
        self,
        category: str,
        parameter_name: str,
        sub_parameter_name: str,
        default: Any = None,
    ) -> Any:
        """Get nested payload parameter value."""
        parent_elem = self.get_payload_parameter(category, parameter_name)
        if parent_elem is not None and hasattr(parent_elem, "find"):
            child_elem = parent_elem.find(sub_parameter_name)
            if child_elem is not None and child_elem.text:
                return child_elem.text.strip()
        return default

    def set_payload_parameter(
        self, category: str, parameter_name: str, value: Any
    ) -> bool:
        """Set payload parameter value in XML structure."""
        if category not in self.payload_params:
            return False

        elem = self.payload_params[category].find(parameter_name)
        if elem is not None:
            elem.text = str(value)
            return True
        return False

    def get_all_payload_parameters(self) -> Dict[str, Any]:
        """Return all payload parameters as a nested dictionary.

        The returned structure converts XML elements into Python-native types
        (strings, dicts, lists) making it easier to inspect payload values in
        tests and analysis code.
        """
        all_params = {}

        for category, element in self.payload_params.items():
            if isinstance(element, ET.Element):
                all_params[category] = self._xml_element_to_clean_dict(element)
            else:
                all_params[category] = element

        return all_params

    def _xml_element_to_clean_dict(self, element: ET.Element) -> Any:
        """Convert XML element to dictionary with clean values (no _text wrappers)."""
        result = {}

        # If element has only text and no children, return the text directly
        if element.text and element.text.strip() and len(element) == 0:
            return element.text.strip()

        # Add attributes directly if they exist
        if element.attrib:
            for attr_name, attr_value in element.attrib.items():
                result[f"@{attr_name}"] = attr_value

        # Add child elements
        for child in element:
            child_name = child.tag
            child_value = self._xml_element_to_clean_dict(child)

            # Handle multiple children with same tag
            if child_name in result:
                if not isinstance(result[child_name], list):
                    result[child_name] = [result[child_name]]
                result[child_name].append(child_value)
            else:
                result[child_name] = child_value

        # If element has text content and children, add the text as a special key
        if element.text and element.text.strip() and len(element) > 0:
            result["_content"] = element.text.strip()

        return result

    def get_flat_payload_parameters(self) -> Dict[str, Any]:
        """Return a flattened mapping of payload parameters.

        Keys use dot-notation to represent nesting (e.g.
        'AcquireVisCamScienceData.ExposureTime_us'). Attributes are
        represented using the '@' prefix (e.g. 'Payload.@attr').
        """
        flat_params = {}

        for category, element in self.payload_params.items():
            if isinstance(element, ET.Element):
                self._flatten_xml_element(element, category, flat_params)
            else:
                flat_params[category] = element

        return flat_params

    def _flatten_xml_element(
        self, element: ET.Element, prefix: str, result_dict: Dict[str, Any]
    ) -> None:
        """Flatten XML element to dot-notation dictionary."""
        # Add element text if it exists
        if element.text and element.text.strip():
            result_dict[prefix] = element.text.strip()

        # Add attributes
        for attr_name, attr_value in element.attrib.items():
            result_dict[f"{prefix}.@{attr_name}"] = attr_value

        # Add child elements
        for child in element:
            child_key = f"{prefix}.{child.tag}"
            self._flatten_xml_element(child, child_key, result_dict)

    def __repr__(self):
        return f"<ObservationSequence {self.id}: {self.target} (P{self.priority}, {self.duration.sec / 60:.1f}min)>"


class Visit:
    """Represents a visit in the science calendar."""

    def __init__(self, id, sequences):
        self.id: str = id
        self.sequences: List["ObservationSequence"] = sequences

    @property
    def total_duration(self):
        """Total duration of all sequences in this visit."""
        return np.sum([seq.duration for seq in self.sequences])

    @property
    def total_duration_minutes(self):
        """Total duration of all sequences in this visit (minutes)."""
        return self.total_duration.sec / 60.0

    @property
    def start_time(self):
        """Start time of the first sequence in this visit."""
        if not self.sequences:
            return None
        return np.min([seq.start_time for seq in self.sequences])

    @property
    def end_time(self):
        """End time of the last sequence in this visit."""
        if not self.sequences:
            return None
        return np.max([seq.stop_time for seq in self.sequences])

    def copy(
        self, sequences: Optional[List["ObservationSequence"]] = None
    ) -> "Visit":
        """Create a copy of this visit, optionally with different sequences."""
        if sequences is None:
            sequences = [seq.copy() for seq in self.sequences]
        return Visit(id=self.id, sequences=sequences)

    def __repr__(self):
        return f"<Visit {self.id}: {len(self.sequences)} sequences, {self.total_duration_minutes:.1f}min>"


class ScienceCalendar:
    """Represents a complete Science Calendar."""

    def __init__(
        self,
        metadata: Optional[Dict[str, Any]],
        visits: List[Visit],
        visibility: Any = None,
    ):
        self.metadata: Dict[str, Any] = metadata or {}
        self.visits: List[Visit] = visits
        self.visibility = visibility

    def set_visibility_calculator(self, visibility: Any) -> None:
        """Set or update the visibility calculator."""
        self.visibility = visibility

    @property
    def total_sequences(self):
        """Total number of observation sequences."""
        return np.sum([len(visit.sequences) for visit in self.visits])

    @property
    def total_duration(self):
        """Total duration of all observations in minutes."""
        return np.sum([visit.total_duration for visit in self.visits])

    @property
    def total_duration_minutes(self) -> float:
        """Total duration of all observations in minutes."""
        return self.total_duration.to(u.s).value / 60

    @property
    def total_duration_hours(self) -> float:
        """Total duration of all observations in hours."""
        return self.total_duration_minutes / 60

    @property
    def total_duration_days(self) -> float:
        """Total duration of all observations in days."""
        return self.total_duration_minutes / (60 * 24)

    @property
    def calendar_span(self):
        """Total span of the calendar from first to last observation."""
        if not self.visits:
            return 0

        all_start_times = []
        all_end_times = []

        for visit in self.visits:
            if visit.start_time:
                all_start_times.append(visit.start_time)
            if visit.end_time:
                all_end_times.append(visit.end_time)

        if not all_start_times or not all_end_times:
            return 0

        earliest = np.min(all_start_times)
        latest = np.max(all_end_times)

        return latest - earliest

    @property
    def calendar_span_days(self) -> float:
        """Total span of the calendar from first to last observation in days."""
        return self.calendar_span.to(u.d).value

    @property
    def duty_cycle_percent(self) -> float:
        """Percentage of calendar span that is actually observing."""
        span_days = self.calendar_span_days
        if span_days == 0:
            return 0
        return (self.total_duration_days / span_days) * 100

    @property
    def priority_breakdown(self):
        """Breakdown of sequences by priority level."""
        breakdown = {}
        for visit in self.visits:
            for seq in visit.sequences:
                priority = seq.priority
                if priority not in breakdown:
                    breakdown[priority] = {"count": 0, "duration_minutes": 0}
                breakdown[priority]["count"] += 1
                breakdown[priority]["duration_minutes"] += (
                    seq.duration.sec / 60
                )

        return breakdown

    @property
    def date_range(self):
        """Start and end dates of the calendar."""
        if not self.visits:
            return None, None

        all_times = []
        for visit in self.visits:
            for seq in visit.sequences:
                all_times.extend([seq.start_time, seq.stop_time])

        if not all_times:
            return None, None

        return np.min(all_times), np.max(all_times)

    def get_summary_stats(self) -> Dict[str, Any]:
        """Get comprehensive summary statistics."""
        start_date, end_date = self.date_range
        priority_stats = self.priority_breakdown

        stats = {
            "total_visits": len(self.visits),
            "total_sequences": self.total_sequences,
            "total_duration_minutes": self.total_duration_minutes,
            "total_duration_hours": self.total_duration_hours,
            "total_duration_days": self.total_duration_days,
            "calendar_span_days": self.calendar_span_days,
            "duty_cycle_percent": self.duty_cycle_percent,
            "start_date": start_date.isot if start_date else None,
            "end_date": end_date.isot if end_date else None,
            "priority_breakdown": {},
        }

        # Format priority breakdown
        for priority, data in priority_stats.items():
            stats["priority_breakdown"][f"priority_{priority}"] = {
                "count": data["count"],
                "duration_hours": data["duration_minutes"] / 60,
                "duration_days": data["duration_minutes"] / (60 * 24),
            }

        return stats

    def copy(self) -> "ScienceCalendar":
        """Create a complete deep copy of this calendar."""
        copied_visits = [visit.copy() for visit in self.visits]
        return ScienceCalendar(
            metadata=self.metadata.copy() if self.metadata else {},
            visits=copied_visits,
            visibility=self.visibility,
        )

    def get_sequence(
        self, visit_id: str, sequence_id: str
    ) -> Optional["ObservationSequence"]:
        """Get observation sequence by visit ID and sequence ID."""
        for visit in self.visits:
            if visit.id == visit_id:
                for sequence in visit.sequences:
                    if sequence.id == sequence_id:
                        return sequence
        return None

    def replace_sequence(
        self,
        visit_id: str,
        sequence_id: str,
        new_sequence: "ObservationSequence",
    ) -> bool:
        """
        Replace an existing sequence with a new one.

        Args:
            visit_id: ID of the visit containing the sequence
            sequence_id: ID of the sequence to replace
            new_sequence: New ObservationSequence object to replace with

        Returns:
            bool: True if replacement was successful, False if sequence not found
        """
        for visit in self.visits:
            if visit.id == visit_id:
                for i, sequence in enumerate(visit.sequences):
                    if sequence.id == sequence_id:
                        visit.sequences[i] = new_sequence
                        return True
        return False

    def __repr__(self):
        return (
            f"<ScienceCalendar: {len(self.visits)} visits, "
            f"{self.total_sequences} sequences, "
            f"{self.total_duration_days:.2f} days>"
        )
