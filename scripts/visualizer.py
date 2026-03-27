"""Visualization helpers for schedule analysis.

This module provides `ScheduleVisualizer`, a class to produce multiple
comparison plots (Gantt timelines, visit summaries, priority statistics)
based on `ScienceCalendar` objects and the scheduler's gap report.

The visualizer attempts to be robust to large time spans by windowing
plots automatically and provides multiple plotting entry points used by
the documentation and analysis notebooks.
"""

# Standard library
import argparse
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import timedelta
from pathlib import Path
from typing import Any, Optional, Tuple

# Third-party
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle


@dataclass(frozen=True)
class Sequence:
    id: str
    target: str
    priority: int
    start_time: Time
    stop_time: Time
    ra: float = float("nan")
    dec: float = float("nan")

    @property
    def duration(self):
        return self.stop_time - self.start_time


@dataclass(frozen=True)
class Visit:
    id: str
    sequences: list[Sequence] = field(default_factory=list)

    @property
    def start_time(self) -> Optional[Time]:
        if not self.sequences:
            return None
        return min(seq.start_time for seq in self.sequences)

    @property
    def end_time(self) -> Optional[Time]:
        if not self.sequences:
            return None
        return max(seq.stop_time for seq in self.sequences)

    @property
    def total_duration_minutes(self) -> float:
        return sum(seq.duration.sec for seq in self.sequences) / 60.0


@dataclass(frozen=True)
class ScienceCalendar:
    metadata: dict[str, str]
    visits: list[Visit] = field(default_factory=list)


class _DummyScheduler:
    def get_gap_report(self):
        return {}


def _child_text(parent: ET.Element, path: str, default: str = "") -> str:
    node = parent.find(path)
    if node is None or node.text is None:
        return default
    return node.text.strip()


def parse_science_calendar(xml_path: Path) -> ScienceCalendar:
    tree = ET.parse(xml_path)
    root = tree.getroot()

    metadata: dict[str, str] = {}
    for child in root:
        tag = child.tag.rsplit("}", 1)[-1]
        if tag != "Visit" and child.text is not None:
            metadata[tag] = child.text.strip()

    visits: list[Visit] = []
    visit_elements = root.findall(".//{*}Visit")
    for visit_index, visit_el in enumerate(visit_elements, start=1):
        visit_id = _child_text(visit_el, "./{*}ID", str(visit_index))
        sequences: list[Sequence] = []
        for seq_el in visit_el.findall("./{*}Observation_Sequence"):
            seq_id = _child_text(seq_el, "./{*}ID", f"{visit_id}-1")
            target = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Target",
                "UNKNOWN",
            )
            priority_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Priority",
                "0",
            )
            start_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Timing/{*}Start",
            )
            stop_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Timing/{*}Stop",
            )
            ra_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Boresight/{*}RA",
                "nan",
            )
            dec_text = _child_text(
                seq_el,
                "./{*}Observational_Parameters/{*}Boresight/{*}DEC",
                "nan",
            )

            if not start_text or not stop_text:
                continue

            try:
                priority = int(float(priority_text))
            except ValueError:
                priority = 0

            try:
                ra = float(ra_text)
            except ValueError:
                ra = float("nan")

            try:
                dec = float(dec_text)
            except ValueError:
                dec = float("nan")

            sequences.append(
                Sequence(
                    id=seq_id,
                    target=target,
                    priority=priority,
                    start_time=Time(start_text, format="isot", scale="utc"),
                    stop_time=Time(stop_text, format="isot", scale="utc"),
                    ra=ra,
                    dec=dec,
                )
            )

        if sequences:
            visits.append(Visit(id=visit_id, sequences=sequences))

    return ScienceCalendar(metadata=metadata, visits=visits)


def _save_figures(figures: list[Figure], output_path: Path) -> list[Path]:
    if len(figures) == 1:
        figures[0].savefig(output_path, dpi=150, bbox_inches="tight")
        return [output_path]

    saved: list[Path] = []
    stem = output_path.stem
    suffix = output_path.suffix or ".png"
    for idx, fig in enumerate(figures, start=1):
        candidate = output_path.with_name(f"{stem}_window_{idx}{suffix}")
        fig.savefig(candidate, dpi=150, bbox_inches="tight")
        saved.append(candidate)
    return saved


class ScheduleVisualizer:
    """Class for creating visualizations of schedule analysis.

    Primary methods include:
    - plot_gantt_timeline(original_calendar, processed_calendar)
    - plot_gantt_timeline_with_visits(...)
    - plot_visit_summary_timeline(...)

    The visualizer reads the scheduler's `gap_report` to annotate plots
    and produce comparison summaries.
    """

    def __init__(self, scheduler: Any) -> None:
        self.scheduler = scheduler
        self.gap_report = scheduler.get_gap_report()

    def plot_gantt_timeline(
        self,
        original_calendar: ScienceCalendar,
        processed_calendar: ScienceCalendar,
        figsize: Tuple[int, int] = (16, 10),
    ) -> Figure:
        """Create a proper Gantt chart style timeline."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

        # Plot both calendars
        self._plot_gantt_chart(original_calendar, ax1, "Original Calendar")
        self._plot_gantt_chart(processed_calendar, ax2, "Processed Calendar")

        # Get overall time range
        all_times = []
        for calendar in [original_calendar, processed_calendar]:
            for visit in calendar.visits:
                for seq in visit.sequences:
                    all_times.extend(
                        [seq.start_time.datetime, seq.stop_time.datetime]
                    )

        if all_times:
            min_time = min(all_times)
            max_time = max(all_times)

            # Add some padding
            time_span = max_time - min_time
            padding = time_span * 0.02
            ax1.set_xlim(min_time - padding, max_time + padding)

            # Format x-axis
            if time_span.total_seconds() < 86400:  # Less than 1 day
                ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax2.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            elif time_span.total_seconds() < 86400 * 7:  # Less than 1 week
                ax2.xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d %H:%M")
                )
                ax2.xaxis.set_major_locator(mdates.HourLocator(interval=6))
            else:
                ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
                ax2.xaxis.set_major_locator(mdates.DayLocator())

        plt.xticks(rotation=45)
        fig.suptitle(
            "Schedule Comparison: Original vs Processed", fontsize=16, y=0.95
        )
        fig.set_constrained_layout(True)

        return fig

    def plot_gantt_timeline_with_visits(
        self,
        original_calendar: ScienceCalendar,
        processed_calendar: ScienceCalendar,
        figsize: Tuple[int, int] = (20, 14),
        show_sequence_labels: bool = False,
        processed_only: bool = False,
    ) -> Figure:
        """Create a Gantt chart organized by visits with cleaner labeling."""

        if processed_only:
            fig, ax = plt.subplots(1, 1, figsize=(figsize[0], figsize[1] / 2))
            changes = self._analyze_sequence_changes(
                original_calendar, processed_calendar
            )
            self._plot_gantt_by_visits(
                processed_calendar,
                ax,
                "Processed Calendar",
                changes,
                is_original=False,
                show_sequence_labels=show_sequence_labels,
            )
            axes = [ax]
        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
            changes = self._analyze_sequence_changes(
                original_calendar, processed_calendar
            )
            self._plot_gantt_by_visits(
                original_calendar,
                ax1,
                "Original Calendar",
                changes,
                is_original=True,
                show_sequence_labels=show_sequence_labels,
            )
            self._plot_gantt_by_visits(
                processed_calendar,
                ax2,
                "Processed Calendar",
                changes,
                is_original=False,
                show_sequence_labels=show_sequence_labels,
            )
            axes = [ax1, ax2]

        # Format time axis
        self._format_time_axis_gantt(
            axes, original_calendar, processed_calendar
        )

        # Add legend
        self._add_change_legend(fig)

        # Title is handled by caller or figure annotations; no local variable needed
        # fig.suptitle("Schedule Comparison by Visits: Original vs Processed", fontsize=16, y=0.96)
        fig.set_constrained_layout(True)

        return fig

    def _plot_gantt_by_visits(
        self,
        calendar,
        ax,
        title,
        changes,
        is_original=False,
        show_sequence_labels=True,
    ):
        """Plot Gantt chart organized by visits with optional sequence labels."""
        if not calendar.visits:
            ax.text(
                0.5,
                0.5,
                "No sequences found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(title)
            return

        # Sort visits by start time
        sorted_visits = sorted(
            calendar.visits,
            key=lambda v: v.start_time if v.start_time else Time("2000-01-01"),
        )

        # Change highlight colors
        change_colors = {
            "extended": "lightgreen",
            "shortened": "orange",
            "unchanged": None,  # Use base color
            "new": "yellow",
            "removed": "lightgray",
        }

        # Color scheme for different targets
        all_targets = set()
        for visit in sorted_visits:
            for seq in visit.sequences:
                all_targets.add(seq.target)

        target_colors = {}
        colors = plt.cm.Set3(np.linspace(0, 1, len(all_targets)))
        for i, target in enumerate(sorted(all_targets)):
            target_colors[target] = colors[i]

        y_pos = 0
        y_labels = []
        y_positions = []
        visit_boundaries = []

        for visit_idx, visit in enumerate(sorted_visits):
            visit_start_y = y_pos

            # Group sequences by target within this visit
            sequences_by_target = {}
            for seq in visit.sequences:
                target = seq.target
                if target not in sequences_by_target:
                    sequences_by_target[target] = []
                sequences_by_target[target].append(seq)

            # Sort sequences within each target by start time
            for target in sequences_by_target:
                sequences_by_target[target].sort(key=lambda s: s.start_time)

            # Sort targets by their first sequence start time
            sorted_targets = sorted(
                sequences_by_target.keys(),
                key=lambda t: sequences_by_target[t][0].start_time,
            )

            # Plot each target within this visit
            for target in sorted_targets:
                sequences = sequences_by_target[target]

                # Plot all sequences for this target
                for seq_idx, seq in enumerate(sequences):
                    start_time = seq.start_time.datetime
                    duration = seq.duration.sec / 3600  # Convert to hours

                    # Create sequence key for change detection
                    seq_key = f"{visit.id}_{seq.id}_{seq.start_time.isot}_{seq.target}"

                    # Determine color based on changes
                    base_color = target_colors[target]
                    change_type = None
                    for change, seq_set in changes.items():
                        if seq_key in seq_set:
                            change_type = change
                            break

                    if change_type and change_colors[change_type]:
                        face_color = change_colors[change_type]
                        alpha = 0.9
                        edge_color = (
                            "red"
                            if change_type in ["extended", "shortened"]
                            else "black"
                        )
                        linewidth = (
                            2
                            if change_type in ["extended", "shortened"]
                            else 1
                        )
                    else:
                        face_color = base_color
                        alpha = 0.7
                        edge_color = "black"
                        linewidth = 1

                    # Create rectangle for sequence
                    rect = Rectangle(
                        (mdates.date2num(start_time), y_pos - 0.35),
                        duration / 24,  # Convert hours to days for matplotlib
                        0.7,
                        facecolor=face_color,
                        alpha=alpha,
                        edgecolor=edge_color,
                        linewidth=linewidth,
                    )
                    ax.add_patch(rect)

                    # Add sequence label only if requested and sequence is long enough
                    if (
                        show_sequence_labels and duration > 0.5
                    ):  # More than 30 minutes
                        mid_time = mdates.date2num(start_time) + (
                            duration / 24 / 2
                        )
                        ax.text(
                            mid_time,
                            y_pos,
                            seq.id,
                            ha="center",
                            va="center",
                            fontsize=6,
                            fontweight="bold",
                            bbox=dict(
                                boxstyle="round,pad=0.1",
                                facecolor="white",
                                alpha=0.8,
                            ),
                        )

                # Create label for this target within the visit
                y_labels.append(f"  {target} ({len(sequences)})")
                y_positions.append(y_pos)
                y_pos += 1

            visit_end_y = y_pos - 1
            visit_boundaries.append((visit_start_y, visit_end_y))

            # Add visit separator line
            if visit_idx < len(sorted_visits) - 1:
                ax.axhline(
                    y=y_pos - 0.5, color="black", linewidth=2, alpha=0.7
                )

            y_pos += 0.2  # Small gap between visits

        # Set up y-axis
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.set_ylim(-0.5, y_pos - 0.7)

        # Add visit labels with better positioning to avoid overlap
        for visit_idx, (visit, (start_y, end_y)) in enumerate(
            zip(sorted_visits, visit_boundaries)
        ):
            mid_y = (start_y + end_y) / 2
            # Normalize y position and offset more to avoid overlap
            norm_y = mid_y / (y_pos - 0.7)
            visit_info = (
                f"Visit {visit.id}\n({len([s for s in visit.sequences])} seq)"
            )

            # Position labels further to the right and stagger them
            x_offset = (
                1.05 + (visit_idx % 2) * 0.08
            )  # Stagger alternating visits
            ax.text(
                x_offset,
                norm_y,
                visit_info,
                transform=ax.transAxes,
                ha="left",
                va="center",
                fontsize=9,
                fontweight="bold",
                bbox=dict(
                    boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7
                ),
            )

        # Grid and labels
        ax.grid(True, alpha=0.3, axis="x")
        ax.set_ylabel("Targets (grouped by visit)", fontsize=12)
        ax.set_title(title, fontsize=14, pad=10)

    def plot_visit_summary_timeline(
        self,
        original_calendar: ScienceCalendar,
        processed_calendar: ScienceCalendar,
        figsize: Tuple[int, int] = (16, 8),
        processed_only: bool = False,
    ) -> Figure:
        """Create a high-level timeline showing just visits and their overall spans."""

        if processed_only:
            fig, ax = plt.subplots(1, 1, figsize=(figsize[0], figsize[1] / 2))
            self._plot_visit_bars(
                processed_calendar, ax, "Processed Calendar - Visit Overview"
            )
            axes = [ax]
        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
            self._plot_visit_bars(
                original_calendar, ax1, "Original Calendar - Visit Overview"
            )
            self._plot_visit_bars(
                processed_calendar, ax2, "Processed Calendar - Visit Overview"
            )
            axes = [ax1, ax2]

        self._format_time_axis_visits(
            axes, original_calendar, processed_calendar
        )

        # Title is handled by caller or figure annotations; no local variable needed
        # fig.suptitle("Visit-Level Schedule Comparison", fontsize=16, y=0.95)
        fig.set_constrained_layout(True)

        return fig

    def _plot_visit_bars(self, calendar, ax, title):
        """Plot visit-level summary bars with better label positioning."""
        if not calendar.visits:
            ax.text(
                0.5,
                0.5,
                "No visits found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(title)
            return

        sorted_visits = sorted(
            calendar.visits,
            key=lambda v: v.start_time if v.start_time else Time("2000-01-01"),
        )
        colors = plt.cm.Set1(np.linspace(0, 1, len(sorted_visits)))

        for i, visit in enumerate(sorted_visits):
            if not visit.sequences:
                continue

            start_time = visit.start_time.datetime
            end_time = visit.end_time.datetime
            duration_hours = (end_time - start_time).total_seconds() / 3600

            # Create bar for entire visit span
            ax.barh(
                i,
                duration_hours / 24,
                left=mdates.date2num(start_time),
                height=0.6,
                color=colors[i],
                alpha=0.7,
                edgecolor="black",
                linewidth=1,
            )

            # Add visit info with better positioning
            mid_time = mdates.date2num(start_time) + (duration_hours / 24 / 2)

            # Only show labels if the visit bar is wide enough
            if duration_hours > 2:  # More than 2 hours
                visit_info = f"Visit {visit.id}\n{len(visit.sequences)} seq\n{visit.total_duration_minutes:.0f} min"
                ax.text(
                    mid_time,
                    i,
                    visit_info,
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                    bbox=dict(
                        boxstyle="round,pad=0.2", facecolor="white", alpha=0.9
                    ),
                )
            else:
                # For narrow visits, put label to the right
                visit_info = f"V{visit.id}"
                right_edge = mdates.date2num(start_time) + (
                    duration_hours / 24
                )
                ax.text(
                    right_edge + 0.001,
                    i,
                    visit_info,
                    ha="left",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                )

        ax.set_yticks(range(len(sorted_visits)))
        ax.set_yticklabels([f"Visit {v.id}" for v in sorted_visits])
        ax.set_ylabel("Visits")
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis="x")

    def _format_time_axis_gantt(
        self, axes, original_calendar, processed_calendar
    ):
        """Format time axis for gantt charts."""
        # Get overall time range
        all_times = []
        for calendar in [original_calendar, processed_calendar]:
            for visit in calendar.visits:
                for seq in visit.sequences:
                    all_times.extend(
                        [seq.start_time.datetime, seq.stop_time.datetime]
                    )

        if all_times:
            min_time = min(all_times)
            max_time = max(all_times)

            # Add some padding
            time_span = max_time - min_time
            padding = time_span * 0.01
            for ax in axes:
                ax.set_xlim(min_time - padding, max_time + padding)

            # Format x-axis based on time span
            if time_span.total_seconds() < 86400:  # Less than 1 day
                axes[-1].xaxis.set_major_formatter(
                    mdates.DateFormatter("%H:%M")
                )
                axes[-1].xaxis.set_major_locator(
                    mdates.HourLocator(interval=3)
                )
            elif time_span.total_seconds() < 86400 * 7:  # Less than 1 week
                axes[-1].xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d %H:%M")
                )
                axes[-1].xaxis.set_major_locator(
                    mdates.HourLocator(interval=12)
                )
            else:
                axes[-1].xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d")
                )
                axes[-1].xaxis.set_major_locator(mdates.DayLocator())

        plt.xticks(rotation=45)

    def _format_time_axis_visits(
        self, axes, original_calendar, processed_calendar
    ):
        """Format time axis for visit overview charts."""
        self._format_time_axis_gantt(
            axes, original_calendar, processed_calendar
        )

    # Updated generate_full_report method
    def generate_full_report(
        self,
        original_calendar: ScienceCalendar,
        processed_calendar: ScienceCalendar,
        save_path: Optional[str] = None,
        show_sequence_labels: bool = False,
        processed_only: bool = False,
    ) -> Tuple[Figure, Figure, Figure, pd.DataFrame]:
        """Generate complete visualization report with options."""

        # Create visit-organized Gantt chart
        fig1 = self.plot_gantt_timeline_with_visits(
            original_calendar,
            processed_calendar,
            show_sequence_labels=show_sequence_labels,
            processed_only=processed_only,
        )

        # Create visit-level summary
        fig2 = self.plot_visit_summary_timeline(
            original_calendar,
            processed_calendar,
            processed_only=processed_only,
        )

        # Create comprehensive summary (always show comparison)
        fig3 = self.plot_schedule_comparison_summary()

        # Create comparison table
        comparison_df = self.create_sequence_comparison_table(
            original_calendar, processed_calendar
        )

        if save_path:
            try:
                suffix = "_processed_only" if processed_only else "_comparison"
                fig1.savefig(
                    f"{save_path}_gantt{suffix}.png",
                    dpi=150,
                    bbox_inches="tight",
                )
                fig2.savefig(
                    f"{save_path}_visits{suffix}.png",
                    dpi=150,
                    bbox_inches="tight",
                )
                fig3.savefig(
                    f"{save_path}_summary.png", dpi=150, bbox_inches="tight"
                )
                comparison_df.to_csv(
                    f"{save_path}_sequence_comparison.csv", index=False
                )
                print(f"Plots and data saved to {save_path}_*")
            except Exception as e:
                print(f"Warning: Could not save plots: {e}")

        return fig1, fig2, fig3, comparison_df

    def _plot_gantt_chart(
        self, calendar: ScienceCalendar, ax: Any, title: str
    ) -> None:
        """Plot calendar as a proper Gantt chart with targets as rows."""
        if not calendar.visits:
            ax.text(
                0.5,
                0.5,
                "No sequences found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(title)
            return

        # Get all sequences and organize by target
        sequences_by_target = {}
        for visit in calendar.visits:
            for seq in visit.sequences:
                target = seq.target
                if target not in sequences_by_target:
                    sequences_by_target[target] = []
                sequences_by_target[target].append(seq)

        # Sort sequences within each target by start time
        for target in sequences_by_target:
            sequences_by_target[target].sort(key=lambda s: s.start_time)

        # Sort targets by first sequence start time
        sorted_targets = sorted(
            sequences_by_target.keys(),
            key=lambda t: sequences_by_target[t][0].start_time,
        )

        # Color scheme
        colors = plt.cm.Set3(np.linspace(0, 1, len(sorted_targets)))
        target_colors = {
            target: colors[i] for i, target in enumerate(sorted_targets)
        }

        # Plot each target's sequences
        y_pos = 0
        y_labels = []
        y_positions = []

        for target in sorted_targets:
            sequences = sequences_by_target[target]

            for seq in sequences:
                start_time = seq.start_time.datetime
                duration = seq.duration.sec / 3600  # Convert to hours

                # Create rectangle for sequence
                rect = Rectangle(
                    (mdates.date2num(start_time), y_pos - 0.4),
                    duration / 24,  # Convert hours to days for matplotlib
                    0.8,
                    facecolor=target_colors[target],
                    alpha=0.8,
                    edgecolor="black",
                    linewidth=1,
                )
                ax.add_patch(rect)

                # Add sequence ID label
                mid_time = mdates.date2num(start_time) + (duration / 24 / 2)
                ax.text(
                    mid_time,
                    y_pos,
                    seq.id,
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                )

            y_labels.append(target)
            y_positions.append(y_pos)
            y_pos += 1

        # Set up y-axis
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels)
        ax.set_ylim(-0.5, len(sorted_targets) - 0.5)

        # Grid and labels
        ax.grid(True, alpha=0.3, axis="x")
        ax.set_ylabel("Targets")
        ax.set_title(title, fontsize=14, pad=10)

        # Invert y-axis so first target is at top
        ax.invert_yaxis()

    def plot_schedule_comparison_summary(self, figsize=(15, 8)):
        """Create a comprehensive summary comparison plot."""
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

        # 1. Duration comparison
        ax1 = fig.add_subplot(gs[0, 0])
        self._plot_duration_comparison(ax1)

        # 2. Duty cycle comparison
        ax2 = fig.add_subplot(gs[0, 1])
        self._plot_duty_cycle_comparison(ax2)

        # 3. Sequence count comparison
        ax3 = fig.add_subplot(gs[0, 2])
        self._plot_sequence_count_comparison(ax3)

        # 4. Gap analysis
        ax4 = fig.add_subplot(gs[1, :2])
        self._plot_gap_analysis(ax4)

        # 5. Processing summary
        ax5 = fig.add_subplot(gs[1, 2])
        self._plot_processing_summary(ax5)

        fig.suptitle("Schedule Processing Analysis Summary", fontsize=16)
        return fig

    def _plot_duration_comparison(self, ax):
        """Plot duration comparison."""
        original = self.gap_report.get("original_calendar_stats", {})
        processed = self.gap_report.get("processed_calendar_stats", {})

        categories = ["Total Duration"]
        original_hours = original.get("total_duration_hours", 0)
        processed_hours = processed.get("total_duration_hours", 0)

        x = np.arange(len(categories))
        width = 0.35

        bars1 = ax.bar(
            x - width / 2,
            [original_hours],
            width,
            label="Original",
            alpha=0.7,
            color="lightblue",
        )
        bars2 = ax.bar(
            x + width / 2,
            [processed_hours],
            width,
            label="Processed",
            alpha=0.7,
            color="lightgreen",
        )

        ax.set_ylabel("Hours")
        ax.set_title("Total Duration Comparison")
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()

        # Add value labels
        for bar, value in zip(bars1, [original_hours]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.1,
                f"{value:.1f}h",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        for bar, value in zip(bars2, [processed_hours]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.1,
                f"{value:.1f}h",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    def _plot_duty_cycle_comparison(self, ax):
        """Plot duty cycle comparison."""
        original = self.gap_report.get("original_calendar_stats", {})
        processed = self.gap_report.get("processed_calendar_stats", {})

        original_dc = original.get("duty_cycle_percent", 0)
        processed_dc = processed.get("duty_cycle_percent", 0)

        categories = ["Duty Cycle"]
        x = np.arange(len(categories))
        width = 0.35

        bars1 = ax.bar(
            x - width / 2,
            [original_dc],
            width,
            label="Original",
            alpha=0.7,
            color="lightcoral",
        )
        bars2 = ax.bar(
            x + width / 2,
            [processed_dc],
            width,
            label="Processed",
            alpha=0.7,
            color="lightgreen",
        )

        ax.set_ylabel("Percentage (%)")
        ax.set_title("Duty Cycle Comparison")
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()

        # Add value labels
        for bar, value in zip(bars1, [original_dc]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{value:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        for bar, value in zip(bars2, [processed_dc]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{value:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    def _plot_sequence_count_comparison(self, ax):
        """Plot sequence count comparison."""
        original = self.gap_report.get("original_calendar_stats", {})
        processed = self.gap_report.get("processed_calendar_stats", {})

        original_count = original.get("total_sequences", 0)
        processed_count = processed.get("total_sequences", 0)

        categories = ["Sequence Count"]
        x = np.arange(len(categories))
        width = 0.35

        bars1 = ax.bar(
            x - width / 2,
            [original_count],
            width,
            label="Original",
            alpha=0.7,
            color="lightblue",
        )
        bars2 = ax.bar(
            x + width / 2,
            [processed_count],
            width,
            label="Processed",
            alpha=0.7,
            color="lightgreen",
        )

        ax.set_ylabel("Count")
        ax.set_title("Sequence Count Comparison")
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()

        # Add value labels
        for bar, value in zip(bars1, [original_count]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{int(value)}",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        for bar, value in zip(bars2, [processed_count]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                f"{int(value)}",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    def _plot_gap_analysis(self, ax):
        """Plot gap analysis."""
        summary = self.gap_report.get("processing_summary", {})

        gaps_found = summary.get("total_gaps_found", 0)
        gaps_filled = summary.get("gaps_filled", 0)
        gaps_remaining = summary.get("gaps_remaining", 0)

        if gaps_found > 0:
            labels = ["Filled", "Remaining"]
            sizes = [gaps_filled, gaps_remaining]
            colors = ["lightgreen", "lightcoral"]

            wedges, texts, autotexts = ax.pie(
                sizes,
                labels=labels,
                colors=colors,
                autopct="%1.1f%%",
                startangle=90,
            )
            ax.set_title(f"Gap Resolution\n({gaps_found} total gaps found)")
        else:
            ax.text(
                0.5,
                0.5,
                "No gap data available",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title("Gap Resolution")

    def _plot_processing_summary(self, ax):
        """Plot processing summary as text."""
        summary = self.gap_report.get("processing_summary", {})

        text_lines = [
            "Processing Summary:",
            "",
            f"Duration Improvement: {summary.get('duration_improvement_hours', 0):.1f} hrs",
            f"Duty Cycle Improvement: {summary.get('duty_cycle_improvement_percent', 0):.1f}%",
            f"Sequences Modified: {summary.get('sequences_modified', 0)}",
            f"Gaps Filled: {summary.get('gaps_filled', 0)}/{summary.get('gaps_filled', 0) + summary.get('gaps_remaining', 0)}",
        ]

        ax.text(
            0.05,
            0.95,
            "\n".join(text_lines),
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment="top",
            fontfamily="monospace",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    def plot_gantt_timeline_with_changes(
        self, original_calendar, processed_calendar, figsize=(18, 12)
    ):
        """Create a Gantt chart that highlights changes between original and processed."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

        # Get sequence changes for highlighting
        changes = self._analyze_sequence_changes(
            original_calendar, processed_calendar
        )

        # Plot both calendars with change highlighting
        self._plot_gantt_with_changes(
            original_calendar,
            ax1,
            "Original Calendar",
            changes,
            is_original=True,
        )
        self._plot_gantt_with_changes(
            processed_calendar,
            ax2,
            "Processed Calendar",
            changes,
            is_original=False,
        )

        # Get overall time range
        all_times = []
        for calendar in [original_calendar, processed_calendar]:
            for visit in calendar.visits:
                for seq in visit.sequences:
                    all_times.extend(
                        [seq.start_time.datetime, seq.stop_time.datetime]
                    )

        if all_times:
            min_time = min(all_times)
            max_time = max(all_times)

            # Add some padding
            time_span = max_time - min_time
            padding = time_span * 0.01
            ax1.set_xlim(min_time - padding, max_time + padding)

            # Format x-axis based on time span
            if time_span.total_seconds() < 86400:  # Less than 1 day
                ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax2.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            elif time_span.total_seconds() < 86400 * 7:  # Less than 1 week
                ax2.xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d %H:%M")
                )
                ax2.xaxis.set_major_locator(mdates.HourLocator(interval=6))
            else:
                ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
                ax2.xaxis.set_major_locator(mdates.DayLocator())

        plt.xticks(rotation=45)

        # Add legend for changes
        # Third-party
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(
                facecolor="lightgreen", alpha=0.9, label="Extended Sequences"
            ),
            Patch(facecolor="orange", alpha=0.9, label="Shortened Sequences"),
            Patch(facecolor="yellow", alpha=0.9, label="New Sequences"),
            Patch(facecolor="lightgray", alpha=0.9, label="Removed Sequences"),
            Patch(
                facecolor="lightblue", alpha=0.7, label="Unchanged Sequences"
            ),
        ]

        fig.legend(
            handles=legend_elements,
            loc="upper right",
            bbox_to_anchor=(0.98, 0.98),
        )

        fig.suptitle(
            "Schedule Comparison: Original vs Processed (Changes Highlighted)",
            fontsize=16,
            y=0.95,
        )
        fig.set_constrained_layout(True)

        return fig

    def _get_sequence_display_key(self, visit_id, seq):
        """Create a display-friendly key for sequences."""
        return f"{visit_id}_{seq.id}_{seq.target}_{seq.start_time.datetime.strftime('%H:%M')}"

    def _add_change_legend(self, fig):
        """Add legend explaining the change highlighting."""
        # Third-party
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(
                facecolor="lightgreen", alpha=0.9, label="Extended Sequences"
            ),
            Patch(facecolor="orange", alpha=0.9, label="Shortened Sequences"),
            Patch(facecolor="yellow", alpha=0.9, label="New Sequences"),
            Patch(facecolor="lightgray", alpha=0.9, label="Removed Sequences"),
            Patch(
                facecolor="lightblue", alpha=0.7, label="Unchanged Sequences"
            ),
        ]

        fig.legend(
            handles=legend_elements,
            loc="upper right",
            bbox_to_anchor=(0.98, 0.98),
        )

    def _analyze_sequence_changes(self, original_calendar, processed_calendar):
        """Analyze what changed between calendars using robust sequence matching."""
        changes = {
            "extended": set(),  # Sequences that were extended
            "shortened": set(),  # Sequences that were shortened
            "unchanged": set(),  # Sequences that stayed the same
            "new": set(),  # New sequences (from splitting/gap filling)
            "removed": set(),  # Sequences that were removed
        }

        # Create robust identifiers for sequences (visit_id + seq_id + start_time + target)
        def create_sequence_key(visit_id, seq):
            return f"{visit_id}_{seq.id}_{seq.start_time.isot}_{seq.target}"

        def create_base_sequence_key(visit_id, seq):
            """Create key for matching split sequences (ignore sub-sequence numbers)."""
            base_id = seq.id.split(".")[
                0
            ]  # Handle split sequences like 001.1 → 001
            return f"{visit_id}_{base_id}_{seq.target}"

        # Create lookup for original sequences
        original_seqs = {}
        original_base_seqs = {}

        for visit in original_calendar.visits:
            for seq in visit.sequences:
                key = create_sequence_key(visit.id, seq)
                base_key = create_base_sequence_key(visit.id, seq)
                original_seqs[key] = seq
                if base_key not in original_base_seqs:
                    original_base_seqs[base_key] = []
                original_base_seqs[base_key].append(seq)

        # Analyze processed sequences
        processed_keys = set()

        for visit in processed_calendar.visits:
            for seq in visit.sequences:
                key = create_sequence_key(visit.id, seq)
                base_key = create_base_sequence_key(visit.id, seq)
                processed_keys.add(key)

                # Check for exact match first
                if key in original_seqs:
                    orig_seq = original_seqs[key]

                    # Check for time changes (extensions/shortenings)
                    orig_duration = orig_seq.duration.sec
                    new_duration = seq.duration.sec
                    time_diff = abs((seq.start_time - orig_seq.start_time).sec)

                    if (
                        time_diff < 60
                        and abs(new_duration - orig_duration) < 60
                    ):
                        changes["unchanged"].add(key)
                    else:
                        # Times changed - determine if extended or shortened
                        if new_duration > orig_duration or time_diff > 60:
                            changes["extended"].add(key)
                        else:
                            changes["shortened"].add(key)

                # Check for base sequence match (could be split sequence)
                elif base_key in original_base_seqs:
                    # This could be a split sequence or modified sequence
                    changes["new"].add(key)  # Mark as new for now

                else:
                    # Completely new sequence
                    changes["new"].add(key)

        # Find removed sequences
        for key in original_seqs:
            if key not in processed_keys:
                changes["removed"].add(key)

        return changes

    def _plot_gantt_with_changes(
        self, calendar, ax, title, changes, is_original=False
    ):
        """Plot Gantt chart with change highlighting using robust sequence identification."""
        if not calendar.visits:
            ax.text(
                0.5,
                0.5,
                "No sequences found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(title)
            return

        # Get all sequences and organize by target
        sequences_by_target = {}
        sequence_keys = {}  # Store the keys for each sequence

        for visit in calendar.visits:
            for seq in visit.sequences:
                target = seq.target
                if target not in sequences_by_target:
                    sequences_by_target[target] = []

                sequences_by_target[target].append(seq)

                # Create the key for this sequence
                seq_key = (
                    f"{visit.id}_{seq.id}_{seq.start_time.isot}_{seq.target}"
                )
                sequence_keys[id(seq)] = seq_key  # Use object id as lookup

        # Sort sequences within each target by start time
        for target in sequences_by_target:
            sequences_by_target[target].sort(key=lambda s: s.start_time)

        # Sort targets by first sequence start time
        sorted_targets = sorted(
            sequences_by_target.keys(),
            key=lambda t: sequences_by_target[t][0].start_time,
        )

        # Base colors for targets
        base_colors = plt.cm.Set3(np.linspace(0, 1, len(sorted_targets)))
        target_colors = {
            target: base_colors[i] for i, target in enumerate(sorted_targets)
        }

        # Change highlight colors
        change_colors = {
            "extended": "lightgreen",
            "shortened": "orange",
            "unchanged": None,  # Use base color
            "new": "yellow",
            "removed": "lightgray",
        }

        # Plot each target's sequences
        y_pos = 0
        y_labels = []
        y_positions = []

        for target in sorted_targets:
            sequences = sequences_by_target[target]

            for seq in sequences:
                start_time = seq.start_time.datetime
                duration = seq.duration.sec / 3600  # Convert to hours
                seq_key = sequence_keys.get(id(seq), "unknown")

                # Determine color based on changes
                base_color = target_colors[target]

                # Check what type of change this sequence represents
                change_type = None
                for change, seq_set in changes.items():
                    if seq_key in seq_set:
                        change_type = change
                        break

                if change_type and change_colors[change_type]:
                    face_color = change_colors[change_type]
                    alpha = 0.9
                    edge_color = (
                        "red"
                        if change_type in ["extended", "shortened"]
                        else "black"
                    )
                    linewidth = (
                        2 if change_type in ["extended", "shortened"] else 1
                    )
                else:
                    face_color = base_color
                    alpha = 0.7
                    edge_color = "black"
                    linewidth = 1

                # Create rectangle for sequence
                rect = Rectangle(
                    (mdates.date2num(start_time), y_pos - 0.4),
                    duration / 24,  # Convert hours to days for matplotlib
                    0.8,
                    facecolor=face_color,
                    alpha=alpha,
                    edgecolor=edge_color,
                    linewidth=linewidth,
                )
                ax.add_patch(rect)

                # Add sequence ID label with better positioning
                mid_time = mdates.date2num(start_time) + (duration / 24 / 2)

                # Create display text (seq_id + start time for uniqueness)
                display_text = (
                    f"{seq.id}\n{seq.start_time.datetime.strftime('%H:%M')}"
                )

                # Adjust text size based on sequence duration
                if duration > 1:  # More than 1 hour
                    fontsize = 8
                    rotation = 0
                else:
                    fontsize = 6
                    rotation = 90
                    display_text = seq.id  # Shorter text for small sequences

                ax.text(
                    mid_time,
                    y_pos,
                    display_text,
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    fontweight="bold",
                    rotation=rotation,
                    bbox=dict(
                        boxstyle="round,pad=0.2", facecolor="white", alpha=0.8
                    ),
                )

            y_labels.append(f"{target} ({len(sequences)} seq)")
            y_positions.append(y_pos)
            y_pos += 1

        # Set up y-axis
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=10)
        ax.set_ylim(-0.5, len(sorted_targets) - 0.5)

        # Grid and labels
        ax.grid(True, alpha=0.3, axis="x")
        ax.set_ylabel("Targets", fontsize=12)
        ax.set_title(title, fontsize=14, pad=10)

        # Invert y-axis so first target is at top
        ax.invert_yaxis()

    def create_sequence_comparison_table(
        self, original_calendar, processed_calendar
    ):
        """Create a detailed table comparing sequences between original and processed calendars."""
        # Third-party

        comparison_data = []

        # Create lookup for processed sequences
        processed_seqs = {}
        for visit in processed_calendar.visits:
            for seq in visit.sequences:
                key = f"{visit.id}_{seq.id}_{seq.target}"
                if key not in processed_seqs:
                    processed_seqs[key] = []
                processed_seqs[key].append(
                    {
                        "visit_id": visit.id,
                        "seq": seq,
                        "start_time": seq.start_time,
                        "stop_time": seq.stop_time,
                        "duration_min": seq.duration.sec / 60,
                    }
                )

        # Compare original sequences
        for visit in original_calendar.visits:
            for seq in visit.sequences:
                key = f"{visit.id}_{seq.id}_{seq.target}"

                original_data = {
                    "visit_id": visit.id,
                    "sequence_id": seq.id,
                    "target": seq.target,
                    "original_start": seq.start_time.datetime,
                    "original_stop": seq.stop_time.datetime,
                    "original_duration_min": seq.duration.sec / 60,
                    "processed_start": None,
                    "processed_stop": None,
                    "processed_duration_min": None,
                    "change_type": "removed",
                    "time_change_min": None,
                    "duration_change_min": None,
                }

                # Look for matching processed sequence(s)
                if key in processed_seqs:
                    processed_matches = processed_seqs[key]

                    # For now, take the first match (could be improved with better matching logic)
                    if processed_matches:
                        proc_seq_info = processed_matches[0]
                        proc_seq = proc_seq_info["seq"]

                        original_data.update(
                            {
                                "processed_start": proc_seq.start_time.datetime,
                                "processed_stop": proc_seq.stop_time.datetime,
                                "processed_duration_min": proc_seq.duration.sec
                                / 60,
                                "time_change_min": (
                                    proc_seq.start_time - seq.start_time
                                ).sec
                                / 60,
                                "duration_change_min": (
                                    proc_seq.duration.sec - seq.duration.sec
                                )
                                / 60,
                            }
                        )

                        # Determine change type
                        time_change = abs(
                            (proc_seq.start_time - seq.start_time).sec
                        )
                        duration_change = (
                            proc_seq.duration.sec - seq.duration.sec
                        )

                        if time_change < 60 and abs(duration_change) < 60:
                            original_data["change_type"] = "unchanged"
                        elif duration_change > 60:
                            original_data["change_type"] = "extended"
                        elif duration_change < -60:
                            original_data["change_type"] = "shortened"
                        else:
                            original_data["change_type"] = "modified"

                comparison_data.append(original_data)

        # Add completely new sequences from processed calendar
        original_keys = set()
        for visit in original_calendar.visits:
            for seq in visit.sequences:
                original_keys.add(f"{visit.id}_{seq.id}_{seq.target}")

        for key, proc_seq_list in processed_seqs.items():
            if key not in original_keys:
                for proc_seq_info in proc_seq_list:
                    proc_seq = proc_seq_info["seq"]
                    new_data = {
                        "visit_id": proc_seq_info["visit_id"],
                        "sequence_id": proc_seq.id,
                        "target": proc_seq.target,
                        "original_start": None,
                        "original_stop": None,
                        "original_duration_min": None,
                        "processed_start": proc_seq.start_time.datetime,
                        "processed_stop": proc_seq.stop_time.datetime,
                        "processed_duration_min": proc_seq.duration.sec / 60,
                        "change_type": "new",
                        "time_change_min": None,
                        "duration_change_min": None,
                    }
                    comparison_data.append(new_data)

        df = pd.DataFrame(comparison_data)
        return df.sort_values(["visit_id", "sequence_id", "target"])

    def _format_time_axis(
        self, ax1, ax2, original_calendar, processed_calendar
    ):
        """Format the time axis appropriately."""
        # Get overall time range
        all_times = []
        for calendar in [original_calendar, processed_calendar]:
            for visit in calendar.visits:
                for seq in visit.sequences:
                    all_times.extend(
                        [seq.start_time.datetime, seq.stop_time.datetime]
                    )

        if all_times:
            min_time = min(all_times)
            max_time = max(all_times)

            # Add some padding
            time_span = max_time - min_time
            padding = time_span * 0.01
            ax1.set_xlim(min_time - padding, max_time + padding)

            # Format x-axis based on time span
            if time_span.total_seconds() < 86400:  # Less than 1 day
                ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax2.xaxis.set_major_locator(mdates.HourLocator(interval=3))
                ax2.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
            elif time_span.total_seconds() < 86400 * 7:  # Less than 1 week
                ax2.xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d %H:%M")
                )
                ax2.xaxis.set_major_locator(mdates.HourLocator(interval=12))
            else:
                ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
                ax2.xaxis.set_major_locator(mdates.DayLocator())

        plt.xticks(rotation=45)

    def plot_priority_statistics(self, calendar, figsize=(12, 8)):
        """Create priority statistics visualization."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)

        # Collect priority data
        priority_data = {}
        for visit in calendar.visits:
            for seq in visit.sequences:
                p = seq.priority
                if p not in priority_data:
                    priority_data[p] = {
                        "count": 0,
                        "duration_minutes": 0,
                        "sequences": [],
                    }
                priority_data[p]["count"] += 1
                priority_data[p]["duration_minutes"] += seq.duration.sec / 60.0
                priority_data[p]["sequences"].append(seq)

        priorities = sorted(priority_data.keys())
        priority_colors = self._get_priority_colors(priorities)
        colors = [priority_colors[p] for p in priorities]

        # 1. Sequence count by priority
        counts = [priority_data[p]["count"] for p in priorities]
        ax1.bar(priorities, counts, color=colors, alpha=0.8)
        ax1.set_xlabel("Priority Level")
        ax1.set_ylabel("Number of Sequences")
        ax1.set_title("Sequence Count by Priority")
        ax1.grid(True, alpha=0.3)

        # 2. Total duration by priority
        durations = [
            priority_data[p]["duration_minutes"] / 60 for p in priorities
        ]  # Convert to hours
        ax2.bar(priorities, durations, color=colors, alpha=0.8)
        ax2.set_xlabel("Priority Level")
        ax2.set_ylabel("Total Duration (hours)")
        ax2.set_title("Total Duration by Priority")
        ax2.grid(True, alpha=0.3)

        # 3. Priority distribution pie chart
        ax3.pie(
            counts,
            labels=[f"P{p}" for p in priorities],
            colors=colors,
            autopct="%1.1f%%",
            startangle=90,
        )
        ax3.set_title("Priority Distribution")

        # 4. Average sequence duration by priority
        avg_durations = [
            priority_data[p]["duration_minutes"] / priority_data[p]["count"]
            for p in priorities
        ]
        ax4.bar(priorities, avg_durations, color=colors, alpha=0.8)
        ax4.set_xlabel("Priority Level")
        ax4.set_ylabel("Average Duration (minutes)")
        ax4.set_title("Average Sequence Duration by Priority")
        ax4.grid(True, alpha=0.3)

        fig.suptitle("Priority Analysis Summary", fontsize=16)
        fig.set_constrained_layout(True)

        return fig

    def _format_time_axis_single(self, ax, calendar):
        """Format time axis for single calendar plot with size limits."""
        # Get time range
        all_times = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_times.extend(
                    [seq.start_time.datetime, seq.stop_time.datetime]
                )

        if all_times:
            min_time = min(all_times)
            max_time = max(all_times)

            # Add padding
            time_span = max_time - min_time
            padding = time_span * 0.01
            ax.set_xlim(min_time - padding, max_time + padding)

            # Format x-axis with reasonable intervals to prevent huge images
            total_seconds = time_span.total_seconds()

            if total_seconds < 3600:  # Less than 1 hour
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=10))
            elif total_seconds < 86400:  # Less than 1 day
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax.xaxis.set_major_locator(
                    mdates.HourLocator(
                        interval=max(1, int(total_seconds / 3600 / 10))
                    )
                )
            elif total_seconds < 86400 * 7:  # Less than 1 week
                ax.xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d %H:%M")
                )
                ax.xaxis.set_major_locator(
                    mdates.HourLocator(
                        interval=max(6, int(total_seconds / 3600 / 20))
                    )
                )
            elif total_seconds < 86400 * 30:  # Less than 1 month
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
                ax.xaxis.set_major_locator(
                    mdates.DayLocator(
                        interval=max(1, int(total_seconds / 86400 / 15))
                    )
                )
            else:  # More than 1 month
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
                ax.xaxis.set_major_locator(
                    mdates.DayLocator(
                        interval=max(7, int(total_seconds / 86400 / 10))
                    )
                )

        plt.xticks(rotation=45)

    # Alternative: Create a time-windowed version
    def plot_gantt_timeline_by_priority_windowed(
        self,
        calendar,
        figsize=(20, 10),
        show_sequence_labels=False,
        title="Schedule by Priority",
        window_days=None,
    ):
        """Create priority Gantt chart with optional time windowing."""

        # If window_days is specified, create multiple plots
        if window_days is not None:
            all_times = []
            for visit in calendar.visits:
                for seq in visit.sequences:
                    all_times.extend(
                        [seq.start_time.datetime, seq.stop_time.datetime]
                    )

            if all_times:
                total_span = max(all_times) - min(all_times)
                total_days = total_span.total_seconds() / 86400

                if total_days > window_days:
                    print(
                        f"Large time span ({total_days:.1f} days). Creating multiple plots..."
                    )
                    return self._create_windowed_plots(
                        calendar,
                        figsize,
                        show_sequence_labels,
                        title,
                        window_days,
                    )

        # Single plot version
        return self.plot_gantt_timeline_by_priority(
            calendar, figsize, show_sequence_labels, title
        )

    def _create_windowed_plots(
        self, calendar, figsize, show_sequence_labels, title, window_days
    ):
        """Create multiple plots for large time spans."""
        # Standard library

        all_times = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_times.extend(
                    [seq.start_time.datetime, seq.stop_time.datetime]
                )

        start_date = min(all_times)
        end_date = max(all_times)

        figures = []
        current_date = start_date
        window_num = 1

        while current_date < end_date:
            window_end = min(
                current_date + timedelta(days=window_days), end_date
            )

            # Create filtered calendar for this window
            windowed_calendar = self._filter_calendar_by_time(
                calendar, current_date, window_end
            )

            if (
                windowed_calendar.visits
            ):  # Only create plot if there are sequences
                window_title = f"{title} - Window {window_num} ({current_date.strftime('%m/%d')} to {window_end.strftime('%m/%d')})"

                fig = self.plot_gantt_timeline_by_priority(
                    windowed_calendar,
                    figsize=figsize,
                    show_sequence_labels=show_sequence_labels,
                    title=window_title,
                )
                figures.append(fig)
                window_num += 1

            current_date = window_end

        return figures

    def _filter_calendar_by_time(self, calendar, start_date, end_date):
        """Filter calendar to only include sequences within time window."""
        # Standard library

        filtered_visits = []

        for visit in calendar.visits:
            filtered_sequences = []

            for seq in visit.sequences:
                seq_start = seq.start_time.datetime
                seq_end = seq.stop_time.datetime

                # Include sequence if it overlaps with the window
                if seq_start < end_date and seq_end > start_date:
                    filtered_sequences.append(seq)

            if filtered_sequences:
                filtered_visit = Visit(
                    id=visit.id, sequences=filtered_sequences
                )
                filtered_visits.append(filtered_visit)

        return ScienceCalendar(
            metadata=calendar.metadata, visits=filtered_visits
        )

    def plot_gantt_timeline_by_priority(
        self,
        calendar,
        figsize=(16, 8),
        show_sequence_labels=False,
        title="Schedule by Priority",
    ):
        """Create a Gantt chart with proper size limits to prevent oversized images."""

        # Check time span and warn if too large
        all_times = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_times.extend(
                    [seq.start_time.datetime, seq.stop_time.datetime]
                )

        if not all_times:
            print("No sequences found in calendar")
            return plt.figure(figsize=figsize)

        time_span = max(all_times) - min(all_times)
        time_span_hours = time_span.total_seconds() / 3600

        # If time span is too large, automatically window it
        if time_span_hours > 168:  # More than 1 week
            print(
                f"Time span too large ({time_span_hours/24:.1f} days). Creating windowed plots..."
            )
            return self._create_windowed_priority_plots(
                calendar, figsize, show_sequence_labels, title
            )

        # Create figure with controlled size
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Set reasonable DPI limit
        fig.set_dpi(100)  # Lower DPI to prevent huge images

        try:
            self._plot_gantt_by_priority(
                calendar, ax, title, show_sequence_labels=show_sequence_labels
            )
            self._format_time_axis_safe(ax, calendar)
            self._add_priority_legend(fig, calendar)

            # fig.suptitle(title, fontsize=14, y=0.95)
            fig.set_constrained_layout(True)

        except Exception as e:
            print(f"Error creating plot: {e}")
            # Return a simple error plot
            ax.text(
                0.5,
                0.5,
                f"Plot too large to render\nTime span: {time_span_hours/24:.1f} days",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
            )
            ax.set_title(title)

        return fig

    def _format_time_axis_safe(self, ax, calendar):
        """Format time axis with strict size limits."""
        all_times = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_times.extend(
                    [seq.start_time.datetime, seq.stop_time.datetime]
                )

        if not all_times:
            return

        min_time = min(all_times)
        max_time = max(all_times)
        time_span = max_time - min_time

        # Set limits with minimal padding
        padding = time_span * 0.001  # Very small padding
        ax.set_xlim(min_time - padding, max_time + padding)

        # Use very conservative tick spacing
        total_hours = time_span.total_seconds() / 3600

        if total_hours < 6:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        elif total_hours < 48:
            ax.xaxis.set_major_locator(
                mdates.HourLocator(interval=max(2, int(total_hours / 10)))
            )
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d %H:%M"))
        else:
            # For longer periods, use day intervals
            interval = max(1, int(total_hours / 24 / 10))
            ax.xaxis.set_major_locator(mdates.DayLocator(interval=interval))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))

        # Limit number of ticks to prevent overcrowding
        # ax.locator_params(axis="x", nbins=10)  # Removed: not supported for date locators
        plt.xticks(rotation=45)

    def _create_windowed_priority_plots(
        self, calendar, figsize, show_sequence_labels, title
    ):
        """Create multiple smaller plots for large time spans."""
        # Standard library

        all_times = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_times.extend(
                    [seq.start_time.datetime, seq.stop_time.datetime]
                )

        start_date = min(all_times)
        end_date = max(all_times)
        total_days = (end_date - start_date).total_seconds() / 86400

        # Create 3-day windows
        window_days = 3
        figures = []
        current_date = start_date
        window_num = 1

        print(
            f"Creating {int(np.ceil(total_days/window_days))} plots for {total_days:.1f} day span..."
        )

        while current_date < end_date:
            window_end = min(
                current_date + timedelta(days=window_days), end_date
            )

            # Filter calendar for this window
            windowed_sequences = []
            for visit in calendar.visits:
                windowed_visit_sequences = []
                for seq in visit.sequences:
                    if (
                        seq.start_time.datetime < window_end
                        and seq.stop_time.datetime > current_date
                    ):
                        windowed_visit_sequences.append(seq)

                if windowed_visit_sequences:
                    windowed_visit = Visit(
                        id=visit.id, sequences=windowed_visit_sequences
                    )
                    windowed_sequences.append(windowed_visit)

            if windowed_sequences:
                # Create windowed calendar
                windowed_calendar = ScienceCalendar(
                    metadata=calendar.metadata, visits=windowed_sequences
                )

                window_title = f"{title} - Days {window_num}-{window_num+2}"

                # Create smaller figure for window
                fig, ax = plt.subplots(1, 1, figsize=(12, 6))
                fig.set_dpi(100)

                try:
                    self._plot_gantt_by_priority(
                        windowed_calendar,
                        ax,
                        window_title,
                        show_sequence_labels=show_sequence_labels,
                    )
                    self._format_time_axis_safe(ax, windowed_calendar)
                    self._add_priority_legend(fig, windowed_calendar)

                    # fig.suptitle(window_title, fontsize=12, y=0.95)
                    fig.set_constrained_layout(True)

                    figures.append(fig)

                except Exception as e:
                    print(f"Error creating window {window_num}: {e}")
                    plt.close(fig)

            current_date = window_end
            window_num += 3

        return figures

    # Alternative: Create a summary plot instead of detailed timeline
    def plot_priority_summary_only(self, calendar, figsize=(12, 8)):
        """Create just the priority statistics without the timeline."""
        return self.plot_priority_statistics(calendar, figsize)

    def plot_simplified_timeline(self, calendar, figsize=(16, 6)):
        """Create a very simplified timeline that won't have size issues."""
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Get all sequences and sort by start time
        all_sequences = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_sequences.append((visit.id, seq))

        all_sequences.sort(key=lambda x: x[1].start_time)

        if not all_sequences:
            ax.text(
                0.5,
                0.5,
                "No sequences found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            return fig

        # Plot as simple horizontal bars
        priority_colors = self._get_priority_colors([1, 2, 3, 4, 5, 6, 7, 8])

        for i, (visit_id, seq) in enumerate(
            all_sequences[:50]
        ):  # Limit to first 50 sequences
            start_num = mdates.date2num(seq.start_time.datetime)
            duration_days = seq.duration.sec / 86400
            color = priority_colors.get(seq.priority, "lightgray")

            ax.barh(
                i,
                duration_days,
                left=start_num,
                color=color,
                alpha=0.8,
                height=0.8,
            )

            # Add simple label
            ax.text(
                start_num + duration_days / 2,
                i,
                f"{seq.id}(P{seq.priority})",
                ha="center",
                va="center",
                fontsize=6,
            )

        ax.set_ylim(-0.5, min(len(all_sequences), 50) - 0.5)
        ax.set_ylabel("Sequences")
        ax.set_title("Simplified Schedule Timeline (First 50 sequences)")

        # Simple time formatting
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        plt.xticks(rotation=45)

        fig.set_constrained_layout(True)
        return fig

    def plot_simple_priority_timeline(
        self,
        calendar,
        figsize=(20, 12),
        show_sequence_labels=True,
        title="Schedule by Priority",
    ):
        """Create a simple timeline colored by priority, similar to the existing change plots."""
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        # Get all unique priorities and create simple color scheme
        all_priorities = set()
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_priorities.add(seq.priority)

        # Simple priority colors
        priority_colors = {
            1: "darkred",
            2: "red",
            3: "orange",
            4: "gold",
            5: "yellow",
            6: "lightgreen",
            7: "lightblue",
            8: "lightgray",
        }

        # For any priorities not in the map, use shades of gray
        sorted_priorities = sorted(all_priorities)
        for i, p in enumerate(sorted_priorities):
            if p not in priority_colors:
                gray_val = (
                    0.3 + (i * 0.1) % 0.5
                )  # Generate grays between 0.3-0.8
                priority_colors[p] = (gray_val, gray_val, gray_val)

        # Get all sequences and organize by target (same as your existing plot)
        sequences_by_target = {}
        for visit in calendar.visits:
            for seq in visit.sequences:
                target = seq.target
                if target not in sequences_by_target:
                    sequences_by_target[target] = []
                sequences_by_target[target].append(seq)

        # Sort sequences within each target by start time
        for target in sequences_by_target:
            sequences_by_target[target].sort(key=lambda s: s.start_time)

        # Sort targets by their first sequence start time
        sorted_targets = sorted(
            sequences_by_target.keys(),
            key=lambda t: sequences_by_target[t][0].start_time,
        )

        # Plot each target's sequences (same layout as existing plot)
        y_pos = 0
        y_labels = []
        y_positions = []

        for target in sorted_targets:
            sequences = sequences_by_target[target]

            for seq in sequences:
                start_time = seq.start_time.datetime
                duration = seq.duration.sec / 3600  # Convert to hours

                # Get color based on priority
                color = priority_colors.get(seq.priority, "lightgray")

                # Create rectangle (same as existing code)
                rect = Rectangle(
                    (mdates.date2num(start_time), y_pos - 0.4),
                    duration / 24,  # Convert hours to days for matplotlib
                    0.8,
                    facecolor=color,
                    alpha=0.8,
                    edgecolor="black",
                    linewidth=1,
                )
                ax.add_patch(rect)

                # Add sequence ID label (same as existing)
                if show_sequence_labels:
                    mid_time = mdates.date2num(start_time) + (
                        duration / 24 / 2
                    )

                    if duration > 1:  # More than 1 hour
                        fontsize = 8
                        rotation = 0
                    else:
                        fontsize = 6
                        rotation = 90

                    ax.text(
                        mid_time,
                        y_pos,
                        seq.id,
                        ha="center",
                        va="center",
                        fontsize=fontsize,
                        fontweight="bold",
                        rotation=rotation,
                        bbox=dict(
                            boxstyle="round,pad=0.2",
                            facecolor="white",
                            alpha=0.8,
                        ),
                    )

            # Target label with sequence count (same as existing)
            y_labels.append(f"{target} ({len(sequences)})")
            y_positions.append(y_pos)
            y_pos += 1

        # Set up y-axis (same as existing)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=10)
        ax.set_ylim(-0.5, len(sorted_targets) - 0.5)

        # Add visit labels on the right (same as existing)
        visit_positions = {}
        for visit in calendar.visits:
            visit_sequences = visit.sequences
            if visit_sequences:
                # Find which y positions this visit occupies
                visit_targets = set(seq.target for seq in visit_sequences)
                visit_y_positions = [
                    i
                    for i, target in enumerate(sorted_targets)
                    if target in visit_targets
                ]
                if visit_y_positions:
                    avg_y = sum(visit_y_positions) / len(visit_y_positions)
                    visit_positions[visit.id] = avg_y

        for i, (visit_id, avg_y) in enumerate(visit_positions.items()):
            norm_y = (
                avg_y / (len(sorted_targets) - 1)
                if len(sorted_targets) > 1
                else 0.5
            )
            visit_info = f"Visit {visit_id}\n({len([s for s in calendar.visits if s.id == visit_id][0].sequences)} seq)"
            ax.text(
                1.02,
                norm_y,
                visit_info,
                transform=ax.transAxes,
                ha="left",
                va="center",
                fontsize=10,
                fontweight="bold",
                bbox=dict(
                    boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7
                ),
            )

        # Format time axis (same as existing)
        all_times = []
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_times.extend(
                    [seq.start_time.datetime, seq.stop_time.datetime]
                )

        if all_times:
            min_time = min(all_times)
            max_time = max(all_times)
            time_span = max_time - min_time
            padding = time_span * 0.01
            ax.set_xlim(min_time - padding, max_time + padding)

            if time_span.total_seconds() < 86400:  # Less than 1 day
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            elif time_span.total_seconds() < 86400 * 7:  # Less than 1 week
                ax.xaxis.set_major_formatter(
                    mdates.DateFormatter("%m/%d %H:%M")
                )
                ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
            else:
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
                ax.xaxis.set_major_locator(mdates.DayLocator())

        plt.xticks(rotation=45)

        # Add simple priority legend
        # Third-party
        from matplotlib.patches import Patch

        legend_elements = []
        for priority in sorted(all_priorities):
            color = priority_colors.get(priority, "lightgray")
            legend_elements.append(
                Patch(facecolor=color, alpha=0.8, label=f"Priority {priority}")
            )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            bbox_to_anchor=(0.98, 0.98),
        )

        # Grid and labels (same as existing)
        ax.grid(True, alpha=0.3, axis="x")
        ax.set_ylabel("Targets (grouped by visit)", fontsize=12)
        ax.set_title(title, fontsize=14, pad=10)
        ax.invert_yaxis()

        fig.set_constrained_layout(True)
        return fig

    def plot_gantt_with_visibility(
        self,
        calendar,
        data_dir: Optional[Path] = None,
        figsize=(14, 8),
        show_sequence_labels=False,
        title="Schedule by Priority — Visibility Overlay",
    ):
        """Gantt chart colored by priority with non-visible minutes
        overlaid in red.

        Queries the scheduler's ``Visibility`` object per-sequence to
        compute minute-by-minute visibility and draws red blocks over
        any minutes that fail a keepout constraint.

        Parameters
        ----------
        calendar : ScienceCalendar
            Processed calendar to plot.
        figsize : tuple
            Figure size (width, height) in inches.
        show_sequence_labels : bool
            Annotate each bar with the sequence ID.
        title : str
            Plot title.

        Returns
        -------
        matplotlib.figure.Figure
        """
        if data_dir is not None:
            return self._plot_gantt_with_visibility_from_data_dir(
                calendar=calendar,
                data_dir=Path(data_dir),
                figsize=figsize,
                show_sequence_labels=show_sequence_labels,
                title=title,
            )

        from astropy import units as u
        from astropy.coordinates import SkyCoord
        from matplotlib.patches import Patch

        scheduler = self.scheduler
        priority_colors = self._get_priority_colors(
            [1, 2, 3, 4, 5, 6, 7, 8]
        )

        fig, ax = plt.subplots(figsize=figsize)

        # Collect all sequences with per-minute visibility
        rows = []
        for visit in sorted(
            calendar.visits,
            key=lambda v: (
                v.sequences[0].start_time
                if v.sequences
                else Time("2000-01-01")
            ),
        ):
            visit_rolls = scheduler._computed_target_rolls.get(
                visit.id, {}
            )
            for seq in visit.sequences:
                n_mins = int(np.rint(seq.duration.sec / 60.0))
                if n_mins <= 0:
                    continue
                coord = SkyCoord(
                    seq.ra, seq.dec, frame="icrs", unit="deg"
                )
                deltas = np.arange(n_mins) * u.min
                times = seq.start_time + deltas

                roll = visit_rolls.get(seq.target)
                if (
                    scheduler._roll_sweep_enabled
                    and roll is not None
                ):
                    vis = scheduler.visibility.get_visibility(
                        coord, times, roll=roll * u.deg
                    )
                else:
                    vis = scheduler.visibility.get_visibility(
                        coord, times
                    )
                rows.append((visit.id, seq, np.asarray(vis)))

        if not rows:
            ax.text(
                0.5, 0.5, "No sequences",
                ha="center", transform=ax.transAxes,
            )
            return fig

        # Build y-axis: one row per (visit, target)
        seen = {}
        y_labels = []
        for vid, seq, _ in rows:
            key = (vid, seq.target)
            if key not in seen:
                seen[key] = len(y_labels)
                y_labels.append(f"{vid} / {seq.target}")

        # Track x-extent for axis limits
        x_min = float("inf")
        x_max = float("-inf")
        non_vis_total = 0
        total_mins = 0

        for vid, seq, vis_arr in rows:
            y = seen[(vid, seq.target)]
            start_num = float(mdates.date2num(seq.start_time.datetime))
            stop_num = float(mdates.date2num(seq.stop_time.datetime))
            dur_days = stop_num - start_num

            x_min = min(x_min, start_num)
            x_max = max(x_max, stop_num)

            # Background bar: priority color
            color = priority_colors.get(seq.priority, "lightgray")
            rect = Rectangle(
                (start_num, y - 0.35), dur_days, 0.7,
                facecolor=color, edgecolor="none",
                alpha=1.0, linewidth=0,
            )
            ax.add_patch(rect)

            # Red overlay for non-visible minutes
            n_mins = len(vis_arr)
            total_mins += n_mins
            if n_mins > 0 and not np.all(vis_arr):
                min_dur = dur_days / n_mins
                in_block = False
                block_start = 0
                for i in range(n_mins + 1):
                    if i < n_mins and not vis_arr[i]:
                        if not in_block:
                            block_start = i
                            in_block = True
                    else:
                        if in_block:
                            x0 = start_num + block_start * min_dur
                            width = (i - block_start) * min_dur
                            r = Rectangle(
                                (x0, y - 0.35), width, 0.7,
                                facecolor="black", edgecolor="none",
                                alpha=1.00, linewidth=0, zorder=1000
                            )
                            ax.add_patch(r)
                            non_vis_total += i - block_start
                            in_block = False

            if show_sequence_labels:
                mid = start_num + dur_days / 2
                ax.text(
                    mid, y, seq.id, ha="center", va="center",
                    fontsize=5, clip_on=True,
                )

        # Axis limits  (patches don't trigger autoscale)
        padding = (x_max - x_min) * 0.005
        ax.set_xlim(x_min - padding, x_max + padding)

        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=7)
        ax.set_ylim(-0.5, len(y_labels) - 0.5)
        ax.invert_yaxis()

        # Time-axis formatting
        self._format_time_axis_safe(ax, calendar)

        vis_pct = (
            100.0 * (total_mins - non_vis_total) / total_mins
            if total_mins > 0
            else 100.0
        )
        ax.set_title(
            f"{title}\n"
            f"({non_vis_total} non-visible min / "
            f"{total_mins} total — {vis_pct:.1f}% visible)",
            fontsize=12, pad=10,
        )
        ax.set_xlabel("Time (UTC)")
        ax.grid(True, axis="x", alpha=0.3)

        # Legend
        legend_items = [
            Patch(facecolor="red", alpha=0.75, label="Non-visible"),
        ]
        used_priorities = sorted(
            set(s.priority for _, s, _ in rows)
        )
        for p in used_priorities:
            c = priority_colors.get(p, "silver")
            legend_items.append(
                Patch(facecolor=c, label=f"Priority {p}")
            )
        ax.legend(handles=legend_items, loc="upper right", fontsize=7)

        fig.tight_layout()
        return fig

    def _plot_gantt_with_visibility_from_data_dir(
        self,
        calendar,
        data_dir: Path,
        figsize=(14, 8),
        show_sequence_labels=False,
        title="Schedule by Priority — Visibility Overlay",
    ):
        fig, ax = plt.subplots(figsize=figsize)
        priority_colors = self._get_priority_colors([1, 2, 3, 4, 5, 6, 7, 8])

        exoplanet_csv = data_dir / "exoplanet_targets.csv"
        planet_to_star: dict[str, str] = {}
        if exoplanet_csv.exists():
            try:
                exoplanet_df = pd.read_csv(exoplanet_csv, usecols=["Planet Name", "Star Name"])
                for _, row in exoplanet_df.iterrows():
                    planet = str(row.get("Planet Name", "")).strip()
                    star = str(row.get("Star Name", "")).strip()
                    if planet and star:
                        planet_to_star[planet] = star
            except Exception:
                pass

        visibility_cache: dict[str, Optional[pd.DataFrame]] = {}

        def _load_visibility(target_name: str) -> Optional[pd.DataFrame]:
            if target_name in visibility_cache:
                return visibility_cache[target_name]

            candidates = [
                data_dir / "targets" / target_name / f"Visibility for {target_name}.parquet",
                data_dir / "aux_targets" / target_name / f"Visibility for {target_name}.parquet",
            ]
            star_name = planet_to_star.get(target_name)
            if star_name:
                candidates.extend(
                    [
                        data_dir / "targets" / star_name / f"Visibility for {star_name}.parquet",
                        data_dir / "aux_targets" / star_name / f"Visibility for {star_name}.parquet",
                    ]
                )

            for candidate in candidates:
                if candidate.exists():
                    try:
                        df = pd.read_parquet(candidate, columns=["Time_UTC", "Visible"])
                    except Exception:
                        try:
                            df = pd.read_parquet(candidate, columns=["Time(MJD_UTC)", "Visible"])
                        except Exception:
                            continue
                    if "Time_UTC" in df.columns:
                        times = pd.to_datetime(df["Time_UTC"], errors="coerce")
                    else:
                        times = Time(
                            df["Time(MJD_UTC)"].to_numpy(dtype=float),
                            format="mjd",
                            scale="utc",
                        ).to_datetime()
                        times = pd.to_datetime(times)
                    index = pd.DatetimeIndex(times)
                    if getattr(index, "tz", None) is not None:
                        index = index.tz_localize(None)
                    # GMAT-derived timestamps can be a few microseconds off each
                    # exact minute boundary. Round to the nearest minute so they
                    # align with the minute-grid XML sequences.
                    index = index.round("min")
                    prepared = pd.DataFrame(
                        {"Visible": (df["Visible"].to_numpy(dtype=float) > 0.5)},
                        index=index,
                    )
                    prepared = prepared.groupby(level=0)["Visible"].max().to_frame()
                    visibility_cache[target_name] = prepared
                    return prepared

            visibility_cache[target_name] = None
            return None

        rows = []
        missing_targets: set[str] = set()
        for visit in sorted(
            calendar.visits,
            key=lambda v: (v.sequences[0].start_time if v.sequences else Time("2000-01-01")),
        ):
            for seq in visit.sequences:
                n_mins = int(np.rint(seq.duration.sec / 60.0))
                if n_mins <= 0:
                    continue
                visibility_df = _load_visibility(seq.target)
                if visibility_df is None:
                    missing_targets.add(seq.target)
                    vis_arr = np.ones(n_mins, dtype=bool)
                else:
                    expected = pd.date_range(
                        seq.start_time.datetime,
                        periods=n_mins,
                        freq="min",
                    )
                    vis_arr = (
                        visibility_df["Visible"]
                        .reindex(expected, fill_value=False)
                        .to_numpy(dtype=bool)
                    )
                rows.append((visit.id, seq, vis_arr))

        if not rows:
            ax.text(0.5, 0.5, "No sequences", ha="center", transform=ax.transAxes)
            return fig

        seen = {}
        y_labels = []
        for vid, seq, _ in rows:
            key = (vid, seq.target)
            if key not in seen:
                seen[key] = len(y_labels)
                y_labels.append(f"{vid} / {seq.target}")

        x_min = float("inf")
        x_max = float("-inf")
        non_vis_total = 0
        total_mins = 0

        for vid, seq, vis_arr in rows:
            y = seen[(vid, seq.target)]
            start_num = float(mdates.date2num(seq.start_time.datetime))
            stop_num = float(mdates.date2num(seq.stop_time.datetime))
            dur_days = stop_num - start_num

            x_min = min(x_min, start_num)
            x_max = max(x_max, stop_num)

            color = priority_colors.get(seq.priority, "lightgray")
            rect = Rectangle(
                (start_num, y - 0.35),
                dur_days,
                0.7,
                facecolor=color,
                edgecolor="none",
                alpha=1.0,
                linewidth=0,
            )
            ax.add_patch(rect)

            n_mins = len(vis_arr)
            total_mins += n_mins
            if n_mins > 0 and not np.all(vis_arr):
                min_dur = dur_days / n_mins
                in_block = False
                block_start = 0
                for i in range(n_mins + 1):
                    if i < n_mins and not vis_arr[i]:
                        if not in_block:
                            block_start = i
                            in_block = True
                    elif in_block:
                        x0 = start_num + block_start * min_dur
                        width = (i - block_start) * min_dur
                        ax.add_patch(
                            Rectangle(
                                (x0, y - 0.35),
                                width,
                                0.7,
                                facecolor="black",
                                edgecolor="none",
                                alpha=1.0,
                                linewidth=0,
                                zorder=1000,
                            )
                        )
                        non_vis_total += i - block_start
                        in_block = False

            if show_sequence_labels:
                mid = start_num + dur_days / 2
                ax.text(mid, y, seq.id, ha="center", va="center", fontsize=5, clip_on=True)

        padding = (x_max - x_min) * 0.005 if x_max > x_min else 0.01
        ax.set_xlim(x_min - padding, x_max + padding)
        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=7)
        ax.set_ylim(-0.5, len(y_labels) - 0.5)
        ax.invert_yaxis()

        self._format_time_axis_safe(ax, calendar)

        vis_pct = (
            100.0 * (total_mins - non_vis_total) / total_mins
            if total_mins > 0
            else 100.0
        )
        title_suffix = ""
        if missing_targets:
            title_suffix = f"\n{len(missing_targets)} target(s) missing visibility parquet"
        ax.set_title(
            f"{title}\n"
            f"({non_vis_total} non-visible min / {total_mins} total — {vis_pct:.1f}% visible)"
            f"{title_suffix}",
            fontsize=12,
            pad=10,
        )
        ax.set_xlabel("Time (UTC)")
        ax.grid(True, axis="x", alpha=0.3)

        from matplotlib.patches import Patch

        legend_items = [Patch(facecolor="black", alpha=1.0, label="Non-visible")]
        used_priorities = sorted(set(s.priority for _, s, _ in rows))
        for p in used_priorities:
            c = priority_colors.get(p, "silver")
            legend_items.append(Patch(facecolor=c, label=f"Priority {p}"))
        ax.legend(handles=legend_items, loc="upper right", fontsize=7)

        fig.tight_layout()
        return fig

    def _get_priority_colors(self, priorities):
        """Generate better color scheme for priorities with more distinct colors."""
        priority_colors = {}

        # Define more distinct colors for priorities
        color_map = {
            0: "lightgray",  # Priority 0 - light gray
            1: "crimson",  # Highest priority - bright red (more distinct from red)
            2: "darkorange",  # High priority - dark orange (more distinct from red)
            3: "gold",  # Medium-high priority - gold
            4: "yellow",  # Medium priority - yellow
            5: "lightgreen",  # Medium-low priority - light green
            6: "lightblue",  # Low priority - light blue
            7: "plum",  # Lower priority - light purple
            8: "lightgray",  # Lowest priority - light gray
        }

        # Use predefined colors for known priorities, generate others
        for priority in priorities:
            if priority in color_map:
                priority_colors[priority] = color_map[priority]
            else:
                # For priorities > 8, use a grayscale that's definitely in 0-1 range
                gray_value = max(0.3, min(0.8, 0.8 - (priority - 8) * 0.05))
                priority_colors[priority] = (
                    gray_value,
                    gray_value,
                    gray_value,
                )

        return priority_colors

    def _plot_gantt_by_priority(
        self, calendar, ax, title, show_sequence_labels=True
    ):
        """Plot Gantt chart with solid colors and no borders."""
        if not calendar.visits:
            ax.text(
                0.5,
                0.5,
                "No sequences found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_title(title)
            return

        # Sort visits by start time
        sorted_visits = sorted(
            calendar.visits,
            key=lambda v: v.start_time if v.start_time else Time("2000-01-01"),
        )

        # Get all unique priorities for color mapping
        all_priorities = set()
        for visit in sorted_visits:
            for seq in visit.sequences:
                all_priorities.add(seq.priority)

        # Create priority color scheme
        priority_colors = self._get_priority_colors(sorted(all_priorities))

        y_pos = 0
        y_labels = []
        y_positions = []
        visit_boundaries = []

        for visit_idx, visit in enumerate(sorted_visits):
            visit_start_y = y_pos

            # Group sequences by target within this visit
            sequences_by_target = {}
            for seq in visit.sequences:
                target = seq.target
                if target not in sequences_by_target:
                    sequences_by_target[target] = []
                sequences_by_target[target].append(seq)

            # Sort sequences within each target by start time
            for target in sequences_by_target:
                sequences_by_target[target].sort(key=lambda s: s.start_time)

            # Sort targets by their first sequence start time
            sorted_targets = sorted(
                sequences_by_target.keys(),
                key=lambda t: sequences_by_target[t][0].start_time,
            )

            # Plot each target within this visit
            for target in sorted_targets:
                sequences = sequences_by_target[target]

                # Plot all sequences for this target
                for seq_idx, seq in enumerate(sequences):
                    start_time = seq.start_time.datetime
                    duration = seq.duration.sec / 3600  # Convert to hours

                    # Get priority-based color
                    priority_color = priority_colors.get(
                        seq.priority, "lightgray"
                    )

                    # Create rectangle with NO BORDER for solid color
                    # Ensure start_time is a scalar datetime
                    start_num = float(mdates.date2num(start_time))
                    rect = Rectangle(
                        (start_num, y_pos - 0.35),
                        duration / 24,  # Convert hours to days for matplotlib
                        0.7,
                        facecolor=priority_color,
                        alpha=1.0,  # Full opacity for solid colors
                        edgecolor="none",  # NO BORDER
                        linewidth=0,  # NO BORDER
                    )
                    ax.add_patch(rect)

                    # Add sequence label with priority info
                    if (
                        show_sequence_labels and duration > 0.5
                    ):  # More than 30 minutes
                        mid_time = mdates.date2num(start_time) + (
                            duration / 24 / 2
                        )
                        label_text = f"{seq.id}\nP{seq.priority}"
                        ax.text(
                            mid_time,
                            y_pos,
                            label_text,
                            ha="center",
                            va="center",
                            fontsize=6,
                            fontweight="bold",
                            bbox=dict(
                                boxstyle="round,pad=0.1",
                                facecolor="white",
                                alpha=0.9,
                            ),
                        )

                # Create label for this target within the visit with priority info
                priorities_in_target = sorted(
                    set(seq.priority for seq in sequences)
                )
                priority_str = f"P{','.join(map(str, priorities_in_target))}"
                y_labels.append(
                    f"  {target} ({len(sequences)} seq, {priority_str})"
                )
                y_positions.append(y_pos)
                y_pos += 1

            visit_end_y = y_pos - 1
            visit_boundaries.append((visit_start_y, visit_end_y))

            # Add visit separator line
            if visit_idx < len(sorted_visits) - 1:
                ax.axhline(
                    y=y_pos - 0.5, color="black", linewidth=2, alpha=0.7
                )

            y_pos += 0.2  # Small gap between visits

        # Set up y-axis
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=9)
        ax.set_ylim(-0.5, y_pos - 0.7)

        # Add visit labels with priority summary
        for visit_idx, (visit, (start_y, end_y)) in enumerate(
            zip(sorted_visits, visit_boundaries)
        ):
            mid_y = (start_y + end_y) / 2
            norm_y = mid_y / (y_pos - 0.7) if (y_pos - 0.7) > 0 else 0.5

            # Calculate priority distribution for this visit
            visit_priorities = [seq.priority for seq in visit.sequences]
            priority_counts = {}
            for p in visit_priorities:
                priority_counts[p] = priority_counts.get(p, 0) + 1

            priority_summary = ", ".join(
                [f"P{p}:{c}" for p, c in sorted(priority_counts.items())]
            )
            visit_info = f"Visit {visit.id}\n({len(visit.sequences)} seq)\n{priority_summary}"

            # Position labels to avoid overlap
            x_offset = 1.05 + (visit_idx % 2) * 0.12
            ax.text(
                x_offset,
                norm_y,
                visit_info,
                transform=ax.transAxes,
                ha="left",
                va="center",
                fontsize=8,
                fontweight="bold",
                bbox=dict(
                    boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.7
                ),
            )

        # Grid and labels
        ax.grid(True, alpha=0.3, axis="x")
        ax.set_ylabel("Targets (grouped by visit)", fontsize=12)
        ax.set_title(title, fontsize=14, pad=10)

    def _add_priority_legend(self, fig, calendar):
        """Add legend with better priority labels."""
        # Third-party
        from matplotlib.patches import Patch

        # Get all unique priorities
        all_priorities = set()
        for visit in calendar.visits:
            for seq in visit.sequences:
                all_priorities.add(seq.priority)

        sorted_priorities = sorted(all_priorities)
        priority_colors = self._get_priority_colors(sorted_priorities)

        # Create legend elements
        legend_elements = []
        for priority in sorted_priorities:
            color = priority_colors[priority]

            # Better labels
            if priority == 0:
                label = f"Priority {priority}"
            elif priority == 1:
                label = f"Priority {priority} "
            elif priority == max(sorted_priorities) and priority > 1:
                label = f"Priority {priority}"
            else:
                label = f"Priority {priority}"

            legend_elements.append(
                Patch(facecolor=color, alpha=1.0, label=label)
            )

        if legend_elements:
            fig.legend(
                handles=legend_elements,
                loc="upper right",
                bbox_to_anchor=(0.98, 0.98),
            )

    def plot_timeline(self, calendar, figsize=(10, 4), show_visits=False):
        """Plot the timeline of visits and their sequences."""
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        for visit in calendar.visits:
            if show_visits:
                ax.vlines(
                    visit.start_time.to_datetime(),
                    -0.45,
                    0.45,
                    lw=2,
                    ls=":",
                    color="k",
                    alpha=1,
                    zorder=100,
                )
            for seq in visit.sequences:

                start = mdates.date2num(seq.start_time.to_datetime())
                duration = (
                    seq.stop_time.to_datetime() - seq.start_time.to_datetime()
                ).total_seconds() / 86400  # days

                color = self._get_priority_colors([seq.priority])[seq.priority]
                ax.barh(0, duration, left=start, color=color, alpha=0.7)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d %H:%M"))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
        ax.set_yticks([])
        ax.set_title("Calendar Timeline", fontsize=14)
        fig.tight_layout()

        return fig, ax

    def plot_target_time(self, calendar, figsize=(10, 4)):
        """Plot the a histogram of the time each target is observed."""

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_title("Target Observation Time")

        # Calculate observation time per target
        proc_targets = self._calculate_target_times(calendar)

        # Get all unique targets
        all_targets = set(proc_targets.keys())

        targets = list(all_targets)
        proc_times = [proc_targets.get(t, 0) for t in targets]

        x = np.arange(len(targets))

        width = 0.7

        ax.bar(x, proc_times, width, alpha=1.0)

        ax.set_xlabel("Target")
        ax.set_ylabel("Total Observation Time (hours)")
        ax.set_xticks(x)
        ax.set_xticklabels(targets, rotation=45, ha="right")

        return fig, ax

    def _calculate_target_times(self, calendar):
        """Calculate the total observation time for each target."""

        target_times = {}
        for visit in calendar.visits:
            for seq in visit.sequences:
                duration_hours = (
                    seq.stop_time.to_datetime() - seq.start_time.to_datetime()
                ).total_seconds() / 3600
                target_times[seq.target] = (
                target_times.get(seq.target, 0) + duration_hours
                )
        return target_times


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize a Pandora science calendar XML.",
    )
    parser.add_argument(
        "xml",
        type=Path,
        help="Path to Pandora_science_calendar.xml",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Output image path. If omitted, show the plot interactively.",
    )
    parser.add_argument(
        "--mode",
        choices=["priority", "timeline", "target-time", "simple", "visibility"],
        default="priority",
        help=(
            "Plot type: 'priority' (default gantt by priority), "
            "'timeline', 'target-time', 'simple', or 'visibility'."
        ),
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Run data directory containing visibility parquet files (required for --mode visibility).",
    )
    parser.add_argument(
        "--show-sequence-labels",
        action="store_true",
        help="Annotate sequence IDs on the plot when supported.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    calendar = parse_science_calendar(args.xml)
    visualizer = ScheduleVisualizer(_DummyScheduler())

    if args.mode == "priority":
        result = visualizer.plot_gantt_timeline_by_priority(
            calendar,
            show_sequence_labels=args.show_sequence_labels,
            title=args.xml.stem,
        )
    elif args.mode == "timeline":
        result = visualizer.plot_timeline(calendar)[0]
    elif args.mode == "target-time":
        result = visualizer.plot_target_time(calendar)[0]
    elif args.mode == "visibility":
        if args.data_dir is None:
            raise SystemExit("--data-dir is required when --mode visibility is used")
        result = visualizer.plot_gantt_with_visibility(
            calendar,
            data_dir=args.data_dir,
            show_sequence_labels=args.show_sequence_labels,
            title=f"{args.xml.stem} — Visibility Check",
        )
    else:
        result = visualizer.plot_simple_priority_timeline(
            calendar,
            show_sequence_labels=args.show_sequence_labels,
            title=args.xml.stem,
        )

    figures = result if isinstance(result, list) else [result]

    if args.out:
        saved_paths = _save_figures(figures, args.out)
        for path in saved_paths:
            print(path)
        plt.close("all")
    else:
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
