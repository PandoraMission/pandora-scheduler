#!/usr/bin/env python3
"""
Analyze gaps, overlaps, and short sequences in Pandora science calendar XML.
"""
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path
import sys


def parse_time(time_str: str) -> datetime:
    """Parse a time string in ISO format."""
    return datetime.fromisoformat(time_str.replace('Z', '+00:00'))


def _fmt_dur(td: timedelta) -> str:
    """Format a timedelta as a human-readable string."""
    total_sec = abs(td.total_seconds())
    if total_sec < 60:
        return f"{total_sec:.0f}s"
    if total_sec < 3600:
        return f"{total_sec / 60:.1f}m"
    return f"{total_sec / 3600:.2f}h"


def _parse_sequences(root: ET.Element) -> list[dict]:
    """Extract all observation sequences from the XML."""
    sequences = []
    for visit in root.findall('.//{/pandora/calendar/}Visit'):
        visit_id = visit.find('.//{/pandora/calendar/}ID')
        visit_id_str = visit_id.text if visit_id is not None else "Unknown"

        for obs_seq in visit.findall('.//{/pandora/calendar/}Observation_Sequence'):
            seq_id = obs_seq.find('.//{/pandora/calendar/}ID')
            seq_id_str = seq_id.text if seq_id is not None else "Unknown"

            obs_params = obs_seq.find('.//{/pandora/calendar/}Observational_Parameters')
            if obs_params is None:
                continue

            target = obs_params.find('.//{/pandora/calendar/}Target')
            target_name = target.text if target is not None else "Unknown"

            timing = obs_params.find('.//{/pandora/calendar/}Timing')
            if timing is None:
                continue

            start_elem = timing.find('.//{/pandora/calendar/}Start')
            stop_elem = timing.find('.//{/pandora/calendar/}Stop')
            if start_elem is None or stop_elem is None:
                continue

            start_time = parse_time(start_elem.text)
            end_time = parse_time(stop_elem.text)

            sequences.append({
                'target': target_name,
                'visit_id': visit_id_str,
                'seq_id': seq_id_str,
                'start': start_time,
                'end': end_time,
                'duration': end_time - start_time,
            })
    sequences.sort(key=lambda x: x['start'])
    return sequences


def _find_gaps_and_overlaps(sequences: list[dict]):
    """Return (gaps, overlaps) lists."""
    gaps, overlaps = [], []
    for i in range(len(sequences) - 1):
        delta = sequences[i + 1]['start'] - sequences[i]['end']
        delta_sec = delta.total_seconds()
        entry = {
            'idx': i,
            'after': sequences[i],
            'before': sequences[i + 1],
            'delta': delta,
            'delta_sec': delta_sec,
        }
        if delta_sec > 0:
            gaps.append(entry)
        elif delta_sec < 0:
            overlaps.append(entry)
    return gaps, overlaps


def _find_short_sequences(sequences: list[dict], threshold_minutes: float = 5.0):
    """Return sequences shorter than threshold."""
    threshold = timedelta(minutes=threshold_minutes)
    return [
        (i, seq) for i, seq in enumerate(sequences)
        if seq['duration'] < threshold
    ]


def _find_zero_or_negative(sequences: list[dict]):
    """Return sequences with zero or negative duration."""
    return [
        (i, seq) for i, seq in enumerate(sequences)
        if seq['duration'].total_seconds() <= 0
    ]


def _print_top_items(items, key_func, label, n=20):
    """Print top N items sorted by key_func descending."""
    if not items:
        print(f"\n  (none)")
        return
    top = sorted(items, key=key_func, reverse=True)[:n]
    for rank, item in enumerate(top, 1):
        yield rank, item


def analyze_schedule(xml_path: Path, short_threshold_min: float = 5.0):
    """Full schedule quality analysis: gaps, overlaps, short sequences."""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    sequences = _parse_sequences(root)
    if not sequences:
        print("No sequences found in XML.")
        return sequences, [], [], []

    gaps, overlaps = _find_gaps_and_overlaps(sequences)
    short_seqs = _find_short_sequences(sequences, short_threshold_min)
    bad_dur_seqs = _find_zero_or_negative(sequences)

    # ── Summary ──────────────────────────────────────────────────────────────
    print("=" * 80)
    print(f"SCHEDULE QUALITY REPORT: {xml_path.name}")
    print("=" * 80)

    total_scheduled = sum(s['duration'].total_seconds() for s in sequences) / 3600
    total_gap = sum(g['delta_sec'] for g in gaps) / 3600
    total_overlap = sum(abs(o['delta_sec']) for o in overlaps) / 3600
    schedule_span = (sequences[-1]['end'] - sequences[0]['start']).total_seconds() / 3600

    print(f"\n  Total sequences:     {len(sequences)}")
    print(f"  Schedule span:       {schedule_span:.2f} hours")
    print(f"  Scheduled time:      {total_scheduled:.2f} hours")
    print(f"  Gap time:            {total_gap:.4f} hours ({len(gaps)} gaps)")
    print(f"  Overlap time:        {total_overlap:.4f} hours ({len(overlaps)} overlaps)")
    print(f"  Efficiency:          {100 * total_scheduled / schedule_span:.2f}%")
    print(f"  Schedule start:      {sequences[0]['start']}")
    print(f"  Schedule end:        {sequences[-1]['end']}")

    # ── Duration statistics ──────────────────────────────────────────────────
    durations_min = [s['duration'].total_seconds() / 60 for s in sequences]
    durations_min.sort()
    n = len(durations_min)
    print(f"\n  Sequence durations:")
    print(f"    Min:    {durations_min[0]:.1f} min")
    print(f"    p5:     {durations_min[int(n*0.05)]:.1f} min")
    print(f"    p25:    {durations_min[int(n*0.25)]:.1f} min")
    print(f"    Median: {durations_min[n//2]:.1f} min")
    print(f"    p75:    {durations_min[int(n*0.75)]:.1f} min")
    print(f"    p95:    {durations_min[int(n*0.95)]:.1f} min")
    print(f"    Max:    {durations_min[-1]:.1f} min")
    print(f"    Mean:   {sum(durations_min)/n:.1f} min")

    # ── Zero/negative duration sequences ─────────────────────────────────────
    if bad_dur_seqs:
        print(f"\n{'=' * 80}")
        print(f"ZERO/NEGATIVE DURATION SEQUENCES: {len(bad_dur_seqs)}")
        print("=" * 80)
        for idx, seq in bad_dur_seqs[:20]:
            print(f"  #{idx+1}: {seq['target']} visit={seq['visit_id']} "
                  f"seq={seq['seq_id']} {seq['start']} -> {seq['end']} "
                  f"({_fmt_dur(seq['duration'])})")
    else:
        print(f"\n  No zero/negative duration sequences.")

    # ── Short sequences ──────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print(f"SHORT SEQUENCES (< {short_threshold_min:.0f} min): {len(short_seqs)}")
    print("=" * 80)
    if short_seqs:
        # By-target summary
        target_counts = Counter(seq['target'] for _, seq in short_seqs)
        print(f"\n  By target:")
        for tgt, cnt in target_counts.most_common(20):
            print(f"    {tgt:40s} {cnt}")

        # Duration distribution
        short_durs = [seq['duration'].total_seconds() / 60 for _, seq in short_seqs]
        short_durs.sort()
        print(f"\n  Duration distribution of short sequences:")
        for bucket_max in [0.5, 1, 2, 3, 5]:
            count = sum(1 for d in short_durs if d < bucket_max)
            print(f"    < {bucket_max} min:  {count}")

        print(f"\n  Shortest 20:")
        for idx, seq in sorted(short_seqs, key=lambda x: x[1]['duration'])[:20]:
            dur_sec = seq['duration'].total_seconds()
            print(f"    #{idx+1}: {seq['target']:35s} {_fmt_dur(seq['duration']):>6s}  "
                  f"{seq['start']} visit={seq['visit_id']}")
    else:
        print(f"\n  No sequences shorter than {short_threshold_min:.0f} min.")

    # ── Overlaps ─────────────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print(f"OVERLAPS: {len(overlaps)}")
    print("=" * 80)
    if overlaps:
        overlap_dur = [abs(o['delta_sec']) / 60 for o in overlaps]
        print(f"\n  Total overlap time:  {total_overlap:.4f} hours")
        print(f"  Min overlap:         {min(overlap_dur):.2f} min")
        print(f"  Max overlap:         {max(overlap_dur):.2f} min")
        print(f"  Mean overlap:        {sum(overlap_dur)/len(overlap_dur):.2f} min")

        print(f"\n  Top 20 largest overlaps:")
        for o in sorted(overlaps, key=lambda x: abs(x['delta_sec']), reverse=True)[:20]:
            print(f"    Overlap of {_fmt_dur(abs(o['delta']))}")
            print(f"      Seq #{o['idx']+1}: {o['after']['target']:35s} ends   {o['after']['end']}")
            print(f"      Seq #{o['idx']+2}: {o['before']['target']:35s} starts {o['before']['start']}")
    else:
        print(f"\n  No overlaps found.")

    # ── Gaps ─────────────────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print(f"GAPS: {len(gaps)}")
    print("=" * 80)
    if gaps:
        gap_sec = [g['delta_sec'] for g in gaps]
        gap_min = [s / 60 for s in gap_sec]
        print(f"\n  Total gap time:  {total_gap:.4f} hours")
        print(f"  Min gap:         {min(gap_min):.2f} min")
        print(f"  Max gap:         {max(gap_min):.2f} min")
        print(f"  Mean gap:        {sum(gap_min)/len(gap_min):.2f} min")
        print(f"  Median gap:      {sorted(gap_min)[len(gap_min)//2]:.2f} min")

        # Distribution
        print(f"\n  Gap distribution:")
        for label, lo, hi in [
            ("< 1 min", 0, 1), ("1-5 min", 1, 5), ("5-30 min", 5, 30),
            ("30-60 min", 30, 60), ("1-6 hr", 60, 360), (">= 6 hr", 360, 1e9),
        ]:
            count = sum(1 for m in gap_min if lo <= m < hi)
            if count:
                print(f"    {label:12s}  {count}")

        print(f"\n  Top 20 largest gaps:")
        for g in sorted(gaps, key=lambda x: x['delta_sec'], reverse=True)[:20]:
            print(f"    Gap of {_fmt_dur(g['delta'])}")
            print(f"      Seq #{g['idx']+1}: {g['after']['target']:35s} ends   {g['after']['end']}")
            print(f"      Seq #{g['idx']+2}: {g['before']['target']:35s} starts {g['before']['start']}")
    else:
        print(f"\n  No gaps found - sequences are perfectly contiguous!")

    # ── Export CSV ───────────────────────────────────────────────────────────
    csv_path = xml_path.parent / "schedule_quality.csv"
    with open(csv_path, 'w') as f:
        f.write("issue_type,seq_num,target,visit_id,seq_id,start,end,duration_min,"
                "other_target,other_start,other_end,delta_min\n")
        for g in gaps:
            f.write(f"gap,{g['idx']+1},{g['after']['target']},{g['after']['visit_id']},"
                    f"{g['after']['seq_id']},{g['after']['start']},{g['after']['end']},"
                    f"{g['after']['duration'].total_seconds()/60:.2f},"
                    f"{g['before']['target']},{g['before']['start']},{g['before']['end']},"
                    f"{g['delta_sec']/60:.4f}\n")
        for o in overlaps:
            f.write(f"overlap,{o['idx']+1},{o['after']['target']},{o['after']['visit_id']},"
                    f"{o['after']['seq_id']},{o['after']['start']},{o['after']['end']},"
                    f"{o['after']['duration'].total_seconds()/60:.2f},"
                    f"{o['before']['target']},{o['before']['start']},{o['before']['end']},"
                    f"{o['delta_sec']/60:.4f}\n")
        for idx, seq in short_seqs:
            f.write(f"short,{idx+1},{seq['target']},{seq['visit_id']},"
                    f"{seq['seq_id']},{seq['start']},{seq['end']},"
                    f"{seq['duration'].total_seconds()/60:.2f},"
                    f",,,"
                    f"\n")
    print(f"\nIssue data exported to: {csv_path}")

    return sequences, gaps, overlaps, short_seqs


if __name__ == "__main__":
    threshold = 5.0
    xml_path = None

    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--short-threshold" and i + 1 < len(args):
            threshold = float(args[i + 1])
            i += 2
        elif not args[i].startswith("-"):
            xml_path = Path(args[i])
            i += 1
        else:
            i += 1

    if xml_path is None:
        xml_path = Path("output_standalone_test2/data_91_20_110/Pandora_science_calendar.xml")

    if not xml_path.exists():
        print(f"Error: XML file not found at {xml_path}")
        sys.exit(1)

    analyze_schedule(xml_path, short_threshold_min=threshold)
