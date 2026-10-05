#!/usr/bin/env python3
"""Detailed XML comparison to identify all differences between Legacy and Rework."""

import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


def parse_xml(path):
    """Parse XML and return root."""
    tree = ET.parse(path)
    return tree.getroot()

def extract_visit_data(visit):
    """Extract all data from a visit element."""
    data = {}
    id_elem = visit.find('.//{*}ID')
    data['id'] = id_elem.text if id_elem is not None else None
    
    # Extract all observation sequences (handle namespaces)
    sequences = []
    for seq in visit.findall('.//{*}Observation_Sequence'):
        seq_data = {}
        seq_id = seq.find('.//{*}ID')
        seq_data['id'] = seq_id.text if seq_id is not None else None
        
        # Observational Parameters
        obs_params = seq.find('.//{*}Observational_Parameters')
        if obs_params is not None:
            for child in obs_params:
                # Strip namespace from tag
                tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                seq_data[f'obs_{tag}'] = child.text
        
        # Payload Parameters
        payload = seq.find('.//{*}Payload_Parameters')
        if payload is not None:
            # NIRDA
            nirda = payload.find('.//{*}AcquireInfCamImages')
            if nirda is not None:
                for child in nirda:
                    tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                    seq_data[f'NIRDA_{tag}'] = child.text
            
            # VDA
            vda = payload.find('.//{*}AcquireVisCamScienceData')
            if vda is not None:
                for child in vda:
                    tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                    seq_data[f'VDA_{tag}'] = child.text
        
        sequences.append(seq_data)
    
    data['sequences'] = sequences
    return data

def compare_xmls(legacy_path, rework_path):
    """Compare two XML files and report differences."""
    print("=" * 100)
    print("DETAILED XML COMPARISON")
    print("=" * 100)
    
    legacy_root = parse_xml(legacy_path)
    rework_root = parse_xml(rework_path)
    
    # Compare Meta (handle namespaces)
    print("\n1. META COMPARISON")
    print("-" * 100)
    legacy_meta = legacy_root.find('.//{*}Meta')
    rework_meta = rework_root.find('.//{*}Meta')
    
    if legacy_meta is not None and rework_meta is not None:
        legacy_attrs = legacy_meta.attrib
        rework_attrs = rework_meta.attrib
        
        for key in sorted(set(legacy_attrs.keys()) | set(rework_attrs.keys())):
            legacy_val = legacy_attrs.get(key, 'MISSING')
            rework_val = rework_attrs.get(key, 'MISSING')
            if legacy_val != rework_val:
                print(f"  {key}:")
                print(f"    Legacy: {legacy_val}")
                print(f"    Rework: {rework_val}")
    
    # Extract all visits (handle namespaces with universal syntax)
    legacy_visits = {}
    for v in legacy_root.findall('.//{*}Visit'):
        id_elem = v.find('.//{*}ID')
        if id_elem is not None:
            legacy_visits[id_elem.text] = extract_visit_data(v)
    
    rework_visits = {}
    for v in rework_root.findall('.//{*}Visit'):
        id_elem = v.find('.//{*}ID')
        if id_elem is not None:
            rework_visits[id_elem.text] = extract_visit_data(v)
    
    print("\n2. VISIT COUNT")
    print("-" * 100)
    print(f"  Legacy: {len(legacy_visits)} visits")
    print(f"  Rework: {len(rework_visits)} visits")
    
    # Find common visits
    common_ids = set(legacy_visits.keys()) & set(rework_visits.keys())
    legacy_only = set(legacy_visits.keys()) - set(rework_visits.keys())
    rework_only = set(rework_visits.keys()) - set(legacy_visits.keys())
    
    if legacy_only:
        print(f"\n  Legacy-only visits: {sorted(legacy_only)}")
    if rework_only:
        print(f"  Rework-only visits: {sorted(rework_only)}")
    
    # Compare common visits
    print(f"\n3. DETAILED PARAMETER COMPARISON (Common visits: {len(common_ids)})")
    print("-" * 100)
    
    diff_categories = defaultdict(list)
    
    for visit_id in sorted(common_ids):
        legacy_visit = legacy_visits[visit_id]
        rework_visit = rework_visits[visit_id]
        
        # Compare sequences
        legacy_seqs = legacy_visit['sequences']
        rework_seqs = rework_visit['sequences']
        
        if len(legacy_seqs) != len(rework_seqs):
            diff_categories['sequence_count'].append(f"Visit {visit_id}: Legacy has {len(legacy_seqs)}, Rework has {len(rework_seqs)}")
            continue
        
        for i, (legacy_seq, rework_seq) in enumerate(zip(legacy_seqs, rework_seqs)):
            # Get all keys
            all_keys = set(legacy_seq.keys()) | set(rework_seq.keys())
            
            for key in all_keys:
                if key == 'id':
                    continue
                    
                legacy_val = legacy_seq.get(key, 'MISSING')
                rework_val = rework_seq.get(key, 'MISSING')
                
                if legacy_val != rework_val:
                    # Skip known differences
                    if 'Created' in key:
                        diff_categories['created_timestamp'].append(f"Visit {visit_id}, Seq {i}: {key}")
                    elif 'RA' in key or 'DEC' in key:
                        try:
                            legacy_float = float(legacy_val)
                            rework_float = float(rework_val)
                            diff = abs(legacy_float - rework_float)
                            if diff > 0.0001:  # Significant difference
                                diff_categories['coordinates'].append(
                                    f"Visit {visit_id}, Seq {i}, {key}: Legacy={legacy_val}, Rework={rework_val}, Diff={diff:.6f}"
                                )
                        except (TypeError, ValueError):
                            diff_categories['coordinates'].append(
                                f"Visit {visit_id}, Seq {i}, {key}: Legacy={legacy_val}, Rework={rework_val}"
                            )
                    elif 'MaxNumStarRois' in key:
                        diff_categories['max_num_star_rois'].append(
                            f"Visit {visit_id}, Seq {i}: Legacy={legacy_val}, Rework={rework_val}"
                        )
                    else:
                        diff_categories['other'].append(
                            f"Visit {visit_id}, Seq {i}, {key}: Legacy={legacy_val}, Rework={rework_val}"
                        )
    
    # Print categorized differences
    print("\nDIFFERENCE CATEGORIES:")
    print()
    
    for category, diffs in sorted(diff_categories.items()):
        print(f"\n{category.upper().replace('_', ' ')} ({len(diffs)} differences):")
        if len(diffs) <= 10:
            for diff in diffs:
                print(f"  {diff}")
        else:
            print(f"  Showing first 10 of {len(diffs)}:")
            for diff in diffs[:10]:
                print(f"  {diff}")
            print(f"  ... and {len(diffs) - 10} more")
    
    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    for category, diffs in sorted(diff_categories.items()):
        print(f"  {category}: {len(diffs)} differences")
    print()

if __name__ == "__main__":
    legacy_xml = Path("src/pandorascheduler/data/Pandora_science_calendar.xml")
    rework_xml = Path("output_standalone/data/Pandora_science_calendar.xml")
    
    compare_xmls(legacy_xml, rework_xml)
