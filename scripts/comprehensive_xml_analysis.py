#!/usr/bin/env python3
"""
Comprehensive XML Analysis Report

Analyzes differences between legacy and rework XMLs, focusing on:
1. Overall statistics
2. Occultation differences
3. Non-occultation differences (auxiliary targets)
4. Missing visibility data issues
"""

import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


def parse_xml(xml_path):
    """Parse XML and extract visit information."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    ns = "{/pandora/calendar/}"
    
    visits = []
    for visit in root.findall(f"{ns}Visit"):
        visit_id = visit.find(f"{ns}ID")
        if visit_id is None:
            continue
            
        visit_id_text = visit_id.text
        
        obs_seq = visit.find(f"{ns}Observation_Sequence")
        if obs_seq is None:
            visits.append({
                'id': visit_id_text,
                'target': None,
                'start': None,
                'stop': None,
                'empty': True
            })
            continue
        
        obs_params = obs_seq.find(f"{ns}Observational_Parameters")
        if obs_params is None:
            continue
            
        target = obs_params.find(f"{ns}Target")
        timing = obs_params.find(f"{ns}Timing")
        
        target_text = target.text if target is not None else None
        start_text = timing.find(f"{ns}Start").text if timing is not None else None
        stop_text = timing.find(f"{ns}Stop").text if timing is not None else None
        
        visits.append({
            'id': visit_id_text,
            'target': target_text,
            'start': start_text,
            'stop': stop_text,
            'empty': False
        })
    
    return visits

def categorize_target(target_name):
    """Categorize target by type."""
    if target_name is None:
        return 'Empty'
    if 'Free Time' in target_name:
        return 'Free Time'
    if 'STD' in target_name:
        return 'Standard'
    if target_name.startswith('G') and len(target_name) > 15:
        return 'Occultation (Gaia)'
    if target_name.startswith('TIC'):
        return 'Occultation (TIC)'
    if 'WARNING' in target_name:
        return 'Warning'
    if any(x in target_name for x in ['b', 'c', 'd', 'e']) and target_name[-1] in 'bcde':
        return 'Exoplanet'
    if any(target_name.startswith(prefix) for prefix in ['HAT-P-', 'WASP-', 'TOI-', 'K2-', 'KELT-', 'Qatar-', 'XO-', 'HD', 'GJ']):
        return 'Exoplanet/Star'
    # BD targets are auxiliary standards
    if target_name.startswith('BD'):
        return 'Auxiliary Standard'
    if target_name.startswith('GD'):
        return 'Auxiliary Standard'
    return 'Other'

def main():
    legacy_xml = Path("src/pandorascheduler/data/Pandora_science_calendar.xml")
    rework_xml = Path("output_standalone/Pandora_science_calendar.xml")
    
    print("=" * 100)
    print("COMPREHENSIVE XML COMPARISON ANALYSIS")
    print("=" * 100)
    print()
    
    # Parse both XMLs
    legacy_visits = parse_xml(legacy_xml)
    rework_visits = parse_xml(rework_xml)
    
    # Create lookup by (target, start)
    legacy_by_key = {}
    rework_by_key = {}
    
    for v in legacy_visits:
        if not v['empty'] and v['target'] and v['start']:
            key = (v['target'], v['start'])
            legacy_by_key[key] = v
    
    for v in rework_visits:
        if not v['empty'] and v['target'] and v['start']:
            key = (v['target'], v['start'])
            rework_by_key[key] = v
    
    # Find differences
    legacy_keys = set(legacy_by_key.keys())
    rework_keys = set(rework_by_key.keys())
    
    matching = legacy_keys & rework_keys
    legacy_only = legacy_keys - rework_keys
    rework_only = rework_keys - legacy_keys
    
    print("OVERALL STATISTICS")
    print("-" * 100)
    print(f"Total legacy visits: {len(legacy_visits)} ({len(legacy_by_key)} non-empty)")
    print(f"Total rework visits: {len(rework_visits)} ({len(rework_by_key)} non-empty)")
    print(f"Matching visits: {len(matching)}")
    print(f"Legacy-only visits: {len(legacy_only)}")
    print(f"Rework-only visits: {len(rework_only)}")
    print(f"Match rate: {100 * len(matching) / (len(matching) + len(legacy_only) + len(rework_only)):.1f}%")
    print()
    
    # Categorize differences
    legacy_by_cat = defaultdict(list)
    rework_by_cat = defaultdict(list)
    
    for key in legacy_only:
        cat = categorize_target(key[0])
        legacy_by_cat[cat].append(key)
    
    for key in rework_only:
        cat = categorize_target(key[0])
        rework_by_cat[cat].append(key)
    
    print("=" * 100)
    print("SUMMARY BY CATEGORY")
    print("=" * 100)
    print()
    
    all_cats = sorted(set(legacy_by_cat.keys()) | set(rework_by_cat.keys()))
    
    print(f"{'Category':<30} {'Legacy-Only':>15} {'Rework-Only':>15} {'Net Difference':>15}")
    print("-" * 100)
    for cat in all_cats:
        legacy_count = len(legacy_by_cat[cat])
        rework_count = len(rework_by_cat[cat])
        net = rework_count - legacy_count
        sign = "+" if net > 0 else ""
        print(f"{cat:<30} {legacy_count:>15} {rework_count:>15} {sign}{net:>14}")
    
    print()
    print("=" * 100)
    print("KEY FINDINGS")
    print("=" * 100)
    print()
    
    # Occultation analysis
    occ_legacy = len(legacy_by_cat.get('Occultation (Gaia)', []))
    occ_rework = len(rework_by_cat.get('Occultation (Gaia)', []))
    
    print("1. OCCULTATION TARGETS:")
    print(f"   - Legacy has {occ_legacy} occultations not in rework")
    print(f"   - Rework has {occ_rework} occultations not in legacy")
    print(f"   - Net change: +{occ_rework - occ_legacy} occultations in rework")
    print()
    
    # Auxiliary standard analysis
    aux_legacy = len(legacy_by_cat.get('Auxiliary Standard', []))
    aux_rework = len(rework_by_cat.get('Auxiliary Standard', []))
    
    print("2. AUXILIARY STANDARDS (BD/GD targets):")
    print(f"   - Legacy has {aux_legacy} auxiliary standards not in rework")
    print(f"   - Rework has {aux_rework} auxiliary standards not in legacy")
    print(f"   - Net change: {aux_rework - aux_legacy} auxiliary standards")
    print()
    
    if aux_legacy > 0:
        print("   Missing auxiliary standards in rework:")
        aux_targets = set(k[0] for k in legacy_by_cat.get('Auxiliary Standard', []))
        for target in sorted(aux_targets):
            count = sum(1 for k in legacy_by_cat['Auxiliary Standard'] if k[0] == target)
            print(f"     - {target}: {count} visits")
    print()
    
    # Exoplanet analysis
    exo_legacy = len(legacy_by_cat.get('Exoplanet/Star', []))
    exo_rework = len(rework_by_cat.get('Exoplanet/Star', []))
    
    print("3. EXOPLANETS/STARS:")
    print(f"   - Legacy has {exo_legacy} exoplanet visits not in rework")
    print(f"   - Rework has {exo_rework} exoplanet visits not in legacy")
    print()
    
    if exo_legacy > 0:
        print("   Exoplanet targets missing in rework:")
        exo_targets = defaultdict(int)
        for k in legacy_by_cat.get('Exoplanet/Star', []):
            exo_targets[k[0]] += 1
        for target, count in sorted(exo_targets.items(), key=lambda x: -x[1])[:15]:
            print(f"     - {target}: {count} visits")
    print()
    
    print("=" * 100)
    print("CRITICAL ISSUES IDENTIFIED")
    print("=" * 100)
    print()
    
    # These are the targets that had visibility errors
    error_targets = ['BD+21_0607', 'HD_233511', 'HD_116405', 'HD_31128', 'GD_153', 'HD_111980']
    
    print("Targets with 'No visibility data' errors during XML generation:")
    for target in error_targets:
        legacy_count = sum(1 for k in legacy_only if k[0] == target)
        rework_count = sum(1 for k in rework_only if k[0] == target)
        if legacy_count > 0 or rework_count > 0:
            print(f"  - {target}: {legacy_count} in legacy-only, {rework_count} in rework-only")
    
    print()
    print("These errors indicate that the rework's science calendar generation is trying")
    print("to schedule these targets but cannot find their visibility files, causing the")
    print("schedule build to abort for those specific time windows.")
    print()
    
    print("=" * 100)
    print("CONCLUSION")
    print("=" * 100)
    print()
    print("The differences between legacy and rework XMLs are NOT limited to occultations.")
    print()
    print("Main differences:")
    print(f"  1. Occultations: +{occ_rework - occ_legacy} net (different targets selected)")
    print(f"  2. Auxiliary standards: -{aux_legacy} (missing BD/GD targets in rework)")
    print(f"  3. Exoplanets: -{exo_legacy - exo_rework} net (some targets missing)")
    print()
    print("Root causes:")
    print("  - Missing visibility files for auxiliary standard targets (BD/GD)")
    print("  - Different occultation selection logic")
    print("  - Possible differences in gap-filling auxiliary target selection")

if __name__ == "__main__":
    main()
