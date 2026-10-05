#!/usr/bin/env python3
"""
Validate integrity of visibility data files.
Checks that all expected planet visibility files exist and contain required data.
"""
import sys
from pathlib import Path

import pandas as pd


def main():
    targets_dir = Path("src/pandorascheduler/data/targets")
    manifest_path = Path("src/pandorascheduler/data/exoplanet_targets.csv")
    
    if not manifest_path.exists():
        print(f"ERROR: Manifest not found: {manifest_path}")
        return 1
    
    manifest = pd.read_csv(manifest_path)
    
    print(f"Validating visibility data for {len(manifest)} planets...")
    print("=" * 80)
    
    missing_files = []
    empty_files = []
    missing_saa = []
    missing_overlap = []
    valid_files = []
    
    for _, row in manifest.iterrows():
        star = row["Star Name"]
        planet = row["Planet Name"]
        vis_file = targets_dir / star / planet / f"Visibility for {planet}.parquet"
        
        if not vis_file.exists():
            missing_files.append(str(vis_file))
            print(f"❌ MISSING: {planet}")
            continue
        
        try:
            df = pd.read_parquet(vis_file)
        except Exception as e:
            print(f"❌ ERROR reading {planet}: {e}")
            continue
        
        if len(df) == 0:
            empty_files.append(str(vis_file))
            print(f"⚠️  EMPTY (no transits): {planet}")
            continue
        
        issues = []
        if "SAA_Overlap" not in df.columns:
            missing_saa.append(str(vis_file))
            issues.append("no SAA_Overlap")
        
        if "Transit_Overlap" not in df.columns:
            missing_overlap.append(str(vis_file))
            issues.append("no Transit_Overlap")
        
        if issues:
            print(f"⚠️  INCOMPLETE: {planet} ({', '.join(issues)}, {len(df)} transits)")
        else:
            valid_files.append(str(vis_file))
            print(f"✅ OK: {planet} ({len(df)} transits)")
    
    print("=" * 80)
    print("\nSUMMARY:")
    print(f"  ✅ Valid files:           {len(valid_files)}")
    print(f"  ⚠️  Empty files:           {len(empty_files)}")
    print(f"  ⚠️  Missing SAA_Overlap:   {len(missing_saa)}")
    print(f"  ⚠️  Missing Transit_Overlap: {len(missing_overlap)}")
    print(f"  ❌ Missing files:         {len(missing_files)}")
    
    
    if empty_files or missing_files or missing_saa:
        print("\n⚠️  VISIBILITY DATA IS INCOMPLETE OR CORRUPTED")
        print("   Comparison script will fail until data is restored/regenerated.")
        return 1
    elif missing_overlap:
        print("\n⚠️  Some files missing Transit_Overlap column")
        print("   This is OK - schedulers will calculate it on-the-fly if needed.")
        return 0
    else:
        print("\n✅ All visibility data is valid and complete!")
        return 0

if __name__ == "__main__":
    sys.exit(main())
