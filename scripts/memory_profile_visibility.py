#!/usr/bin/env python3
"""
Memory profiling script for visibility data caching analysis.

Analyzes the size and count of visibility files to determine safe cache limits.
"""

import sys
from pathlib import Path

import pandas as pd


def analyze_visibility_files(data_dir: Path):
    """Analyze all visibility files in the data directory."""
    
    results = {
        'target_visibility': [],
        'planet_visibility': [],
        'aux_visibility': []
    }
    
    # Analyze primary target visibility
    targets_dir = data_dir / "targets"
    if targets_dir.exists():
        for star_dir in targets_dir.iterdir():
            if not star_dir.is_dir():
                continue
            
            # Star visibility file
            star_vis = star_dir / f"Visibility for {star_dir.name}.parquet"
            if star_vis.exists():
                size = star_vis.stat().st_size
                results['target_visibility'].append({
                    'file': star_vis.name,
                    'size_bytes': size,
                    'size_mb': size / (1024 * 1024)
                })
            
            # Planet visibility files
            for planet_dir in star_dir.iterdir():
                if not planet_dir.is_dir():
                    continue
                planet_vis = planet_dir / f"Visibility for {planet_dir.name}.parquet"
                if planet_vis.exists():
                    size = planet_vis.stat().st_size
                    results['planet_visibility'].append({
                        'file': planet_vis.name,
                        'size_bytes': size,
                        'size_mb': size / (1024 * 1024)
                    })
    
    # Analyze auxiliary target visibility
    aux_dir = data_dir / "aux_targets"
    if aux_dir.exists():
        for star_dir in aux_dir.iterdir():
            if not star_dir.is_dir():
                continue
            star_vis = star_dir / f"Visibility for {star_dir.name}.parquet"
            if star_vis.exists():
                size = star_vis.stat().st_size
                results['aux_visibility'].append({
                    'file': star_vis.name,
                    'size_bytes': size,
                    'size_mb': size / (1024 * 1024)
                })
    
    return results

def estimate_memory_usage(results):
    """Estimate memory usage if all files are cached."""
    
    # DataFrame overhead: roughly 2-3x the CSV file size when loaded into memory
    # (includes index, column names, numpy arrays, etc.)
    overhead_multiplier = 2.5
    
    print("=" * 80)
    print("VISIBILITY FILE ANALYSIS")
    print("=" * 80)
    
    total_size_mb = 0
    
    for category, files in results.items():
        if not files:
            continue
            
        count = len(files)
        sizes = [f['size_mb'] for f in files]
        total_cat_mb = sum(sizes)
        total_size_mb += total_cat_mb
        
        print(f"\n{category.upper().replace('_', ' ')}:")
        print(f"  Count: {count} files")
        print(f"  Total size: {total_cat_mb:.2f} MB")
        print(f"  Average size: {total_cat_mb/count:.2f} MB" if count > 0 else "  N/A")
        print(f"  Min size: {min(sizes):.2f} MB" if sizes else "  N/A")
        print(f"  Max size: {max(sizes):.2f} MB" if sizes else "  N/A")
    
    print("\n" + "=" * 80)
    print("MEMORY ESTIMATION")
    print("=" * 80)
    print(f"Total Parquet size on disk: {total_size_mb:.2f} MB")
    print(f"Estimated memory if all cached: {total_size_mb * overhead_multiplier:.2f} MB")
    print(f"  (using {overhead_multiplier}x multiplier for DataFrame overhead)")
    
    # Recommendations
    print("\n" + "=" * 80)
    print("CACHING RECOMMENDATIONS")
    print("=" * 80)
    
    estimated_mem = total_size_mb * overhead_multiplier
    
    if estimated_mem < 500:  # Less than 500 MB
        print("✅ SAFE: Total memory usage is low.")
        print(f"   Recommended: Unlimited caching (or set maxsize={sum(len(f) for f in results.values())})")
    elif estimated_mem < 2000:  # Less than 2 GB
        print("⚠️  MODERATE: Memory usage is acceptable but should be limited.")
        print("   Recommended: Set maxsize=256 to limit cache")
    else:  # > 2 GB
        print("❌ HIGH: Memory usage could be problematic.")
        print("   Recommended: Set maxsize=128 or use LRU eviction")
        print("   Consider: Lazy loading or streaming for large files")
    
    print()

def sample_load_test(data_dir: Path, num_samples=5):
    """Load a few sample files to measure actual memory overhead."""
    print("=" * 80)
    print("SAMPLE LOAD TEST")
    print("=" * 80)
    
    # Find some sample files
    targets_dir = data_dir / "targets"
    sample_files = []
    
    if targets_dir.exists():
        for star_dir in list(targets_dir.iterdir())[:num_samples]:
            if not star_dir.is_dir():
                continue
            star_vis = star_dir / f"Visibility for {star_dir.name}.parquet"
            if star_vis.exists():
                sample_files.append(star_vis)
    
    if not sample_files:
        print("No sample files found.")
        return
    
    print(f"\nLoading {len(sample_files)} sample files...\n")
    
    for file_path in sample_files:
        file_size_mb = file_path.stat().st_size / (1024 * 1024)
        
        try:
            df = pd.read_parquet(file_path)
            # Rough memory estimate
            mem_usage_mb = df.memory_usage(deep=True).sum() / (1024 * 1024)
            overhead = mem_usage_mb / file_size_mb if file_size_mb > 0 else 0
            
            print(f"  {file_path.name}:")
            print(f"    File size: {file_size_mb:.3f} MB")
            print(f"    Memory usage: {mem_usage_mb:.3f} MB")
            print(f"    Overhead: {overhead:.2f}x")
        except Exception as e:
            print(f"  {file_path.name}: Error loading - {e}")
    
    print()

if __name__ == "__main__":
    # Default to output_standalone/data
    if len(sys.argv) > 1:
        data_dir = Path(sys.argv[1])
    else:
        data_dir = Path("output_standalone/data")
    
    if not data_dir.exists():
        print(f"Error: Data directory not found: {data_dir}")
        print("Usage: python memory_profile_visibility.py [data_dir]")
        sys.exit(1)
    
    print(f"Analyzing visibility files in: {data_dir}\n")
    
    results = analyze_visibility_files(data_dir)
    estimate_memory_usage(results)
    sample_load_test(data_dir)
