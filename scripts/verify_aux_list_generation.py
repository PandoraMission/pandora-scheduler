
import sys
from pathlib import Path

import pandas as pd

from pandorascheduler_rework.pipeline import _generate_target_manifests


def test_aux_list_generation():
    # Setup temporary directory structure
    tmp_dir = Path("output_verification")
    tmp_dir.mkdir(exist_ok=True)
    data_dir = tmp_dir / "data"
    data_dir.mkdir(exist_ok=True, parents=True)
    
    # Create dummy target definition files
    def_dir = Path("comparison_outputs/target_definition_files_limited")
    if not def_dir.exists():
        print(f"Error: Target definition directory {def_dir} not found.")
        sys.exit(1)

    # Define output paths
    primary_csv = data_dir / "exoplanet_targets.csv"
    aux_csv = data_dir / "auxiliary-standard_targets.csv"
    mon_csv = data_dir / "monitoring-standard_targets.csv"
    occ_csv = data_dir / "occultation-standard_targets.csv"
    
    print("Running _generate_target_manifests...")
    try:
        _generate_target_manifests(
            target_definition_files=["exoplanet", "auxiliary-standard", "monitoring-standard", "occultation-standard"],
            base_dir=def_dir,
            primary_target_csv=primary_csv,
            auxiliary_target_csv=aux_csv,
            monitoring_target_csv=mon_csv,
            occultation_target_csv=occ_csv
        )
    except Exception as e:
        print(f"Error during manifest generation: {e}")
        sys.exit(1)

    # Check if all_targets.csv was created
    all_targets = data_dir / "all_targets.csv"
    if all_targets.exists():
        print(f"SUCCESS: {all_targets} was generated.")
        df = pd.read_csv(all_targets)
        print(f"File contains {len(df)} rows.")
        # Verify content - should contain combined targets
        if len(df) > 0:
             print("Content verification passed (file is not empty).")
        else:
             print("WARNING: File is empty.")
    else:
        print(f"FAILURE: {all_targets} was NOT generated.")
        sys.exit(1)

if __name__ == "__main__":
    test_aux_list_generation()
