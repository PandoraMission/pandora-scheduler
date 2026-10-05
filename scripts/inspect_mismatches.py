import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


def inspect_mismatch():
    legacy_dir = Path("src/pandorascheduler/data/baseline")
    rework_dir = Path("output_standalone")
    legacy_data_dir = Path("src/pandorascheduler/data")
    
    legacy_xml = legacy_data_dir / "Pandora_science_calendar.xml"
    rework_xml = rework_dir / "Pandora_science_calendar.xml"
    
    print(f"Parsing Legacy XML: {legacy_xml}")
    tree_leg = ET.parse(legacy_xml)
    root_leg = tree_leg.getroot()
    
    print(f"Parsing Rework XML: {rework_xml}")
    tree_rew = ET.parse(rework_xml)
    root_rew = tree_rew.getroot()
    
    def get_visit_target(root, visit_id):
        ns = "{/pandora/calendar/}"
        # Find visit with ID
        for visit in root.findall(f"{ns}Visit"):
            vid = visit.find(f"{ns}ID").text
            if vid == visit_id:
                obs = visit.find(f"{ns}Observation_Sequence")
                if obs is not None:
                    params = obs.find(f"{ns}Observational_Parameters")
                    if params is not None:
                        return params.find(f"{ns}Target").text
        return "NOT FOUND"

    target_leg = get_visit_target(root_leg, "0211")
    target_rew = get_visit_target(root_rew, "0211")
    
    print("\nXML Visit 0211 Target:")
    print(f"  Legacy: {target_leg}")
    print(f"  Rework: {target_rew}")
    
    # Check CSV for HAT-P-18 and Gaia ID
    legacy_sched_path = list(legacy_dir.glob("Pandora_Schedule_*.csv"))[0]
    df_leg = pd.read_csv(legacy_sched_path)
    print("Legacy CSV Row 210:")
    print(df_leg.iloc[210])

    rework_sched_path = list(rework_dir.glob("Pandora_Schedule_*.csv"))[0]
    print(f"Reading CSV: {rework_sched_path}")
    df = pd.read_csv(rework_sched_path)
    print(f"Columns: {df.columns.tolist()}")
    
    print(f"Unique Targets: {df['Target'].unique().tolist()}")
    
    print("\nSearching CSV for HAT-P-18:")
    matches = df[df['Target'] == 'HAT-P-18']
    if not matches.empty:
        print(matches[['Observation Start']].to_string())
    else:
        print("HAT-P-18 NOT FOUND in CSV")

    print("\nSearching CSV for G3738601634818313472:")
    matches_gaia = df[df['Target'] == 'G3738601634818313472']
    if not matches_gaia.empty:
        print(matches_gaia.to_string())
        print(f"Index of Gaia ID: {matches_gaia.index.tolist()}")
    else:
        print("Gaia ID NOT FOUND in CSV")
    
    # Find the row corresponding to Visit 211 (skipping Free Time)
    visit_count = 0
    target_at_211 = None
    row_index_211 = -1
    
    for idx, row in df.iterrows():
        target = str(row['Target'])
        if not target or target == 'nan' or target == 'Free Time' or target.startswith('WARNING'):
            continue
        
        visit_count += 1
        if visit_count == 211:
            target_at_211 = target
            row_index_211 = idx
            break
            
    print("\nCalculated CSV Row for Visit 211:")
    print(f"  Row Index: {row_index_211}")
    print(f"  Target: {target_at_211}")
    
    if row_index_211 != -1:
        print(df.iloc[row_index_211])
        
    print("\nCSV Rows 200-240:")
    print(df.iloc[200:240][['Target', 'Observation Start']].to_string())

if __name__ == "__main__":
    inspect_mismatch()
