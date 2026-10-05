import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


def compare_csvs(file1, file2, index_col=None):
    print(f"Comparing {file1} vs {file2}...")
    try:
        df1 = pd.read_csv(file1, index_col=index_col)
        df2 = pd.read_csv(file2, index_col=index_col)
        
        if df1.shape != df2.shape:
            print(f"  SHAPE MISMATCH: {df1.shape} vs {df2.shape}")
            return False
        
        # Sort by index if present, or by first column
        if index_col is None:
            df1 = df1.sort_values(by=df1.columns[0]).reset_index(drop=True)
            df2 = df2.sort_values(by=df2.columns[0]).reset_index(drop=True)
            
        pd.testing.assert_frame_equal(df1, df2, check_dtype=False, rtol=1e-5)
        print("  MATCH: DataFrames are identical (within tolerance).")
        return True
    except AssertionError as e:
        print(f"  MISMATCH: {e}")
        return False
    except Exception as e:
        print(f"  ERROR: {e}")
        return False

def compare_xml_deep(file1, file2):
    print(f"Deep comparing XML {file1} vs {file2}...")
    try:
        tree1 = ET.parse(file1)
        tree2 = ET.parse(file2)
        root1 = tree1.getroot()
        root2 = tree2.getroot()

        # Compare Meta
        meta1 = root1.find("{/pandora/calendar/}Meta")
        if meta1 is None:
            # Try without namespace if failed
            meta1 = root1.find("Meta")

        meta2 = root2.find("{/pandora/calendar/}Meta")
        if meta2 is None:
            meta2 = root2.find("Meta")

        if meta1 is not None and meta2 is not None:
            attr1 = dict(meta1.attrib)
            attr2 = dict(meta2.attrib)
            attr1.pop("Created", None)
            attr2.pop("Created", None)
            
            if attr1 != attr2:
                print("  MISMATCH: Meta attributes differ.")
                for k in set(attr1) | set(attr2):
                    v1 = attr1.get(k)
                    v2 = attr2.get(k)
                    if v1 != v2:
                        print(f"    {k}: {v1} vs {v2}")
        
        # Helper to find visits regardless of namespace
        def get_visits(root):
            visits = root.findall("{/pandora/calendar/}Visit")
            if not visits:
                visits = root.findall("Visit")
            return visits

        visits1 = get_visits(root1)
        visits2 = get_visits(root2)

        print(f"  Visit Count: {len(visits1)} vs {len(visits2)}")
        if len(visits1) != len(visits2):
            print("  MISMATCH: Visit count differs.")
            # Don't return, try to compare what we can

        mismatches = 0
        for i, (v1, v2) in enumerate(zip(visits1, visits2)):
            # Helper to find child text
            def get_text(elem, tag):
                # Try with namespace
                child = elem.find(f"{{/pandora/calendar/}}{tag}")
                if child is None:
                    child = elem.find(tag)
                return child.text if child is not None else None

            def get_child(elem, tag):
                 child = elem.find(f"{{/pandora/calendar/}}{tag}")
                 if child is None:
                     child = elem.find(tag)
                 return child

            id1 = get_text(v1, "ID")
            id2 = get_text(v2, "ID")
            
            if id1 != id2:
                print(f"  MISMATCH at index {i}: ID {id1} vs {id2}")
                mismatches += 1
                if mismatches > 10:
                    break
                continue

            obs1 = get_child(v1, "Observation_Sequence")
            obs2 = get_child(v2, "Observation_Sequence")
            
            if obs1 is None and obs2 is None:
                continue
            if (obs1 is None) != (obs2 is None):
                print(f"  MISMATCH at Visit {id1}: One has ObsSeq, other doesn't.")
                mismatches += 1
                continue

            p1 = get_child(obs1, "Observational_Parameters")
            p2 = get_child(obs2, "Observational_Parameters")
            
            t1 = get_text(p1, "Target")
            t2 = get_text(p2, "Target")
            
            if t1 != t2:
                print(f"  MISMATCH at Visit {id1}: Target {t1} vs {t2}")
                mismatches += 1
            
            time1 = get_child(p1, "Timing")
            time2 = get_child(p2, "Timing")
            start1 = get_text(time1, "Start")
            start2 = get_text(time2, "Start")
            stop1 = get_text(time1, "Stop")
            stop2 = get_text(time2, "Stop")
            
            if start1 != start2 or stop1 != stop2:
                 print(f"  MISMATCH at Visit {id1}: Timing ({start1}-{stop1}) vs ({start2}-{stop2})")
                 mismatches += 1

        if mismatches == 0:
            print("  MATCH: All visits match in key parameters (ID, Target, Timing).")
            return True
        else:
            print(f"  Found {mismatches} mismatches in Visits.")
            return False

    except Exception as e:
        print(f"  ERROR comparing XML: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    legacy_dir = Path("src/pandorascheduler/data/baseline")
    rework_dir = Path("output_standalone")
    legacy_data_dir = Path("src/pandorascheduler/data")

    # 1. Compare Schedule
    legacy_schedule = list(legacy_dir.glob("Pandora_Schedule_*.csv"))[0]
    rework_schedule = list(rework_dir.glob("Pandora_Schedule_*.csv"))[0]
    compare_csvs(legacy_schedule, rework_schedule)

    # 2. Compare Tracker
    compare_csvs(legacy_dir / "tracker.csv", rework_dir / "tracker.csv")

    # 3. Compare Science Calendar
    legacy_xml = legacy_data_dir / "Pandora_science_calendar.xml"
    rework_xml = rework_dir / "Pandora_science_calendar.xml"
    
    if legacy_xml.exists() and rework_xml.exists():
        compare_xml_deep(legacy_xml, rework_xml)
    else:
        print(f"Skipping XML comparison: Files not found ({legacy_xml}, {rework_xml})")
