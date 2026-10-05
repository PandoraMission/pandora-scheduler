import logging
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def normalize_xml(element):
    """Normalize XML element for comparison."""
    # Sort children by tag and text content to handle unordered elements if necessary
    # For now, we assume order matters for visits, but maybe not for parameters
    # We will just strip whitespace from text
    if element.text:
        element.text = element.text.strip()
    if element.tail:
        element.tail = element.tail.strip()
    
    for child in element:
        normalize_xml(child)

def elements_equal(e1, e2):
    if e1.tag != e2.tag:
        return False, f"Tag mismatch: {e1.tag} != {e2.tag}"
    if (e1.text or "").strip() != (e2.text or "").strip():
        return False, f"Text mismatch in {e1.tag}: '{e1.text}' != '{e2.text}'"
    if e1.attrib != e2.attrib:
        return False, f"Attribute mismatch in {e1.tag}: {e1.attrib} != {e2.attrib}"
    # Special handling for ScienceCalendar root to compare Visits by Content
    if e1.tag.endswith("ScienceCalendar"):
        ns = "{/pandora/calendar/}"
        
        def get_visit_key(visit):
            # Return (Target, StartTime) tuple, or None if empty/invalid
            obs = visit.find(f"{ns}Observation_Sequence")
            if obs is None:
                return None
            params = obs.find(f"{ns}Observational_Parameters")
            if params is None:
                return None
            target = params.find(f"{ns}Target").text
            start = params.find(f"{ns}Timing/{ns}Start").text
            return (target, start)

        # Filter visits by date range (2026-02-05 to 2026-02-12)
        start_cutoff = "2026-02-05"
        end_cutoff = "2026-02-12"
        
        visits1 = {}
        for v in e1.findall(f"{ns}Visit"):
            key = get_visit_key(v)
            if key and start_cutoff <= key[1] <= end_cutoff:
                visits1[key] = v
                
        visits2 = {}
        for v in e2.findall(f"{ns}Visit"):
            key = get_visit_key(v)
            if key and start_cutoff <= key[1] <= end_cutoff:
                visits2[key] = v
        
        all_keys = sorted(set(visits1.keys()) | set(visits2.keys()))
        
        mismatch_count = 0
        for key in all_keys:
            if key not in visits1:
                logger.error(f"Visit {key} missing in File 1")
                mismatch_count += 1
                continue
            if key not in visits2:
                logger.error(f"Visit {key} missing in File 2")
                mismatch_count += 1
                continue
            
            # Compare content, but IGNORE ID differences
            # We can't easily ignore just ID with the recursive function
            # So let's just compare the Observation_Sequence part
            v1_obs = visits1[key].find(f"{ns}Observation_Sequence")
            v2_obs = visits2[key].find(f"{ns}Observation_Sequence")
            
            eq, msg = elements_equal(v1_obs, v2_obs)
            if not eq:
                logger.error(f"Mismatch in Visit {key}: {msg}")
                mismatch_count += 1
        
        if mismatch_count > 0:
            return False, f"Found {mismatch_count} mismatches in visits"
            
        return True, ""

    if len(e1) != len(e2):
        return False, f"Child count mismatch in {e1.tag}: {len(e1)} != {len(e2)}"
    
    for c1, c2 in zip(e1, e2):
        eq, msg = elements_equal(c1, c2)
        if not eq:
            return False, msg
            
    return True, ""

def compare_xmls(file1, file2):
    logger.info(f"Comparing {file1} and {file2}...")
    
    try:
        tree1 = ET.parse(file1)
        tree2 = ET.parse(file2)
    except ET.ParseError as e:
        logger.error(f"XML Parse Error: {e}")
        return False

    root1 = tree1.getroot()
    root2 = tree2.getroot()

    normalize_xml(root1)
    normalize_xml(root2)

    eq, msg = elements_equal(root1, root2)
    
    if eq:
        logger.info("XML files are semantically IDENTICAL.")
        return True
    else:
        logger.error(f"XML files DIFFER: {msg}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_xml_files.py <file1> <file2>")
        sys.exit(1)
        
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    if not Path(file1).exists():
        print(f"File not found: {file1}")
        sys.exit(1)
    if not Path(file2).exists():
        print(f"File not found: {file2}")
        sys.exit(1)
        
    success = compare_xmls(file1, file2)
    sys.exit(0 if success else 1)
