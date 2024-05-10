import xml.etree.ElementTree as ET
from xml.dom import minidom
import os
from tqdm import tqdm

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

calendar = ET.Element('Calendar', xmlns="/pandora/calendar/")
# meta=ET.SubElement(calendar, 'Meta',
#                    Calendar_Weights='0.0, 0.0, 1.0',
#                    Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
#                    Keepout_Angles='90.0, 40.0, 63.0',
#                    Author="P Bonney",
#                    Delivery_Id='',
#                    )


t_name_ = ['T1', "T2", "T3"]
priority_ = [0, 0, 0]
START_ = ["0", "0", "0"]
STOP_ = ["999", "999", "999"]

for ii in tqdm(range(3)):
    t_name = t_name_[ii]
    pr = priority_[ii]
    START, STOP = START_[ii], STOP_[ii]


    params = {
        "ID_1": str(ii),
        "ID_2": "N2_value",
        "Target": t_name,
        "Priority": f'{pr}',
        "Start": START,
        "Stop": STOP,
        "ra_1": "ra_1_value",
        "ra_2": "ra_2_value",
        "number_3": "number_3_value",
        "Target_Name": "Name Value ZZZ"
    }
    params_NIRDA = {
        "AverageGroups": "1", 
        "ROI_StartX": "0", 
        "ROI_StartY": "824", 
        "ROI_SizeX": "80", 
        "ROI_SizeY": "400", 
        "SC_Resets1": "1", 
        "SC_Resets2": "1", 
        "SC_DropFrames1": "0", 
        "SC_DropFrames2": "16", 
        "SC_DropFrames3": "0", 
        "SC_ReadFrames": "4", 
        "TargetID": t_name, 
        "SC_Groups": "2", 
        "SC_Integrations": "525", 
    }

    params_VDA = {
        "StartRoiDetMethod": 0,
        "FramesPerCoadd": 50,
        "NumTotalFramesRequested": 9000,
        "TargetRA": 60.000000,
        "TargetDEC": 30.000000,
        "IncludeFieldSolnsInResp": 1,
        "StarRoiDimension": 50,
        "MaxNumStarRois": 0,
        "numPredefinedStarRois": 5,
        "PredefinedStarRoiRa": [60.1, 60.2, 60.3, 60.4, 60.5], 
        "PredefinedStarRoiDec": [-30.1, -30.2, -30.3, -30.4, -30.5],
        "TargetID": t_name,
        "NumExposuresMax": 1,
        "ExposureTime_us": 200000,
    }

    # Add child elements, populating with parameters from the 'params' dictionary
    visit = ET.SubElement(calendar, "Visit")

    id0 = ET.SubElement(visit, "ID")
    id0.text = params["ID_1"]

    observation_sequence = ET.SubElement(visit, "Observation_Sequence")

    id1 = ET.SubElement(observation_sequence, "ID")
    id1.text = params["ID_2"]

    target = ET.SubElement(observation_sequence, "Target")
    target.text = params["Target"]

    priority = ET.SubElement(observation_sequence, "Priority")
    priority.text = params["Priority"]

    timing = ET.SubElement(observation_sequence, "Timing")
    start = ET.SubElement(timing, "Start")
    start.text = params["Start"]
    stop = ET.SubElement(timing, "Stop")
    stop.text = params["Stop"]

    boresight = ET.SubElement(observation_sequence, "Boresight")
    ra = ET.SubElement(boresight, "RA")
    ra.text = params["ra_1"]
    dec = ET.SubElement(boresight, "Dec")
    dec.text = params["ra_2"]

    payload_parameters = ET.SubElement(observation_sequence, "Payload_Parameters")

    ### NIRDA Parameters
    nirda = ET.SubElement(payload_parameters, "NIRDA")

    average_groups = ET.SubElement(nirda, "AverageGroups")
    average_groups.text = str(params_NIRDA["AverageGroups"])

    ROI_StartX = ET.SubElement(nirda, "ROI_StartX")
    ROI_StartX.text = str(params_NIRDA["ROI_StartX"])

    ROI_StartY = ET.SubElement(nirda, "ROI_StartY")
    ROI_StartY.text = str(params_NIRDA["ROI_StartY"])

    ROI_SizeX = ET.SubElement(nirda, "ROI_SizeX")
    ROI_SizeX.text = str(params_NIRDA["ROI_SizeX"])

    ROI_SizeY = ET.SubElement(nirda, "ROI_SizeY")
    ROI_SizeY.text = str(params_NIRDA["ROI_SizeY"])

    SC_Resets1 = ET.SubElement(nirda, "SC_Resets1")
    SC_Resets1.text = str(params_NIRDA["SC_Resets1"])

    SC_Resets2 = ET.SubElement(nirda, "SC_Resets2")
    SC_Resets2.text = str(params_NIRDA["SC_Resets2"])

    SC_DropFrames1 = ET.SubElement(nirda, "SC_DropFrames1")
    SC_DropFrames1.text = str(params_NIRDA["SC_DropFrames1"])

    SC_DropFrames2 = ET.SubElement(nirda, "SC_DropFrames2")
    SC_DropFrames2.text = str(params_NIRDA["SC_DropFrames2"])

    SC_DropFrames3 = ET.SubElement(nirda, "SC_DropFrames3")
    SC_DropFrames3.text = str(params_NIRDA["SC_DropFrames3"])

    SC_ReadFrames = ET.SubElement(nirda, "SC_ReadFrames")
    SC_ReadFrames.text = str(params_NIRDA["SC_ReadFrames"])

    TargetID = ET.SubElement(nirda, "TargetID")
    TargetID.text = str(params_NIRDA["TargetID"])

    SC_Groups = ET.SubElement(nirda, "SC_Groups")
    SC_Groups.text = str(params_NIRDA["SC_Groups"])

    SC_Integrations = ET.SubElement(nirda, "SC_Integrations")
    SC_Integrations.text = str(params_NIRDA["SC_Integrations"])

    ### VDA Parameters:
    vda = ET.SubElement(payload_parameters, "VDA")
    StartRoiDetMethod = ET.SubElement(vda, "StartRoiDetMethod")
    StartRoiDetMethod.text = str(params_VDA["StartRoiDetMethod"])

    FramesPerCoadd = ET.SubElement(vda, "FramesPerCoadd")
    FramesPerCoadd.text = str(params_VDA["FramesPerCoadd"])

    NumTotalFramesRequested = ET.SubElement(vda, "NumTotalFramesRequested")
    NumTotalFramesRequested.text = str(params_VDA["NumTotalFramesRequested"])

    TargetRA = ET.SubElement(vda, "TargetRA")
    TargetRA.text = str(params_VDA["TargetRA"])

    TargetDEC = ET.SubElement(vda, "TargetDEC")
    TargetDEC.text = str(params_VDA["TargetDEC"])

    IncludeFieldSolnsInResp = ET.SubElement(vda, "IncludeFieldSolnsInResp")
    IncludeFieldSolnsInResp.text = str(params_VDA["IncludeFieldSolnsInResp"])

    StarRoiDimension = ET.SubElement(vda, "StarRoiDimension")
    StarRoiDimension.text = str(params_VDA["StarRoiDimension"])

    MaxNumStarRois = ET.SubElement(vda, "MaxNumStarRois")
    MaxNumStarRois.text = str(params_VDA["MaxNumStarRois"])

    numPredefinedStarRois = ET.SubElement(vda, "numPredefinedStarRois")
    numPredefinedStarRois.text = str(params_VDA["numPredefinedStarRois"])

    PredefinedStarRoiRa = ET.SubElement(vda, "PredefinedStarRoiRa")
    PredefinedStarRoiRa.text = str(params_VDA["PredefinedStarRoiRa"])

    PredefinedStarRoiDec = ET.SubElement(vda, "PredefinedStarRoiDec")
    PredefinedStarRoiDec.text = str(params_VDA["PredefinedStarRoiDec"])

    TargetID = ET.SubElement(vda, "TargetID")
    TargetID.text = str(params_VDA["TargetID"])

    NumExposuresMax = ET.SubElement(vda, "NumExposuresMax")
    NumExposuresMax.text = str(params_VDA["NumExposuresMax"])

    ExposureTime_us = ET.SubElement(vda, "ExposureTime_us")
    ExposureTime_us.text = str(params_VDA["ExposureTime_us"])

etstr=ET.tostring(calendar, xml_declaration=True)

from xml.dom import minidom
dom = minidom.parseString(etstr)

#dom = xml.dom.minidom.parseString(etstr)
pretty_xml_as_string = dom.toprettyxml()
f=open(f'{PACKAGEDIR}/data/test_cal_pretty.xml', 'w+')
f.write(pretty_xml_as_string)
f.close()


