import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

def general_parameters():
    observation_sequence_duration = 90 # minutes
    return observation_sequence_duration

def observation_sequence(visit, obs_seq_ID, t_name, priority, start, stop, ra, dec):

    o_seq = ET.SubElement(visit,'Observation_Sequence')
    obs_seq_id = ET.SubElement(o_seq, "ID")
    obs_seq_id.text = obs_seq_ID

    observational_parameters, params_NIRDA, params_VDA = params_obs_NIRDA_VDA(t_name, priority, start, stop, ra, dec)

    obs_parameters = ET.SubElement(o_seq, "Observational_Parameters")
    ### Observational Parameters
    for ii, jj in zip(observational_parameters.keys(), observational_parameters.values()):
        if (ii != "Timing") & (ii != "Boresight"):
            obs_param_element = ET.SubElement(obs_parameters, ii)
            obs_param_element.text = str(jj)
        else:
            obs_param_element = ET.SubElement(obs_parameters, ii)
            for kk in range(2):
                sub_element_tmp = ET.SubElement(obs_param_element, jj[kk])
                sub_element_tmp.text = jj[kk+2]

    ### Payload Parameters
    payload_parameters = ET.SubElement(o_seq, "Payload_Parameters")
    ### NIRDA Parameters
    nirda = ET.SubElement(payload_parameters, "NIRDA")
    for nirda_key, nirda_values in zip(params_NIRDA.keys(), params_NIRDA.values()):
        nirda_subelement_ = ET.SubElement(nirda, nirda_key)
        nirda_subelement_.text = nirda_values
    ### VDA Parameters:
    vda = ET.SubElement(payload_parameters, "VDA")
    for vda_key, vda_values in zip(params_VDA.keys(), params_VDA.values()):
        vda_subelement_ = ET.SubElement(vda, vda_key)
        vda_subelement_.text = str(vda_values)

    return

def params_obs_NIRDA_VDA(t_name, priority, start, stop, ra, dec):

    try:
        start_format, stop_format = f'{datetime.strftime(start, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(stop, "%Y-%m-%dT%H:%M:%SZ")}'
    except:
        start_format, stop_format = start, stop

    try:
        ra_tmp, dec_tmp = f'{float(ra)}', f'{float(dec)}'
    except:
        ra_tmp, dec_tmp = f'{float(-999)}', f'{float(-999)}'

    observational_parameters = {
        "Target": t_name,
        "Priority": f'{priority}',
        "Timing": ["Start", "Stop", start_format, stop_format],#f'{datetime.strftime(start, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(stop, "%Y-%m-%dT%H:%M:%SZ")}'], #f'{start}', f'{stop}'],#
        "Boresight": ["RA", "DEC", ra_tmp, dec_tmp],#f'{float(ra)}', f'{float(dec)}'], 
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
        "TargetRA": ra_tmp,#f'{float(ra)}',
        "TargetDEC": dec_tmp,#f'{float(dec)}',
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

    return observational_parameters, params_NIRDA, params_VDA

def remove_short_sequences(array, sequence_too_short):
    A = array#[1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]

    # Initialize variables to keep track of the start index of a subarray and the positions of qualifying subarrays
    start_index = None
    positions = []

    # Iterate through the array, keeping track of the start and end of each subarray of 1s
    for i in range(len(A)):
        # Check if we are at the start of a subarray of 1s
        if A[i] == 1 and start_index is None:
            start_index = i
        # Check if we are at the end of a subarray of 1s
        elif A[i] == 0 and start_index is not None:
            # Check if the subarray is shorter than 2 elements
            if i - start_index < sequence_too_short:
                positions.append((start_index, i - 1))
            start_index = None
    # Check if the last element of A is part of a qualifying subarray
    if start_index is not None and len(A) - start_index < sequence_too_short:
        positions.append((start_index, len(A) - 1))

    A_new = A.copy()
    if len(positions) != 0:
        for ii in range(len(positions)):
            A_new[positions[ii][0]:positions[ii][1]+1] = 0.

    return A_new, positions
#
#
def break_long_sequences(start, end, step):
    ranges = []
    current = start
    while current < end:
        next_val = min(current + step, end)
        ranges.append([current, next_val])
        current += step
    return ranges
