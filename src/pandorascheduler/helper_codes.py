import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import json
import pandas as pd
from scheduler import PACKAGEDIR

def general_parameters(obs_sequence_duration = 90, occ_sequence_limit = 30):
    observation_sequence_duration = obs_sequence_duration # minutes
    occultation_sequence_limit = occ_sequence_limit # minutes
    return observation_sequence_duration, occultation_sequence_limit

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

# def no_phase_event():
#     return

def read_json_files(targ_list, fn_tmp):
    import pandas as pd
    import numpy as np
    target_list_copy = targ_list.copy()
    with open(fn_tmp, 'r') as file:
        data = json.load(file)

    # Iterate through the key-value pairs in the JSON data
        for key, value in data.items():
            # Convert lists or arrays to strings
            if isinstance(value, (list, np.ndarray)):
                value = str(value)

            # Check if the column exists in the DataFrame
            if key not in target_list_copy.columns:
                # If it doesn't exist, add it as a new column
                target_list_copy[key] = np.nan
            
            # Check if the value already exists in the column
            if value not in target_list_copy[key].values:
                # Find the first NaN value in the column and replace it
                nan_indices = target_list_copy[key].isna()
                if nan_indices.any():
                    nan_index = nan_indices.idxmax()
                    target_list_copy.at[nan_index, key] = value
                else:
                    # If no NaN values, append a new row
                    new_row = pd.DataFrame({key: [value]})
                    target_list_copy = pd.concat([target_list_copy, new_row], ignore_index=True)

        old_column_name = "Transit Epoch (BJD_TDB-ZZZZZ)"
        column_index = target_list_copy.columns.get_loc(old_column_name)
        new_column_name = "Transit Epoch (BJD_TDB-2400000.5)"
        if old_column_name in target_list_copy.columns:
            target_list_copy[old_column_name] = target_list_copy[old_column_name] - 2400000.5
            # print(f"Column '{old_column_name}' has been updated.")
            target_list = target_list_copy.rename(columns={old_column_name: new_column_name})
        else:
            target_list = target_list_copy

        # targ_list_copy.loc[0, "Transit Duration (hrs)"] = data["pl_trandur (hrs)"]
    return target_list

def update_target_list(targ_list, pl_names):
    import os
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    filtered_targ_list = targ_list[targ_list["Planet Name"].isin(pl_names)]
    updated_targ_list = filtered_targ_list.copy()

    for pl_name in pl_names:#targ_list["Planet Name"]:
        fn_tmp = f'{PACKAGEDIR}/data/target_json_files/' + pl_name + '.json'
        if os.path.exists(fn_tmp):

            tmp_arr = read_json_files(filtered_targ_list[filtered_targ_list["Planet Name"] == pl_name], fn_tmp)
            # filtered_targ_list_copy = filtered_targ_list.copy()
            # Ensure filtered_targ_list has all columns from tmp_arr
            for col in tmp_arr.columns:
                if col not in updated_targ_list.columns:
                    updated_targ_list[col] = None
            # Update filtered_targ_list with tmp_arr data
            updated_targ_list.update(tmp_arr)
        # else:
        #     print(f"The JSON file for '{pl_name}' does not exist.")
    return updated_targ_list

# def load_visibility_data(name):
#     try:
#         vis = pd.read_csv(f"{PACKAGEDIR}/data/aux_targets/{name}/Visibility for {name}.csv")
#         vis['Time(MJD_UTC)'] = pd.to_numeric(vis['Time(MJD_UTC)'])
#         return name, vis
#     except FileNotFoundError:
#         return name, None

# def process_target(data, start, stop):
#     from astropy.time import Time
#     name, vis = data
#     if vis is None:
#         return None
    
#     mask = (vis['Time(MJD_UTC)'] >= Time(start).mjd) & (vis['Time(MJD_UTC)'] <= Time(stop).mjd)
#     vis_filtered = vis.loc[mask]
    
#     if vis_filtered['Visible'].all():
#         return name, 100, True
#     elif vis_filtered['Visible'].any():
#         visibility_percentage = 100 * vis_filtered['Visible'].mean()
#         return name, visibility_percentage, False
#     return None

def find_first_visible_target(start, stop, names):
    from tqdm import tqdm
    from astropy.time import Time
    # for n, name in tqdm(enumerate(names), desc=f"Finding visible aux target for {start} to {stop}", total=len(names)):
    for n in tqdm(range(len(names)), desc="Finding visible aux target for " + str(start) + ' to ' + str(stop)):
        try:
            vis_file = f"{PACKAGEDIR}/data/aux_targets/{names[n]}/Visibility for {names[n]}.csv"
            
            # Read only the necessary columns
            vis = pd.read_csv(vis_file, usecols=["Time(MJD_UTC)", "Visible"])
            
            # Convert to Time object and filter in one step
            time_mask = (Time(vis["Time(MJD_UTC)"], format='mjd', scale='utc') >= start) & \
                        (Time(vis["Time(MJD_UTC)"], format='mjd', scale='utc') <= stop)
            
            vis_filtered = vis[time_mask]
            
            if not vis_filtered.empty and vis_filtered['Visible'].all():
                print(f'VK use the 1st aux target that is 100% visible; start = {start}; stop = {stop}, {n}')
                return n, 100.0  # Return index and visibility percentage
            
            elif not vis_filtered.empty and vis_filtered['Visible'].any():
                visibility_percentage = 100 * (vis_filtered['Visible'].sum() / len(vis_filtered))
                return n, visibility_percentage
            
        except FileNotFoundError:
            continue
    
    return None, 0.0  # If no suitable target found


def find_visible_targets(names, start, stop):
    # Load all visibility data
    with Pool() as pool:
        all_data = dict(pool.map(helper_codes.load_visibility_data, names))
    
    # Process targets
    process_func = partial(helper_codes.process_target, start=start, stop=stop)
    with Pool() as pool:
        results = pool.map(process_func, all_data.items())
    
    vis_all_targs = []
    vis_any_targs = []
    targ_vis = []
    
    for result in results:
        if result:
            name, visibility, fully_visible = result
            if fully_visible:
                vis_all_targs.append(name)
                break
            else:
                vis_any_targs.append(name)
                targ_vis.append(visibility)
    
    return vis_all_targs, vis_any_targs, targ_vis
