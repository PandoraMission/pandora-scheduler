import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import json
import pandas as pd
import numpy as np
from astropy.time import Time
from tqdm import tqdm
import os

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

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

    return o_seq

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

def process_visibility(v_names, starts, stops, path):

    results = []

    for n, v_name in enumerate(tqdm(v_names, desc="Processing targets")):

        vis = pd.read_csv(f"{path}/{v_name}/Visibility for {v_name}.csv")
        vis_times = vis['Time(MJD_UTC)']
        visibility = vis['Visible']
        
        valid_intervals = []
        
        for s, (start, stop) in enumerate(zip(starts, stops)):
            # Check if the target is visible for the entire interval
            interval_mask = (vis_times >= start) & (vis_times <= stop)
            if np.all(visibility[interval_mask] == 1):
                valid_intervals.append((start, stop, 1))
            else:
                valid_intervals.append((start, stop, 0))
        
        if valid_intervals:
            results.append({
                'v_name': v_name,
                'valid_intervals': valid_intervals,
                'visibility': visibility[interval_mask].tolist()
            })
    
    return results

def create_visibility_dataframe(results):
    # Create a dictionary to store data for the DataFrame
    data = {}
    # Iterate through all results
    for result in results:
        target_name = result['v_name']
        visibilities = []
        start_times = []
        stop_times = []
        for start, stop, visibility in result['valid_intervals']:
            # Convert to ISO format and remove milliseconds
            start_iso = Time(start, format='mjd', scale='utc').iso.split('.')[0]
            stop_iso = Time(stop, format='mjd', scale='utc').iso.split('.')[0]
            visibilities.append(visibility)
            start_times.append(start_iso)
            stop_times.append(stop_iso)
        # Add this target's data to the dictionary
        data[target_name] = {
            'Start': start_times,
            'Stop': stop_times,
            'Visibility': visibilities
        }
    
    # Create the DataFrame
    df = pd.DataFrame([(target, start, stop, vis) 
                       for target, d in data.items() 
                       for start, stop, vis in zip(d['Start'], d['Stop'], d['Visibility'])],
                      columns=['Target', 'Start', 'Stop', 'Visibility'])
    
    # Sort by Start time
    # df = df.sort_values('Target')
    
    # Set Target as index
    df = df.set_index('Target')
    
    return df

def schedule_observations(visibility_data, max_targets=3):
    if not isinstance(visibility_data, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame")

    # Create a pivot table to get the visibility matrix
    visibility_matrix = visibility_data.pivot_table(index='Target', 
                                                    columns='Start', 
                                                    values='Visibility', 
                                                    aggfunc='first')
    visibility_matrix = visibility_matrix.fillna(0)  # Fill NaN with 0 (not visible)
    
    visibility_array = visibility_matrix.values
    target_names = visibility_matrix.index
    time_index = visibility_matrix.columns

    num_targets, num_intervals = visibility_array.shape

    # Calculate total visibility for each target
    total_visibility = np.sum(visibility_array, axis=1)

    # Sort targets by total visibility
    sorted_indices = np.argsort(total_visibility)[::-1]

    schedule = pd.DataFrame(index=time_index, columns=['Stop', 'Target'], dtype='object')
    used_targets = 0

    for target_idx in sorted_indices:
        if used_targets >= max_targets:
            break

        target_visibility = visibility_array[target_idx]

        # Schedule this target where it's visible and not already scheduled
        new_schedules = schedule['Target'].isnull() & (target_visibility == 1)
        if new_schedules.any():
            schedule.loc[new_schedules, 'Target'] = target_names[target_idx]
            used_targets += 1

        # Check if all intervals are now scheduled
        if not schedule['Target'].isnull().any():
            break

    # Fill in the Stop times
    stop_times = visibility_data.groupby('Start')['Stop'].first()
    schedule['Stop'] = schedule.index.map(stop_times.to_dict())

    # schedule['Target'] = schedule['Target'].fillna('No target')
    schedule['Target'] = schedule['Target'].fillna(np.nan)

    schedule_reset = schedule.reset_index()

    return schedule_reset, used_targets

def add_all_targets_visibility(schedule_array, vis_df):
    def get_target_visibility(start_time, target):
        if target in vis_df.index:
            mask = (vis_df.loc[target]['Start'] == start_time)
            if mask.any():
                return vis_df.loc[target][mask]['Visibility'].values[0]
        return np.nan  # Return NaN if no matching visibility found

    # Convert schedule array to DataFrame
    df_schedule = pd.DataFrame(schedule_array, columns=['Start', 'Stop', 'Target'])
    
    # Get unique targets from the schedule
    unique_targets = df_schedule['Target'].unique()
    
    # Remove 'No target' if present
    unique_targets = [target for target in unique_targets if target != 'No target']
    
    # Add Visibility columns for each target
    for target in unique_targets:
        column_name = f'{target} Visibility'
        df_schedule[column_name] = df_schedule['Start'].apply(lambda x: get_target_visibility(x, target))
    
    # Convert back to numpy array
    new_schedule_array = df_schedule.values
    
    return new_schedule_array, list(df_schedule.columns)

def add_random_targets_visibility(schedule_array, vis_df, target_names):
    def get_target_visibility(start_time, target):
        if target in vis_df.index:
            mask = (vis_df.loc[target]['Start'] == start_time)
            if mask.any():
                return vis_df.loc[target][mask]['Visibility'].values[0]
        return np.nan  # Return NaN if no matching visibility found

    # Convert schedule array to DataFrame
    df_schedule = pd.DataFrame(schedule_array, columns=['Start', 'Stop', 'Target'])
    
    # Add Visibility columns for each target
    for target in target_names:
        column_name = f'{target} Visibility'
        df_schedule[column_name] = df_schedule['Start'].apply(lambda x: get_target_visibility(x, target))
    
    # Convert back to numpy array
    new_schedule_array = df_schedule.values
    
    return new_schedule_array

def schedule_occultation_targets_all(v_names, starts, stops, path, o_df, o_list):

    results = process_visibility(v_names, starts, stops, path)

    vis_df = create_visibility_dataframe(results)
    # print("Targets:", vis_df.index.tolist())
    # print("Start Times:", vis_df.columns.tolist())
    # print(vis_df.loc["GJ 1214"])
    # visibility_array = vis_df.values

    schedule, num_targets_used = schedule_observations(vis_df, max_targets=3)

    # scheduled_observations = schedule[schedule['Target'].notna()] 
    unscheduled_times = schedule[schedule['Target'].isna()]

    d_flag = False
    if len(unscheduled_times) == 0:
        # print('All occultation times covered')
        for s in range(len(starts)):
            idx, = np.where(o_list["Star Name"] == schedule['Target'][s])
            if len(idx) > 0:
                o_df['Target'][s] = schedule['Target'][s]
                o_df['RA'][s] = o_list["RA"][idx].values[0]
                o_df['DEC'][s] = o_list["DEC"][idx].values[0]
        d_flag = True
    # else:
    #     print('Missing occultation times')
    #     print(unscheduled_times)

    # new_schedule, column_names = add_all_targets_visibility(schedule, vis_df)
    # print("Column names:", column_names)
    # print(new_schedule)
    # print()

    # targets = ["GJ 1214", "HIP 65 A", "WASP-107"]  # Add all target names you want to include
    # new_schedule = add_random_targets_visibility(schedule, vis_df, targets)
    # print(new_schedule)

    return o_df, d_flag

def schedule_occultation_targets(v_names, starts, stops, path, o_df, o_list, try_occ_targets):#, position):
    schedule = pd.DataFrame(index=starts, columns=['Stop', 'Target', 'Visibility'], dtype='object')
    schedule['Stop'] = stops
    schedule['Target'] = np.nan
    schedule['Visibility'] = np.nan

    # Add 'Visibility' column to o_df if it doesn't exist
    if 'Visibility' not in o_df.columns:
        o_df['Visibility'] = np.nan

    for v_name in tqdm(v_names, desc=f"Finding visible occultation target from {try_occ_targets}", leave = False):#, position=position):#, leave=leave):#, leave=(position != 0)):#desc="Processing targets"):
    # for v_name in v_names:
        # Process visibility for this target
        vis = pd.read_csv(f"{path}/{v_name}/Visibility for {v_name}.csv")
        vis_times = vis['Time(MJD_UTC)']
        visibility = vis['Visible']

        for s, (start, stop) in enumerate(zip(starts, stops)):
            if pd.isna(schedule.loc[start, 'Target']):
                # Check if the target is visible for the entire interval
                interval_mask = (vis_times >= start) & (vis_times <= stop)

                # print(Time(start, format='mjd').datetime.strftime("%Y-%m-%dT%H:%M:%SZ"), \
                #     Time(stop, format='mjd').datetime.strftime("%Y-%m-%dT%H:%M:%SZ"), v_name, visibility[interval_mask].values.astype(int))
                
                if np.all(visibility[interval_mask] == 1):
                    schedule.loc[start, 'Target'] = v_name
                    schedule.loc[start, 'Visibility'] = 1
                    
                    # Update o_df
                    idx, = np.where(o_list["Star Name"] == v_name)
                    if len(idx) > 0:
                        o_df.loc[s, 'Target'] = v_name
                        o_df.loc[s, 'RA'] = o_list.loc[idx[0], "RA"]
                        o_df.loc[s, 'DEC'] = o_list.loc[idx[0], "DEC"]
                        o_df.loc[s, 'Visibility'] = 1
                else:
                    # If the target is not visible for the entire interval, mark it as 0
                    if pd.isna(schedule.loc[start, 'Visibility']):
                        schedule.loc[start, 'Visibility'] = 0
                        o_df.loc[s, 'Visibility'] = 0

        # print(v_name, o_df)

        # Check if schedule is completely filled
        if not schedule['Target'].isna().any():
            return o_df, True

    # If we've gone through all targets and still have empty slots
    # Fill remaining slots with 'No target' and Visibility 0
    mask = schedule['Target'].isna()
    schedule.loc[mask, 'Target'] = 'No target'
    schedule.loc[mask, 'Visibility'] = 0
    
    o_df.loc[o_df['Target'].isna(), 'Target'] = 'No target'
    o_df.loc[o_df['Visibility'].isna(), 'Visibility'] = 0

    return o_df, False

def break_long_visibility_changes(arr, max_sequence = 90):
    result = [arr[0]]
    for i in range(1, len(arr)):
        diff = arr[i] - arr[i-1]
        if diff > max_sequence:
            # Calculate how many elements to insert
            num_to_insert = diff // max_sequence
            step = diff / (num_to_insert + 1)
            for j in range(1, num_to_insert + 1):
                result.append(int(arr[i-1] + j * step))
        result.append(arr[i])
    return np.array(result)

def print_element_from_xml(elem, level=0):
    print("  " * level + f"{elem.tag}: {elem.text.strip() if elem.text else ''}")
    for child in elem:
        print_element_from_xml(child, level + 1)
