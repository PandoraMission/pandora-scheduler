import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import json
import pandas as pd
import numpy as np
from astropy.time import Time
from tqdm import tqdm
import os
import re

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

def general_parameters(obs_sequence_duration = 90, occ_sequence_limit = 30):
    observation_sequence_duration = obs_sequence_duration # minutes
    occultation_sequence_limit = occ_sequence_limit # minutes
    return observation_sequence_duration, occultation_sequence_limit

def observation_sequence(visit, obs_seq_ID, t_name, priority, start, stop, ra, dec, targ_info):

    import logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    o_seq = ET.SubElement(visit,'Observation_Sequence')
    obs_seq_id = ET.SubElement(o_seq, "ID")
    obs_seq_id.text = obs_seq_ID

    # observational_parameters, params_NIRDA, params_VDA = params_obs_NIRDA_VDA(t_name, priority, start, stop, ra, dec)
    observational_parameters = params_obs_NIRDA_VDA(t_name, priority, start, stop, ra, dec)

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
    ### Calculate integration duration in seconds, from "Start" to "Stop"
    if isinstance(stop, str):
        diff_in_sec = (datetime.strptime(stop, '%Y-%m-%dT%H:%M:%SZ') - datetime.strptime(start, '%Y-%m-%dT%H:%M:%SZ')).total_seconds()
    elif isinstance(stop, datetime):
        diff_in_sec = (stop - start).total_seconds()

    ### Payload Parameters
    payload_parameters = ET.SubElement(o_seq, "Payload_Parameters")
    ### NIRDA Parameters
    nirda = ET.SubElement(payload_parameters, "AcquireInfCamImages")
    nirda_columns = targ_info.columns[targ_info.columns.str.startswith('NIRDA_')]
    columns_to_ignore = ['IncludeFieldSolnsInResp', 'NIRDA_TargetID', 'NIRDA_SC_Integrations', 'NIRDA_FramesPerIntegration', 'NIRDA_IntegrationTime_s']
    for nirda_key, nirda_values in targ_info[nirda_columns].iloc[0].items():
    # for nirda_key, nirda_values in zip(params_NIRDA.keys(), params_NIRDA.values()):
        if pd.notna(nirda_values):  # This condition checks if the value is not NaN
            xml_key = nirda_key.replace('NIRDA_', '')
            if (nirda_key not in columns_to_ignore):# and (nirda_key != 'NIRDA_SC_Integrations'):
                nirda_subelement_ = ET.SubElement(nirda, xml_key)
                nirda_subelement_.text = str(nirda_values)
            elif nirda_key == 'NIRDA_TargetID':
                if pd.isnull(targ_info['Planet Name'].iloc[0]):
                    tmp_t_name = targ_info['Star Name'].iloc[0]
                else:
                    tmp_t_name = targ_info['Planet Name'].str.replace(r'\s+([bcd])$', r'\1', regex=True).iloc[0]
                nirda_subelement_ = ET.SubElement(nirda, xml_key)
                nirda_subelement_.text = tmp_t_name#targ_info['Planet Name'].iloc[0]
            elif nirda_key == 'NIRDA_SC_Integrations':
                nirda_subelement_ = ET.SubElement(nirda, xml_key)
                nirda_subelement_.text = str(np.round(diff_in_sec/targ_info['NIRDA_IntegrationTime_s'].iloc[0]).astype(int))
            pass
        # else:
        #     logging.info(f"Searching for occultation targets from {st} to {sp}")

    ### VDA Parameters:
    vda = ET.SubElement(payload_parameters, "AcquireVisCamScienceData")
    vda_columns = targ_info.columns[targ_info.columns.str.startswith('VDA_')]
    # columns_to_ignore = ['VDA_IntegrationTime']
    columns_to_ignore = ['VDA_NumExposuresMax', 'VDA_NumTotalFramesRequested', 'VDA_TargetID', 'VDA_TargetRA', 'VDA_TargetDEC', \
        'VDA_StarRoiDetMethod', 'VDA_numPredefinedStarRois', 'VDA_PredefinedStarRoiRa', 'VDA_PredefinedStarRoiDec', \
            'VDA_IntegrationTime_s', 'VDA_MaxNumStarRois']
    # for vda_key, vda_values in zip(params_VDA.keys(), params_VDA.values()):
    for vda_key, vda_values in targ_info[vda_columns].iloc[0].items():
        if pd.notna(vda_values):  # This condition checks if the value is not NaN
            xml_key = vda_key.replace('VDA_', '')
            if vda_key not in columns_to_ignore:
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(vda_values)
            elif vda_key == 'VDA_TargetID':
                if pd.isnull(targ_info['Planet Name'].iloc[0]):
                    tmp_t_name = targ_info['Star Name'].iloc[0]
                else:
                    tmp_t_name = targ_info['Planet Name'].str.replace(r'\s+([bcd])$', r'\1', regex=True).iloc[0]
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = tmp_t_name#targ_info['Planet Name'].iloc[0]
            elif vda_key == 'VDA_TargetRA':
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(targ_info['RA'].iloc[0])
            elif vda_key == 'VDA_TargetDEC':
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(targ_info['DEC'].iloc[0])
            elif vda_key == 'VDA_StarRoiDetMethod':
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(targ_info['StarRoiDetMethod'].iloc[0])
            elif vda_key == 'VDA_MaxNumStarRois' and targ_info['StarRoiDetMethod'].iloc[0] == 0:
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(0)
            elif vda_key == 'VDA_MaxNumStarRois' and targ_info['StarRoiDetMethod'].iloc[0] == 1:
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(9)
            elif vda_key == 'VDA_numPredefinedStarRois' and targ_info['StarRoiDetMethod'].iloc[0] != 1:
                vda_subelement_ = ET.SubElement(vda, xml_key)
                try:
                    vda_subelement_.text = str(targ_info['numPredefinedStarRois'].iloc[0])
                except:
                    vda_subelement_.text = str('-999')
            elif vda_key == 'VDA_PredefinedStarRoiRa' and targ_info['StarRoiDetMethod'].iloc[0] != 1:
                roi_coord_columns = [col for col in targ_info.columns if col.startswith('ROI_coord_') and col != 'ROI_coord_epoch']
                roi_coord_values = targ_info[roi_coord_columns].dropna(axis = 1)
                import ast
                all_columns = np.asarray([ast.literal_eval(item) for item in roi_coord_values.values[0]])
                vda_subelement_ = ET.SubElement(vda, xml_key)
                # vda_subelement_.text = str(all_columns[:,0])
                for jj in range(all_columns.shape[0]):
                    vda_subelement_tmp = ET.SubElement(vda_subelement_, f'RA{jj+1}')
                    vda_subelement_tmp.text = f'{all_columns[jj,0]:.6f}'
            elif vda_key == 'VDA_PredefinedStarRoiDec' and targ_info['StarRoiDetMethod'].iloc[0] != 1:
                roi_coord_columns = [col for col in targ_info.columns if col.startswith('ROI_coord_') and col != 'ROI_coord_epoch']
                roi_coord_values = targ_info[roi_coord_columns].dropna(axis = 1)
                import ast
                all_columns = np.asarray([ast.literal_eval(item) for item in roi_coord_values.values[0]])
                vda_subelement_ = ET.SubElement(vda, xml_key)
                # vda_subelement_.text = str(all_columns[:,1])
                for jj in range(all_columns.shape[0]):
                    vda_subelement_tmp = ET.SubElement(vda_subelement_, f'Dec{jj+1}')
                    vda_subelement_tmp.text = f'{all_columns[jj,1]:.6f}'
            elif vda_key == 'VDA_NumTotalFramesRequested':
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(np.round(diff_in_sec/targ_info['VDA_IntegrationTime_s'].iloc[0]).astype(int))
            elif vda_key == 'VDA_NumExposuresMax':
                vda_subelement_ = ET.SubElement(vda, xml_key)
                vda_subelement_.text = str(np.round(diff_in_sec/targ_info['VDA_IntegrationTime_s'].iloc[0]).astype(int))
            
            pass
        # else:
        #     logging.info(f"Searching for occultation targets from {st} to {sp}")

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
        "TargetID": t_name, 
        "SC_Resets1": "1", 
        "SC_Resets2": "1", 
        "SC_DropFrames1": "0", 
        "SC_DropFrames2": "16", 
        "SC_DropFrames3": "0", 
        "SC_ReadFrames": "4", 
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

    # return observational_parameters, params_NIRDA, params_VDA
    return observational_parameters

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
#
#
def read_json_files(targ_list, fn_tmp):
    import pandas as pd
    import numpy as np
    import json
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

        old_column_name = "Transit Epoch (BJD_TDB)"
        # column_index = target_list_copy.columns.get_loc(old_column_name)
        new_column_name = "Transit Epoch (BJD_TDB-2400000.5)"
        if old_column_name in target_list_copy.columns:
            target_list_copy[old_column_name] = target_list_copy[old_column_name] - 2400000.5
            # print(f"Column '{old_column_name}' has been updated.")
            target_list = target_list_copy.rename(columns={old_column_name: new_column_name})
        else:
            target_list = target_list_copy

        # targ_list_copy.loc[0, "Transit Duration (hrs)"] = data["pl_trandur (hrs)"]
    return target_list
#
#
def update_target_list(targ_list, pl_names, which_targets):
    import os
    import warnings
    import glob
    warnings.filterwarnings("ignore", category=FutureWarning)
    # filtered_targ_list = targ_list[targ_list["Planet Name"].isin(pl_names)]
    # updated_targ_list = filtered_targ_list.copy()
    updated_targ_list = targ_list.copy()

    dir_tmp = '/Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files/' + which_targets
    json_files = glob.glob(f'{dir_tmp}/*.json')

    for file in json_files:
        json_data = read_json_file_new(file)
        updated_targ_list = update_dataframe_with_json(updated_targ_list, json_data)

    if 'Transit Epoch (BJD_TDB)' in updated_targ_list.columns:
        # Create the new column
        updated_targ_list['Transit Epoch (BJD_TDB) - 2400000.5'] = updated_targ_list['Transit Epoch (BJD_TDB)'] - 2400000.5
        # Handle any potential NaN values in the original column
        updated_targ_list['Transit Epoch (BJD_TDB) - 2400000.5'] = updated_targ_list['Transit Epoch (BJD_TDB) - 2400000.5'].fillna(-999)

    dir_tmp = '/Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files/'
    with open(dir_tmp + 'nirda_readout_schemes.json', 'r') as file:
        nirda_settings = json.load(file)['data']

    with open(dir_tmp + 'vda_readout_schemes.json', 'r') as file:
        vda_settings = json.load(file)['data']

    # Function to get NIRDA settings
    def get_nirda_settings(setting):
        return nirda_settings.get(setting, {})

    # Function to get VDA settings
    def get_vda_settings(setting):
        return vda_settings.get(setting, {})

        # Update the DataFrame
    for index, row in updated_targ_list.iterrows():
        # Update NIRDA settings
        nirda_setting = row['NIRDA Setting']
        nirda_values = get_nirda_settings(nirda_setting)
        for key, value in nirda_values.items():
            updated_targ_list.at[index, f'NIRDA_{key}'] = value
        
        # Update VDA settings
        vda_setting = row['VDA Setting']
        vda_values = get_vda_settings(vda_setting)
        for key, value in vda_values.items():
            updated_targ_list.at[index, f'VDA_{key}'] = value

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

def sch_occultation_targets_all(v_names, starts, stops, path, o_df, o_list):

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

def schedule_occultation_targets(v_names, starts, stops, st, sp, path, o_df, o_list, try_occ_targets):#, position):
    schedule = pd.DataFrame(index=starts, columns=['Stop', 'Target', 'Visibility'], dtype='object')
    schedule['Stop'] = stops
    schedule['Target'] = np.nan
    schedule['Visibility'] = np.nan

    # Add 'Visibility' column to o_df if it doesn't exist
    if 'Visibility' not in o_df.columns:
        o_df['Visibility'] = np.nan

    for v_name in tqdm(v_names, desc=f"{st} to {sp}: Searching for occultation target from {try_occ_targets}", leave = False):#, position=position):#, leave=leave):#, leave=(position != 0)):#desc="Processing targets"):
    # for v_name in v_names:
        # Process visibility for this target
        vis = pd.read_csv(find_file(v_name))#f"{path}/{v_name}/Visibility for {v_name}.csv")
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

def get_targets_table(which_targets):
    directory = '/Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files/' + which_targets
    df = pd.DataFrame(parse_json_files(directory))
    df = df.sort_values('Planet Name')
    df = df.reset_index(drop=True)
    df = df[['Planet Name', 'Planet Simbad Name', 'Star Name', 'Star Simbad Name', 'Number of Transits to Capture', 'Priority', 'Original Filename']]
    return df

def read_json_file_new(filename):
    with open(filename, 'r') as file:
        return json.load(file)

def update_dataframe_with_json(df, json_data):
    planet_name = f"{json_data['Star Name']} {json_data['Planet Letter']}"
    
    # Update existing columns
    df.loc[df['Planet Name'] == planet_name, 'Number of Transits to Capture'] = json_data['Number of Transits to Capture']
    
    # Add new columns from JSON data
    new_columns = ['RA', 'DEC', 'coord_epoch', 'pm_RA', 'pm_DEC', 'Jmag', 'Gmag', 'Teff (K)', 'logg', 
                   'Period (days)', 'Period Uncertainty (days)', 'Transit Duration (hrs)', 
                   'Transit Epoch (BJD_TDB)', 'Transit Epoch Uncertainty (days)', 'VDA Setting', 'NIRDA Setting']
    
    for column in new_columns:
        if column in json_data:
            df.loc[df['Planet Name'] == planet_name, column] = json_data[column]
    
    # Handle Additional Planets
    if 'Additional Planets' in json_data and json_data['Additional Planets']:
        additional_planets = ', '.join([p['Planet Letter'] for p in json_data['Additional Planets']])
        df.loc[df['Planet Name'] == planet_name, 'Additional Planets'] = additional_planets
    
    return df

def parse_json_files(directory):
    data = []
    for filename in os.listdir(directory):
        if filename.endswith('_target_definition.json'):
            with open(os.path.join(directory, filename), 'r') as file:
                json_data = json.load(file)
                
                planet_name = format_planet_name(filename)
                star_name = planet_name.rsplit(' ', 1)[0]
                if planet_name == 'HIP 65 A b':
                    star_name = 'HIP 65 A'
                
                row = {
                    'Planet Name': planet_name,
                    'Planet Simbad Name': planet_name,
                    'Star Name': star_name,
                    'Star Simbad Name': star_name,
                    'Number of Transits to Capture': json_data.get('Number of Transits to Capture', 10),
                    'Priority': 1,
                    'Original Filename': filename.replace('_target_definition.json', '')
                }
                data.append(row)
    return data


def format_planet_name(filename):
    # Remove '_target_definition.json' from the end
    name = filename.replace('_target_definition.json', '')
    
    # Split the name by underscores
    parts = name.split('_')
    
    if name == 'HIP_65_Ab':
        return f"{parts[0]} {parts[1]} {parts[-1][0]} {parts[-1][-1]}"

    # Handle special cases
    if name.startswith('L_') or name.startswith('GJ_') or name.startswith('HD_') or name.startswith('HIP_'):
        return f"{parts[0]} {parts[1][:-1]} {parts[1][-1]}"
    
    if name.startswith('LTT_'):
        return f"LTT {parts[1]} {parts[2][0]} {parts[2][1]}"
    
    if len(parts) > 2:  # For other cases with multiple underscores
        # Join all parts except the last one with spaces
        base = ' '.join(parts[:-1])
        last_part = parts[-1]
        
        # Separate the last letter from the last part
        if len(last_part) > 1 and last_part[-1].isalpha() and last_part[-2].isdigit():
            return f"{base} {last_part[:-1]} {last_part[-1]}"
        else:
            return f"{base} {last_part}"
    
    elif len(parts) == 2:  # For cases like GJ_1214b
        return f"{parts[0]} {parts[1][:-1]} {parts[1][-1]}"
    
    elif '-' in name:  # For cases like TOI-1416b
        base, letter = name[:-1], name[-1]
        return f"{base} {letter}"
    
    # If none of the above apply, return the name as is
    return name


def read_json_from_github(url):
    # Convert GitHub URL to raw content URL
    raw_url = url.replace("github.com", "raw.githubusercontent.com").replace("/blob/", "/")
    
    # Fetch the content
    response = requests.get(raw_url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse JSON content
        return json.loads(response.text)
    else:
        print(f"Failed to fetch file. Status code: {response.status_code}")
        return None
#
#
#
def load_readout_schemes(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data['data']
#
#
#
def process_target_files(keyword):
    base_dir = '/Users/vkostov/Documents/GitHub/PandoraTargetList/target_definition_files/'  # Set this to your base directory if needed
    directory = os.path.join(base_dir, keyword)
    
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return None

    # Load readout schemes
    nirda_schemes = load_readout_schemes(base_dir + 'nirda_readout_schemes.json')
    vda_schemes = load_readout_schemes(base_dir + 'vda_readout_schemes.json')

    # def flatten_dict(d, parent_key='', sep='_'):
    #     items = []
    #     for k, v in d.items():
    #         new_key = f"{parent_key}{sep}{k}" if parent_key else k
    #         if isinstance(v, dict):
    #             items.extend(flatten_dict(v, new_key, sep=sep).items())
    #         elif isinstance(v, list):
    #             for i, item in enumerate(v):
    #                 items.extend(flatten_dict(item, f"{new_key}{sep}{i}", sep=sep).items())
    #         else:
    #             items.append((new_key, v))
    #     return dict(items)

    def flatten_dict(d, parent_key='', sep='_'):
        items = []
        for k, v in d.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else k
            if isinstance(v, dict):
                items.extend(flatten_dict(v, new_key, sep=sep).items())
            elif isinstance(v, list):
                for i, item in enumerate(v):
                    if isinstance(item, dict):
                        items.extend(flatten_dict(item, f"{new_key}{sep}{i}", sep=sep).items())
                    else:
                        items.append((f"{new_key}{sep}{i}", item))
            else:
                items.append((new_key, v))
        return dict(items)

    # if keyword == 'auxiliary-standard':
    #     aaa = 333.

    data_list = []
    for filename in tqdm(os.listdir(directory)):
        if filename.endswith('_target_definition.json'):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Remove unwanted keys
            for key in ['Time Created', 'Version',  'Author', 'Time Updated']:
                data.pop(key, None)
            
            flat_data = flatten_dict(data)
            
            original_filename = filename.replace('_target_definition.json', '')
            flat_data['Original Filename'] = original_filename

            # if keyword == 'exoplanet' or keyword == 'auxiliary-exoplanet' or 'primary-exoplanet' or 'secondary-exoplanet':#!= "occultation-standard":
            if keyword in ('exoplanet', 'auxiliary-exoplanet', 'primary-exoplanet', 'secondary-exoplanet'):
                priority_fn = os.path.join(directory, f"{keyword}_priorities.csv")
                metadata, data = read_priority_csv(priority_fn)
                priority_ = data[data["target"] == original_filename]["priority"].values[0]
                flat_data['Priority'] = priority_
                flat_data['Number of Transits to Capture'] = int(data[data["target"] == original_filename]["transits_req"].values[0])
            elif keyword in ('auxiliary-standard', 'monitoring-standard'):
                priority_fn = os.path.join(directory, f"{keyword}_priorities.csv")
                metadata, data = read_priority_csv(priority_fn)
                priority_ = data[data["target"] == original_filename]["priority"].values[0]
                flat_data['Priority'] = priority_
                flat_data['Number of Hours Requested'] = int(data[data["target"] == original_filename]["hours_req"].values[0])
            elif keyword in ('occultation-standard'):
                flat_data['Priority'] = 0.1 # Default priority
                # flat_data['Number of Transits to Capture'] = 0

            # if (keyword == "primary-exoplanet") or (keyword == "secondary-exoplanet"):
            #     flat_data['Priority'] = 1  # Default priority
            # else:
            #     non_primary_priority_fn = os.path.join(directory, f"{keyword}_priorities.csv")
            #     metadata, data = read_priority_csv(non_primary_priority_fn)
            #     non_primary_priority = data[data["target"] == original_filename]["priority"].values[0]
            #     flat_data['Priority'] = non_primary_priority
            
            # Add NIRDA and VDA readout scheme data
            nirda_setting = flat_data.get('NIRDA Setting')
            vda_setting = flat_data.get('VDA Setting')

            # NIRDA fixed parameters
            for key, value in nirda_schemes["FixedParameters"].items():
                flat_data[f'NIRDA_{key}'] = value

            # VDA fixed parameters
            for key, value in vda_schemes["FixedParameters"].items():
                flat_data[f'VDA_{key}'] = value
            
            if nirda_setting in nirda_schemes:
                for key, value in nirda_schemes[nirda_setting].items():
                    flat_data[f'NIRDA_{key}'] = value
            
            if vda_setting in vda_schemes:
                for key, value in vda_schemes[vda_setting].items():
                    flat_data[f'VDA_{key}'] = value
            
            # if not filename.startswith('DR3') or keyword != 'monitoring-standard':
            # if keyword != 'monitoring-standard' or keyword !='occultation-standard':
            # if keyword != 'monitoring-standard' and keyword !='occultation-standard':
            if keyword not in ('monitoring-standard', 'occultation-standard'):
                # Separate the last lowercase letter with a space
                planet_name = re.sub(r'([a-z])$', r' \1', flat_data.get('Planet Name', ''))
                flat_data['Planet Name'] = planet_name
                flat_data['Planet Simbad Name'] = planet_name
                flat_data['Star Simbad Name'] = flat_data.get('Star Name', '')
                
                # Add the new column
                if 'Transit Epoch (BJD_TDB)' in flat_data:
                    flat_data['Transit Epoch (BJD_TDB-2400000.5)'] = flat_data['Transit Epoch (BJD_TDB)'] - 2400000.5
            else:
                flat_data['Star Simbad Name'] = flat_data.get('Star Name', '')
                flat_data['Planet Name'] = flat_data.get('Star Name', '')
                flat_data['Planet Simbad Name'] = flat_data.get('Star Name', '')

            flat_data['RA'], flat_data['DEC']= update_coordinates_astropy(flat_data['RA'], flat_data['DEC'], flat_data['pm_RA'], flat_data['pm_DEC'])
            
            data_list.append(flat_data)

    if not data_list:
        print(f"No valid JSON files found in {directory}")
        return None

    df = pd.DataFrame(data_list)
    
    # Determine columns based on whether it's a Gaia DR3 file or not
    # if df['Original Filename'].str.startswith('DR3').any() or keyword in ('auxiliary-standard', 'monitoring-standard'):#keyword == 'monitoring-standard':
    if keyword in ('auxiliary-standard', 'monitoring-standard', 'occultation-standard'):
        columns_order = ['Star Name', 'Star Simbad Name']
        columns_order.extend([col for col in df.columns if col not in columns_order and col != 'Priority'])
        columns_order.append('Priority')
    else:
        columns_order = ['Planet Name', 'Planet Simbad Name', 'Star Name', 'Star Simbad Name', 
                         'Number of Transits to Capture', 'Priority', 'Original Filename']
        if 'Transit Epoch (BJD_TDB-2400000.5)' in df.columns:
            columns_order.append('Transit Epoch (BJD_TDB-2400000.5)')
        columns_order.extend([col for col in df.columns if col not in columns_order])
    
    # Move NIRDA and VDA settings and their related columns to the end
    nirda_vda_columns = [col for col in df.columns if col.startswith(('NIRDA', 'VDA'))]
    for col in nirda_vda_columns:
        columns_order.remove(col)
    columns_order.extend(nirda_vda_columns)
    
    return df[columns_order]


def check_and_update_target(info, flag):
    if flag:
        for i, target in enumerate(info['Target']):
            current_time = occ_target_times.get(target, timedelta())
            new_time = current_time + (oc_stops[i] - oc_starts[i])
            if new_time <= occ_limit:
                occ_target_times[target] = new_time
                return info.iloc[[i]], True
            else:
                logging.info(f"Skipping target {target} as it has reached the time limit.")
    return None, False

def save_observation_time_report(all_target_obs_time, target_list, output_file):
    primary_targets = set(target_list["Planet Name"])
    
    with open(output_file, 'w') as f:
        f.write("Target,Is Primary,Total Observation Time (hours)\n")
        
        for target, time in all_target_obs_time.items():
            is_primary = "Yes" if target in primary_targets else "No"
            hours = time.total_seconds() / 3600  # Convert to hours
            f.write(f"{target},{is_primary},{hours:.2f}\n")


def read_priority_csv(file_path):
    # Read the entire file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract metadata
    metadata = {}
    for line in lines:
        if line.startswith('#'):
            if ':' in line:
                key, value = line.strip('# ').split(':', 1)
                metadata[key.strip()] = value.strip()
        else:
            break  # Stop when we hit non-comment lines

    # Find the index where the actual CSV data starts
    data_start = next(i for i, line in enumerate(lines) if not line.startswith('#'))

    # Read the CSV data
    df = pd.read_csv(file_path, skiprows=data_start)

    return metadata, df

def create_aux_list(target_definition_files, PACKAGEDIR):
    import pandas as pd
    from functools import reduce
    import os
    
    # Create full file paths
    file_paths = [f"{PACKAGEDIR}/data/{file}_targets.csv" for file in target_definition_files]

    # Read all CSV files into a list of DataFrames
    dfs = []
    for file in file_paths:
        if os.path.exists(file):
            df = pd.read_csv(file)
            dfs.append(df)
        else:
            print(f"Warning: File {file} not found. Skipping.")

    if not dfs:
        print("No valid files found. Exiting.")
        exit()

    # Store the column order of the first file
    first_file_columns = dfs[0].columns.tolist()

    # Find common columns
    common_columns = list(reduce(set.intersection, [set(df.columns) for df in dfs]))

    # Select only common columns from each DataFrame
    dfs = [df[common_columns] for df in dfs]

    # Concatenate all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)

    # Remove duplicate rows if any
    combined_df = combined_df.drop_duplicates()

    # Sort the columns of combined_df to match the order in the first file
    sorted_columns = [col for col in first_file_columns if col in common_columns]
    combined_df = combined_df[sorted_columns]

    # Write the result to a new CSV file
    output_file = f"{PACKAGEDIR}/data/aux_list_new.csv"
    combined_df.to_csv(output_file, index=False)

    # print(f"Combined CSV created with {len(common_columns)} common columns.")
    # print(f"Output file: {output_file}")
    # print("Common columns in order:", sorted_columns)

    return

def check_if_transits_in_obs_window(tracker, temp_df, target_list, start, pandora_start, pandora_stop, \
    sched_start, sched_stop, obs_rng, obs_window, sched_wts, transit_coverage_min):
    for i in range(len(tracker)):
        planet_name = tracker["Planet Name"][i]
        ra_tar, dec_tar = tracker["RA"][i], tracker["DEC"][i]

        if (
            tracker.loc[(tracker["Planet Name"] == planet_name), "Transits Needed"][
                i
            ]
            == 0
        ):
            pass

        else:
            star_name = target_list["Star Name"][
                np.where(target_list["Planet Name"] == planet_name)[0][0]
            ]
            planet_data = pd.read_csv(
                f"{PACKAGEDIR}/data/targets/{star_name}/{planet_name}/Visibility for {planet_name}.csv"
            )

            planet_data = planet_data.drop(
                planet_data.index[
                    (planet_data["Transit_Coverage"] < transit_coverage_min)
                ]
            ).reset_index(drop=True)
            planet_data["Transit_Start"] = Time(
                planet_data["Transit_Start"], format="mjd", scale="utc"
            ).to_value("datetime")
            planet_data["Transit_Stop"] = Time(
                planet_data["Transit_Stop"], format="mjd", scale="utc"
            ).to_value("datetime")
            planet_data = planet_data.drop(
                planet_data.index[(planet_data["Transit_Start"] < start)]
            ).reset_index(drop=True)
            start_transits = planet_data["Transit_Start"].copy()
            end_transits = planet_data["Transit_Stop"].copy()
            # start_transits = Time(planet_data["Transit_Start"], format="mjd", scale="utc").to_value("datetime")
            # end_transits = Time(planet_data["Transit_Stop"], format="mjd", scale="utc").to_value("datetime")

            p_trans = planet_data.index[
                (pandora_start <= start_transits) & (end_transits <= pandora_stop)
            ]
            s_trans = planet_data.index[
                (sched_start <= start_transits) & (end_transits <= sched_stop)
            ]
            tracker.loc[
                (tracker["Planet Name"] == planet_name), "Transits Left in Lifetime"
            ] = len(p_trans)
            tracker.loc[
                (tracker["Planet Name"] == planet_name), "Transits Left in Schedule"
            ] = len(s_trans)
            tracker.loc[
                (tracker["Planet Name"] == planet_name), "Transit Priority"
            ] = (
                tracker.loc[
                    (tracker["Planet Name"] == planet_name),
                    "Transits Left in Lifetime",
                ]
                - tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transits Needed"
                ]
            )

            # Remove seconds and below from times
            for j in range(len(start_transits)):
                start_transits.iloc[j] = start_transits.iloc[j] - timedelta(
                    seconds=start_transits.iloc[j].second,
                    microseconds=start_transits.iloc[j].microsecond,
                )
                end_transits.iloc[j] = end_transits.iloc[j] - timedelta(
                    seconds=end_transits.iloc[j].second,
                    microseconds=end_transits.iloc[j].microsecond,
                )

            early_start = end_transits - timedelta(
                hours=20
            )  # Earliest start time to capture transit plus >=4 hours post transit
            late_start = start_transits - timedelta(
                hours=4
            )  # Latest start time to capture transit plus >=4 hours pre transit

            # Check if any transit occurs during observing window
            for j in range(len(early_start)):
                start_rng = pd.date_range(early_start[j], late_start[j], freq="min")
                overlap_times = obs_rng.intersection(start_rng)
                if len(overlap_times) > 0:
                    # Calc a 'transit factor'
                    t_left = tracker.loc[
                        (tracker["Planet Name"] == planet_name),
                        "Transits Left in Lifetime",
                    ].iloc[0]
                    t_need = tracker.loc[
                        (tracker["Planet Name"] == planet_name), "Transits Needed"
                    ].iloc[0]
                    t_factor = t_left / t_need

                    # Calc scheduling efficiency factor
                    obs_start = overlap_times[0]
                    gap_time = obs_start - obs_rng[0] # THIS IS THE time between the start of the current observation window and the start of the potential observation.
                    s_factor = 1 - (gap_time / obs_window)  # maximize

                    # Calc a quality factor (currently based on transit coverage, SAA crossing, scheduling efficiency)
                    trans_cover = planet_data["Transit_Coverage"][j]  # maximize
                    saa_cover = planet_data["SAA_Overlap"][j]
                    q_factor = (
                        (sched_wts[0] * trans_cover)
                        + (sched_wts[1] * (1 - saa_cover))
                        + (sched_wts[2] * s_factor)
                    )

                    temp = [
                        [
                            planet_name,
                            ra_tar, 
                            dec_tar,
                            obs_start,
                            gap_time,
                            planet_data["Transit_Coverage"][j],
                            saa_cover,
                            s_factor,
                            t_factor,
                            q_factor,
                            np.nan,
                        ]
                    ]
                    temp = pd.DataFrame(
                        temp,
                        columns=[
                            "Planet Name",
                            "RA",
                            "DEC",
                            "Obs Start",
                            "Obs Gap Time",
                            "Transit Coverage",
                            "SAA Overlap",
                            "Schedule Factor",
                            "Transit Factor",
                            "Quality Factor",
                            "Comments",
                        ],
                    )
                    # temp_df = pd.concat([temp_df, temp], axis=0) LATEST PYTHON DOESN'T LIKE THIS SO CHANGE IT TO THE LINE BELOW
                    if temp_df.empty:
                        temp_df = temp.copy()
                    else:
                        temp_df = pd.concat([temp_df, temp], axis=0)#, ignore_index=True)

    return temp_df

def remove_suffix(name):
    import re
    return re.sub(r' [a-z]$', '', name)


def find_file(filename):
    possible_paths = [
        f"{PACKAGEDIR}/data/targets/",
        f"{PACKAGEDIR}/data/aux_targets/"
    ]
    
    for path in possible_paths:
        full_path = os.path.join(path, filename, f"Visibility for {filename}.csv")
        if os.path.isfile(full_path):
            return full_path
    
    return None  # File not found in any of the directories

def sch_occ_old_but_working(starts, stops, list_path, sort_key=None, prev_obs = None):#, position = 0):#**kwargs):
    
    #build empty dataframe except for starts and stops
    e_sched = [['',datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"),datetime.strftime(stops[s], "%Y-%m-%dT%H:%M:%SZ"), '', ''] for s in range(len(starts))]
    o_df = pd.DataFrame(e_sched,columns=["Target","start","stop", "RA", "DEC"])
    
    #convert to mjd to compare with visibility data
    starts=Time(starts, format='datetime').to_value('mjd')
    stops=Time(stops, format='datetime').to_value('mjd')
    
    if sort_key == None:
        #No occluded target scheduling, free time
        starts=starts.to_value('datetime')
        stops=stops.to_value('datetime')
        free = [["Free Time",datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"),datetime.strftime(stops[s], "%Y/%m/%d, %H:%M:%S"), '', ''] for s in range(len(starts))]
        o_df = pd.DataFrame(free, columns=["Target", "start", "stop", "RA", "DEC"])
    
    else:       
        o_list = pd.read_csv(list_path)
        ras=o_list['RA']
        decs=o_list['DEC']
        
        if sort_key == 'closest':
            #sort name list based on minimized sky distance
            #prev_obs must be specified as an array consisting of ra and dec (in degrees) for the previous observation
            #e.g. [359.10775132017, -49.28901740485]
            #this seeks to minimize slew distance to an occultation target, though actual slew sims are not performed
            try:
                po_sc=SkyCoord(unit='deg', ra=prev_obs[0], dec=prev_obs[1])
                oc_sc=[SkyCoord(unit='deg', ra=ras[n], dec=decs[n]) for n in range(len(ras))]

                dif=[oc_sc[n].separation(po_sc).deg for n in range(len(oc_sc))]
                o_list['sky_dif'] = dif
                o_list = o_list.sort_values(by='sky_dif').reset_index(drop=True)
                
            except NameError:
                print('No previous observation was specified, defaulting to random auxiliary target.')
                o_list.sample(frac=1).reset_index(drop=True)
        else:
            #default sort is random
            o_list.sample(frac=1).reset_index(drop=True)
        
        v_names = o_list['Star Name']
        v_names=np.array(v_names)
        #For prioritization via flag later
        o_flag = o_list['Priority']
        #Reload these
        ras=o_list['RA']
        decs=o_list['DEC']

        multi_target_occultation = True#False

        d_flag = False
        
        #empty dataframe to hold visibility information for multiple targets
        v_ = np.asarray(np.zeros(len(starts)), dtype=bool)
        v_df = pd.DataFrame([v_])
        vis_df = pd.DataFrame(columns=range(len(starts))).astype(bool)
        d_flag = False
        for n in range(len(v_names)):
            try:
                if (list_path == tar_path) or (list_path == tar_path_ALL):
                    vis=pd.read_csv(f"{PACKAGEDIR}/data/targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
                elif list_path == aux_path:
                    vis=pd.read_csv(f"{PACKAGEDIR}/data/aux_targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
                vis_ = vis['Time(MJD_UTC)']
                
                #array to hold visibility info for this target
                v_ar = np.asarray(np.zeros(len(starts)), dtype=bool)
                v_ar_any = np.asarray(np.zeros(len(starts)), dtype=bool)
                #iterate through occultation periods and see if v_names[n] is visible for all of them
                vis_f = False
                for s in range(len(starts)):
                    #get occultation time window
                    win = vis.index[(vis_ >= starts[s]) & (vis_ <= stops[s])].tolist()
                    if np.all(vis['Visible'][win] == 1):# or vis_ratio >= 0.6: 
                        # UNCOMMENT THE REST OF THE LINE ABOVE FOR TARGETS THAT ARE NOT VISIBLE FOR MANY HOURS AND ALLOW "OCCULTATION"
                        # TARGET TO BE CONSIDERED AS VISIBLE >= 60% OF 90 MINUTES!!!
                        v_ar[s] = True               
                    else:
                        vis_f = True

                    # VK START
                    # CHECK THE FRACTION OF TIME OCCULTATION TARGET IS VISIBILE DURING OCCULTATION 
                    if len(vis['Visible'][win]) > 0.:
                        vis_ratio = len(vis['Visible'][win][vis['Visible'][win] == 1])/len(vis['Visible'][win])
                    # VK END

                    # print(v_names[n], Time(starts[s], format="mjd", scale="utc").to_value("datetime").strftime("%H:%M:%S"), \
                    #     Time(stops[s], format="mjd", scale="utc").to_value("datetime").strftime("%H:%M:%S"), vis_ratio)
                            # len(vis['Visible'][win][vis['Visible'][win] == 1]), len(vis['Visible'][win]))#np.asarray(vis['Visible'][win]))

                # vis_df.loc[len(vis_df)] = v_ar
                # vis_df_sum = np.sum(np.asarray(vis_df), axis = 0)
                # print(n, vis_df_sum)
                # if np.all(vis_df_sum > 0):
                #     stop

                #if not visible for all times, check if any entry in v_df and this one cover the total occultation time
                if vis_f:
                    if not d_flag:

                        # vis_df_sum = np.sum(np.asarray(vis_df), axis = 0)
                        # # print(n, vis_df_sum)
                        # if np.all(vis_df_sum > 0):
                        #     d_flag = True
                        # else:
                        #     vis_df.loc[len(vis_df)] = v_ar

                        v_arr=np.asarray(v_df)
                        overlap=np.where([np.all((v_arr+np.asarray(v_ar, dtype=bool))[i]) for i in range(len(v_arr))])[0]
                        if len(overlap) > 0:
                            #at least one entry has visibility that covers the total time along with the current target
                            #take both and enter them in o_df in their respective times, prefering the closer one
                            m=overlap[0]
                            v1=v_arr[m]
                            for s in range(len(starts)):
                                if v1[s]:
                                    o_df['Target'][s]=v_names[m]
                                    o_df['RA'][s]=str(ras[m])
                                    o_df['DEC'][s]=str(decs[m])
                                else:
                                    o_df['Target'][s]=v_names[n]
                                    o_df['RA'][s]=str(ras[n])
                                    o_df['DEC'][s]=str(decs[n])
                            d_flag=True
                        else:
                            #add the current visibility array to the master list
                            v_df.loc[len(v_df.index)] = v_ar
                            # vis_df.loc[len(vis_df)] = v_ar
                    else:
                        # break
                        continue
                
                else:
                    #If we made it here, the target is visible for the entire occultation time
                    #since the list is already sorted, break and use this one!
                    #add to the list of visible targets (not necessary, but keeps functionality with prev code)
                    o_df['Target'][:]=v_names[n]
                    o_df['RA'][:]=str(ras[n])
                    o_df['DEC'][:]=str(decs[n])
                    d_flag=True
                    # print(st_name, ': ', n, v_names[n], v_names[n], 'Occ target visible ALL')
                    break
                    # THIS BREAK DOESNT SEEM TO WORK!!!! REPLACE WITH return o_df, d_flag
                    # return o_df, d_flag

                # print(f"{st_name}: {n} ({v_names[n]}) not 100% visible, try next on the list")
            
            #If a target(s) on the list don't have visibility data, ignore them!
            except FileNotFoundError:
                continue    
        
        #if there were not <=2 occultation targets that covered the time (cya clause)
        if not d_flag:
            #only if considering aux targets, real last ditch effort here
            if list_path.endswith('aux_list.csv'):# or st_name.startswith('Gaia'):
                #check one last time to make sure there are no gaps
                v_arr=np.asarray(v_df)
                if np.all([np.any(v_arr[i]) for i in range(len(v_arr))]):
                    #iterate and set the nearest that is visible for each window as the occultation target
                    for m in range(len(v_arr)):
                        # VK START: IT LOOKS LIKE THE i IN v_arr[i] IS NOT PART OF THIS FOR LOOP. MAYBE A BUG?
                        #n=np.where(v_arr[i])[0][0]
                        # CHANGE i TO m
                        # VK END
                        n=np.where(v_arr[m])[0][0]
                        o_df['Target'][s]=v_names[n]
                        o_df['RA'][s]=str(ras[n])
                        o_df['DEC'][s]=str(decs[n])
                    d_flag=True
                    print('More than two auxiliary targets were needed to cover the occultation time.')

    return o_df, d_flag

def update_coordinates_astropy(ra0, dec0, pm_ra, pm_dec):
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    import astropy.units as u

    t0 = Time('J2016.0')
    t1 = Time('2025-11-15')
    coord = SkyCoord(ra=ra0*u.degree, dec=dec0*u.degree, 
                     pm_ra_cosdec=pm_ra*u.mas/u.yr, 
                     pm_dec=pm_dec*u.mas/u.yr, 
                     frame='icrs', obstime=t0, distance=None)
    
    new_coord = coord.apply_space_motion(new_obstime=Time(t1))
    
    return new_coord.ra.degree, new_coord.dec.degree

def check_visibility():

    targ = ['TOI-270', 'TOI-2076', 'TOI-942', 'TOI-244', 'HD_73583', 'HAT-P-12']
    fig, axs = plt.subplots(len(targ), 1, figsize=(20, 8))
    # fig.tight_layout(pad=0.5)

    for ii, tt in enumerate(targ):
        tf_vis = pd.read_csv(f"/Users/vkostov/Documents/GitHub/pandora-scheduler/src/pandorascheduler/data/targets/{tt}/Visibility for {tt}.csv")
        time_mjd = tf_vis['Time(MJD_UTC)'].values
        vis = tf_vis['Visible'].values
        
        nbins = 100
        
        axs[ii].plot(np.mean(time_mjd[:-(time_mjd.size % nbins)].reshape(-1, nbins), axis=1),
                    np.mean(vis[:-(vis.size % nbins)].reshape(-1, nbins), axis=1), '.-', ms=0.1, label=tt)
        axs[ii].set_xlim(min(time_mjd), max(time_mjd))
        axs[ii].legend()
        
        if ii < len(targ) - 1:  # Remove x-axis labels for all but the last subplot
            axs[ii].set_xticks([])
        
        # axs[ii].set_yticks([])

        # plt.tight_layout()

    plt.subplots_adjust(hspace=0.1) 
    plt.show()