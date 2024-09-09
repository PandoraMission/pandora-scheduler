import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import json
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
from astropy.time import Time
from tqdm import tqdm
import os

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

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

def remove_short_sequences(v_time, v_flag, sequence_too_short):
    result_flag = v_flag.copy()
    result_time = v_time.copy()
    
    # Find sequences of 1s and 0s
    changes = np.diff(np.concatenate(([0], result_flag, [0])))
    starts = np.where(changes == 1)[0]
    ends = np.where(changes == -1)[0]

    for start, end in zip(starts, ends):
        duration = (result_time[end-1] - result_time[start]).total_seconds() / 60
        if duration < sequence_too_short:
            result_flag[start:end] = 0

    return result_time, result_flag

def break_long_sequences(ranges, step):
    broken_ranges = []
    for start, stop in ranges:
        current = start
        while current < stop:
            next_val = min(current + step, stop)
            broken_ranges.append([current, next_val])
            current += step
    return broken_ranges

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

def print_element_from_xml(elem, level=0):
    print("  " * level + f"{elem.tag}: {elem.text.strip() if elem.text else ''}")
    for child in elem:
        print_element_from_xml(child, level + 1)

def sch_occ(starts, stops, list_path, sort_key=None, prev_obs=None, tar_path=None, tar_path_ALL=None, aux_path=None):
    # Build empty dataframe except for starts and stops
    e_sched = [['', datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"), datetime.strftime(stops[s], "%Y-%m-%dT%H:%M:%SZ"), '', ''] for s in range(len(starts))]
    o_df = pd.DataFrame(e_sched, columns=["Target", "start", "stop", "RA", "DEC"])

    # Convert to mjd to compare with visibility data
    starts = Time(starts, format='datetime').to_value('mjd')
    stops = Time(stops, format='datetime').to_value('mjd')

    if sort_key is None:
        # No occluded target scheduling, free time
        starts = Time(starts, format='mjd').to_value('datetime')
        stops = Time(stops, format='mjd').to_value('datetime')
        free = [["Free Time", datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"), datetime.strftime(stops[s], "%Y-%m-%d %H:%M:%S"), '', ''] for s in range(len(starts))]
        o_df = pd.DataFrame(free, columns=["Target", "start", "stop", "RA", "DEC"])
    else:
        o_list = pd.read_csv(list_path)
        ras = o_list['RA']
        decs = o_list['DEC']
        if sort_key == 'closest':
            # Sort name list based on minimized sky distance
            try:
                po_sc = SkyCoord(unit='deg', ra=prev_obs[0], dec=prev_obs[1])
                oc_sc = [SkyCoord(unit='deg', ra=ras[n], dec=decs[n]) for n in range(len(ras))]
                dif = [oc_sc[n].separation(po_sc).deg for n in range(len(oc_sc))]
                o_list['sky_dif'] = dif
                o_list = o_list.sort_values(by='sky_dif').reset_index(drop=True)
            except NameError:
                print('No previous observation was specified, defaulting to random auxiliary target.')
                o_list = o_list.sample(frac=1).reset_index(drop=True)
        else:
            # Default sort is random
            o_list = o_list.sample(frac=1).reset_index(drop=True)

    v_names = o_list['Star Name']
    v_names = np.array(v_names)

    # For prioritization via flag later
    o_flag = o_list['Flag']

    # Reload these
    ras = o_list['RA']
    decs = o_list['DEC']

    multi_target_occultation = True
    d_flag = False
    if multi_target_occultation:
        if (list_path == tar_path) or (list_path == tar_path_ALL):
            path_ = f"{PACKAGEDIR}/data/targets"
            try_occ_targets = 'target list'
        elif list_path == aux_path:
            path_ = f"{PACKAGEDIR}/data/aux_targets"
            try_occ_targets = 'aux list'
        else:
            raise ValueError(f"Unrecognized list_path: {list_path}")

        o_df, d_flag = schedule_occultation_targets(v_names, starts, stops, path_, o_df, o_list, try_occ_targets)

    return o_df, d_flag

def visualize_schedule(schedule_data, output_dir, filename):
    if not schedule_data:
        print("Warning: schedule_data is empty. No visualization will be created.")
        return

    fig, ax = plt.subplots(figsize=(20, 10))  # Increased figure size
    
    targets = []
    starts = []
    durations = []
    colors = []
    
    for target, start, stop, is_visible in schedule_data:
        targets.append(target)
        starts.append(start)
        duration = (stop - start).total_seconds() / 3600  # Convert to hours
        durations.append(duration)
        colors.append('green' if is_visible else 'red')
    
    # Convert starts to numbers for plotting
    start_nums = mdates.date2num(starts)
    
    # Create the Gantt chart
    unique_targets = list(dict.fromkeys(targets))  # Remove duplicates while preserving order
    y_positions = range(len(unique_targets))
    
    for i, target in enumerate(unique_targets):
        target_indices = [j for j, t in enumerate(targets) if t == target]
        ax.barh([i] * len(target_indices), 
                [durations[j] for j in target_indices], 
                left=[start_nums[j] for j in target_indices], 
                height=0.5, align='center', 
                color=[colors[j] for j in target_indices], 
                alpha=0.8)
    
    # Customize the chart
    ax.set_yticks(y_positions)
    ax.set_yticklabels(unique_targets)
    ax.set_ylabel('Targets')
    ax.set_xlabel('Time')
    ax.set_title('Observation Schedule Visualization')
    
    # Format x-axis to show dates
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    fig.autofmt_xdate()  # Rotate and align the tick labels
    
    # Add a legend
    ax.bar(0, 0, color='green', label='Visible', alpha=0.8)
    ax.bar(0, 0, color='red', label='Not Visible', alpha=0.8)
    ax.legend()
    
    plt.tight_layout()
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the figure
    output_path = os.path.join(output_dir, f"{filename}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)  # Close the figure to free up memory
    
    print(f"Schedule visualization saved as {output_path}")
    print(f"Number of unique targets: {len(unique_targets)}")
    print(f"X-axis range: {ax.get_xlim()}")
    print(f"Y-axis range: {ax.get_ylim()}")


# Add this function to test the visualization with sample data
def test_visualize_schedule():
    now = datetime.now()
    test




