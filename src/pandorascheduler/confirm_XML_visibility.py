import xml.etree.ElementTree as ET
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from astropy.time import Time
import os
import helper_codes_aux as hcc
import tqdm as tqdm

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Parse the XML file
fname = f'{PACKAGEDIR}/data/calendar_Pandora_Schedule_TEST_wip.xml'
tree = ET.parse(fname)
root = tree.getroot()

aux_vis_path=f'{PACKAGEDIR}/data/aux_targets/'
tar_vis_path=f'{PACKAGEDIR}/data/targets/'

namespace = {'ns': '/pandora/calendar/'}

# Function to check visibility
def check_visibility(target, start_time, stop_time):
    # Convert ISO time to MJD
    start_mjd = Time(start_time, format='isot', scale='utc').mjd
    stop_mjd = Time(stop_time, format='isot', scale='utc').mjd

    # Read visibility data

    if target.endswith(('b', 'c', 'd', 'e', 'f')):
        st_name = target[:-2]
    else:
        st_name = target

    # st_name = target[:-2] if target.endswith('b') or target.endswith('c') else target

    if not target.startswith('DR3'):
        v_data = pd.read_csv(tar_vis_path+f'{st_name}/Visibility for {st_name}.csv')
    else:
        v_data = pd.read_csv(aux_vis_path+f'{target}/Visibility for {target}.csv')

    # Filter data for the observation period
    mask = (v_data['Time(MJD_UTC)'] >= start_mjd) & (v_data['Time(MJD_UTC)'] <= stop_mjd)
    period_data = v_data[mask]
    
    return period_data['Visible'].tolist(), period_data['Time(MJD_UTC)'].tolist()

# Function to create and save a figure for a visit
def create_visit_figure(visit_data, visit_id, target):
    import numpy as np
    fig, ax = plt.subplots(figsize=(10, 5))

    unique_targets = list(dict.fromkeys(d['target'] for d in visit_data))
    target_to_y = {target: i for i, target in enumerate(unique_targets)}

    

    # for data in visit_data:
    for i, data in enumerate(visit_data):
        visible_ = np.asarray(data['visibility']) == 1.
        times = np.asarray(data['times'])
        y_value = target_to_y[data['target']]

        vis_times = data['times']
        vis_values = [i+1 if v else None for v in data['visibility']]
        # ax.scatter(vis_times, [y_value] * len(vis_times), color='green', marker='.', s=2)
        ax.plot(times, [y_value] * len(times), color='green', linewidth=1)

        if not(all(visible_)):
            for t in times[~visible_]:
                ax.axvline(t, ymin=y_value/len(target_to_y), 
                    ymax=(y_value+1)/len(target_to_y), 
                    color='r', linewidth=1, alpha=1)
            # vis_times = data['times']
            # vis_values = [i+1 if v else None for v in data['visibility']]
            # ax.scatter(vis_times, [i+1] * len(vis_times), c=vis_values, cmap='RdYlGn', vmin=0, vmax=1, s=3)
        #     ax.scatter(times[visible_], [y_value] * np.sum(visible_), color='green', marker='.', s=2)
        # else:
        # Plot green points for visible times, red for invisible
        # ax.scatter(times[visible_], [y_value] * np.sum(visible_), color='green', marker='.', s=2)
            # ax.scatter(times[~visible_], [y_value] * np.sum(~visible_), color='red', marker='o', s=4)


    # for i, data in enumerate(visit_data):
    #     visible_ = np.asarray(data['visibility']) == 1.
    #     times = np.linspace(data['start'], data['stop'], len(visible_))

    #     # Plot green points for visible times, red for invisible
    #     ax.scatter(times[visible_], [i] * np.sum(visible_), color='green', marker='.', s=2)
    #     ax.scatter(times[~visible_], [i] * np.sum(~visible_), color='red', marker='o', s=4)

    #     # if visible_:
    #     #     color_ = 'green'
    #     # else:
    #     #     color_ = 'red'

    #     # if all(data['visibility']):
    #     #     color_ = 'green'
    #     #     marker_ = 'o' 
    #     #     ms_ = 2
    #     # elif data['visibility']:
    #     #     color_ = 'red'
    #     #     marker_ = 'o' 
    #     #     ms_ = 4
    #     # else:
    #     #     color_ = 'gray'
    #     #     marker_ = 'x' 
    #     #     ms_ = 8

    #     # color_ = 'green' if all(data['visibility']) else 'red' if data['visibility'] else 'gray'

    #     # ax.plot([data['start'], data['stop']], [i, i], color=color_, marker = 'o', ms = 2)
        
    #     # # Plot visibility
    #     # if all(data['visibility']):
    #     #     vis_times = data['times']
    #     #     vis_values = [i+1 if v else None for v in data['visibility']]
    #     #     ax.scatter(vis_times, [i+1] * len(vis_times), c=vis_values, cmap='RdYlGn', vmin=0, vmax=1, s=3)
    #     #     # print(data['start'], data['stop'], all(data['visibility']))

    # Set y-axis labels
    ax.set_yticks(range(len(unique_targets)))
    ax.set_yticklabels(unique_targets)

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    fig.autofmt_xdate()

    plt.title(f'Target Visibility During Visit {visit_id}')
    plt.xlabel('Time')
    plt.ylabel('Target')

    plt.tight_layout()

    # Save the figure as a PNG file
    try:
        output_file = os.path.join(output_dir, f'visit_{visit_id}_visibility_for_{target}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        # print(f"Figure for Visit {visit_id} saved as {output_file}")
    except:
        no_fig = True

# Prepare output directory
output_dir = PACKAGEDIR + '/data/confirm_visibility'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Process each Visit
for visit in tqdm.tqdm(root.findall('ns:Visit', namespace)):
    visit_id = visit.find('ns:ID', namespace).text
    visit_data = []

    for obs_seq in visit.findall('ns:Observation_Sequence', namespace):
        target_elem = obs_seq.find('./ns:Observational_Parameters/ns:Target', namespace)
        target = target_elem.text if target_elem is not None and target_elem.text else "No target"
        
        start_elem = obs_seq.find('./ns:Observational_Parameters/ns:Timing/ns:Start', namespace)
        stop_elem = obs_seq.find('./ns:Observational_Parameters/ns:Timing/ns:Stop', namespace)
        
        if start_elem is None or stop_elem is None:
            print(f"Warning: Missing timing data for sequence in Visit {visit_id}")
            continue
        
        start_time = start_elem.text
        stop_time = stop_elem.text
        
        if target != "No target":
            try:
                visibility, times = check_visibility(target, start_time, stop_time)
            except Exception as e: 
                print(f"Error checking visibility for target {target} in Visit {visit_id}: {str(e)}")
                visibility, times = [], []
        else:
            visibility, times = [], []
        
        visit_data.append({
            'target': target,
            'start': datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%SZ"),
            'stop': datetime.strptime(stop_time, "%Y-%m-%dT%H:%M:%SZ"),
            'visibility': visibility,
            'times': [Time(t, format='mjd').datetime for t in times]
        })

    try:
        target#print(target)
    except:
        target = None#standard_ = True
    # Create and save figure for this visit
    create_visit_figure(visit_data, visit_id, target)



# # Prepare data for plotting
# plot_data = []

# for visit in root.findall('ns:Visit', namespace):
#     for obs_seq in visit.findall('ns:Observation_Sequence', namespace):
#         target_elem = obs_seq.find('./ns:Observational_Parameters/ns:Target', namespace)
#         target = target_elem.text if target_elem is not None and target_elem.text else "No target"
        
#         start_elem = obs_seq.find('./ns:Observational_Parameters/ns:Timing/ns:Start', namespace)
#         stop_elem = obs_seq.find('./ns:Observational_Parameters/ns:Timing/ns:Stop', namespace)

#         print(target, start_elem.text, stop_elem.text)
        
#         if start_elem is None or stop_elem is None:
#             print(f"Warning: Missing timing data for sequence {obs_seq.find('ns:ID', namespace).text}")
#             continue
        
#         start_time = start_elem.text
#         stop_time = stop_elem.text
        
#         if target != "No target":
#             try:
#                 visibility, times = check_visibility(target, start_time, stop_time)
#             except Exception as e:
#                 print(f"Error checking visibility for target {target}: {str(e)}")
#                 visibility, times = [], []
#         else:
#             visibility, times = [], []
        
#         plot_data.append({
#             'target': target,
#             'start': datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%SZ"),
#             'stop': datetime.strptime(stop_time, "%Y-%m-%dT%H:%M:%SZ"),
#             'visibility': visibility,
#             'times': [Time(t, format='mjd').datetime for t in times]
#         })

# # Create the plot
# fig, ax = plt.subplots(figsize=(15, 10))

# for i, data in enumerate(plot_data):
#     color = 'green' if all(data['visibility']) else 'red' if data['visibility'] else 'gray'
#     ax.plot([data['start'], data['stop']], [i, i], color=color, linewidth=5)
    
#     # Plot visibility
#     if data['visibility']:
#         vis_times = data['times']
#         vis_values = [i if v else None for v in data['visibility']]
#         ax.scatter(vis_times, [i] * len(vis_times), c=vis_values, cmap='RdYlGn', vmin=0, vmax=1, s=20)

# # Set y-axis labels
# ax.set_yticks(range(len(plot_data)))
# ax.set_yticklabels([d['target'] for d in plot_data])

# # Format x-axis
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
# fig.autofmt_xdate()

# plt.title('Target Visibility During Observation Sequences') 
# plt.xlabel('Time') 
# plt.ylabel('Targets')
# plt.tight_layout()

# #Save the figure as a PNG file
# output_dir = PACKAGEDIR # Replace with your desired output directory 
# if not os.path.exists(output_dir): 
#     os.makedirs(output_dir)

# output_file = os.path.join(output_dir, 'target_visibility.png') 
# plt.savefig(output_file, dpi=300, bbox_inches='tight') 
# plt.close(fig)

# print(f"Figure saved as {output_file}")
