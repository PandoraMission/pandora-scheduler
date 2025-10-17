import xml.etree.ElementTree as ET
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from astropy.time import Time
import os

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Parse the XML file
fname = f'{PACKAGEDIR}/data/calendar_Pandora_Schedule_27Aug2024_Claude.xml'
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
    st_name = target[:-2] if target.endswith('b') else target
    if not target.startswith('Gaia'):
        v_data = pd.read_csv(tar_vis_path+f'{st_name}/Visibility for {st_name}.csv')
    else:
        v_data = pd.read_csv(aux_vis_path+f'{target}/Visibility for {target}.csv')

    # Filter data for the observation period
    mask = (v_data['Time(MJD_UTC)'] >= start_mjd) & (v_data['Time(MJD_UTC)'] <= stop_mjd)
    period_data = v_data[mask]
    
    return period_data['Visible'].tolist(), period_data['Time(MJD_UTC)'].tolist()

# Prepare data for plotting
plot_data = []

for visit in root.findall('ns:Visit', namespace):
    for obs_seq in visit.findall('ns:Observation_Sequence', namespace):
        target_elem = obs_seq.find('./ns:Observational_Parameters/ns:Target', namespace)
        target = target_elem.text if target_elem is not None and target_elem.text else "No target"
        
        start_elem = obs_seq.find('./ns:Observational_Parameters/ns:Timing/ns:Start', namespace)
        stop_elem = obs_seq.find('./ns:Observational_Parameters/ns:Timing/ns:Stop', namespace)

        print(target, start_elem.text, stop_elem.text)
        
        if start_elem is None or stop_elem is None:
            print(f"Warning: Missing timing data for sequence {obs_seq.find('ns:ID', namespace).text}")
            continue
        
        start_time = start_elem.text
        stop_time = stop_elem.text
        
        if target != "No target":
            try:
                visibility, times = check_visibility(target, start_time, stop_time)
            except Exception as e:
                print(f"Error checking visibility for target {target}: {str(e)}")
                visibility, times = [], []
        else:
            visibility, times = [], []
        
        plot_data.append({
            'target': target,
            'start': datetime.strptime(start_time, "%Y-%m-%dT%H:%M:%SZ"),
            'stop': datetime.strptime(stop_time, "%Y-%m-%dT%H:%M:%SZ"),
            'visibility': visibility,
            'times': [Time(t, format='mjd').datetime for t in times]
        })

# Create the plot
fig, ax = plt.subplots(figsize=(15, 10))

for i, data in enumerate(plot_data):
    color = 'green' if all(data['visibility']) else 'red' if data['visibility'] else 'gray'
    ax.plot([data['start'], data['stop']], [i, i], color=color, linewidth=5)
    
    # Plot visibility
    if data['visibility']:
        vis_times = data['times']
        vis_values = [i if v else None for v in data['visibility']]
        ax.scatter(vis_times, [i] * len(vis_times), c=vis_values, cmap='RdYlGn', vmin=0, vmax=1, s=20)

# Set y-axis labels
ax.set_yticks(range(len(plot_data)))
ax.set_yticklabels([d['target'] for d in plot_data])

# Format x-axis
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
fig.autofmt_xdate()

plt.title('Target Visibility During Observation Sequences') 
plt.xlabel('Time') 
plt.ylabel('Targets')
plt.tight_layout()

#Save the figure as a PNG file
output_dir = PACKAGEDIR # Replace with your desired output directory 
if not os.path.exists(output_dir): 
    os.makedirs(output_dir)

output_file = os.path.join(output_dir, 'target_visibility.png') 
plt.savefig(output_file, dpi=300, bbox_inches='tight') 
plt.close(fig)

print(f"Figure saved as {output_file}")
