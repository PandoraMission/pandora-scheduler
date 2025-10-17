#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import xml.dom.minidom
from astropy.time import Time
from datetime import datetime, timedelta
import os
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import warnings
import helper_codes_claude as hcc

warnings.filterwarnings("ignore")

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Global variables
obs_seq_duration = None
occ_seq_limit = None
dt = None
occultation_sequence_limit = None
too_short_sequences = None
tar_path = None
aux_path = None
tv_st = None
tv_sp = None

def load_data():
    schedule_path = f'{PACKAGEDIR}/data/Pandora_Schedule_2025-08-04_3months_29Aug2024.csv'
    tar_path = f'{PACKAGEDIR}/data/Pandora_Target_List_Top20_14May2024.csv'
    aux_path = f'{PACKAGEDIR}/data/aux_list.csv'
    
    t_list = pd.read_csv(tar_path)
    a_list = pd.read_csv(aux_path)
    sch = pd.read_csv(schedule_path)
    
    return t_list, a_list, sch

def create_calendar(sch):
    cal = ET.Element('ScienceCalendar', xmlns="/pandora/calendar/")
    meta = ET.SubElement(cal, 'Meta',
                         Valid_From=f"{sch['Observation Start'][0]}",
                         Expires=f"{sch['Observation Stop'][len(sch)-1]}",
                         Calendar_Weights='0.0, 0.0, 1.0',
                         Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
                         Keepout_Angles='90.0, 40.0, 63.0',
                         Observation_Sequence_Duration_hrs=f'{obs_seq_duration}',
                         Removed_Sequences_Shorter_Than_min=f'{too_short_sequences}',
                         Created=f'{str(datetime.now())}',
                         Delivery_Id='')
    return cal

def process_visit(cal, sch, i, t_list, a_list, tar_path, tar_path_ALL, aux_path):
    t_name = sch['Target'][i]
    st_name = t_name if t_name.startswith('Gaia') else t_name[:-2]
    
    visit = ET.SubElement(cal, 'Visit')
    id0 = ET.SubElement(visit, "ID")
    id0.text = f'{("0"*(4-len(str(i))))+str(i)}'
    
    start = pd.to_datetime(sch['Observation Start'][i])
    stop = pd.to_datetime(sch['Observation Stop'][i])
    
    v_data, targ_info, i_flag = get_visibility_data(t_name, st_name, t_list, a_list)
    
    ra, dec = get_coordinates(st_name, targ_info)
    
    v_time, v_flag = process_visibility(v_data, start, stop)
    
    # Remove short sequences by setting their visibility to 0
    # print(f"Before removal of too short sequences: {sum(v_flag)} visible periods")
    _, v_flag_updated = hcc.remove_short_sequences(v_time, v_flag, too_short_sequences)
    # print(f"After removal of too short sequences: {sum(v_flag_updated)} visible periods")
    
    # Recalculate v_change based on the updated v_flag
    v_change = np.where(v_flag_updated[:-1] != v_flag_updated[1:])[0]
    
    process_sequences(visit, v_change, v_flag_updated, v_time, t_name, ra, dec, i_flag, tar_path, tar_path_ALL, aux_path)

    return visit


def get_visibility_data(t_name, st_name, t_list, a_list):
    if not t_name.startswith('Gaia'):
        v_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/{st_name}/Visibility for {st_name}.csv')
        targ_info = t_list.loc[t_list['Planet Name'] == t_name]
        i_flag = 1
    else:
        v_data = pd.read_csv(f'{PACKAGEDIR}/data/aux_targets/{t_name}/Visibility for {t_name}.csv')
        targ_info = a_list.loc[a_list['Star Name'] == t_name]
        i_flag = 0
    
    return v_data, targ_info, i_flag

def get_coordinates(st_name, targ_info):
    try:
        star_sc = SkyCoord.from_name(st_name)
        ra, dec = star_sc.ra.deg, star_sc.dec.deg
    except:
        ra, dec = targ_info['RA'].iloc[0], targ_info['DEC'].iloc[0]
    return ra, dec

def process_visibility(v_data, start, stop):
    v_time_all = Time(v_data["Time(MJD_UTC)"], format="mjd", scale="utc").to_value("datetime")
    start_dt = pd.to_datetime(start)
    stop_dt = pd.to_datetime(stop)
    v_time = v_time_all[(v_time_all >= start_dt) & (v_time_all <= stop_dt)]
    v_flag = np.asarray(v_data['Visible'])[(v_time_all >= start_dt) & (v_time_all <= stop_dt)]
    return v_time, v_flag

def process_sequences(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, tar_path, tar_path_ALL, aux_path):
    if len(v_change) == 0:
        process_full_visibility(visit, v_time, t_name, ra, dec, i_flag)
    else:
        process_partial_visibility(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, tar_path, tar_path_ALL, aux_path)

def process_full_visibility(visit, v_time, t_name, ra, dec, i_flag):
    st, sp = v_time[0], v_time[-1]
    n = (sp - st) / dt
    sps = [st + (dt * (i + 1)) for i in range(int(n))]
    if int(n) < n:
        sps.append(sp)
    if sps[-1] == v_time[-1]:
        sps[-1] = v_time[-2]
    
    sps_all = list(np.hstack((st, sps)))
    for s in range(len(sps_all) - 1):
        priority = get_priority(i_flag, sps_all[s], sps_all[s+1])
        hcc.observation_sequence(visit, f'{("0"*(3-len(str(s+1))))+str(s+1)}',
                                 t_name, priority, 
                                 sps_all[s].strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                 sps_all[s+1].strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                 ra, dec)

def process_partial_visibility(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, tar_path, tar_path_ALL, aux_path):
    occultation_ranges = get_occultation_times(v_change, v_flag, v_time)
    oc_starts = [range[0] for range in occultation_ranges]
    oc_stops = [range[1] for range in occultation_ranges]
    
    info, flag = hcc.sch_occ(oc_starts, oc_stops, tar_path, sort_key='closest', prev_obs=[ra, dec], 
                             tar_path=tar_path, tar_path_ALL=tar_path_ALL, aux_path=aux_path)
    
    if not flag:
        info, flag = hcc.sch_occ(oc_starts, oc_stops, aux_path, sort_key='closest', prev_obs=[ra, dec], 
                                 tar_path=tar_path, tar_path_ALL=tar_path_ALL, aux_path=aux_path)
    
    if not flag:
        print(f"\nMore targets are necessary to cover these occultation times. Neither target_list nor aux_list work. {v_time[0]}, {v_time[-1]}")
    
    schedule_observation_sequences(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, info)

def get_occultation_times(v_change, v_flag, v_time):
    oc_ranges = []
    if not v_flag[-1]:
        v_change = np.append(v_change, len(v_time) - 2)
    
    if not v_flag[0]:
        oc_ranges.append([v_time[0], v_time[v_change[0]]])
    
    for v in range(len(v_change) - 1):
        if not v_flag[v_change[v+1]]:
            oc_ranges.append([v_time[v_change[v] + 1], v_time[v_change[v+1]]])
    
    return hcc.break_long_sequences(oc_ranges, occultation_sequence_limit)

def schedule_observation_sequences(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, info):
    oc_tr = 0
    v_change = np.append(v_change, len(v_time) - 1)  # Add the last index
    seq_counter = 1
    schedule_data = []

    for v in range(len(v_change)):
        if v == 0:
            st = v_time[0]
        else:
            st = v_time[v_change[v-1] + 1]
        
        sp = v_time[v_change[v]]

        print(v, st, sp)
        
        # Calculate duration in minutes
        duration = (sp - st).total_seconds() / 60
        
        if v_flag[v_change[v]] == 1:  # Visible period
            if duration > dt.total_seconds() / 60:
                # Break into sequences of length dt
                current = st
                while current < sp:
                    next_val = min(current + dt, sp)
                    priority = get_priority(i_flag, current, next_val)
                    hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                             t_name, priority, 
                                             current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                             next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                             ra, dec)
                    schedule_data.append((t_name, current, next_val, True))
                    seq_counter += 1
                    current = next_val
            else:
                # Execute single observation sequence
                priority = get_priority(i_flag, st, sp)
                hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                         t_name, priority, 
                                         st.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                         sp.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                         ra, dec)
                schedule_data.append((t_name, st, sp, True))
                seq_counter += 1
        else:  # Non-visible period (v_flag[v_change[v]] == 0)
            if duration > occultation_sequence_limit.total_seconds() / 60:
                # Break into sequences of length occultation_sequence_limit
                current = st
                while current < sp:
                    next_val = min(current + occultation_sequence_limit, sp)
                    target_, ra_, dec_ = info['Target'][oc_tr], info['RA'][oc_tr], info['DEC'][oc_tr]
                    priority = '0'
                    hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                             target_, priority, 
                                             current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                             next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                             ra_, dec_)
                    schedule_data.append((t_name, current, next_val, True))
                    seq_counter += 1
                    current = next_val
                oc_tr += 1
            else:
                # Execute single observation sequence for occultation target
                target_, ra_, dec_ = info['Target'][oc_tr], info['RA'][oc_tr], info['DEC'][oc_tr]
                priority = '0'
                hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                         target_, priority, 
                                         st.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                         sp.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                         ra_, dec_)
                schedule_data.append((target_, st, sp, False))
                seq_counter += 1
                oc_tr += 1

    # output_dir = "path/to/your/output/directory"
    # filename = f"schedule_visualization_{t_name}"  # You can customize this as needed
    # hcc.visualize_schedule(schedule_data, PACKAGEDIR, filename)

    return visit






def get_priority(i_flag, start, stop):
    if i_flag:
        return '2' if np.any((tv_st <= stop) & (tv_sp >= start)) else '1'
    else:
        return '0'

def save_xml(cal):
    etstr = ET.tostring(cal, xml_declaration=True)
    dom = xml.dom.minidom.parseString(etstr)
    pretty_xml_as_string = dom.toprettyxml()
    
    with open(f'{PACKAGEDIR}/data/calendar_Pandora_Schedule_27Aug2024.xml', 'w+') as f:
        f.write(pretty_xml_as_string)

def main():
    global obs_seq_duration, occ_seq_limit, dt, occultation_sequence_limit, too_short_sequences, tar_path, aux_path, tv_st, tv_sp

    obs_seq_duration, occ_seq_limit = hcc.general_parameters()
    dt = timedelta(minutes=obs_seq_duration)
    occultation_sequence_limit = timedelta(minutes=occ_seq_limit + 1.)
    too_short_sequences = 5
    
    t_list, a_list, sch = load_data()
    cal = create_calendar(sch)
    
    tar_path = f'{PACKAGEDIR}/data/Pandora_Target_List_Top20_14May2024.csv'
    tar_path_ALL = f'{PACKAGEDIR}/data/Pandora_Target_List_Top40_16Feb2024_Top40_SDM.csv'
    aux_path = f'{PACKAGEDIR}/data/aux_list.csv'
    
    for i in tqdm(range(8,9)):#len(sch))):
        t_name = sch['Target'][i]
        st_name = t_name if t_name.startswith('Gaia') else t_name[:-2]
        
        if not t_name.startswith('Gaia'):
            tv_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/{st_name}/{t_name}/Visibility for {t_name}.csv')
            tv_st = Time(tv_data['Transit_Start'], format='mjd', scale='utc').to_value('datetime')
            tv_sp = Time(tv_data['Transit_Stop'], format='mjd', scale='utc').to_value('datetime')
        
        process_visit(cal, sch, i, t_list, a_list, tar_path, tar_path_ALL, aux_path)
    
    save_xml(cal)

if __name__ == "__main__":
    main()

