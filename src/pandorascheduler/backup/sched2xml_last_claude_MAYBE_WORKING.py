#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import logging
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

def process_visit(cal, sch, i, t_list, a_list, tar_path, tar_path_ALL, aux_path):
    logging.info(f"Processing visit {i}")
    try:
        t_name = sch['Target'][i]
        logging.info(f"Target name: {t_name}")
        st_name = t_name if t_name.startswith('Gaia') else t_name[:-2]
        
        visit = ET.SubElement(cal, 'Visit')
        id0 = ET.SubElement(visit, "ID")
        id0.text = f'{("0"*(4-len(str(i))))+str(i)}'
        
        start = pd.to_datetime(sch['Observation Start'][i])
        stop = pd.to_datetime(sch['Observation Stop'][i])
        logging.info(f"Visit period: {start} to {stop}")
        
        v_data, targ_info, i_flag = get_visibility_data(t_name, st_name, t_list, a_list)
        logging.info(f"Visibility data retrieved, i_flag: {i_flag}")
        
        ra, dec = get_coordinates(st_name, targ_info)
        logging.info(f"Coordinates: RA={ra}, Dec={dec}")
        
        v_time, v_flag = process_visibility(v_data, start, stop)
        logging.info(f"Visibility processed. v_time shape: {v_time.shape}, v_flag shape: {v_flag.shape}")
        
        v_change = np.where(v_flag[:-1] != v_flag[1:])[0]
        logging.info(f"Visibility change points: {v_change}")
        
        if np.all(v_flag == 1):
            logging.info("Processing full visibility")
            process_full_visibility(visit, v_time, t_name, ra, dec, i_flag)
        else:
            logging.info("Processing partial visibility")
            process_partial_visibility(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, tar_path, tar_path_ALL, aux_path)

        return visit
    except Exception as e:
        logging.error(f"Error in process_visit for visit {i}: {str(e)}")
        raise

def get_visibility_data(t_name, st_name, t_list, a_list):
    logging.info(f"Getting visibility data for {t_name}")
    try:
        if not t_name.startswith('Gaia'):
            visibility_file = f'{PACKAGEDIR}/data/targets/{st_name}/Visibility for {st_name}.csv'
            logging.info(f"Attempting to read visibility data from: {visibility_file}")
            v_data = pd.read_csv(visibility_file)
            logging.info(f"Searching for target info in t_list with Planet Name: {t_name}")
            targ_info = t_list.loc[t_list['Planet Name'] == t_name]
            i_flag = 1
        else:
            visibility_file = f'{PACKAGEDIR}/data/aux_targets/{t_name}/Visibility for {t_name}.csv'
            logging.info(f"Attempting to read visibility data from: {visibility_file}")
            v_data = pd.read_csv(visibility_file)
            logging.info(f"Searching for target info in a_list with Star Name: {t_name}")
            targ_info = a_list.loc[a_list['Star Name'] == t_name]
            i_flag = 0
        
        if v_data.empty:
            logging.error(f"Visibility data is empty for {t_name}")
            return None, None, i_flag
        
        if targ_info.empty:
            logging.error(f"Target info is empty for {t_name}")
            return None, None, i_flag
        
        logging.info(f"Visibility data shape: {v_data.shape}, Target info shape: {targ_info.shape}")
        return v_data, targ_info, i_flag
    
    except FileNotFoundError as e:
        logging.error(f"File not found: {str(e)}")
        return None, None, i_flag
    except Exception as e:
        logging.error(f"Unexpected error in get_visibility_data: {str(e)}")
        return None, None, i_flag

def get_coordinates(st_name, targ_info):
    try:
        ra = float(targ_info['RA'].iloc[0])
        dec = float(targ_info['DEC'].iloc[0])
    except (ValueError, IndexError):
        print(f"Warning: Could not get coordinates for {st_name} from target info. Using SkyCoord.")
        try:
            star_sc = SkyCoord.from_name(st_name)
            ra, dec = star_sc.ra.deg, star_sc.dec.deg
        except:
            print(f"Error: Could not resolve coordinates for {st_name}")
            ra, dec = None, None
    return ra, dec

def process_visibility(v_data, start, stop):
    v_time = Time(v_data["Time(MJD_UTC)"], format="mjd", scale="utc").to_value("datetime")
    v_flag = np.asarray(v_data['Visible'])
    
    # Filter to the specified time range
    mask = (v_time >= start) & (v_time <= stop)
    v_time = v_time[mask]
    v_flag = v_flag[mask]
    
    return v_time, v_flag

def process_full_visibility(visit, v_time, t_name, ra, dec, i_flag):
    st, sp = v_time[0], v_time[-1]
    total_duration = sp - st
    n = total_duration / dt
    
    current = st
    seq_counter = 1
    
    while current < sp:
        next_val = min(current + dt, sp)
        duration = next_val - current
        
        priority = get_priority(i_flag, current, next_val)
        
        hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                 t_name, priority, 
                                 current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                 next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                 ra, dec)
        
        seq_counter += 1
        current = next_val

def get_priority(i_flag, start, stop):
    if i_flag:
        return '2' if np.any((tv_st <= stop) & (tv_sp >= start)) else '1'
    else:
        return '0'
#

def process_partial_visibility(visit, v_change, v_flag, v_time, t_name, ra, dec, i_flag, tar_path, tar_path_ALL, aux_path):
    logging.info(f"Entering process_partial_visibility for {t_name}")
    seq_counter = 1
    v_change = np.append(v_change, len(v_time) - 1)  # Add the last index
    occultation_limit = timedelta(minutes=30)  # 30 minutes limit for occultation segments

    for v in range(len(v_change)):
        if v == 0:
            st = v_time[0]
        else:
            st = v_time[v_change[v-1] + 1]
        
        sp = v_time[v_change[v]]
        duration = sp - st

        logging.info(f"Processing segment {v+1}/{len(v_change)}: {st} to {sp}")

        if v_flag[v_change[v]]:  # Visible period
            logging.info(f"Visible period: {st} to {sp}")
            current = st
            while current < sp:
                next_val = min(current + dt, sp)
                priority = get_priority(i_flag, current, next_val)
                try:
                    hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                             t_name, priority, 
                                             current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                             next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                             ra, dec)
                    logging.info(f"Added observation sequence: {current} to {next_val}")
                except Exception as e:
                    logging.error(f"Error adding observation sequence: {str(e)}")
                seq_counter += 1
                current = next_val
        else:  # Non-visible period (occultation)
            logging.info(f"Non-visible period: {st} to {sp}")
            current = st
            while current < sp:
                next_val = min(current + occultation_limit, sp)
                occ_target = find_occultation_target(current, next_val, tar_path, tar_path_ALL, aux_path, ra, dec)
                if occ_target is not None:
                    try:
                        hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                                 occ_target['name'], '0', 
                                                 current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                                 next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                                 occ_target['ra'], occ_target['dec'])
                        logging.info(f"Added occultation sequence: {current} to {next_val}")
                        seq_counter += 1
                    except Exception as e:
                        logging.error(f"Error adding occultation sequence: {str(e)}")
                else:
                    logging.warning(f"No occultation target found for period: {current} to {next_val}")
                current = next_val

    logging.info(f"Completed process_partial_visibility for {t_name}")
    return visit



def find_occultation_target(start, stop, tar_path, tar_path_ALL, aux_path, ra, dec):
    logging.info(f"Searching for occultation target from {start} to {stop}")
    
    # Try to find a target from tar_path
    info, flag = hcc.sch_occ(start, stop, tar_path, sort_key='closest', prev_obs=[ra, dec], 
                             tar_path=tar_path, tar_path_ALL=tar_path_ALL, aux_path=aux_path)
    logging.info(f"Search in tar_path result: flag={flag}")
    
    # if not flag:
    #     # If not found, try tar_path_ALL
    #     info, flag = hcc.sch_occ(start, stop, tar_path_ALL, sort_key='closest', prev_obs=[ra, dec], 
    #                              tar_path=tar_path, tar_path_ALL=tar_path_ALL, aux_path=aux_path)
    #     logging.info(f"Search in tar_path_ALL result: flag={flag}")
    
    if not flag:
        # If still not found, try aux_path
        info, flag = hcc.sch_occ(start, stop, aux_path, sort_key='closest', prev_obs=[ra, dec], 
                                 tar_path=tar_path, tar_path_ALL=tar_path_ALL, aux_path=aux_path)
        logging.info(f"Search in aux_path result: flag={flag}")
    
    if flag:
        target_info = {
            'name': info['Target'][0],
            'ra': info['RA'][0],
            'dec': info['DEC'][0]
        }
        logging.info(f"Occultation target found: {target_info['name']}")
        return target_info
    else:
        logging.warning(f"No suitable occultation target found for period {start} to {stop}")
        return None



def save_xml(cal, output_path):
    etstr = ET.tostring(cal, xml_declaration=True)
    dom = xml.dom.minidom.parseString(etstr)
    pretty_xml_as_string = dom.toprettyxml()
    
    with open(output_path, 'w+') as f:
        f.write(pretty_xml_as_string)
    
    print(f"XML saved to {output_path}")

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        global obs_seq_duration, occ_seq_limit, dt, occultation_sequence_limit, too_short_sequences, tar_path, aux_path, tv_st, tv_sp, tar_path_ALL

        # Set up parameters
        obs_seq_duration, occ_seq_limit = hcc.general_parameters()
        dt = timedelta(minutes=obs_seq_duration)
        occultation_sequence_limit = timedelta(minutes=occ_seq_limit)
        too_short_sequences = 5  # minutes

        # Define file paths
        PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
        schedule_path = f'{PACKAGEDIR}/data/Pandora_Schedule_2025-08-04_3months_29Aug2024.csv'
        tar_path = f'{PACKAGEDIR}/data/Pandora_Target_List_Top20_14May2024.csv'
        aux_path = f'{PACKAGEDIR}/data/aux_list.csv'
        tar_path_ALL = f'{PACKAGEDIR}/data/Pandora_Target_List_Top40_16Feb2024_Top40_SDM.csv'

        # Load data
        t_list = pd.read_csv(tar_path)
        a_list = pd.read_csv(aux_path)
        sch = pd.read_csv(schedule_path)

        # Create XML structure
        cal = ET.Element('ScienceCalendar', xmlns="/pandora/calendar/")
        meta = ET.SubElement(cal, 'Meta', 
                             Valid_From=f"{sch['Observation Start'][0]}",
                             Expires=f"{sch['Observation Stop'][len(sch)-1]}",
                             Calendar_Weights='0.0, 0.0, 1.0',
                             Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
                             Keepout_Angles='90.0, 40.0, 63.0',
                             Observation_Sequence_Duration_hrs = f'{dt}',
                             Removed_Sequences_Shorter_Than_min = f'{too_short_sequences}',
                             Created=f'{str(datetime.now())}',
                             Delivery_Id='',
                             )

        # Process visits
        for i in tqdm(range(8,9)):#len(sch))):
            try:
                process_visit(cal, sch, i, t_list, a_list, tar_path, tar_path_ALL, aux_path)
            except Exception as e:
                logging.error(f"Error processing visit {i}: {str(e)}")

        # Save the XML
        output_path = f'{PACKAGEDIR}/data/calendar_Pandora_Schedule_27Aug2024.xml'
        save_xml(cal, output_path)

        logging.info("Schedule processing completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred during schedule processing: {str(e)}")

if __name__ == "__main__":
    main()




