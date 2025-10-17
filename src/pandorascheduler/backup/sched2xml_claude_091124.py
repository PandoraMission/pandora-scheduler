#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 14:09:15 2023

@author: paul
"""

import numpy as np
import pandas as pd
import xml
import xml.etree.ElementTree as ET
from astropy.time import Time
from datetime import datetime, timedelta
import os
from tqdm import tqdm
from astropy.coordinates import SkyCoord
import random
import json
import logging

# VK BEGIN:
import helper_codes_claude as hcc
import importlib
import warnings
warnings.filterwarnings("ignore")
# VK END

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
schedule_path=f'{PACKAGEDIR}/data/Pandora_Schedule_2025-08-04_3months_29Aug2024.csv'#Pandora_Schedule_2025-08-04_2months.csv'#Pandora_Schedule_2025-08-04.csv'
tar_vis_path=f'{PACKAGEDIR}/data/targets/'
aux_vis_path=f'{PACKAGEDIR}/data/aux_targets/'
tar_path=f'{PACKAGEDIR}/data/Pandora_Target_List_Top20_14May2024.csv'#target_list_top20_16Feb2024.csv'
aux_path=f'{PACKAGEDIR}/data/aux_list.csv'
tar_path_ALL = f'{PACKAGEDIR}/data/Pandora_Target_List_Top40_16Feb2024_Top40_SDM.csv'
t_list=pd.read_csv(tar_path)
a_list=pd.read_csv(aux_path)
sch=pd.read_csv(schedule_path)
# author='Paul Bonney'

#save a stripped version as csv for LLNL
save_csv=False

#function to schedule an target during occultation
#takes a start/stop in datetime and 
#list_path is a path to a csv file with target info in it like aux_list.csv
#visibility data needs to be in oc_targets for now
#if 'closest' is given as sort_key, prev_obs needs to be given as a kwarg
def sch_occ(starts, stops, list_path, sort_key=None, prev_obs = None):#, position = 0):#**kwargs):
    
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
        o_flag = o_list['Flag']
        #Reload these
        ras=o_list['RA']
        decs=o_list['DEC']

        multi_target_occultation = True#False

        d_flag = False

        if multi_target_occultation:
            if (list_path == tar_path) or (list_path == tar_path_ALL):
                path_ = f"{PACKAGEDIR}/data/targets"
                try_occ_targets = 'target list'
            elif list_path == aux_path:
                path_ = f"{PACKAGEDIR}/data/aux_targets"
                try_occ_targets = 'aux list'
            
            # importlib.reload(helper_codes_claude)
            o_df, d_flag = hcc.schedule_occultation_targets(v_names, starts, stops, path_, o_df, o_list, try_occ_targets)#, position)

        return o_df, d_flag

#max time for an observation sequence
obs_sequence_duration = 90 # minutes
occ_sequence_limit = 30 # minutes
obs_seq_duration, occ_seq_limit = hcc.general_parameters(obs_sequence_duration, occ_sequence_limit)
dt = timedelta(minutes = obs_seq_duration)
occultation_sequence_limit = timedelta(minutes = occ_seq_limit + 1.)

#Remove too short sequences
too_short_sequences = 5 # minutes

cal=ET.Element('ScienceCalendar', xmlns="/pandora/calendar/")
meta=ET.SubElement(cal, 'Meta', 
                   Valid_From=f"{sch['Observation Start'][0]}",
                   Expires=f"{sch['Observation Stop'][len(sch)-1]}",
                   Calendar_Weights='0.0, 0.0, 1.0',
                   Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
                   Keepout_Angles='90.0, 40.0, 63.0',
                   Observation_Sequence_Duration_hrs = f'{dt}',
                   Removed_Sequences_Shorter_Than_min = f'{too_short_sequences}',
                   Created=f'{str(datetime.now())}',
#                   Author="P Bonney",
                   Delivery_Id='',
                   )

for i in tqdm(range(7,9), desc=f"Processing visit"):#, position = 0, leave = True):#len(sch))):#3)):#len(18,19)):#

    logging.basicConfig(level=logging.INFO, format='%(message)s')#format='%(asctime)s - %(levelname)s - %(message)s')

    t_name=sch['Target'][i]
    st_name = t_name if t_name.startswith('Gaia') else t_name[:-2]
    
    #set visit number and visit element
    visit=ET.SubElement(cal,'Visit')
    id0 = ET.SubElement(visit, "ID")
    id0.text = f'{("0"*(4-len(str(i))))+str(i)}'
    
    start = datetime.strptime(sch['Observation Start'][i], "%Y-%m-%d %H:%M:%S")
    stop = datetime.strptime(sch['Observation Stop'][i], "%Y-%m-%d %H:%M:%S")
    
    #Get visibility data, replace if then with the flag later
    if not t_name.startswith('Gaia'):# or t_name.startswith('Free'):
        v_data=pd.read_csv(tar_vis_path+f'{st_name}/Visibility for {st_name}.csv')
        targ_info=t_list.loc[(t_list['Planet Name'] == t_name)]
        i_flag=1
        tv_data=pd.read_csv(tar_vis_path+f'{st_name}/{t_name}/Visibility for {t_name}.csv')
        tv_st=Time(tv_data['Transit_Start'], format='mjd', scale='utc').to_value('datetime')
        tv_sp=Time(tv_data['Transit_Stop'], format='mjd', scale='utc').to_value('datetime')
    else:
        v_data=pd.read_csv(aux_vis_path+f'{t_name}/Visibility for {t_name}.csv')
        targ_info=a_list.loc[(a_list['Star Name'] == t_name)]
        i_flag=0
    
    # VK BEGIN: try getting RA & Dec from SkyCoord
    try:
        star_sc = SkyCoord.from_name(st_name)
        ra = star_sc.ra.deg
        dec = star_sc.dec.deg
    except:
        ra=targ_info['RA'].iloc[0]
        dec=targ_info['DEC'].iloc[0]
    # VK END
    
    #get times during this visit
    v_time_all = Time(v_data["Time(MJD_UTC)"], format="mjd", scale="utc").to_value("datetime")
    v_time = v_time_all[(v_time_all >= start) & (v_time_all <= stop)]
    v_time = np.vectorize(hcc.round_to_nearest_second)(v_time)
    v_flag = np.asarray(v_data['Visible'])[(v_time_all >= start) & (v_time_all <= stop)]

    # VK START: REMOVE SEQUENCES THAT ARE TOO SHORT:
    v_flag_update, positions = hcc.remove_short_sequences(v_flag, too_short_sequences)
    # if v_flag_update[-1] == 0.:
    #     v_flag_update[-1] = 1.
    v_flag = v_flag_update.copy()
    # VK END

    #figure out where the visibility changes (gives final element where the visibility is the same)
    v_change = np.where(v_flag[:-1] != v_flag[1:])[0]
    
    st = start
    sp = v_time[-1]
    
    #offset for observation sequence numbering
    obs_off=0
    
    #if full visibility
    def full_visibility():
        # print(f'Target is visible for the entire visit ({st} to {sp})')
        
        #break observation sequence into <= 90 minute blocks
        n = (sp - st)/dt
        
        sps=[st+(dt*(i+1)) for i in range(int(n))]
        if int(n) < n:
           sps.append(sp)

        # if sps[-1] == v_time[-1]:
        #     sps[-1] = v_time[-2]
   
        sps_all = list(np.hstack((st, sps)))
        for s in range(len(sps_all)-1):

            if i_flag:
                pr = 2 if np.any((tv_st <= sps_all[s+1])*(tv_sp >= sps_all[s])) else 1
            else:
                pr=0

            aa = hcc.observation_sequence(visit, f'{("0"*(3-len(str(s+1))))+str(s+1)}', \
                t_name, pr, sps_all[s], sps_all[s+1], ra, dec)
            pass
    if len(v_change) == 0:
        full_visibility_ = full_visibility()
    
    #if NOT full visibility
    else:
        #identify occultation times
        oc_starts=[]
        oc_stops=[]

        if not v_flag[-1]:
            v_change=v_change.tolist()
            v_change.append(len(v_time)-2)
            v_change=np.array(v_change)
        
        if not v_flag[0]:
            #append in first section if occluded
            oc_starts.append(v_time[0])
            oc_stops.append(v_time[v_change[0]])

        for v in range(len(v_change)-1):
            #loop through the rest of v_change
            if not v_flag[v_change][v+1]:
                #only append if v_flag is 0 (i.e. occluded)
                oc_starts.append(v_time[v_change[v]+1])
                oc_stops.append(v_time[v_change[v+1]])


        if v_flag[-1] == 1:
            v_change = np.append(v_change, len(v_time)-2)

        #visibility change tracker (for convenience)
        v_t = v_flag[v_change]

        # VK BEGIN: BREAK OCCULTATION SEQUENCES LONGER THAN 90 MINUTES
        start_tmp, stop_tmp = [], []
        for ii in range(len(oc_stops)):
            ranges = hcc.break_long_sequences(oc_starts[ii], oc_stops[ii], occultation_sequence_limit)
            if len(ranges) > 1:
                for jj in range(len(ranges)):
                    start_tmp.append(ranges[jj][0])
                    stop_tmp.append(ranges[jj][1] - 0.*timedelta(minutes=1))
            else:
                start_tmp.append(oc_starts[ii])
                stop_tmp.append(oc_stops[ii] - 0*timedelta(minutes=1))
        oc_starts_bak, oc_stops_bak = oc_starts.copy(), oc_stops.copy()
        oc_starts, oc_stops = start_tmp, stop_tmp
        # VK END

        #find an occultation target for this visit that will always be visible
        def find_occultation_target(oc_starts, oc_stops, tar_path, tar_path_ALL, aux_path, ra, dec):
            # logging.info(f"Searching for occultation targets from {st} to {sp}")
            # Try to find a target from tar_path
            info, flag = hcc.sch_occ(oc_starts, oc_stops, tar_path, tar_path, tar_path_ALL, aux_path, sort_key = 'closest', prev_obs=[ra,dec])
            # logging.info(f"From target list itself? {flag}")
            
            oc_flag=1
            # if not flag:
            #     # If not found, try tar_path_ALL
            #     info, flag = hcc.sch_occ(start, stop, tar_path_ALL, sort_key='closest', prev_obs=[ra, dec], 
            #                              tar_path=tar_path, tar_path_ALL=tar_path_ALL, aux_path=aux_path)
            #     logging.info(f"Search in tar_path_ALL result: flag={flag}")
            
            if not flag:
                # If still not found, try aux_path
                info, flag = hcc.sch_occ(oc_starts, oc_stops, aux_path, tar_path, tar_path_ALL, aux_path, sort_key = 'closest', prev_obs = [ra,dec])#, position = 2)
                # logging.info(f"From aux list? {flag}")
            
            if flag:
                target_info = {
                    'name': info['Target'][0],
                    'ra': info['RA'][0],
                    'dec': info['DEC'][0]
                }
                # logging.info(f"Occultation targets found: {target_info['name']}")
                return info, flag#target_info
            else:
                logging.warning(f"No suitable occultation targets found for period {st} to {sp}")
                return None

        info, flag = find_occultation_target(oc_starts, oc_stops, tar_path, tar_path_ALL, aux_path, ra, dec)

        oc_flag=0
        #schedule first observation sequence
        #occultation tracker
        oc_tr = 0
        
        st = start
        sp = v_time[v_change[0]]

        seq_counter = 1

        def get_priority(i_flag, start, stop):
            if i_flag:
                return '2' if np.any((tv_st <= stop) & (tv_sp >= start)) else '1'
            else:
                return '0'

        for v in range(len(v_change)):
            if v == 0:
                st = v_time[0]
            else:
                st = v_time[v_change[v-1] + 1]

            sp = v_time[v_change[v]]
            # duration = sp - st
            # logging.info(f"Processing observing sequence {v+1}/{len(v_change)}: {hcc.round_to_nearest_second(st)} to {hcc.round_to_nearest_second(sp)}")
        
            if v_flag[v_change[v]]:  # Visible period
                # logging.info(f"Visible period: {st} to {sp}")
                current = st
                while current < sp:
                    try:
                        next_val = min(current + dt, sp)
                        priority = get_priority(i_flag, current, next_val)
                        aa = hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}', \
                            t_name, priority, 
                            current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                            ra, dec)
                        # logging.info(f"Vis sequence: {hcc.round_to_nearest_second(current)} to {hcc.round_to_nearest_second(next_val)}......DONE!")
                        # print()
                    except Exception as e:
                        logging.error(f"Error adding visible sequence: {str(e)}")
                    seq_counter += 1
                    current = next_val

            else:  # Non-visible period (occultation)
                # logging.info(f"Occultation duration: {hcc.round_to_nearest_second(st)} to {hcc.round_to_nearest_second(sp)}")
                current = st
                while current < sp:
                    next_val = min(current + occultation_sequence_limit, sp)
                    try:
                        aa = hcc.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                            info['Target'][oc_tr], '0', 
                            current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                            next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                            info['RA'][oc_tr], info['DEC'][oc_tr])
                        # logging.info(f"Occ sequence: {hcc.round_to_nearest_second(current)} to {hcc.round_to_nearest_second(next_val)}...DONE")
                        # print()
                        oc_tr += 1
                        seq_counter += 1
                    except Exception as e:
                        logging.error(f"Error adding occultation sequence: {str(e)}")
                    current = next_val
            # logging.info(f"Done with observing sequence {seq_counter}: {hcc.round_to_nearest_second(st)} to {hcc.round_to_nearest_second(sp)}")
  
        #                 # hcc.print_element_from_xml(aa)

    # print(f"Done with visit {i}")
    print()
        
etstr=ET.tostring(cal, xml_declaration=True)

from xml.dom import minidom
dom = minidom.parseString(etstr)

#dom = xml.dom.minidom.parseString(etstr)
pretty_xml_as_string = dom.toprettyxml()
f=open(f'{PACKAGEDIR}/data/calendar_Pandora_Schedule_27Aug2024.xml','w+')#test.xml', 'w+')
f.write(pretty_xml_as_string)
f.close()