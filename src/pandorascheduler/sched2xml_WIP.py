#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
import helper_codes
import helper_codes_aux as hcc
import importlib
import warnings
warnings.filterwarnings("ignore")
# VK END

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
schedule_path = f'{PACKAGEDIR}/data/Pandora_Schedule_2026-02-05_to_2027-02-05_092425.csv'#Pandora_Schedule_2025-08-04_to_2026-08-03_last.csv'#Pandora_Schedule_2025-08-04_3months_29Aug2024.csv'#Pandora_Schedule_2025-08-04_2months.csv'#Pandora_Schedule_2025-08-04.csv'
tar_vis_path = f'{PACKAGEDIR}/data/targets/'
aux_vis_path = f'{PACKAGEDIR}/data/aux_targets/'
tar_path = f'{PACKAGEDIR}/data/exoplanet_targets.csv'#primary-exoplanet-extended_targets.csv'#Pandora_Target_List_Top20_14May2024.csv'#target_list_top20_16Feb2024.csv'
# tar_path_ALL = f'{PACKAGEDIR}/data/primary-exoplanet-extended_targets.csv'#Pandora_Target_List_Top40_16Feb2024_Top40_SDM.csv'
aux_path = f'{PACKAGEDIR}/data/aux_list_new.csv'
t_list = pd.read_csv(tar_path)
a_list = pd.read_csv(aux_path)
sch = pd.read_csv(schedule_path)

# target_definition_files = ['primary-exoplanet-extended', 'auxiliary-exoplanet-reduced', 'auxiliary-standard', 'occultation-standard', \
#     'monitoring-standard', 'secondary-exoplanet']
target_definition_files = ['exoplanet', 'auxiliary-exoplanet', 'auxiliary-standard', 'monitoring-standard', 'secondary-exoplanet', 'occultation-standard']
target_definition_files = ['exoplanet', 'auxiliary-standard', 'monitoring-standard', 'occultation-standard']

t_list = pd.read_csv(f"{PACKAGEDIR}/data/{target_definition_files[0]}_targets.csv")

# author='VK'

#save a stripped version as csv for LLNL
save_csv = False

#function to schedule an target during occultation
#takes a start/stop in datetime and 
#list_path is a path to a csv file with target info in it like aux_list.csv
#visibility data needs to be in oc_targets for now
#if 'closest' is given as sort_key, prev_obs needs to be given as a kwarg
#def sch_occ(starts, stops, list_path, sort_key=None, **kwargs):

#
def sch_occ_new(starts, stops, st, sp, list_path, sort_key=None, prev_obs = None):#, position = 0):#**kwargs):
    
    #build empty dataframe except for starts and stops
    e_sched = [['',datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"),datetime.strftime(stops[s], "%Y-%m-%dT%H:%M:%SZ"), '', ''] for s in range(len(starts))]
    o_df = pd.DataFrame(e_sched,columns=["Target","start","stop", "RA", "DEC"])
    
    #convert to mjd to compare with visibility data
    starts=Time(starts, format='datetime').to_value('mjd')
    stops=Time(stops, format='datetime').to_value('mjd')
    

    # if sort_key == None:
    #     #No occluded target scheduling, free time
    #     starts=starts.to_value('datetime')
    #     stops=stops.to_value('datetime')
    #     free = [["Free Time",datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"),datetime.strftime(stops[s], "%Y/%m/%d, %H:%M:%S"), '', ''] for s in range(len(starts))]
    #     o_df = pd.DataFrame(free, columns=["Target", "start", "stop", "RA", "DEC"])
    
    # else:       
    #     o_list = pd.read_csv(list_path)
    #     ras=o_list['RA']
    #     decs=o_list['DEC']
        
    #     if sort_key == 'closest':
    #         #sort name list based on minimized sky distance
    #         #prev_obs must be specified as an array consisting of ra and dec (in degrees) for the previous observation
    #         #this seeks to minimize slew distance to an occultation target, though actual slew sims are not performed
    #         try:
    #             po_sc=SkyCoord(unit='deg', ra=prev_obs[0], dec=prev_obs[1])
    #             oc_sc=[SkyCoord(unit='deg', ra=ras[n], dec=decs[n]) for n in range(len(ras))]

    #             dif=[oc_sc[n].separation(po_sc).deg for n in range(len(oc_sc))]
    #             o_list['sky_dif'] = dif
    #             o_list = o_list.sort_values(by='sky_dif').reset_index(drop=True)
                
    #         except NameError:
    #             print('No previous observation was specified, defaulting to random auxiliary target.')
    #             o_list.sample(frac=1).reset_index(drop=True)
    #     else:
    #         #default sort is random
    #         o_list.sample(frac=1).reset_index(drop=True)
        
    o_list = pd.read_csv(list_path)
    v_names = o_list['Star Name']
    v_names=np.array(v_names)
    #For prioritization via flag later
    o_flag = o_list['Priority']#['Flag']
    #Reload these
    ras=o_list['RA']
    decs=o_list['DEC']

    multi_target_occultation = True#False

    d_flag = False

    if multi_target_occultation:
        if (list_path == tar_path):# or (list_path == tar_path_ALL):
            path_ = f"{PACKAGEDIR}/data/targets"
            try_occ_targets = 'target list'
        elif list_path == aux_path:
            path_ = f"{PACKAGEDIR}/data/aux_targets"
            try_occ_targets = 'aux list'
        
        # importlib.reload(helper_codes)
        o_df, d_flag = helper_codes.schedule_occultation_targets(v_names, starts, stops, st, sp, path_, o_df, o_list, try_occ_targets)#, position)

    return o_df, d_flag



#max time for an observation sequence
obs_sequence_duration = 90 # minutes
occ_sequence_limit = 30 # minutes
obs_seq_duration, occ_seq_limit = helper_codes.general_parameters(obs_sequence_duration, occ_sequence_limit)
dt = timedelta(minutes = obs_seq_duration)
occultation_sequence_limit = timedelta(minutes = occ_seq_limit + 1.)

#Remove too short sequences
too_short_sequences = 5 # minutes

cal=ET.Element('ScienceCalendar', xmlns="/pandora/calendar/")
meta=ET.SubElement(cal, 'Meta', 
                   Valid_From=f"{sch['Observation Start'][0]}",
                   Expires=f"{sch['Observation Stop'][len(sch)-1]}",
                   Calendar_Weights='0.8, 0.0, 0.2',
                #    Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
                   Keepout_Angles='91.0, 25.0, 63.0',
                   Observation_Sequence_Duration_hrs_max = f'{dt}',
                   Removed_Sequences_Shorter_Than_min = f'{too_short_sequences}',
                   Created=f'{str(hcc.round_to_nearest_second(datetime.now()))}',
#                   Author="P Bonney",
                   Delivery_Id='',
                   )
#
#
#
for i in tqdm(range(10)):#len(sch))):#, position = 0, leave = True):#len(sch))):#3)):#len(18,19)):#

    logging.basicConfig(level=logging.INFO, format='%(message)s')#format='%(asctime)s - %(levelname)s - %(message)s')

    t_name = sch['Target'][i]

    if np.isfinite(sch.loc[i]['Transit Coverage']):
        exoplanet_tdf = True
    else:
        exoplanet_tdf = False
    # st_name = t_name if t_name.startswith('Gaia') else t_name[:-2]
    
    if t_name.endswith(('b', 'c', 'd', 'e', 'f')) and (t_name != 'EV_Lac'):
        st_name = t_name[:-2]
    elif t_name.endswith(('STD')):
        t_name = t_name[:-4]
        st_name = t_name
        # st_name = t_name[:-2]
    else:
        st_name = t_name
    
    #set visit number and visit element
    visit=ET.SubElement(cal,'Visit')
    id0 = ET.SubElement(visit, "ID")
    id0.text = f'{("0"*(4-len(str(i))))+str(i)}'
    
    start = datetime.strptime(sch['Observation Start'][i], "%Y-%m-%d %H:%M:%S")
    stop = datetime.strptime(sch['Observation Stop'][i], "%Y-%m-%d %H:%M:%S")
    
    #Get visibility data, replace if then with the flag later
    # if not t_name.startswith('Gaia'):# or t_name.startswith('Free'):
    # if t_name.endswith(('b', 'c', 'd', 'e', 'f')):
    # del v_data
    if t_name in t_list['Planet Name'].values and exoplanet_tdf:
        v_data = pd.read_csv(tar_vis_path+f'{st_name}/Visibility for {st_name}.csv')
        tmp_idx = t_list.index[t_list['Planet Name'] == t_name].tolist()
        targ_info = t_list.loc[[tmp_idx[0]]]#t_list.loc[(t_list['Planet Name'] == t_name)]
        i_flag = 1
        tv_data = pd.read_csv(tar_vis_path+f'{st_name}/{t_name}/Visibility for {t_name}.csv')
        tv_st = Time(tv_data['Transit_Start'], format='mjd', scale='utc').to_value('datetime')
        tv_sp = Time(tv_data['Transit_Stop'], format='mjd', scale='utc').to_value('datetime')
    elif exoplanet_tdf == False and t_name != 'Free Time' and not t_name.startswith(('WARNING')):#t_name in a_list['Star Name'].values:
        v_data = pd.read_csv(aux_vis_path+f'{st_name}/Visibility for {t_name}.csv')
        tmp_idx = a_list.index[(a_list['Star Name'] == t_name)].tolist()# & (pd.isnull(a_list['Planet Name']))].tolist()
        targ_info = pd.DataFrame(a_list.loc[tmp_idx[0]]).T#a_list.loc[[tmp_idx[0]]]#a_list.loc[(a_list['Star Name'] == t_name) & (a_list['Planet Name'].notna())]
        i_flag = 0
    elif t_name == 'Free Time':
        continue
    elif t_name.startswith(('WARNING')):#t_name == 'STD':
        print(f'-------> WARNING: need visible STD <--------')
        continue
    else:
        print(f"No visibility data for {t_name}. Stop code")
        break
        xxxx
    # elif:
    #     v_data = pd.read_csv(aux_vis_path+f'{t_name}/Visibility for {t_name}.csv')
    #     targ_info = a_list.loc[(a_list['Star Name'] == t_name) & (a_list['Planet Name'].notna())]#a_list.loc[(a_list['Star Name'] == t_name)]
    #     i_flag = 0

    # if t_name == 'TRAPPIST-1 f':
    #     print('TRAPPIST-1 f')

    # VK BEGIN: try getting RA & Dec from SkyCoord
    try:
        ra = targ_info['RA'].iloc[0]
        dec = targ_info['DEC'].iloc[0]
    except:
        star_sc = SkyCoord.from_name(st_name)
        ra = star_sc.ra.deg
        dec = star_sc.dec.deg
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
        # print('Target is visible for the entire visit')
        tqdm.write(f'{st} to {sp}: No occultations needed; {t_name} is visible for the entire visit')
        
        #break observation sequence into <= 90 minute blocks
        n = (sp - st)/dt
        
        sps=[st+(dt*(i+1)) for i in range(int(n))]
        if int(n) < n:
           sps.append(sp)

        if sps[-1] == v_time[-1]:
            sps[-1] = v_time[-2]
   
        sps_all = list(np.hstack((st, sps)))
        for s in range(len(sps_all)-1):

            if i_flag:
                pr = 2 if np.any((tv_st <= sps_all[s+1])*(tv_sp >= sps_all[s])) else 1
            else:
                pr=0

            aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(s+1))))+str(s+1)}', \
                t_name, pr, sps_all[s], sps_all[s+1], ra, dec, targ_info)
            # print('xxx')
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

        # for ii, jj in zip(oc_starts, oc_stops):
        #     print(ii, jj)

        # VK BEGIN: BREAK OCCULTATION SEQUENCES LONGER THAN 90 MINUTES
        break_occ_seq_longer_than_occultation_sequence_limit = True
        if break_occ_seq_longer_than_occultation_sequence_limit:
            start_tmp, stop_tmp = [], []
            for ii in range(len(oc_stops)):
                ranges = helper_codes.break_long_sequences(oc_starts[ii], oc_stops[ii], occultation_sequence_limit)
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
        def find_occultation_target(oc_starts, oc_stops, st, sp, tar_path, aux_path, ra, dec):
            # logging.info(f"Searching for occultation targets from {st} to {sp}")
            # tqdm.write(f"{st} to {sp}: Searching for occultation targets from {st} to {sp}")
            # # Try to find a target from aux_list
            # info, flag = sch_occ_new(oc_starts, oc_stops, st, sp, aux_path, sort_key = 'closest', prev_obs = [ra,dec])
            # if flag:
            #     tqdm.write(f"{st} to {sp}:     Found occultation target from aux list")
            # # logging.info(f"From aux list? {flag}")
            # oc_flag=1
            # if not flag:
            #     # If still not found, try tar_path
            #     info, flag = sch_occ_new(oc_starts, oc_stops, st, sp, tar_path, sort_key = 'closest', prev_obs = [ra,dec])#, position = 2)
            #     if flag:
            #          tqdm.write(f"{st} to {sp}:         Found occultation target from target list itself")
            #     # logging.info(f"From tar list? {flag}")

            # Try to find a target from tar_list
            info, flag = sch_occ_new(oc_starts, oc_stops, st, sp, tar_path, sort_key = 'closest', prev_obs=[ra,dec])
            if flag:
                tqdm.write(f"{st} to {sp}:     Found occultation target from target list itself")
            # logging.info(f"From target list itself? {flag}")
            oc_flag=1
            if not flag:
                # If still not found, try tar_path
                info, flag = sch_occ_new(oc_starts, oc_stops, st, sp, aux_path, sort_key = 'closest', prev_obs = [ra,dec])#, position = 2)
                if flag:
                     tqdm.write(f"{st} to {sp}:         Found occultation target from aux list")
                # logging.info(f"From tar list? {flag}")
            
            if flag:
                target_info = {
                    'name': info['Target'][0],
                    'ra': info['RA'][0],
                    'dec': info['DEC'][0]
                }
                # logging.info(f"Occultation targets found: {target_info['name']}")
                return info, flag#target_info
            else:
                tqdm.write(f"{st} to {sp}:             No suitable occultation targets found!!!!")
                # logging.warning(f"No suitable occultation targets found for period {st} to {sp}")
                return None

        info, flag = find_occultation_target(oc_starts, oc_stops,  st, sp, tar_path, aux_path, ra, dec)

        # #find an occultation target for this visit that will always be visible
        # # VK BEGIN: there is no "nearest" in
        # # info, flag = sch_occ(oc_starts, oc_stops, tar_path, sort_key = 'nearest', prev_obs=[ra,dec])
        # info, flag = sch_occ(oc_starts, oc_stops, tar_path, sort_key = 'closest', prev_obs=[ra,dec])#, position = 1)
        # # print()
        # if flag:
        #     # tqdm.write('Find occultation target from target list itself...DONE!')
        #     print('\nFind occultation target from target list itself...DONE!')
        # # VK END
        # oc_flag=1
        # if not flag:
        #     #: VK BEGIN: there is no "nearest" in
        #     # tqdm.write('Target list doesnt work, try aux list instead...')
        #     print('\nTarget list doesnt work, try aux list instead...') 
        #     info, flag = sch_occ(oc_starts, oc_stops, aux_path, sort_key = 'closest', prev_obs = [ra,dec])#, position = 2)
        #     if flag:
        #         # tqdm.write('Find occultation targets from aux list...DONE!')
        #         print('\nFind occultation targets from aux list...DONE!')
        #     # VK END
        # if not flag:
        #     print("\nMore targets are necessary to cover these occultation times. Neither target_list nor aux_list work.", st, sp)

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
                while current < sp: # break observation sequences longer than 90 min
                    if 1 == 1:
                    # try:
                        next_val = min(current + dt, sp)
                        priority = get_priority(i_flag, current, next_val)
                        aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}', \
                            t_name, priority, 
                            current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                            ra, dec, targ_info)
                        # logging.info(f"Vis sequence: {hcc.round_to_nearest_second(current)} to {hcc.round_to_nearest_second(next_val)}......DONE!")
                        # print()
                    else:
                    # except Exception as e:
                        # logging.info(f"current: {current}, next_val: {next_val}, priority: {priority}")
                        # logging.info(f"ra: {ra}, dec: {dec}")
                        # logging.info(f"t_name: {t_name}, seq_counter: {seq_counter}")
                        # logging.info(f"targ_info: {targ_info}")
                        logging.error(f"Error adding visible sequence: {str(e)}")
                        # print('xxx')
                    seq_counter += 1
                    current = next_val

            else:  # Non-visible period (occultation)
                if break_occ_seq_longer_than_occultation_sequence_limit:
                    # logging.info(f"Occultation duration: {hcc.round_to_nearest_second(st)} to {hcc.round_to_nearest_second(sp)}")
                    current = st
                    while current < sp:
                        next_val = min(current + occultation_sequence_limit, sp)
                        try:
                            occ_targ_info = t_list.loc[(t_list['Star Name'] == info['Target'][oc_tr])]
                            if occ_targ_info.empty:
                                occ_targ_info = a_list.loc[(a_list['Star Name'] == info['Target'][oc_tr])]
                            aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                                info['Target'][oc_tr], '0', 
                                current.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                next_val.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                                info['RA'][oc_tr], info['DEC'][oc_tr], occ_targ_info)
                            # logging.info(f"Occ sequence: {hcc.round_to_nearest_second(current)} to {hcc.round_to_nearest_second(next_val)}...DONE")
                            # print()
                            oc_tr += 1
                            seq_counter += 1
                        except Exception as e:
                            logging.error(f"Error adding occultation sequence: {str(e)}")
                            # print('xxx')
                        current = next_val
                else:
                    try:
                        occ_targ_info = t_list.loc[(t_list['Star Name'] == info['Target'][oc_tr])]
                        if occ_targ_info.empty:
                            occ_targ_info = a_list.loc[(a_list['Star Name'] == info['Target'][oc_tr])]
                        aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(seq_counter))))+str(seq_counter)}',
                            info['Target'][oc_tr], '0', 
                            st.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                            sp.strftime("%Y-%m-%dT%H:%M:%SZ"), 
                            info['RA'][oc_tr], info['DEC'][oc_tr], occ_targ_info)
                        # logging.info(f"Occ sequence: {hcc.round_to_nearest_second(st)} to {hcc.round_to_nearest_second(sp)}...DONE")
                        oc_tr += 1
                        seq_counter += 1
                    except Exception as e:
                        logging.error(f"Error adding occultation sequence: {str(e)}")
                        # print('xxx')
            # logging.info(f"Done with observing sequence {seq_counter}: {hcc.round_to_nearest_second(st)} to {hcc.round_to_nearest_second(sp)}")


        # #case where the main target isn't visible at first
        # if not v_t[0]:
        #     target_, start_, stop_, ra_, dec_ = info['Target'][0], \
        #         info['start'][0], info['stop'][0], \
        #             info['RA'][0], info['DEC'][0]
        #     priority_ = f'{oc_flag}'
        #     oc_tr += 1
            
        # #case where the main target is visible at first
        # else:
        #     target_, start_, stop_, ra_, dec_ = t_name,\
        #         f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(sp, "%Y-%m-%dT%H:%M:%SZ")}',\
        #             f'{float(ra)}', f'{float(dec)}'

        # # VK BEGIN: Create first observation sequence
        # start_format, stop_format = Time(start_).to_value('datetime'), Time(stop_).to_value('datetime')
        # if stop_format - start_format <= dt:
        #     if i_flag:
        #         priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
        #     else:
        #         priority_ = '0'

        #     aa = helper_codes.observation_sequence(visit, "001", target_, priority_, start_format, stop_format, ra_, dec_)
        #     long_sequence = 0
        # # If first sequence longer than dt, break it into sections of length dt:
        # else:
        #     nn = (stop_format - start_format)/dt
        #     sps = [start_format+(dt*(i+1)) for i in range(int(nn))]
        #     oc_starts_dt = list(np.sort(np.hstack((start_format, sps))))
        #     oc_stops_dt = list(np.sort(np.hstack((stop_format, sps))))
        #     long_sequence = 1
        #     for ii, jj in zip(oc_starts_dt, oc_stops_dt):

        #         if i_flag:
        #             priority_ = '2' if np.any((tv_st <= jj)*(tv_sp >= ii)) else '1'
        #         else:
        #             priority_ = '0'

        #         aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(long_sequence))))+str(long_sequence)}', target_, priority_, ii, jj, ra_, dec_)
        #         long_sequence += 1
        
        # # First sequence done! Now loop to schedule consecutive observation sequences
        # for v in range(len(v_change)-1):

        #     st = v_time[v_change[v]+1]
        #     sp = v_time[v_change[v+1]]
        #     start_format, stop_format = f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(sp, "%Y-%m-%dT%H:%M:%SZ")}'
        #     os_i = v + long_sequence if long_sequence > 0 else (v + 2)

        #     #set elements for the target if target is visible for this sequence
        #     if v_t[v+1]: 
        #         target_, ra_, dec_ = t_name, f'{float(ra)}', f'{float(dec)}'
        #         #check for a visible transit if primary science target and visible
        #         #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit 
        #         if i_flag:
        #             priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
        #         else:
        #             priority_ = '0'

        #         aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(os_i))))+str(os_i)}', target_, priority_, start_format, stop_format, ra_, dec_)

        #     # #otherwise, set elements for an occultation target
        #     else:
        #         nn = (sp - st)/occultation_sequence_limit
        #         print('FIX long_sequence, it is slipping!!!!!!!!!')
        #         if nn <= 1:
        #             target_, ra_, dec_ = info['Target'][oc_tr], info['RA'][oc_tr], info['DEC'][oc_tr]
        #             priority_ = f'{oc_flag}'
        #             oc_tr+=1
        #             aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(os_i))))+str(os_i)}', target_, priority_, start_format, stop_format, ra_, dec_)
        #         else:
        #             sps = [st + (occultation_sequence_limit*(i + 1)) for i in range(int(nn))]
        #             oc_starts_dt = list(np.sort(np.hstack((st, sps))))
        #             oc_stops_dt = list(np.sort(np.hstack((sp, sps))))
        #             long_sequence = 1
        #             for ii, jj in zip(oc_starts_dt, oc_stops_dt):
        #                 target_, ra_, dec_ = info['Target'][oc_tr], info['RA'][oc_tr], info['DEC'][oc_tr]
        #                 priority_ = f'{oc_flag}'
        #                 aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(os_i + long_sequence))))+str(os_i + long_sequence)}', target_, priority_, ii, jj, ra_, dec_)
        #                 long_sequence += 1
        #                 oc_tr+=1
        #                 # helper_codes.print_element_from_xml(aa)
        

# for child in visit:
#     print(child.tag, child.attrib)

# def float_to_str(element):
#     for el in element.iter():
#         for k, v in el.attrib.items():
#             if isinstance(v, float):
#                 el.set(k, str(v))
#         if el.text and isinstance(el.text, float):
#             el.text = str(el.text)

# float_to_str(cal)

def convert_to_string(element):
    for el in element.iter():
        for k, v in el.attrib.items():
            if isinstance(v, (int, float, np.int64, np.integer, np.floating)):
                el.set(k, str(v))
        if el.text and isinstance(el.text, (int, float, np.int64, np.integer, np.floating)):
            el.text = str(el.text)


convert_to_string(cal)
etstr=ET.tostring(cal, xml_declaration=True)

from xml.dom import minidom
dom = minidom.parseString(etstr)

#dom = xml.dom.minidom.parseString(etstr)
pretty_xml_as_string = dom.toprettyxml()
f=open(f'{PACKAGEDIR}/data/Pandora_science_calendar.xml','w+')#test.xml', 'w+')
f.write(pretty_xml_as_string)
f.close()

# # After the main loop, add a summary of occultation target times
# logging.info("Summary of occultation target observation times:")
# for target, time in occ_target_times.items():
#     logging.info(f"{target}: {time}")