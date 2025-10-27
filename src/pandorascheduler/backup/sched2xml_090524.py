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

# VK BEGIN:
import helper_codes
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
#def sch_occ(starts, stops, list_path, sort_key=None, **kwargs):
#
# VK BEGIN: adding explicit prev_obs keyword
def sch_occ(starts, stops, list_path, sort_key=None, prev_obs = None):#, position = 0):#**kwargs):
# VK END
    
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
            
            # importlib.reload(helper_codes)
            o_df, d_flag = helper_codes.schedule_occultation_targets(v_names, starts, stops, path_, o_df, o_list, try_occ_targets)#, position)

        return o_df, d_flag

            # if d_flag_copy:
            #     return o_df_copy.drop(columns=['Visibility']), d_flag_copy
        
        # #empty dataframe to hold visibility information for multiple targets
        # v_ = np.asarray(np.zeros(len(starts)), dtype=bool)
        # v_df = pd.DataFrame([v_])
        # vis_df = pd.DataFrame(columns=range(len(starts))).astype(bool)
        # d_flag = False
        # for n in range(len(v_names)):
        #     try:
        #         if (list_path == tar_path) or (list_path == tar_path_ALL):
        #             vis=pd.read_csv(f"{PACKAGEDIR}/data/targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
        #         elif list_path == aux_path:
        #             vis=pd.read_csv(f"{PACKAGEDIR}/data/aux_targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
        #         vis_ = vis['Time(MJD_UTC)']
                
        #         #array to hold visibility info for this target
        #         v_ar = np.asarray(np.zeros(len(starts)), dtype=bool)
        #         v_ar_any = np.asarray(np.zeros(len(starts)), dtype=bool)
        #         #iterate through occultation periods and see if v_names[n] is visible for all of them
        #         vis_f = False
        #         for s in range(len(starts)):
        #             #get occultation time window
        #             win = vis.index[(vis_ >= starts[s]) & (vis_ <= stops[s])].tolist()
        #             if np.all(vis['Visible'][win] == 1):# or vis_ratio >= 0.6: 
        #                 # UNCOMMENT THE REST OF THE LINE ABOVE FOR TARGETS THAT ARE NOT VISIBLE FOR MANY HOURS AND ALLOW "OCCULTATION"
        #                 # TARGET TO BE CONSIDERED AS VISIBLE >= 60% OF 90 MINUTES!!!
        #                 v_ar[s] = True               
        #             else:
        #                 vis_f = True

        #             # VK START
        #             # CHECK THE FRACTION OF TIME OCCULTATION TARGET IS VISIBILE DURING OCCULTATION 
        #             if len(vis['Visible'][win]) > 0.:
        #                 vis_ratio = len(vis['Visible'][win][vis['Visible'][win] == 1])/len(vis['Visible'][win])
        #             # VK END

        #             # print(v_names[n], Time(starts[s], format="mjd", scale="utc").to_value("datetime").strftime("%H:%M:%S"), \
        #             #     Time(stops[s], format="mjd", scale="utc").to_value("datetime").strftime("%H:%M:%S"), vis_ratio)
        #                     # len(vis['Visible'][win][vis['Visible'][win] == 1]), len(vis['Visible'][win]))#np.asarray(vis['Visible'][win]))

        #         # vis_df.loc[len(vis_df)] = v_ar
        #         # vis_df_sum = np.sum(np.asarray(vis_df), axis = 0)
        #         # print(n, vis_df_sum)
        #         # if np.all(vis_df_sum > 0):
        #         #     stop

        #         #if not visible for all times, check if any entry in v_df and this one cover the total occultation time
        #         if vis_f:
        #             if not d_flag:

        #                 # vis_df_sum = np.sum(np.asarray(vis_df), axis = 0)
        #                 # # print(n, vis_df_sum)
        #                 # if np.all(vis_df_sum > 0):
        #                 #     d_flag = True
        #                 # else:
        #                 #     vis_df.loc[len(vis_df)] = v_ar

        #                 v_arr=np.asarray(v_df)
        #                 overlap=np.where([np.all((v_arr+np.asarray(v_ar, dtype=bool))[i]) for i in range(len(v_arr))])[0]
        #                 if len(overlap) > 0:
        #                     #at least one entry has visibility that covers the total time along with the current target
        #                     #take both and enter them in o_df in their respective times, prefering the closer one
        #                     m=overlap[0]
        #                     v1=v_arr[m]
        #                     for s in range(len(starts)):
        #                         if v1[s]:
        #                             o_df['Target'][s]=v_names[m]
        #                             o_df['RA'][s]=str(ras[m])
        #                             o_df['DEC'][s]=str(decs[m])
        #                         else:
        #                             o_df['Target'][s]=v_names[n]
        #                             o_df['RA'][s]=str(ras[n])
        #                             o_df['DEC'][s]=str(decs[n])
        #                     d_flag=True
        #                 else:
        #                     #add the current visibility array to the master list
        #                     v_df.loc[len(v_df.index)] = v_ar
        #                     # vis_df.loc[len(vis_df)] = v_ar
        #             else:
        #                 # break
        #                 continue
                
        #         else:
        #             #If we made it here, the target is visible for the entire occultation time
        #             #since the list is already sorted, break and use this one!
        #             #add to the list of visible targets (not necessary, but keeps functionality with prev code)
        #             o_df['Target'][:]=v_names[n]
        #             o_df['RA'][:]=str(ras[n])
        #             o_df['DEC'][:]=str(decs[n])
        #             d_flag=True
        #             # print(st_name, ': ', n, v_names[n], v_names[n], 'Occ target visible ALL')
        #             break
        #             # THIS BREAK DOESNT SEEM TO WORK!!!! REPLACE WITH return o_df, d_flag
        #             # return o_df, d_flag

        #         # print(f"{st_name}: {n} ({v_names[n]}) not 100% visible, try next on the list")
            
        #     #If a target(s) on the list don't have visibility data, ignore them!
        #     except FileNotFoundError:
        #         continue    
        
        # #if there were not <=2 occultation targets that covered the time (cya clause)
        # if not d_flag:
        #     #only if considering aux targets, real last ditch effort here
        #     if list_path.endswith('aux_list.csv'):# or st_name.startswith('Gaia'):
        #         #check one last time to make sure there are no gaps
        #         v_arr=np.asarray(v_df)
        #         if np.all([np.any(v_arr[i]) for i in range(len(v_arr))]):
        #             #iterate and set the nearest that is visible for each window as the occultation target
        #             for m in range(len(v_arr)):
        #                 # VK START: IT LOOKS LIKE THE i IN v_arr[i] IS NOT PART OF THIS FOR LOOP. MAYBE A BUG?
        #                 #n=np.where(v_arr[i])[0][0]
        #                 # CHANGE i TO m
        #                 # VK END
        #                 n=np.where(v_arr[m])[0][0]
        #                 o_df['Target'][s]=v_names[n]
        #                 o_df['RA'][s]=str(ras[n])
        #                 o_df['DEC'][s]=str(decs[n])
        #             d_flag=True
        #             print('More than two auxiliary targets were needed to cover the occultation time.')

    # return o_df, d_flag

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
                   Calendar_Weights='0.0, 0.0, 1.0',
                   Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
                   Keepout_Angles='90.0, 40.0, 63.0',
                   Observation_Sequence_Duration_hrs = f'{dt}',
                   Removed_Sequences_Shorter_Than_min = f'{too_short_sequences}',
                   Created=f'{str(datetime.now())}',
#                   Author="P Bonney",
                   Delivery_Id='',
                   )

for i in tqdm(range(8, 9)):#, position = 0, leave = True):#len(sch))):#3)):#len(18,19)):#
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
    # VK BEGIN: simplify the following 4 lines
    # times=np.where(v_time >= start)
    # times_=np.where(v_time[times] <= stop)
    # v_time=v_time[times][times_]
    #visibility flag during those times
    #v_flag=np.asarray(v_data['Visible'])[times][times_]
    v_time = v_time_all[(v_time_all >= start) & (v_time_all <= stop)]
    v_flag = np.asarray(v_data['Visible'])[(v_time_all >= start) & (v_time_all <= stop)]
    # VK END


    # VK START: REMOVE SEQUENCES THAT ARE TOO SHORT:
    v_flag_update, positions = helper_codes.remove_short_sequences(v_flag, too_short_sequences)
    # if v_flag_update[-1] == 0.:
    #     v_flag_update[-1] = 1.
    v_flag = v_flag_update.copy()
    # VK END

    #figure out where the visibility changes (gives final element where the visibility is the same)
    v_change = np.where(v_flag[:-1] != v_flag[1:])[0]


    # for v in range(len(v_change)-1):
    #     # print(v_time[v_change[v]+1], v_time[v_change[v+1]])
    #     start_, end_ = v_time[v_change[v]], v_time[v_change[v+1]]
    #     print(start_, end_)


    # new_v_change = helper_codes.break_long_visibility_changes(v_change, max_sequence = obs_sequence_duration)

    # print(t_name, v_change[-1], len(v_time)-1, v_flag[-1])
    
    st=start
    sp=v_time[-1]
    
    #offset for observation sequence numbering
    obs_off=0
    
    #if full visibility
    if len(v_change) == 0:

        print('Target is visible for the entire visit')
        
        #break observation sequence into <= 90 minute blocks
        n=(sp-st)/dt
        
        sps=[st+(dt*(i+1)) for i in range(int(n))]
        if int(n) < n:
           sps.append(sp)

        if sps[-1] == v_time[-1]:
            sps[-1] = v_time[-2]

        #get priority
        #check for a visible transit if primary science target and visible
        #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit
        # if i_flag:
        #     pr = 2 if np.any((tv_st <= v_time[-1])*(tv_sp >= st)) else 1
        # else:
        #     pr=0
   
        sps_all = list(np.hstack((st, sps)))
        for s in range(len(sps_all)-1):

            if i_flag:
                pr = 2 if np.any((tv_st <= sps_all[s+1])*(tv_sp >= sps_all[s])) else 1
            else:
                pr=0

            aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(s+1))))+str(s+1)}', \
                t_name, pr, sps_all[s], sps_all[s+1], ra, dec)
            pass
    
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

        #add final index of v_time for iteration
        # v_change.tolist().append(len(v_time)-1)
        # v_change=np.array(v_change)
        #
        #
        # VK BEGIN:
        # print(t_name, v_change[-1], len(v_time)-1, v_flag[-1])
        if v_flag[-1] == 1:
            v_change = np.append(v_change, len(v_time)-2)
        # v_change = np.append(v_change, len(v_time)-2)
        # if v_change[-1] != len(v_time)-1:
        #     v_change = np.append(v_change, len(v_time)-2)
        # else:
        #     v_change[-1] = len(v_time)-2
        #
        #
        #visibility change tracker (for convenience)
        v_t = v_flag[v_change]
        # VK END

        #oc_times=[pd.date_range(oc_starts[o],oc_stops[o], freq='min') for o in range(len(oc_starts))]

        # VK BEGIN: BREAK OCCULTATION SEQUENCES LONGER THAN 90 MINUTES
        start_tmp, stop_tmp = [], []
        for ii in range(len(oc_stops)):
            ranges = helper_codes.break_long_sequences(oc_starts[ii], oc_stops[ii], occultation_sequence_limit)
            if len(ranges) > 1:
                for jj in range(len(ranges)):
                    start_tmp.append(ranges[jj][0])
                    stop_tmp.append(ranges[jj][1])
            else:
                start_tmp.append(oc_starts[ii])
                stop_tmp.append(oc_stops[ii])
        oc_starts_bak, oc_stops_bak = oc_starts.copy(), oc_stops.copy()
        oc_starts, oc_stops = start_tmp, stop_tmp
        # VK END

        #find an occultation target for this visit that will always be visible
        # VK BEGIN: there is no "nearest" in
        # info, flag = sch_occ(oc_starts, oc_stops, tar_path, sort_key = 'nearest', prev_obs=[ra,dec])
        info, flag = sch_occ(oc_starts, oc_stops, tar_path, sort_key = 'closest', prev_obs=[ra,dec])#, position = 1)
        # print()
        if flag:
            # tqdm.write('Find occultation target from target list itself...DONE!')
            print('\nFind occultation target from target list itself...DONE!')
        # VK END
        oc_flag=1
        if not flag:
            #: VK BEGIN: there is no "nearest" in
            # tqdm.write('Target list doesnt work, try aux list instead...')
            print('\nTarget list doesnt work, try aux list instead...') 
            info, flag = sch_occ(oc_starts, oc_stops, aux_path, sort_key = 'closest', prev_obs = [ra,dec])#, position = 2)
            if flag:
                # tqdm.write('Find occultation targets from aux list...DONE!')
                print('\nFind occultation targets from aux list...DONE!')
            # VK END
        if not flag:
            print("\nMore targets are necessary to cover these occultation times. Neither target_list nor aux_list work.", st, sp)
        oc_flag=0
        
        #schedule first observation sequence
        #occultation tracker
        oc_tr = 0
        
        st = start
        sp = v_time[v_change[0]]

        #case where the main target isn't visible at first
        if not v_t[0]:
            target_, start_, stop_, ra_, dec_ = info['Target'][0], \
                info['start'][0], info['stop'][0], \
                    info['RA'][0], info['DEC'][0]
        # VK START
        # THIS IS WIP, FOR THE CASE WHEN THE TARGET IS NOT VISIBLE FOR MANY HOURS IN THE BEGINNING. 
        # IN THAT CASE, sch_occ WILL ALLOW "OCCULTATION" TARGET TO BE VISIBLE FOR ONLY >=60% OF 90 MIN 
        # if not v_t[0]:
        #     if info['stop'][0] != info['start'][1]:
        #         target_, start_, stop_, ra_, dec_ = info['Target'][0], \
        #             info['start'][0], info['stop'][0], \
        #                 info['RA'][0], info['DEC'][0]
        #     elif info['stop'][0] == info['start'][1]:
        #         target_, start_, stop_, ra_, dec_ = info['Target'][0], \
        #             info['start'][0], oc_stops_bak[0], \
        #                 info['RA'][0], info['DEC'][0]
        # VK END
            priority_ = f'{oc_flag}'
            oc_tr += 1
            

        #case where the main target is visible at first
        else:
            target_, start_, stop_, ra_, dec_ = t_name,\
                f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(sp, "%Y-%m-%dT%H:%M:%SZ")}',\
                    f'{float(ra)}', f'{float(dec)}'
            #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit 
            # if i_flag:
            #         priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
            # else:
            #     priority_ = '0'

        # VK BEGIN: Create first observation sequence
        start_format, stop_format = Time(start_).to_value('datetime'), Time(stop_).to_value('datetime')
        if stop_format - start_format <= dt:

            if i_flag:
                priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
            else:
                priority_ = '0'

            aa = helper_codes.observation_sequence(visit, "001", target_, priority_, start_format, stop_format, ra_, dec_)
            long_sequence = 0
         # If first sequence longer than dt, break it into sections of length dt:
        else:
            nn = (stop_format - start_format)/dt
            sps = [start_format+(dt*(i+1)) for i in range(int(nn))]
            oc_starts_dt = list(np.sort(np.hstack((start_format, sps))))
            oc_stops_dt = list(np.sort(np.hstack((stop_format, sps))))
            long_sequence = 1
            for ii, jj in zip(oc_starts_dt, oc_stops_dt):

                if i_flag:
                    priority_ = '2' if np.any((tv_st <= jj)*(tv_sp >= ii)) else '1'
                else:
                    priority_ = '0'

                aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(long_sequence))))+str(long_sequence)}', target_, priority_, ii, jj, ra_, dec_)
                long_sequence += 1

        # print(target_, start_, stop_)
        
        #loop to schedule consecutive observation sequences
        for v in range(len(v_change)-1):

            st = v_time[v_change[v]+1]
            sp = v_time[v_change[v+1]]

            # nn = (sp - st)/dt
            # sps = [st+(dt*(i+1)) for i in range(int(nn))]
            # oc_starts_dt = list(np.sort(np.hstack((st, sps))))
            # oc_stops_dt = list(np.sort(np.hstack((sp, sps))))
            # long_sequence = 1
            # for ii, jj in zip(oc_starts_dt, oc_stops_dt):
            #     aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(long_sequence))))+str(long_sequence)}', target_, priority_, ii, jj, ra_, dec_)
            #     long_sequence += 1

            # if nn > 1:
            #     print(t_name, st, sp)

            #set elements for the target if target is visible for this sequence
            if v_t[v+1]: 
                target_, ra_, dec_ = t_name, f'{float(ra)}', f'{float(dec)}'
                #check for a visible transit if primary science target and visible
                #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit 
                if i_flag:
                    priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
                else:
                    priority_ = '0'

            #otherwise, set elements for an occultation target
            else:
                target_, ra_, dec_ = info['Target'][oc_tr], info['RA'][oc_tr], info['DEC'][oc_tr]
                priority_ = f'{oc_flag}'
                oc_tr+=1

            # Create the rest of the observation sequences
            os_i = v + long_sequence if long_sequence > 0 else (v + 2)
            start_format, stop_format = f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(sp, "%Y-%m-%dT%H:%M:%SZ")}'#Time(start_).to_value('datetime'), Time(stop_).to_value('datetime')
            aa = helper_codes.observation_sequence(visit, f'{("0"*(3-len(str(os_i))))+str(os_i)}', target_, priority_, start_format, stop_format, ra_, dec_)

    #     print(target_, start_, stop_)
    # print()

etstr=ET.tostring(cal, xml_declaration=True)


from xml.dom import minidom
dom = minidom.parseString(etstr)

#dom = xml.dom.minidom.parseString(etstr)
pretty_xml_as_string = dom.toprettyxml()
f=open(f'{PACKAGEDIR}/data/calendar_Pandora_Schedule_27Aug2024.xml','w+')#test.xml', 'w+')
f.write(pretty_xml_as_string)
f.close()