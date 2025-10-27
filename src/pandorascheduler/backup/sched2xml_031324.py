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

# VK BEGIN: remove warnings:
import warnings
warnings.filterwarnings("ignore")
# VK END

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
schedule_path=f'{PACKAGEDIR}/data/Pandora_Schedule_0.0_0.0_1.0_2025-05-25.csv'
tar_vis_path=f'{PACKAGEDIR}/data/targets/'
aux_vis_path=f'{PACKAGEDIR}/data/aux_targets/'
tar_path=f'{PACKAGEDIR}/data/target_list.csv'
aux_path=f'{PACKAGEDIR}/data/aux_list.csv'
t_list=pd.read_csv(tar_path)
a_list=pd.read_csv(aux_path)
sch=pd.read_csv(schedule_path)
author='Paul Bonney'

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
def sch_occ(starts, stops, list_path, sort_key=None, prev_obs = None):#**kwargs):
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
                o_list['sk_dif']=dif
                o_list.sort_values(by='sk_dif').reset_index(drop=True)
                
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
        
        #empty dataframe to hold visibility information for multiple targets
        v_=np.asarray(np.zeros(len(starts)), dtype=bool)
        v_df = pd.DataFrame([v_])
        d_flag=False
        for n in range(len(v_names)):
            try:
                #vis=pd.read_csv(f"{PACKAGEDIR}/data/oc_targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
                # VK BEGIN: there is no oc_targets directory, trying with targets or aux_targets directories:
                if list_path == tar_path:
                    vis=pd.read_csv(f"{PACKAGEDIR}/data/targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
                elif list_path == aux_path:
                    vis=pd.read_csv(f"{PACKAGEDIR}/data/aux_targets/{v_names[n]}/Visibility for {v_names[n]}.csv")
                # VK END
                vis_=vis['Time(MJD_UTC)']
                
                #array to hold visibility info for this target
                v_ar=np.asarray(np.zeros(len(starts)), dtype=bool)
                #iterate through occultation periods and see if v_names[n] is visible for all of them
                vis_f=False
                for s in range(len(starts)):
                    
                    #get occultation time window
                    ast=set(vis.index[vis_ >= starts[s]])
                    bsp=set(vis.index[vis_ <= stops[s]])
                    win=list(ast.intersection(bsp))
                                      
                    if np.all(vis['Visible'][win] == 1):
                        v_ar[s]=True
                        
                    else:
                        vis_f=True
                
                #if not visible for all times, check if any entry in v_df and this one cover the 
                #   total occultation time
                if vis_f:
                    if not d_flag:
                        v_arr=np.asarray(v_df)
                        overlap=np.where([np.all((v_arr+np.asarray(v_ar, dtype=bool))[i]) for i in range(len(v_arr))])[0]
                        if len(overlap) > 0:
                            #at least one entry has visibility that covers the total time along with the
                            #   current target
                            #take both and enter them in o_df in their respective times, prefering the 
                            #   closer one
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
                    else:
                        continue
                
                else:
                    #If we made it here, the target is visible for the entire occultation time
                    #since the list is already sorted, break and use this one!
                    #add to the list of visible targets (not necessary, but keeps functionality with prev code)
                    o_df['Target'][:]=v_names[n]
                    o_df['RA'][:]=str(ras[n])
                    o_df['DEC'][:]=str(decs[n])
                    d_flag=True
                    break
                
                #print(st_name, ': ', n, v_names[n], ' not visible, try next on the list')
            
            #If a target(s) on the list don't have visibility data, ignore them!
            except FileNotFoundError:
                continue    
        
        #if there were not <=2 occultation targets that covered the time (cya clause)
        if not d_flag:
            #only if considering aux targets, real last ditch effort here
            if list_path.endswith('aux_list.csv'):
                #check one last time to make sure there are no gaps
                v_arr=np.asarray(v_df)
                if np.all([np.any(v_arr[i]) for i in range(len(v_arr))]):
                    #iterate and set the nearest that is visible for each window as the occultation target
                    for m in range(len(v_arr)):
                        n=np.where(v_arr[i])[0][0]
                        o_df['Target'][s]=v_names[n]
                        o_df['RA'][s]=str(ras[n])
                        o_df['DEC'][s]=str(decs[n])
                    d_flag=True
                    print('More than two auxiliary targets were needed to cover the occultation time.')
                    
        
    return o_df, d_flag



#max time for an observation sequence
observation_sequence_duration = 90 # minutes
dt = timedelta(minutes = observation_sequence_duration)

cal=ET.Element('ScienceCalendar', xmlns="/pandora/calendar/")
meta=ET.SubElement(cal, 'Meta', 
                   Valid_From=f"{sch['Observation Start'][0]}",
                   Expires=f"{sch['Observation Stop'][len(sch)-1]}",
                   Calendar_Weights='0.0, 0.0, 1.0',
                   Ephemeris='sma=6828.14, ecc=0.0, inc=97.2188, aop=0.0, raan=303.263, ta=0.0',
                   Keepout_Angles='90.0, 40.0, 63.0',
                   Created=f'{str(datetime.now())}',
#                   Author="P Bonney",
                   Delivery_Id='',
                   )

for i in tqdm(range(1,2)):#len(sch))):
    t_name=sch['Target'][i]
    st_name=t_name[:-2]
    
    #set visit number and visit element
    visit=ET.SubElement(cal,'Visit')
    id0 = ET.SubElement(visit, "ID")
    id0.text = f'{("0"*(4-len(str(i))))+str(i)}'
    
    start=datetime.strptime(sch['Observation Start'][i], "%Y-%m-%d %H:%M:%S")
    stop=datetime.strptime(sch['Observation Stop'][i], "%Y-%m-%d %H:%M:%S")
    
    #Get visibility data, replace if then with the flag later
    if not t_name.startswith('Gaia'):
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
    
    ra=targ_info['RA'].iloc[0]
    dec=targ_info['DEC'].iloc[0]

    # VK BEGIN: the original target_list.csv has incorrect RA & Dec. 
    # Instead of manually fixing these, I'll just get them from SkyCoord
    try:
        star_sc = SkyCoord.from_name(st_name)
        ra = star_sc.ra.deg
        dec = star_sc.dec.deg
    except:
        pass
    # VK END

    ####
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
        "TargetRA": f'{float(ra)}',
        "TargetDEC": f'{float(dec)}',
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
    ####
    
    #get times during this visit
    v_time_all = Time(v_data["Time(MJD_UTC)"], format="mjd", scale="utc").to_value("datetime")
    # VK BEGIN: simplify the following 4 lines
    # times=np.where(v_time >= start)
    # times_=np.where(v_time[times] <= stop)
    # v_time=v_time[times][times_]
    #visibility flag during those times
    #v_flag=np.asarray(v_data['Visible'])[times][times_]
    v_time = v_time_all[(v_time_all >= start) & (v_time_all <= stop)]
    v_time_mjd = v_data["Time(MJD_UTC)"][(v_time_all >= start) & (v_time_all <= stop)].values
    v_flag = np.asarray(v_data['Visible'])[(v_time_all >= start) & (v_time_all <= stop)]
    # VK END

    #figure out where the visibility changes (gives final element where the vis is the same)
    v_change = np.where(v_flag[:-1] != v_flag[1:])[0]
    
    st=start
    sp=v_time[-1]
    
    #offset for observation sequence numbering
    obs_off=0
    
    #if full visibility
    if len(v_change) == 0:
        
        #break observation sequence into <= 90 minute blocks
        n=(sp-st)/dt
        
        sps=[st+(dt*(i+1)) for i in range(int(n))]
        if int(n) < n:
            sps.append(sp)

        for s in range(-1, len(sps)-1):
            if s == -1:
                st = start
                sp = v_time[-1]
                stop_tmp_ = sps[0]
                idx1_tmp = "001"
            else:
                st = sps[s]
                sp = sps[s+1]
                stop_tmp_ = sp
                idx1_tmp = f'{("0"*(3-len(str(s+2))))+str(s+2)}'


            #get priority
            #check for a visible transit if primary science target and visible
            #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit
            if i_flag:
                pr = 2 if np.any((tv_st <= sp)*(tv_sp >= st)) else 1
            else:
                pr=0
            ###
            ###
            o_seq = ET.SubElement(visit,'Observation_Sequence')
            obs_seq_id = ET.SubElement(o_seq, "ID")
            obs_seq_id.text = idx1_tmp
            ###
            ### Observational Parameters
            obs_parameters = ET.SubElement(o_seq, "Observational_Parameters")
            observational_parameters = {
                #"ID": idx1_tmp,
                "Target": t_name,
                "Priority": f'{pr}',
                "Timing": ["Start", "Stop", f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(stop_tmp_, "%Y-%m-%dT%H:%M:%SZ")}'], 
                "Boresight": ["RA", "DEC", f'{float(ra)}', f'{float(dec)}'], 
            }
            for ii, jj in zip(observational_parameters.keys(), observational_parameters.values()):
                if (ii != "Timing") & (ii != "Boresight"):
                    obs_param_element = ET.SubElement(obs_parameters, ii)
                    obs_param_element.text = str(jj)
                else:
                    obs_param_element = ET.SubElement(obs_parameters, ii)
                    for kk in range(2):
                        sub_element_tmp = ET.SubElement(obs_param_element, jj[kk])
                        sub_element_tmp.text = jj[kk+2]
            ###
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

            pass
    
    else:
        #identify occultation times
        oc_starts=[]
        oc_stops=[]

        if not v_flag[-1]:
            v_change=v_change.tolist()
            v_change.append(len(v_time)-1)
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

        oc_times=[pd.date_range(oc_starts[o],oc_stops[o], freq='min') for o in range(len(oc_starts))]
        #visibility change tracker (for convenience)
        v_t = v_flag[v_change]
        
        # VK BEGIN: Break up sequences longer than dt into sections of length dt:
        oc_starts_extend, oc_stops_extend = [], []
        #oc_starts_bak, oc_stops_bak = oc_starts, oc_stops
        for ii in range(len(oc_stops)):
            if oc_stops[ii] - oc_starts[ii] > dt:
                nn = (oc_stops[ii] - oc_starts[ii])/dt
                sps = [oc_starts[ii]+(dt*(i+1)) for i in range(int(nn))]
                oc_starts_extend = list(np.sort(np.hstack((oc_starts, sps))))
                oc_stops_extend = list(np.sort(np.hstack((oc_stops, sps))))
        # VK END

        #find an occultation target for this visit that will always be visible
        # VK BEGIN: there is no "nearest" in
        # info, flag = sch_occ(oc_starts, oc_stops, tar_path, sort_key = 'nearest', prev_obs=[ra,dec])
        info, flag = sch_occ(oc_starts, oc_stops, tar_path, sort_key = 'closest', prev_obs=[ra,dec])
        if flag:
            print('Find occultation targets from target_list itself...DONE')
        # VK END
        oc_flag=1
        if not flag:
            #: VK BEGIN: there is no "nearest" in
            # info, flag = sch_occ(oc_starts, oc_stops, aux_path, sort_key = 'nearest', prev_obs=[ra,dec])
            info, flag = sch_occ(oc_starts, oc_stops, aux_path, sort_key = 'closest', prev_obs=[ra,dec])
            print('target_list doesnt work, find occultation targets from aux_list instead...DONE')
            # VK END
        if not flag:
            print("More targets are necessary to cover these occultation times. Neither target_list nor aux_list work.")
        oc_flag=0
        
        # #add final index of v_time for iteration
        # v_change.tolist().append(len(v_time)-1)
        # v_change=np.array(v_change)
        
        #schedule first observation sequence
        
        #occultation tracker
        oc_tr=0
        
        st=start
        sp=v_time[v_change[0]]

        #case where the main target isn't visible at first
        if not v_t[0]:       

            #VBK BEGIN: break sequences longer than dt
            long_sequence = Time(info['stop'][0]).to_value('datetime') - Time(info['start'][0]).to_value('datetime')
            if long_sequence > observation_sequence_duration:

            (Time(info['stop'][0]).to_value('datetime') - Time(info['start'][0]).to_value('datetime'))*24*60




            target_, start_, stop_, ra_, dec_ = info['Target'][0], \
                info['start'][0], info['stop'][0], \
                    info['RA'][0], info['DEC'][0]
            priority_ = f'{oc_flag}'
            oc_tr += 1
        #case where the main target is visible at first
        else:
            target_, start_, stop_, ra_, dec_ = t_name,\
                f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(sp, "%Y-%m-%dT%H:%M:%SZ")}',\
                    f'{float(ra)}', f'{float(dec)}'
            #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit 
            if i_flag:
                    priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
            else:
                priority_ = '0'
        ###
        ###
        o_seq = ET.SubElement(visit,'Observation_Sequence')
        obs_seq_id = ET.SubElement(o_seq, "ID")
        obs_seq_id.text = "001"
        ###
        ### Observational Parameters
        obs_parameters = ET.SubElement(o_seq, "Observational_Parameters")
        observational_parameters = {
            "Target": target_,
            "Priority": priority_,
            "Timing": ["Start", "Stop", start_, stop_], 
            "Boresight": ["RA", "DEC", ra_, dec_], 
        }
        for ii, jj in zip(observational_parameters.keys(), observational_parameters.values()):
            if (ii != "Timing") & (ii != "Boresight"):
                obs_param_element = ET.SubElement(obs_parameters, ii)
                obs_param_element.text = str(jj)
            else:
                obs_param_element = ET.SubElement(obs_parameters, ii)
                for kk in range(2):
                    sub_element_tmp = ET.SubElement(obs_param_element, jj[kk])
                    sub_element_tmp.text = jj[kk+2]
        ###
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
        
        #loop to schedule consecutive observation sequences
        for v in range(len(v_change)-1):
            st=v_time[v_change[v]+1]
            sp=v_time[v_change[v+1]]

            os_i=v+2
            
            #set elements for the target if target is visible for this sequence
            if v_t[v+1]:
                target_, start_, stop_, ra_, dec_ = t_name, \
                    f'{datetime.strftime(st, "%Y-%m-%dT%H:%M:%SZ")}', f'{datetime.strftime(sp, "%Y-%m-%dT%H:%M:%SZ")}', \
                        f'{float(ra)}', f'{float(dec)}'
                #check for a visible transit if primary science target and visible
                #set flag 0 = non-primary target; 1 = primary target; 2 = in-transit 
                if i_flag:
                    priority_ = '2' if np.any((tv_st <= sp)*(tv_sp >= st)) else '1'
                else:
                    priority_ = '0'

            #otherwise, set elements for an occultation target
            else:
                st=start
                sp=v_time[v_change[0]]

                target_, start_, stop_, ra_, dec_ = info['Target'][oc_tr], \
                    info['start'][oc_tr], info['stop'][oc_tr], \
                        info['RA'][oc_tr], info['DEC'][oc_tr]
                priority_ = f'{oc_flag}'
                oc_tr+=1
            
            ###
            ###
            o_seq = ET.SubElement(visit,'Observation_Sequence')
            obs_seq_id = ET.SubElement(o_seq, "ID")
            obs_seq_id.text = f'{("0"*(3-len(str(os_i))))+str(os_i)}'
            ###
            ### Observational Parameters
            obs_parameters = ET.SubElement(o_seq, "Observational_Parameters")
            observational_parameters = {
                "Target": target_,
                "Priority": priority_,
                "Timing": ["Start", "Stop", start_, stop_], 
                "Boresight": ["RA", "DEC", ra_, dec_], 
            }

            for ii, jj in zip(observational_parameters.keys(), observational_parameters.values()):
                if (ii != "Timing") & (ii != "Boresight"):
                    obs_param_element = ET.SubElement(obs_parameters, ii)
                    obs_param_element.text = str(jj)
                else:
                    obs_param_element = ET.SubElement(obs_parameters, ii)
                    for kk in range(2):
                        sub_element_tmp = ET.SubElement(obs_param_element, jj[kk])
                        sub_element_tmp.text = jj[kk+2]

            ###
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
                
etstr=ET.tostring(cal, xml_declaration=True)

from xml.dom import minidom
dom = minidom.parseString(etstr)

#dom = xml.dom.minidom.parseString(etstr)
pretty_xml_as_string = dom.toprettyxml()
f=open(f'{PACKAGEDIR}/data/cal_pretty.xml', 'w+')
f.write(pretty_xml_as_string)
f.close()








