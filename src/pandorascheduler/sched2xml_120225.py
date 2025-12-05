#!/usr/bin/env python3
"""Generate a Pandora science calendar XML from schedule CSV inputs."""

import logging
import os
from datetime import datetime, timedelta
from numbers import Number

import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from astropy.coordinates import SkyCoord
from astropy.time import Time
from xml.dom import minidom
from tqdm import tqdm
import warnings

import helper_codes
# import helper_codes_aux as hcc

warnings.filterwarnings("ignore")

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
# schedule_path = f'{PACKAGEDIR}/data/Pandora_Schedule_0.8_0.0_0.2_2026-02-05_to_2026-03-05.csv'
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(PACKAGEDIR, "data")

schedule_path = os.path.join(
    DATA_DIR, "baseline", "Pandora_Schedule_0.8_0.0_0.2_2026-02-05_to_2027-02-05.csv"
)
tar_vis_path = os.path.join(DATA_DIR, "targets")
aux_vis_path = os.path.join(DATA_DIR, "aux_targets")
tar_path = os.path.join(DATA_DIR, "exoplanet_targets.csv")
aux_path = os.path.join(DATA_DIR, "all_targets.csv")
occ_path = os.path.join(DATA_DIR, "occultation-standard_targets.csv")

a_list = pd.read_csv(aux_path)
sch = pd.read_csv(schedule_path)
occ_list = pd.read_csv(occ_path)

target_definition_files = ['exoplanet', 'auxiliary-standard', 'monitoring-standard', 'occultation-standard']
t_list = pd.read_csv(os.path.join(DATA_DIR, f"{target_definition_files[0]}_targets.csv"))


def sch_occ_new(starts, stops, visit_start, visit_stop, list_path, sort_key=None, prev_obs=None):
    """Build an occultation schedule for the provided intervals."""

    e_sched = [['',datetime.strftime(starts[s], "%Y-%m-%dT%H:%M:%SZ"), datetime.strftime(stops[s], "%Y-%m-%dT%H:%M:%SZ"), '', ''] for s in range(len(starts))]
    o_df = pd.DataFrame(e_sched, columns=["Target", "start", "stop", "RA", "DEC"])

    starts_mjd = Time(starts, format="datetime").to_value("mjd")
    stops_mjd = Time(stops, format="datetime").to_value("mjd")

    o_list = pd.read_csv(list_path)
    v_names = o_list["Star Name"].to_numpy()

    path_ = aux_vis_path
    try_occ_targets = "aux list"
    if list_path == tar_path:
        path_ = tar_vis_path
        try_occ_targets = "target list"
    elif list_path == aux_path:
        path_ = aux_vis_path
        try_occ_targets = "aux list"
    elif list_path == occ_path:
        path_ = aux_vis_path
        try_occ_targets = "occ list"

    o_df, d_flag = helper_codes.schedule_occultation_targets_new(
        v_names,
        starts_mjd,
        stops_mjd,
        visit_start,
        visit_stop,
        path_,
        o_df,
        o_list,
        try_occ_targets,
    )

    return o_df, d_flag

#max time for an observation sequence
obs_sequence_duration = 90 # minutes
occ_sequence_limit = 50 # minutes
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
                   Keepout_Angles='91.0, 25.0, 63.0',
                   Observation_Sequence_Duration_hrs_max = f'{dt}',
                   Removed_Sequences_Shorter_Than_min = f'{too_short_sequences}',
                   Created=f'{str(helper_codes.round_to_nearest_second(datetime.now()))}',
                   Delivery_Id='',
                   )
#
#
#
for i in tqdm(range(10)):#len(sch))):

    t_name = sch['Target'][i]

    exoplanet_tdf = bool(np.isfinite(sch.loc[i]['Transit Coverage']))

    if t_name.endswith(('b', 'c', 'd', 'e', 'f')) and (t_name != 'EV_Lac'):
        st_name = t_name[:-1]
    elif t_name.endswith(('STD')):
        t_name = t_name[:-4]
        st_name = t_name
    else:
        st_name = t_name
    
    #set visit number and visit element
    visit=ET.SubElement(cal,'Visit')
    id0 = ET.SubElement(visit, "ID")
    id0.text = f'{("0"*(4-len(str(i))))+str(i+1)}'
    
    start = datetime.strptime(sch['Observation Start'][i], "%Y-%m-%d %H:%M:%S")
    stop = datetime.strptime(sch['Observation Stop'][i], "%Y-%m-%d %H:%M:%S")
    
    if t_name in t_list['Planet Name'].values and exoplanet_tdf:
        v_data = pd.read_csv(
            os.path.join(tar_vis_path, st_name, f"Visibility for {st_name}.csv")
        )
        tmp_idx = t_list.index[t_list['Planet Name'] == t_name].tolist()
        targ_info = t_list.loc[[tmp_idx[0]]]
        i_flag = 1
        tv_data = pd.read_csv(
            os.path.join(
                tar_vis_path, st_name, t_name, f"Visibility for {t_name}.csv"
            )
        )

        tv_st = Time(tv_data['Transit_Start'], format='mjd', scale='utc').to_value('datetime')
        tv_sp = Time(tv_data['Transit_Stop'], format='mjd', scale='utc').to_value('datetime')

        # tv_st = pd.to_datetime(tv_data['Transit_Start_UTC']).dt.to_pydatetime()
        # tv_sp = pd.to_datetime(tv_data['Transit_Stop_UTC']).dt.to_pydatetime()

    elif not exoplanet_tdf and t_name != 'Free Time' and not t_name.startswith(('WARNING')):
        v_data = pd.read_csv(
            os.path.join(aux_vis_path, st_name, f"Visibility for {t_name}.csv")
        )
        tmp_idx = a_list.index[(a_list['Star Name'] == t_name)].tolist()
        if len(tmp_idx) == 1:
            targ_info = pd.DataFrame(a_list.loc[tmp_idx[0]]).T
        else:
            targ_info = pd.DataFrame(a_list.loc[tmp_idx][a_list.loc[tmp_idx]['numPredefinedStarRois'] == 0].iloc[0]).T
        i_flag = 0
    elif t_name == 'Free Time':
        continue
    elif t_name.startswith(('WARNING')):
        print('Need visible STD during %s', t_name)
        logger.warning('Need visible STD during %s', t_name)
        continue
    else:
        logger.error("No visibility data for %s. Aborting schedule build.", t_name)
        break

    try:
        ra = targ_info['RA'].iloc[0]
        dec = targ_info['DEC'].iloc[0]
    except (KeyError, IndexError, AttributeError):
        try:
            star_sc = SkyCoord.from_name(st_name)
            ra = star_sc.ra.deg
            dec = star_sc.dec.deg
        except Exception as exc:  # SkyCoord lookup failed
            logger.error("Unable to resolve coordinates for %s: %s", st_name, exc)
            continue

    if "Time_UTC" in v_data.columns:
        v_time_all = pd.to_datetime(v_data["Time_UTC"]).dt.to_pydatetime()
    else:
        v_time_all = Time(v_data["Time(MJD_UTC)"], format="mjd", scale="utc").to_value("datetime")

    v_time = v_time_all[(v_time_all >= start) & (v_time_all <= stop)]
    v_time = np.vectorize(helper_codes.round_to_nearest_second)(v_time)
    v_flag = np.asarray(v_data['Visible'])[(v_time_all >= start) & (v_time_all <= stop)]

    v_flag, _ = helper_codes.remove_short_sequences(v_flag, too_short_sequences)

    #figure out where the visibility changes (gives final element where the visibility is the same)
    v_change = np.where(v_flag[:-1] != v_flag[1:])[0]
    
    st = start
    sp = v_time[-1]
    
    #if full visibility
    def full_visibility():
        tqdm.write(
            f'{st} to {sp}: No occultations needed; {t_name} is visible for the entire visit'
        )

        n = (sp - st) / dt

        sps = [st + (dt * (i + 1)) for i in range(int(n))]
        if int(n) < n:
            sps.append(sp)

        # if sps[-1] == v_time[-1]:
        #     sps[-1] = v_time[-2]

        sps_all = [st, *sps]
        for s in range(len(sps_all) - 1):
            if i_flag:
                pr = "2" if np.any((tv_st <= sps_all[s + 1]) * (tv_sp >= sps_all[s])) else "1"
            else:
                pr = "0"

            helper_codes.observation_sequence(
                visit,
                f'{("0" * (3 - len(str(s + 1)))) + str(s + 1)}',
                t_name,
                pr,
                sps_all[s],
                sps_all[s + 1],
                ra,
                dec,
                targ_info,
            )

    if len(v_change) == 0:
        full_visibility()
    
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

        break_occ_seq_longer_than_occultation_sequence_limit = True
        if break_occ_seq_longer_than_occultation_sequence_limit:
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
            oc_starts, oc_stops = start_tmp, stop_tmp


        def find_occultation_target(
            oc_starts,
            oc_stops,
            visit_start,
            visit_stop,
            primary_path,
            fallback_path,
            ra,
            dec,
            prefer_primary,
        ):
            """Loosely mimic prior logic: try one list, then fallback to another."""

            candidate_paths = (
                (primary_path, "target list"),
                (fallback_path, os.path.relpath(fallback_path, DATA_DIR)),
            )
            if not prefer_primary:
                candidate_paths = tuple(reversed(candidate_paths))

            for path, label in candidate_paths:
                info, flag = sch_occ_new(
                    oc_starts,
                    oc_stops,
                    visit_start,
                    visit_stop,
                    path,
                    sort_key="closest",
                    prev_obs=[ra, dec],
                )
                if flag:
                    tqdm.write(f"{visit_start} to {visit_stop}: Found occultation target from {label}")
                    return info, True

            tqdm.write(f"{visit_start} to {visit_stop}: No suitable occultation targets found")
            return None, False

        use_tar_list_for_occultations = False

        info, flag = find_occultation_target(
            oc_starts,
            oc_stops,
            st,
            sp,
            tar_path,
            occ_path,
            ra,
            dec,
            use_tar_list_for_occultations,
        )
        if not flag or info is None:
            logger.warning(
                "Unable to schedule occultation target for %s between %s and %s",
                t_name,
                st,
                sp,
            )
            continue

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
            if v_flag[v_change[v]]:  # Visible period
                current = st
                while current < sp:  # break observation sequences longer than 90 min
                    next_val = min(current + dt, sp)
                    priority = get_priority(i_flag, current, next_val)
                    try:
                        helper_codes.observation_sequence(
                            visit,
                            f'{("0" * (3 - len(str(seq_counter)))) + str(seq_counter)}',
                            t_name,
                            priority,
                            current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            next_val.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            ra,
                            dec,
                            targ_info,
                        )
                    except Exception as exc:  # helper downstream may fail
                        logger.error(
                            "Error adding visible sequence for %s (%s to %s): %s",
                            t_name,
                            current,
                            next_val,
                            exc,
                        )
                    seq_counter += 1
                    current = next_val

            else:  # Non-visible period (occultation)
                if break_occ_seq_longer_than_occultation_sequence_limit:
                    current = st
                    while current < sp:
                        next_val = min(current + occultation_sequence_limit, sp)
                        try:
                            if oc_tr >= len(info):
                                logger.warning(
                                    "Ran out of occultation targets for %s between %s and %s",
                                    t_name,
                                    current,
                                    next_val,
                                )
                                break
                            target_name = info['Target'].iloc[oc_tr]
                            occ_targ_info = t_list.loc[(t_list['Star Name'] == target_name)]
                            if occ_targ_info.empty:
                                occ_targ_info = a_list.loc[(a_list['Star Name'] == target_name)]
                            ra_occ = info['RA'].iloc[oc_tr]
                            dec_occ = info['DEC'].iloc[oc_tr]
                            helper_codes.observation_sequence(
                                visit,
                                f'{("0" * (3 - len(str(seq_counter)))) + str(seq_counter)}',
                                target_name,
                                '0',
                                current.strftime("%Y-%m-%dT%H:%M:%SZ"),
                                next_val.strftime("%Y-%m-%dT%H:%M:%SZ"),
                                ra_occ,
                                dec_occ,
                                occ_targ_info,
                            )
                            oc_tr += 1
                            seq_counter += 1
                        except Exception as exc:
                            target_label = (
                                info['Target'].iloc[oc_tr]
                                if oc_tr < len(info)
                                else "unknown"
                            )
                            logger.error(
                                "Error adding occultation sequence for %s (%s to %s): %s",
                                target_label,
                                current,
                                next_val,
                                exc,
                            )
                        current = next_val
                else:
                    try:
                        if oc_tr >= len(info):
                            logger.warning(
                                "Ran out of occultation targets for %s between %s and %s",
                                t_name,
                                st,
                                sp,
                            )
                            break
                        target_name = info['Target'].iloc[oc_tr]
                        occ_targ_info = t_list.loc[(t_list['Star Name'] == target_name)]
                        if occ_targ_info.empty:
                            occ_targ_info = a_list.loc[(a_list['Star Name'] == target_name)]
                        ra_occ = info['RA'].iloc[oc_tr]
                        dec_occ = info['DEC'].iloc[oc_tr]
                        helper_codes.observation_sequence(
                            visit,
                            f'{("0" * (3 - len(str(seq_counter)))) + str(seq_counter)}',
                            target_name,
                            '0',
                            st.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            sp.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            ra_occ,
                            dec_occ,
                            occ_targ_info,
                        )
                        oc_tr += 1
                        seq_counter += 1
                    except Exception as exc:
                        target_label = (
                            info['Target'].iloc[oc_tr]
                            if oc_tr < len(info)
                            else "unknown"
                        )
                        logger.error(
                            "Error adding occultation sequence for %s (%s to %s): %s",
                            target_label,
                            st,
                            sp,
                            exc,
                        )
def convert_to_string(element):
    for el in element.iter():
        for k, v in el.attrib.items():
            if isinstance(v, Number):
                el.set(k, str(v))
        if el.text and isinstance(el.text, Number):
            el.text = str(el.text)


convert_to_string(cal)
etstr = ET.tostring(cal, xml_declaration=True)
dom = minidom.parseString(etstr)
pretty_xml_as_string = dom.toprettyxml()

output_path = os.path.join(DATA_DIR, "Pandora_science_calendar.xml")
with open(output_path, "w+", encoding="utf-8") as handle:
    handle.write(pretty_xml_as_string)