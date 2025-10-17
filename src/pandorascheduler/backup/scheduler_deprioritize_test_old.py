import os
import numpy as np
import pandas as pd
import pickle
from astropy.time import Time
from astropy.coordinates import SkyCoord
from datetime import datetime, timedelta
import logging
import transits
from tqdm import tqdm
import helper_codes
import json
from multiprocessing import Pool
from functools import partial

# from . import PACKAGEDIR
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))


def Schedule(
    pandora_start: str,
    pandora_stop: str,
    target_list: str,
    obs_window: timedelta,
    transit_coverage_min: float,
    sched_wts: list,
    aux_key: str,#='random',
    aux_list:str,#="aux_list.csv",
    fname_tracker: str,
    commissioning_time: int = 30,
    sched_start: str = None,
    sched_stop: str = None,
    ):
    """Determine visibility for target(s) host star with Pandora given avoidance angles
    for Sun, Moon, and Earth limb.

    Parameters
    ----------
    pandora_start:          string
                            Date and time of start of Pandora science observing
                            ex. '2025-04-25 00:00:00'
    pandora_stop:           string
                            Date and time of end of Pandora science observing
                            ex. '2026-04-25 00:00:00'
    obs_window:             timedelta
                            Time for Pandora to stay on given target
    transit_coverage_min:   float
                            Minimum coverage required to capture a planetary transit
                            (scale 0-1)
    sched_wts:              list
                            weights given to transit coverage, saa overlap, and schedule factor (should sum to 1)
    commissioning_time:     int
                            Number of days for commissioning
                            ex. 30
    sched_start:            string
                            Date and time to start scheduling
                            ex. '2025-04-25 00:00:00'
    sched_stop:             string
                            Date and time to stop scheduling
                            ex. '2026-04-25 00:00:00'
    aux_key:                string or None
                            Key to determine if or how to schedule auxiliary observations
                            'closest' or None right now, add priority and others? later
    aux_list:               string
                            Path to the auxiliary target list. If there isn't one, specify None in aux_key
                            e.x. f"{PACKAGEDIR}/data/aux_list.csv" [default]

    Returns
    -------
    csv file
        file containing schedule for Pandora
    """
    assert sum(sched_wts) == 1., "Sum of weights should equal 1"

    # Convert times to datetime
    pandora_start = datetime.strptime(pandora_start, "%Y-%m-%d %H:%M:%S")
    # Add in commissioning time
    # pandora_start = pandora_start+timedelta(days=commissioning_time)
    pandora_stop = datetime.strptime(pandora_stop, "%Y-%m-%d %H:%M:%S")

    if sched_start == None:
        sched_start = pandora_start
    #VK BEGIN: the following elif crashes so I removed it
    #elif datetime.strptime(sched_start) < pandora_start:
    #    sched_start = pandora_start
    else:
        sched_start = datetime.strptime(sched_start, "%Y-%m-%d %H:%M:%S")# + timedelta(days=commissioning_time)
    if sched_stop == None:
        sched_stop = pandora_stop
    else:
        sched_stop = datetime.strptime(sched_stop, "%Y-%m-%d %H:%M:%S")

    # Import target list

    # if 'updated_targ_list' in globals():
    #     target_list = updated_targ_list.reset_index(drop=True)
    # else:
    #     target_list = pd.read_csv(f"{PACKAGEDIR}/data/" + target_list, sep=",")
    target_list = pd.read_csv(f"{PACKAGEDIR}/data/" + target_list, sep=",")
    #target_list = pd.read_csv(f"{PACKAGEDIR}/data/target_list.csv", sep=",")

    # Import no phase events
    if os.path.exists(f"{PACKAGEDIR}/data/no_phase_list.csv") == True:
        nophase_list = pd.read_csv(f"{PACKAGEDIR}/data/no_phase_list.csv", sep=",")
        nophase_targets = nophase_list["Target"]
        nophase_starts = [
            datetime.strptime(d, "%Y-%m-%d %H:%M:%S")
            for d in nophase_list["Obs Window Start"]
        ]
        nophase_stops = [
            datetime.strptime(d, "%Y-%m-%d %H:%M:%S")
            for d in nophase_list["Obs Window Stop"]
        ]

    ### Initialize schedule tracker
    planet_names = pd.DataFrame(
        np.array(target_list["Planet Name"]), columns=["Planet Name"]
    )
    transit_need = pd.DataFrame(
        np.array(target_list["Number of Transits to Capture"]), columns=["Transits Needed"]
    )
    transit_have = pd.DataFrame(
        np.zeros(len(planet_names)), columns=["Transits Acquired"]
    )
    tracker = pd.concat([planet_names, transit_need, transit_have], axis=1)

    ### Check if previous observations already exist and if so update tracker
    if os.path.exists(f"{PACKAGEDIR}/data/Pandora_archive.csv") == True:
        archive = pd.read_csv(f"{PACKAGEDIR}/data/Pandora_archive.csv", sep=",")
        ### ADD: Warning if conflict between scheduling start and previous observations
        if len(archive) > 0:
            for i in range(len(archive)):
                tracker.loc[
                    (tracker["Planet Name"] == archive["Target"][i]), "Transits Needed"
                ] = (
                    tracker.loc[(tracker["Planet Name"] == archive["Target"][i])][
                        "Transits Needed"
                    ]
                    - 1
                )
                tracker.loc[
                    (tracker["Planet Name"] == archive["Target"][i]),
                    "Transits Acquired",
                ] = (
                    tracker.loc[(tracker["Planet Name"] == archive["Target"][i])][
                        "Transits Acquired"
                    ]
                    + 1
                )

    ### Determine how many transits exist for each target within Pandora's lifetime
    ### and the specified observation window if different
    pandora_transits_left = []
    schedule_transits_left = []
    for i in range(len(target_list)):
        planet_name = target_list["Planet Name"][i]
        star_name = target_list["Star Name"][i]
        planet_data = pd.read_csv(
            f"{PACKAGEDIR}/data/targets/{star_name}/{planet_name}/Visibility for {planet_name}.csv"
        )
        planet_data = planet_data.drop(
            planet_data.index[(planet_data["Transit_Coverage"] < transit_coverage_min)]
        ).reset_index(drop=True)

        logging.info(planet_name, "has", len(planet_data), "transits with greater transit coverage than", transit_coverage_min,)

        start_transits = Time(
            planet_data["Transit_Start"], format="mjd", scale="utc").to_value("datetime")
        end_transits = Time(
            planet_data["Transit_Stop"], format="mjd", scale="utc").to_value("datetime")
        p_trans = planet_data.index[
            (pandora_start <= start_transits) & (end_transits <= pandora_stop)
        ]
        s_trans = planet_data.index[
            (sched_start <= start_transits) & (end_transits <= sched_stop)
        ]
        pandora_transits_left.append(len(p_trans))
        schedule_transits_left.append(len(s_trans))
        
    pandora_transits_left = pd.DataFrame(
        pandora_transits_left, columns=["Transits Left in Lifetime"]
    )
    schedule_transits_left = pd.DataFrame(
        schedule_transits_left, columns=["Transits Left in Schedule"]
    )
    tracker = pd.concat(
        [tracker, pandora_transits_left, schedule_transits_left], axis=1
    )

    ### Calculate an initial priority number for each planet
    trans_priority = []
    for i in range(len(tracker)):
        trans_priority.append(
            tracker["Transits Left in Lifetime"][i] - tracker["Transits Needed"][i]
        )

    trans_priority = pd.DataFrame(trans_priority, columns=["Transit Priority"])
    tracker = pd.concat([tracker, trans_priority], axis=1)

    ### Begin scheduling
    start = sched_start
    stop = start + obs_window
    sched_df = pd.DataFrame(
        [],
        columns=[
            "Target",
            "Observation Start",
            "Observation Stop",
            "Transit Coverage",
            "SAA Overlap",
            "Schedule Factor",
            "Transit Factor",
            "Quality Factor",
        ],
    )
    
    #code to check if losing 10% of days causes us to fail
    # cut_arr=np.random.rand(37)*10
    # cut_dt=[timedelta(days=np.round(cut_arr[0]))]
    # for i in range(len(cut_arr)-1):
    #     cut_dt.append(timedelta(days=np.round(np.sum(cut_arr[:i+1]))))
    # dates_to_cut=[pandora_start+t for t in cut_dt]
    # cut_i=0

    non_primary_obs_time = {}
    deprioritized_targets = set()
    all_target_obs_time = {}
    
    while stop <= sched_stop:

        #code to check if losing 10% of days causes us to fail
        # if np.abs((start - dates_to_cut[cut_i]).total_seconds()) < 86400:
        #     #skip the day, update start and stop
        #     start = stop
        #     stop = start + obs_window
        #     cut_i=np.mod(cut_i+1, 37)
        #     continue
        
        tracker = tracker.sort_values(by=["Transit Priority"]).reset_index(drop=True)
        obs_rng = pd.date_range(start, stop, freq="min")
        temp_df = pd.DataFrame(
            [],
            columns=[
                "Planet Name",
                "Obs Start",
                "Obs Gap Time",
                "Transit Coverage",
                "SAA Overlap",
                "Schedule Factor",
                "Transit Factor",
                "Quality Factor",
            ],
        )

        ### First check if a no phase event is within observing window
        overlap_nophase = obs_rng.intersection(nophase_starts)
        if len(overlap_nophase) > 0:
            print('here nophase')
            obs_start = nophase_starts[nophase_starts.index(overlap_nophase[0])]
            obs_stop = nophase_stops[nophase_starts.index(overlap_nophase[0])]
            target = nophase_targets[nophase_starts.index(overlap_nophase[0])]

            if obs_rng[0] < obs_start:
                free = [["Free Time", obs_rng[0], obs_start]]
                free = pd.DataFrame(
                    free, columns=["Target", "Observation Start", "Observation Stop"]
                )
                sched_df = pd.concat([sched_df, free], axis=0)

            sched = [[target, obs_start, obs_stop]]
            sched = pd.DataFrame(
                sched, columns=["Target", "Observation Start", "Observation Stop"]
            )
            sched_df = pd.concat([sched_df, sched], axis=0)

            logging.info("Scheduled no phase event", target)
            start = obs_stop
            stop = start + obs_window
            
            #Check MTR for each planet
            for i in range(len(tracker)):
                tf=tracker["Transit Factor"][i]
                if tf == 1.:
                    #Force compliance with minimum targeting requirements if the 
                    #final transit occurs within the observation window
                    planet_name = tracker["Planet Name"][i]
                    planet_data = pd.read_csv(
                    f"{PACKAGEDIR}/data/targets/{star_name}/{planet_name}/Visibility for {planet_name}.csv"
                    )
                    
                    start_transit = planet_data["Transit_Start"][-1]
                    end_transit = planet_data["Transit_Stop"][-1]    
                    
                    start_transit = start_transit - timedelta(
                        seconds=start_transit.second,
                        microseconds=start_transit.microsecond,
                    )
                    end_transit = end_transit - timedelta(
                        seconds=end_transit.second,
                        microseconds=end_transit.microsecond,
                    )

                    early_start = end_transit - timedelta(
                        hours=20
                    )  # Earliest start time to capture transit plus >=4 hours post transit
                    late_start = start_transit - timedelta(
                        hours=4
                    )  # Latest start time to capture transit plus >=4 hours pre transit
    
                    # Check if any transit occurs during observing window
                    
                    
                    start_rng = pd.date_range(early_start, late_start, freq="min")
                    overlap_times = obs_rng.intersection(start_rng)
                    
                    if len(overlap_times) > 0:

                        # Calc a 'transit factor'
                        t_left = tracker.loc[
                            (tracker["Planet Name"] == planet_name),
                            "Transits Left in Lifetime",
                        ].iloc[0]
                        t_need = tracker.loc[
                            (tracker["Planet Name"] == planet_name), "Transits Needed"
                        ].iloc[0]
                        t_factor = t_left / t_need

                        # Calc scheduling efficiency factor
                        obs_start = overlap_times[0]
                        gap_time = obs_start - obs_rng[0]
                        s_factor = 1 - (gap_time / obs_window)  # maximize

                        # Calc a quality factor (currently based on transit coverage, SAA crossing, scheduling efficiency)
                        trans_cover = planet_data["Transit_Coverage"][-1]  # maximize
                        saa_cover = planet_data["SAA_Overlap"][-1]
                        q_factor = (
                            (sched_wts[0] * trans_cover)
                            + (sched_wts[1] * (1 - saa_cover))
                            + (sched_wts[2] * s_factor)
                        )
                        
                        # Schedule observation with warning
            
                        if obs_rng[0] < obs_start:
                            free = [["Free Time", obs_rng[0], obs_start]]
                            free = pd.DataFrame(
                                free, columns=["Target", "Observation Start", "Observation Stop"]
                            )
                            sched_df = pd.concat([sched_df, free], axis=0)
            
                        sched = [
                            [
                                planet_name,
                                obs_start,
                                obs_stop,
                                trans_cover,
                                saa_cover,
                                s_factor,
                                q_factor,
                            ]
                        ]
                        sched = pd.DataFrame(
                            sched,
                            columns=[
                                "Target",
                                "Observation Start",
                                "Observation Stop",
                                "Transit Coverage",
                                "SAA Overlap",
                                "Schedule Factor",
                                "Quality Factor",
                            ],
                        )
                        sched_df = pd.concat([sched_df, sched], axis=0)
            
                        # update tracker info
                        tracker.loc[(tracker["Planet Name"] == planet_name), "Transits Needed"] = (
                            tracker.loc[(tracker["Planet Name"] == planet_name)]["Transits Needed"]
                            - 1
                        )
                        tracker.loc[
                            (tracker["Planet Name"] == planet_name), "Transits Acquired"
                        ] = (
                            tracker.loc[(tracker["Planet Name"] == planet_name)][
                                "Transits Acquired"
                            ]
                            + 1
                        )
                        tracker.loc[(tracker["Planet Name"] == planet_name), "Transit Priority"] = (
                            tracker.loc[
                                (tracker["Planet Name"] == planet_name), "Transits Left in Lifetime"
                            ]
                            - tracker.loc[
                                (tracker["Planet Name"] == planet_name), "Transits Needed"
                            ]
                        )
                        # logging.warning(planet_name, "{planet_name} Observation Forced Over No Phase Event")
                        logging.info(
                            "Scheduled the ",
                            tracker.loc[
                                (tracker["Planet Name"] == planet_name), "Transits Acquired"
                            ].iloc[0],
                            " transit of ",
                            planet_name,
                            "transit coverage: ",
                            trans_cover,
                        )
                        start = obs_stop
                        stop = start + obs_window
                        continue
                    
                elif tf <= 2.:
                    #Flag that the MTRM is getting low for the planet, but do
                    #not force compliance
                    planet_name = tracker["Planet Name"][i]
                    if obs_rng[0] < obs_start:
                        free = [["Free Time", obs_rng[0], obs_start]]
                        free = pd.DataFrame(
                            free, columns=["Target", "Observation Start", "Observation Stop"]
                        )
                        sched_df = pd.concat([sched_df, free], axis=0)
        
                    sched = [[target, obs_start, obs_stop]]
                    sched = pd.DataFrame(
                        sched, columns=["Target", "Observation Start", "Observation Stop"]
                    )
                    sched_df = pd.concat([sched_df, sched], axis=0)
        
                    logging.info("Scheduled no phase event", target)
                    # logging.warning(target, "Warning: {planet_name} has MTRM < 2 and is transiting during the event")
                    start = obs_stop
                    stop = start + obs_window
                    #recalculate MTRM?
                    continue
                
            
            
            continue
        
        
            #******ADD: functionality to force MTR prioritization and flag if nophase event causes the minimum schedule to not be met
       

        ### Next look at each planet and determine if transit occurs in window
        for i in range(len(tracker)):
            planet_name = tracker["Planet Name"][i]

            if (
                tracker.loc[(tracker["Planet Name"] == planet_name), "Transits Needed"][
                    i
                ]
                == 0
            ):
                pass

            else:
                star_name = target_list["Star Name"][
                    np.where(target_list["Planet Name"] == planet_name)[0][0]
                ]
                planet_data = pd.read_csv(
                    f"{PACKAGEDIR}/data/targets/{star_name}/{planet_name}/Visibility for {planet_name}.csv"
                )
                planet_data = planet_data.drop(
                    planet_data.index[
                        (planet_data["Transit_Coverage"] < transit_coverage_min)
                    ]
                ).reset_index(drop=True)
                planet_data["Transit_Start"] = Time(
                    planet_data["Transit_Start"], format="mjd", scale="utc"
                ).to_value("datetime")
                planet_data["Transit_Stop"] = Time(
                    planet_data["Transit_Stop"], format="mjd", scale="utc"
                ).to_value("datetime")
                planet_data = planet_data.drop(
                    planet_data.index[(planet_data["Transit_Start"] < start)]
                ).reset_index(drop=True)
                start_transits = planet_data["Transit_Start"].copy()
                end_transits = planet_data["Transit_Stop"].copy()
                # start_transits = Time(planet_data["Transit_Start"], format="mjd", scale="utc").to_value("datetime")
                # end_transits = Time(planet_data["Transit_Stop"], format="mjd", scale="utc").to_value("datetime")

                p_trans = planet_data.index[
                    (pandora_start <= start_transits) & (end_transits <= pandora_stop)
                ]
                s_trans = planet_data.index[
                    (sched_start <= start_transits) & (end_transits <= sched_stop)
                ]
                tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transits Left in Lifetime"
                ] = len(p_trans)
                tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transits Left in Schedule"
                ] = len(s_trans)
                tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transit Priority"
                ] = (
                    tracker.loc[
                        (tracker["Planet Name"] == planet_name),
                        "Transits Left in Lifetime",
                    ]
                    - tracker.loc[
                        (tracker["Planet Name"] == planet_name), "Transits Needed"
                    ]
                )

                # Remove seconds and below from times
                for j in range(len(start_transits)):
                    start_transits.iloc[j] = start_transits.iloc[j] - timedelta(
                        seconds=start_transits.iloc[j].second,
                        microseconds=start_transits.iloc[j].microsecond,
                    )
                    end_transits.iloc[j] = end_transits.iloc[j] - timedelta(
                        seconds=end_transits.iloc[j].second,
                        microseconds=end_transits.iloc[j].microsecond,
                    )

                early_start = end_transits - timedelta(
                    hours=20
                )  # Earliest start time to capture transit plus >=4 hours post transit
                late_start = start_transits - timedelta(
                    hours=4
                )  # Latest start time to capture transit plus >=4 hours pre transit

                # Check if any transit occurs during observing window
                for j in range(len(early_start)):
                    start_rng = pd.date_range(early_start[j], late_start[j], freq="min")
                    overlap_times = obs_rng.intersection(start_rng)
                    if len(overlap_times) > 0:
                        # Calc a 'transit factor'
                        t_left = tracker.loc[
                            (tracker["Planet Name"] == planet_name),
                            "Transits Left in Lifetime",
                        ].iloc[0]
                        t_need = tracker.loc[
                            (tracker["Planet Name"] == planet_name), "Transits Needed"
                        ].iloc[0]
                        t_factor = t_left / t_need

                        # Calc scheduling efficiency factor
                        obs_start = overlap_times[0]
                        gap_time = obs_start - obs_rng[0] # THIS IS THE time between the start of the current observation window and the start of the potential observation.
                        s_factor = 1 - (gap_time / obs_window)  # maximize

                        # Calc a quality factor (currently based on transit coverage, SAA crossing, scheduling efficiency)
                        trans_cover = planet_data["Transit_Coverage"][j]  # maximize
                        saa_cover = planet_data["SAA_Overlap"][j]
                        q_factor = (
                            (sched_wts[0] * trans_cover)
                            + (sched_wts[1] * (1 - saa_cover))
                            + (sched_wts[2] * s_factor)
                        )

                        temp = [
                            [
                                planet_name,
                                obs_start,
                                gap_time,
                                planet_data["Transit_Coverage"][j],
                                saa_cover,
                                s_factor,
                                t_factor,
                                q_factor,
                            ]
                        ]
                        temp = pd.DataFrame(
                            temp,
                            columns=[
                                "Planet Name",
                                "Obs Start",
                                "Obs Gap Time",
                                "Transit Coverage",
                                "SAA Overlap",
                                "Schedule Factor",
                                "Transit Factor",
                                "Quality Factor",
                            ],
                        )
                        temp_df = pd.concat([temp_df, temp], axis=0)

        ### Check if there's no transits occuring during the observing window
        ### Schedule auxiliary observation if possible
        if len(temp_df) == 0:
            # print(f'NO TRANSITS IN {start}, {stop}, WHAT TO USE FOR STAR NAME?!')

            if sched_df['Target'].iloc[-1].endswith(('b', 'c', 'd')):
                star_name_tmp = sched_df['Target'].iloc[-1][0:-2]
            else:
                star_name_tmp = sched_df['Target'].iloc[-1]
                
            star_sc = SkyCoord.from_name(star_name_tmp)#star_name)
            ra = star_sc.ra.deg
            dec = star_sc.dec.deg

            aux_df, log_info, non_primary_obs_time, deprioritized_targets = Schedule_aux(start, stop, aux_key, aux_list, \
                prev_obs=[ra, dec], non_primary_obs_time=non_primary_obs_time, deprioritized_targets=deprioritized_targets); print()
            sched_df = pd.concat([sched_df, aux_df], axis=0)

            all_target_obs_time[planet_name] = all_target_obs_time.get(planet_name, timedelta()) + (obs_stop - obs_start)
            # Update observation time for auxiliary targets
            if len(aux_df) > 0:
                for _, row in aux_df.iterrows():
                    if row['Target'] != "Free Time":
                        aux_target = row['Target']
                        aux_duration = row['Observation Stop'] - row['Observation Start']
                        all_target_obs_time[aux_target] = all_target_obs_time.get(aux_target, timedelta()) + aux_duration

            logging.info(log_info, start, stop)
            start = stop
            stop = start + obs_window
            continue                

        else:
            # VK BEGIN
            # COMMENT OUT THE FOLLOWING, REPLACE WITH CODE BELOW
            # # First, check Transit Factor for planets running out of available transits
            # # if (temp_df["Transit Factor"] < 1.5).any():
            # if (temp_df["Transit Factor"] < 2).any():
            #     temp_df = temp_df.sort_values(by=["Transit Factor"]).reset_index(
            #         drop=True
            #     )
            #     # logging.warning(temp_df["Planet Name"][0], "Transit Factor Warning", temp_df["Transit Factor"][0])

            # # Otherwise sort by Quality Factor and schedule best target
            # else:
            #     temp_df = temp_df.sort_values(
            #         by=["Quality Factor"], ascending=False
            #     ).reset_index(drop=True)
            #######
            # REPLACE WITH THE FOLLOWING SORTING ORDER:
            # 1) Transit Factor; 2) Quality Factor first and then Transit Factor
            #######
            # First, check Transit Factor for planets running out of available transits
            # if (temp_df["Transit Factor"] < 1.5).any():
            if (temp_df["Transit Factor"] < 2).any():
                temp_df = temp_df.sort_values(by=["Transit Factor"]).reset_index(
                    drop=True
                )
                # logging.warning(temp_df["Planet Name"][0], "Transit Factor Warning", temp_df["Transit Factor"][0])
            else:
                # Otherwise, sort by "Quality Factor" (descending) and then by "Transit Factor" (ascending)
                temp_df = temp_df.sort_values(
                    by=["Quality Factor", "Transit Factor"],
                    ascending=[False, True]
                ).reset_index(drop=True)

                # # 3) if more than one planet has Quality Factor = 1, sort by both Quality Factor and Transit Factor.
                # # Check if there are multiple planets with Quality Factor = 1
                # if (temp_df["Quality Factor"] == 1).sum() > 1:
                #     # If yes, re-sort those planets by both "Quality Factor" and "Transit Factor"
                #     mask = temp_df["Quality Factor"] == 1
                #     temp_df.loc[mask] = temp_df.loc[mask].sort_values(
                #         by=["Quality Factor", "Transit Factor"],
                #         ascending=[False, True]
                #     )

            # Reset the index after all sorting operations
            # temp_df = temp_df.reset_index(drop=True)
            #######
            #######
            # VK END

            planet_name = temp_df["Planet Name"][0]
            star_name = target_list["Star Name"][np.where(target_list["Planet Name"] == planet_name)[0][0]]
            obs_start = temp_df["Obs Start"][0]
            obs_stop = obs_start + timedelta(hours=24)
            trans_cover = temp_df["Transit Coverage"][0]
            saa_cover = temp_df["SAA Overlap"][0]
            s_factor = temp_df["Schedule Factor"][0]
            q_factor = temp_df["Quality Factor"][0]

            # VK BEGIN
            print_every_day = False#True
            if print_every_day:
                num_sig_digits_ = 2
                if obs_start == pandora_start:
                    print(temp_df.head(1).round(num_sig_digits_).to_string(index=False))
                else:
                    df_string = temp_df.head(1).round(num_sig_digits_).to_string(index=False)# Get the string representation of the DataFrame with headers
                    lines = df_string.split('\n')# Split the string into lines
                    if len(lines) > 1:# Print only the data line (second line, index 1)
                        print(lines[1])
            # VK END

            if obs_rng[0] < obs_start:
                star_sc = SkyCoord.from_name(star_name)
                ra = star_sc.ra.deg
                dec = star_sc.dec.deg

                aux_df, log_info, non_primary_obs_time, deprioritized_targets = Schedule_aux(start, obs_start, aux_key, aux_list, \
                    prev_obs=[ra, dec], non_primary_obs_time=non_primary_obs_time, deprioritized_targets=deprioritized_targets); print()
                sched_df = pd.concat([sched_df, aux_df], axis=0)

                all_target_obs_time[planet_name] = all_target_obs_time.get(planet_name, timedelta()) + (obs_stop - obs_start)
                # Update observation time for auxiliary targets
                if len(aux_df) > 0:
                    for _, row in aux_df.iterrows():
                        if row['Target'] != "Free Time":
                            aux_target = row['Target']
                            aux_duration = row['Observation Stop'] - row['Observation Start']
                            all_target_obs_time[aux_target] = all_target_obs_time.get(aux_target, timedelta()) + aux_duration


                logging.info(log_info, start, obs_start)

            sched = [
                [
                    planet_name,
                    obs_start,
                    obs_stop,
                    trans_cover,
                    saa_cover,
                    s_factor,
                    q_factor,
                ]
            ]
            sched = pd.DataFrame(
                sched,
                columns=[
                    "Target",
                    "Observation Start",
                    "Observation Stop",
                    "Transit Coverage",
                    "SAA Overlap",
                    "Schedule Factor",
                    "Quality Factor",
                ],
            )
            sched_df = pd.concat([sched_df, sched], axis=0)

            # update tracker info
            tracker.loc[(tracker["Planet Name"] == planet_name), "Transits Needed"] = (
                tracker.loc[(tracker["Planet Name"] == planet_name)]["Transits Needed"]
                - 1
            )
            tracker.loc[
                (tracker["Planet Name"] == planet_name), "Transits Acquired"
            ] = (
                tracker.loc[(tracker["Planet Name"] == planet_name)][
                    "Transits Acquired"
                ]
                + 1
            )
            tracker.loc[(tracker["Planet Name"] == planet_name), "Transit Priority"] = (
                tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transits Left in Lifetime"
                ]
                - tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transits Needed"
                ]
            )

            logging.info(
                "Scheduled the ",
                tracker.loc[
                    (tracker["Planet Name"] == planet_name), "Transits Acquired"
                ].iloc[0],
                " transit of ",
                planet_name,
                "transit coverage: ",
                trans_cover,
            )
            start = obs_stop
            stop = start + obs_window
            continue

    ### Save results
    sched_df = sched_df.sort_values(by=["Observation Start"]).reset_index(drop=True)
    save_fname = f"Pandora_Schedule_{sched_wts[0]}_{sched_wts[1]}_{sched_wts[2]}_{pandora_start}.csv"
    sched_df.to_csv((f"{PACKAGEDIR}/data/" + save_fname), sep=",", index=False)

    # Save tracker
    # save_tracker_fname = f"{PACKAGEDIR}/data/schedule_tracker.pkl"
    with open(fname_tracker, 'wb') as file:
        pickle.dump(tracker, file)

            # Save non_primary_obs_time
    non_primary_obs_time_fname = f"Non_Primary_Obs_Time_{pandora_start}.json"
    with open(f"{PACKAGEDIR}/data/" + non_primary_obs_time_fname, 'w') as file:
        json.dump({k: str(v) for k, v in non_primary_obs_time.items()}, file)

    # Save deprioritized_targets
    deprioritized_targets_fname = f"Deprioritized_Targets_{pandora_start}.json"
    with open(f"{PACKAGEDIR}/data/" + deprioritized_targets_fname, 'w') as file:
        json.dump(list(deprioritized_targets), file)

    # Call this function at the end of the Schedule function
    helper_codes.save_observation_time_report(all_target_obs_time, target_list, f"{PACKAGEDIR}/data/Observation_Time_Report_{pandora_start}.csv")

    return tracker

def Schedule_aux(start, stop, aux_key, aux_list, prev_obs, \
    non_primary_obs_time,deprioritized_targets, **kwargs):

    obs_rng = pd.date_range(start, stop, freq="min")

    if aux_key == None:
        # No aux target scheduling, free time
        free = [["Free Time", start, stop]]
        aux_df = pd.DataFrame(
            free, columns=["Target", "Observation Start", "Observation Stop"]
        )
        log_info = 'Free time, no observation scheduled.'
        return aux_df, log_info, non_primary_obs_time, deprioritized_targets
        
    for tdf_idx in range(1,4):
        aux_fn_tmp = f"{PACKAGEDIR}/data/{target_definition_files[tdf_idx]}_targets.csv"
        aux_targs = pd.read_csv(aux_fn_tmp)
        mask = ~aux_targs['Star Name'].isin(deprioritized_targets)
        aux_targs = aux_targs[mask].reset_index(drop=True)
        names = aux_targs['Star Name']
        ras = aux_targs['RA']
        decs = aux_targs['DEC']

        if aux_key == 'closest':
            po_sc = SkyCoord(unit='deg', ra=prev_obs[0], dec=prev_obs[1])
            aux_sc = SkyCoord(unit='deg', ra=ras, dec=decs)
            dif = aux_sc.separation(po_sc).deg
            aux_sc_sorted_by_distance = dif.argsort()
            names = names.iloc[aux_sc_sorted_by_distance].reset_index(drop=True)
            ras = ras.iloc[aux_sc_sorted_by_distance].reset_index(drop=True)
            decs = decs.iloc[aux_sc_sorted_by_distance].reset_index(drop=True)

        vis_all_targs = []
        vis_any_targs = []
        targ_vis = []

        for n in tqdm(range(len(names)), desc=f"Finding visible non-primary target for {str(start)} to {str(stop)} from '{target_definition_files[tdf_idx]}'"):

            current_obs_time = non_primary_obs_time.get(names[n], timedelta())
            if current_obs_time >= timedelta(hours=12):
                print(names[n])
                continue

            try:
                vis_file = f"{PACKAGEDIR}/data/aux_targets/{names[n]}/Visibility for {names[n]}.csv"
                vis = pd.read_csv(vis_file, usecols=["Time(MJD_UTC)", "Visible"])
                time_mask = (Time(vis["Time(MJD_UTC)"], format='mjd', scale='utc') >= start) & \
                            (Time(vis["Time(MJD_UTC)"], format='mjd', scale='utc') <= stop)
                vis_filtered = vis[time_mask]

                if not vis_filtered.empty and vis_filtered['Visible'].all():
                    vis_all_targs.append(n)
                    break
                elif not vis_filtered.empty and vis_filtered['Visible'].any():
                    vis_any_targs.append(n)
                    targ_vis.append(100*(np.sum(vis_filtered['Visible'])/len(vis_filtered)))
            except FileNotFoundError:
                pass

            # vis_all_targs, vis_any_targs, targ_vis = helper_codes.find_visible_targets(names, start, stop)

        # Select target
        if len(vis_all_targs) > 0:
            idx = 0 if aux_key == 'closest' else np.random.randint(0,len(vis_all_targs))
            name = names[vis_all_targs[idx]]
            print(f"{name} scheduled for observations with full visibility from {target_definition_files[tdf_idx]}.")
            # print()
            log_info = f"{name} scheduled for observations with full visibility from {target_definition_files[tdf_idx]}."
        # elif len(vis_any_targs) > 0:
        #     if aux_key == 'closest':
        #         idx = 0
        #     elif aux_key == 'max_visibility_any':
        #         idx = np.asarray(targ_vis).argmax()
        #     else:
        #         idx = np.random.randint(0, len(vis_any_targs))
        #     name = names[vis_any_targs[idx]]
        #     vis_perc = targ_vis[idx]
        #     print(f"{name} scheduled for observations with {vis_perc:.2f}% visibility from {target_definition_files[tdf_idx]}.")
        #     # print()
        #     log_info = f"{name} scheduled for observations with {vis_perc:.2f}% visibility from {target_definition_files[tdf_idx]}."
        else:
            print(f"No fully visible targets in {target_definition_files[tdf_idx]}, try next...")
            # print()
            log_info = f"No fully visible targets in {target_definition_files[tdf_idx]}, try next..."
            continue  # No visible targets in this file, move to next

        # Update observation time and check for deprioritization
        obs_duration = stop - start
        new_obs_time = non_primary_obs_time.get(name, timedelta()) + obs_duration

        if new_obs_time > timedelta(hours=12):
            deprioritized_targets.add(name)
            print(f"Removing {name} due to deprioritization")
            # print()
            # continue  # Skip this target and try the next file
        
        # If we've reached here, we have a valid target
        non_primary_obs_time[name] = new_obs_time
        
        sched = [[name, start, stop]]
        aux_df = pd.DataFrame(sched, columns=["Target", "Observation Start", "Observation Stop"])

        return aux_df, log_info, non_primary_obs_time, deprioritized_targets

    # If we've gone through all files and found no suitable target
    free = [["Free Time", start, stop]]
    aux_df = pd.DataFrame(free, columns=["Target", "Observation Start", "Observation Stop"])
    log_info = 'Free time, no observation scheduled.'

    return aux_df, log_info, non_primary_obs_time, deprioritized_targets

def Schedule_all_scratch(
        blocks:list,
        pandora_start:str,
        pandora_stop:str,
        target_definition_files:list,
        # target_list:str,
        # target_partner_list:str,
        obs_window: timedelta,
        transit_coverage_min: float,
        sched_wts: list,
        aux_key: str,
        # aux_list: str,
        fname_tracker: str,
        commissioning_time: int,
        sched_start=None,
        sched_stop=None,
        ):
    
    """ Make a full run of the scheduler from scratch including all targets.
    This should only be needed for the first run of the scheduler for a given
    target list or parter target list, though the function checks if individual
    steps have been completed in situ. If these are updated, the visibility
    and SAA anomoly crossing will need to be updated to cover any changes.
    This also requires a GMAT file for the Pandora mission, available at:
    https://github.com/PandoraMission/pandora-scheduler/blob/main/src/pandorascheduler/data/GMAT_Pandora.txt
    stored with git LFS
    
    Parameters
    ----------
    blocks:         list
                    Avoidance angle for the Sun, Moon, and Earth limb in order
    pandora_start:              string
                            Date and time of start of Pandora science observing 
                            ex. '2025-04-25 00:00:00'
    pandora_stop:               string
                            Date and time of end of Pandora science observing 
                            ex. '2026-04-25 00:00:00'
    target_list:                string
                            Name of csv file with list of targets and their parameters
                            
    obs_window:                 timedelta
                            Time for Pandora to stay on given target
    transit_coverage_min:       float
                            Minimum coverage required to capture a planetary transit
                            (scale 0-1)
    sched_wts:                  list
                            weights given to transit coverage, saa overlap, and schedule factor (should sum to 1)
    commissioning_time:     int
                            Number of days for commissioning
                            ex. 30
    sched_start:                string
                            Date and time to start scheduling, defaults to pandora_start
                            ex. '2025-04-25 00:00:00'
    sched_stop:                 string
                            Date and time to stop scheduling, defaults to pandora_stop
                            ex. '2026-04-25 00:00:00'
                    
    Returns
    -------
    csv files
        files containing targets' host star visibility by Pandora (one per star)
    csv files
        files containing targets' transits during Pandora's lifetime (one per target and partner)
        this includes transit overlaps, if applicaple, and SAA anomaly overlap
    csv file
        file containing a schedule for Pandora observation
    """
    
    #Get planet and star names from target_list and target_partner_list for iterables
    # t_csv = pd.read_csv(f'{PACKAGEDIR}/data/' + target_list, sep=',')


    #Create visibity information for target host stars given avoidance angles 
    #for Sun, Moon, and Earth limb during Pandora's science observation lifetime.

    tdf = np.asarray(target_definition_files)
    for ii1 in range(3):
        fn_ = tdf[tdf == tdf[ii1]][0]
        tmp_csv = pd.read_csv(f'{PACKAGEDIR}/data/{fn_}_targets.csv', sep=',')
        tmp_planet_name_lst = tmp_csv['Planet Name']
        tmp_star_name_lst = tmp_csv['Star Name']

        vis_done=[]
        for s in range(len(tmp_star_name_lst)):
            vis_done.append(os.path.exists(f'{PACKAGEDIR}/data/targets/{tmp_star_name_lst[s]}/Visibility for {tmp_star_name_lst[s]}.csv'))

        if np.all(vis_done) == False:
            transits.star_vis(blocks[0], blocks[1], blocks[2], pandora_start, pandora_stop, gmat_file, obs_name, \
                save_pth = f'{PACKAGEDIR}/data/targets/', targ_list = f'{PACKAGEDIR}/data/{fn_}_targets.csv')

        # #Determine primary transits for targets and partners during Pandora's 
        # #science observation lifetime.
        if fn_ == 'primary-exoplanet':
            for s in range(len(tmp_star_name_lst)):
                if not os.path.exists(f'{PACKAGEDIR}/data/targets/{tmp_star_name_lst[s]}/{tmp_planet_name_lst[s]}'):
                    transits.transit_timing(f'{fn_}_targets.csv', tmp_planet_name_lst[s], tmp_star_name_lst[s])

    # vis_done=[]
    # for s in range(len(prim_star_name_list)):
    #     vis_done.append(os.path.exists(f'{PACKAGEDIR}/data/targets/{prim_star_name_list[s]}/Visibility for {prim_star_name_list[s]}.csv'))
    # if np.all(vis_done) == False:
    #     transits.star_vis(blocks[0], blocks[1], blocks[2], pandora_start, pandora_stop) 

    # #Determine primary transits for targets and partners during Pandora's 
    # #science observation lifetime.
    # for t in range(len(prim_planet_name_list)):
    #     if not os.path.exists(f'{PACKAGEDIR}/data/targets/{prim_star_name_list[t]}/{prim_planet_name_list[t]}'):
    #         # transits.transit_timing(target_list, t_list[t], ts_list[t])
    #         transits.transit_timing(f'{target_definition_files[0]}_targets.csv', prim_planet_name_list[t], prim_star_name_list[t])

    # for t in range(len(tp_list)):
    #     if not os.path.exists(f'{PACKAGEDIR}/data/targets/{ts_list[t]}/{tp_list[t]}'):
    #         transits.transit_timing(target_partner_list, tp_list[t], ps_list[t])
    
    t_csv=pd.read_csv(f'{PACKAGEDIR}/data/primary-exoplanet_targets.csv', sep=',')
    tp_csv=pd.read_csv(f'{PACKAGEDIR}/data/auxiliary-exoplanet_targets.csv', sep=',')
    t_list=t_csv['Planet Name']
    tp_list=tp_csv['Planet Name']
    ts_list=t_csv['Star Name']
    ps_list=tp_csv['Star Name']
    for t in tqdm(range(len(t_list))):
        #Determine if there is overlap between target planets' transits and any 
        #companion planets
        try:
            vis=pd.read_csv(f'{PACKAGEDIR}/data/targets/{ts_list[t]}/{t_list[t]}/Visibility for {t_list[t]}.csv')
            t_over=vis['Transit_Overlap']
        except KeyError:
            transits.Transit_overlap('primary-exoplanet_targets.csv','auxiliary-exoplanet_targets.csv',ts_list[t])
        
        #Determine if there is overlap between target planets' transits
        #and South Atlantic Anomaly crossing
        try:
            vis=pd.read_csv(f'{PACKAGEDIR}/data/targets/{ts_list[t]}/{t_list[t]}/Visibility for {t_list[t]}.csv')
            saa=vis['SAA_Overlap']
        except KeyError:    
            transits.SAA_overlap(t_list[t], ts_list[t])
    
    
    #Schedule observations for the scheduling period
    # tracker = Schedule(pandora_start, pandora_stop, target_list, obs_window, transit_coverage_min, \
    #      aux_key, aux_list, fname_tracker, sched_wts, commissioning_time); print(tracker)#, aux_key=None)

    tracker = Schedule(pandora_start, pandora_stop, target_list, obs_window, transit_coverage_min, sched_wts, \
        aux_key='closest', aux_list=aux_list, fname_tracker = fname_tracker, commissioning_time=commissioning_time, \
            sched_start=sched_start, sched_stop=sched_stop)
    
#Default values for parameters used
if __name__ == "__main__":

    # Specify observing parameters
    obs_window = timedelta(hours=24.0)
    pandora_start = "2025-10-15 00:00:00"#"2025-09-01 00:00:00"
    pandora_stop = "2025-11-15 00:00:00"#"2026-10-01 00:00:00"
    sched_start= "2025-10-15 00:00:00"#"2025-09-01 00:00:00"
    sched_stop= "2025-11-15 00:00:00"#"2026-10-01 00:00:00"

    commissioning_time_ = 0

    # sched_wts[transit coverage, saa overlap, schedule factor]
    # sched_wts = [0.5, 0.25, 0.25]
    # sched_wts = [0.0, 0.0, 1.0]
    sched_wts = [0.8, 0.0, 0.2]
    # transit_coverage_min = 0.25
    transit_coverage_min = 0.5
    deprioritization_limit = 12 # hours
    
    #Mission requirements: >= 91 deg avoidance for Sun, >= 20 deg avoidance for Moon and Earth limbs
    # blocks = [91., 40., 63.]
    blocks = [91., 25., 63.]
    # target_list_name = 'Pandora_Target_List_Top20_29Aug2024_names_only'#'Pandora_Target_List_Top20_14May2024'
    # target_list =  target_list_name + '.csv'#'target_list_top20_16Feb2024.csv'
    # target_partner_list='target_partner_list.csv'#'target_list_top5_16Feb2024.csv'#
    gmat_file = 'GMAT_pandora_600_20240512.txt'#'GMAT_pandora_450_20230713.csv'#
    obs_name = 'Pandora_600km_20240518'#'Pandora_450km_20230713'#

    update_target_list_as_per_json_files = True
    if update_target_list_as_per_json_files:
        # targ_list = pd.read_csv(f"{PACKAGEDIR}/data/" + target_list, sep=",")
        # pl_names = ['GJ 9827 b', 'L 98-59 d', 'WASP-69 b']
        # which_targets = 'primary-exoplanet'#'occultation-standard'#'auxiliary-standard'#'auxiliary-exoplanet'#'secondary-exoplanet'#
        # Example usage:
        target_definition_files = ['primary-exoplanet', 'auxiliary-exoplanet', 'auxiliary-standard', 'occultation-standard', \
            'monitoring-standard', 'secondary-exoplanet']

        # keyword_ = 'primary-exoplanet'
        for keyword_ in target_definition_files:
            updated_targ_list_df = helper_codes.process_target_files(keyword_)#; print(updated_targ_list[updated_targ_list.columns[:10]].head())
            save_df_as_csv = True
            if save_df_as_csv:
                updated_targ_list_df.to_csv(f"{PACKAGEDIR}/data/" + keyword_ + '_targets.csv', index=False)
        
        updated_targ_list = keyword_ + '_targets.csv'
        primary_targ_list = target_definition_files[0] + '_targets.csv'

        occ_std_df = helper_codes.process_target_files(target_definition_files[3])
        occ_std_df.to_csv(f"{PACKAGEDIR}/data/aux_list_new.csv", index=False)

    #     targ_list = helper_codes.get_targets_table(which_targets)
    #     pl_names = targ_list['Planet Name'].values
    #     updated_targ_list = helper_codes.update_target_list(targ_list, pl_names, which_targets)
    #     # updated_targ_list.reset_index(drop=True)
    #     target_list_name = which_targets#'Pandora_Target_List_Top20_29Aug2024_updated_targ_list'
    #     # updated_targ_list.to_csv(PACKAGEDIR + "/data/Pandora_Target_List_Top20_29Aug2024_updated_targ_list.csv", index=False)

    # # url = "https://github.com/PandoraMission/PandoraTargetList/blob/main/target_definition_files/primary-exoplanet/GJ_1214b_target_definition.json"
    # # data = read_json_from_github(url)
    # print(updated_targ_list)
    fname_tracker = f"{PACKAGEDIR}/data/Tracker_{pandora_start[0:10]}_to_{pandora_stop[0:10]}.pkl"#f"{PACKAGEDIR}/data/Tracker_" + target_list_name + ".pkl"

    # Schedule_all_scratch(blocks, pandora_start, pandora_stop, target_definition_files, \
    #     #updated_targ_list, target_partner_list, \
    #     obs_window, transit_coverage_min, sched_wts = sched_wts, aux_key='max_visibility_any', \
    #         # aux_list=f"{PACKAGEDIR}/data/aux_list_new.csv", 
    #         fname_tracker = fname_tracker, commissioning_time=commissioning_time_, \
    #             sched_start = sched_start, sched_stop = sched_stop)
    #         # aux_key='closest', aux_list=f"{PACKAGEDIR}/data/aux_list.csv", commissioning_time=30)
    #
    # transits.star_vis(blocks[0], blocks[1], blocks[2], pandora_start, pandora_stop, gmat_file, obs_name, \
    # #     # save_pth = f'{PACKAGEDIR}/data/aux_targets/', targ_list = f"{PACKAGEDIR}/data/aux_list_new.csv")
    # #     save_pth = f'{PACKAGEDIR}/data/targets/', targ_list = f'{PACKAGEDIR}/data/{primary_targ_list}')
    #     save_pth = f'{PACKAGEDIR}/data/aux_targets/', targ_list = f'{PACKAGEDIR}/data/{target_definition_files[2]}_targets.csv')
        # save_pth = f'{PACKAGEDIR}/data/targets/', targ_list = f'{PACKAGEDIR}/data/target_partner_list.csv')
        # save_pth = f'{PACKAGEDIR}/data/aux_targets/', targ_list = f'{PACKAGEDIR}/data/aux_list.csv')
    # #                  
    #
    Schedule(pandora_start, pandora_stop, primary_targ_list, obs_window, transit_coverage_min, sched_wts, \
        aux_key='closest', aux_list=f"{PACKAGEDIR}/data/aux_list_new.csv", fname_tracker = fname_tracker, commissioning_time=30, \
            sched_start=sched_start, sched_stop=sched_stop)

    # Schedule(pandora_start, pandora_stop, obs_window, transit_coverage_min, sched_wts, \
    #          commissioning_time=30, sched_start=sched_start, sched_stop=sched_stop,
    #          aux_key='closest', aux_list=f"{PACKAGEDIR}/data/aux_list.csv")