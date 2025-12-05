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
from multiprocessing import Pool
from functools import partial
from typing import Optional

import warnings
from erfa import ErfaWarning

# Suppress only ERFA warnings
warnings.filterwarnings('ignore', category=ErfaWarning)

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
OUTPUT_DIR = os.path.join(PACKAGEDIR, "data", "baseline")
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def Schedule(
    pandora_start: str,
    pandora_stop: str,
    target_list: str,
    obs_window: timedelta,
    transit_coverage_min: float,
    sched_wts: list,
    min_visibility: float,
    deprioritization_limit: float, 
    aux_key: str,
    aux_list:str,
    fname_tracker: str,
    commissioning_time: int,
    sched_start: str = None,
    sched_stop: str = None,
    output_dir: Optional[str] = None,
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

    # Convert times to datetime and add in commissioning time
    pandora_start = datetime.strptime(pandora_start, "%Y-%m-%d %H:%M:%S") + timedelta(days = commissioning_time)
    pandora_stop = datetime.strptime(pandora_stop, "%Y-%m-%d %H:%M:%S")

    if sched_start == None:
        sched_start = pandora_start
    else:
        sched_start = datetime.strptime(sched_start, "%Y-%m-%d %H:%M:%S") + timedelta(days = commissioning_time)
    if sched_stop == None:
        sched_stop = pandora_stop
    else:
        sched_stop = datetime.strptime(sched_stop, "%Y-%m-%d %H:%M:%S")

    target_list = pd.read_csv(f"{target_list}", sep=",")

    output_dir = output_dir or OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    # Import no phase events if they exist; keep empty lists otherwise
    nophase_targets = []
    nophase_starts = []
    nophase_stops = []
    too_list_path = f"{PACKAGEDIR}/data/ToO_list.csv"
    if os.path.exists(too_list_path):
        nophase_list = pd.read_csv(too_list_path, sep=",")
        nophase_targets = nophase_list["Target"].tolist()
        nophase_starts = [
            datetime.strptime(d, "%Y-%m-%d %H:%M:%S")
            for d in nophase_list["Obs Window Start"]
        ]
        nophase_stops = [
            datetime.strptime(d, "%Y-%m-%d %H:%M:%S")
            for d in nophase_list["Obs Window Stop"]
        ]
    else:
        logger.info("No ToO_list.csv found; skipping Target of Opportunity scheduling.")

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
    primary_target = pd.DataFrame(
        np.array(target_list["Primary Target"]), columns=["Primary Target"]
    )

    ras_tmp = pd.DataFrame(
        np.array(target_list["RA"]), columns=["RA"]
    )

    decs_tmp = pd.DataFrame(
        np.array(target_list["DEC"]), columns=["DEC"]
    )

    tracker = pd.concat([planet_names, primary_target, ras_tmp, decs_tmp, transit_need, transit_have], axis=1)

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

        logger.info(
            "%s has %s transits with coverage above %.2f",
            planet_name,
            len(planet_data),
            transit_coverage_min,
        )

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
            "RA",
            "Dec",
            "Observation Start",
            "Observation Stop",
            "Transit Coverage",
            "SAA Overlap",
            "Schedule Factor",
            "Transit Factor",
            "Quality Factor",
            "Comments",
        ],
    )
    non_primary_obs_time = {}
    all_target_obs_time = {}
    last_std_obs = datetime(2025, 12, 1)
    logger.warning(
        "Update last standard observation reference if required (currently %s)",
        last_std_obs,
    )
    
    while stop <= sched_stop:

        print(f"Start/Stop: {start}, {stop}")

        logger.debug("Evaluating window %s to %s", start, stop)
        tracker = tracker.sort_values(by=['Primary Target', 'Transit Priority'], ascending=[False, True]).reset_index(drop=True)
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
                "Comments",
            ],
        )
        ### First check if a Target of Opportunity is within observing window
        overlap_nophase = obs_rng.intersection(nophase_starts)
        if len(overlap_nophase) > 0:
            logger.info("Attempting to schedule Target of Opportunity")
            obs_start = nophase_starts[nophase_starts.index(overlap_nophase[0])]
            obs_stop = nophase_stops[nophase_starts.index(overlap_nophase[0])]
            ToO = nophase_targets[nophase_starts.index(overlap_nophase[0])]

            # 1) Check for planets with tf <= 1
            critical_planets = tracker[(tracker["Transits Needed"] > 0) & (tracker["Transits Left in Lifetime"] / tracker["Transits Needed"] <= 1)]

            forced_observation = False
            for _, planet in critical_planets.iterrows():
                planet_name = planet["Planet Name"]
                
                # 2) Check if a transit occurs during observing window
                star_name_tmp = pd.Series(planet_name).apply(helper_codes.remove_suffix)
                planet_data = pd.read_csv(f"{PACKAGEDIR}/data/targets/{star_name_tmp.iloc[0]}/{planet_name}/Visibility for {planet_name}.csv")

                start_transit = Time(planet_data["Transit_Start"].iloc[-1], format='mjd', scale='utc')
                end_transit = Time(planet_data["Transit_Stop"].iloc[-1], format='mjd', scale='utc')  
                start_transit = start_transit.datetime.replace(second=0, microsecond=0)
                end_transit = end_transit.datetime.replace(second=0, microsecond=0)

                early_start = end_transit - timedelta(hours=20)
                late_start = start_transit - timedelta(hours=4)
                
                start_rng = pd.date_range(early_start, late_start, freq="min")
                overlap_times = obs_rng.intersection(start_rng)
                
                if len(overlap_times) > 0:
                    # 4) Schedule the planet observation
                    forced_observation = True
                    obs_start = overlap_times[0]

                    if obs_rng[0] < obs_start:
                        free = [[f"FREE PRE-TOO, REPLACE WITH AUX", obs_rng[0], obs_start]]
                        free = pd.DataFrame(free, columns=["Target", "Observation Start", "Observation Stop"])
                        sched_df = pd.concat([sched_df, free], axis=0).reset_index(drop=True)
                    
                    sched = [[planet_name, obs_start, obs_stop, f"{planet_name} forced over ToO due to transit factor <=1"]]
                    sched = pd.DataFrame(sched, columns=["Target", "Observation Start", "Observation Stop", "Comments"])
                    
                    if sched_df.empty:
                        sched_df = sched.copy()
                    else:
                        sched_df = pd.concat([sched_df, sched], axis=0).reset_index(drop=True)
                    
                    logger.info(
                        "Forced observation of %s over ToO due to critical transit factor",
                        planet_name,
                    )
                    break  # Stop after scheduling the first critical planet
            
            if not forced_observation:
                # 3) or 5) Schedule the ToO
                tf_warning = ""
                for _, planet in tracker[tracker["Transits Needed"] > 0].iterrows():
                    planet_name = planet["Planet Name"]
                    tf = planet["Transits Left in Lifetime"] / planet["Transits Needed"]
                    
                    # Check if planet is transiting during ToO
                    star_name_tmp = pd.Series(planet_name).apply(helper_codes.remove_suffix)
                    planet_data = pd.read_csv(f"{PACKAGEDIR}/data/targets/{star_name_tmp.iloc[0]}/{planet_name}/Visibility for {planet_name}.csv")

                    start_transit = Time(planet_data["Transit_Start"].iloc[-1], format='mjd', scale='utc')
                    end_transit = Time(planet_data["Transit_Stop"].iloc[-1], format='mjd', scale='utc')  
                    start_transit = start_transit.datetime.replace(second=0, microsecond=0)
                    end_transit = end_transit.datetime.replace(second=0, microsecond=0)

                    early_start = end_transit - timedelta(hours=20)
                    late_start = start_transit - timedelta(hours=4)

                    start_rng = pd.date_range(early_start, late_start, freq="min")
                    overlap_times = obs_rng.intersection(start_rng)
            
                    if len(overlap_times) > 0 and tf > 1:
                        tf_warning += f"Warning: {planet_name} has MTRM > 1 and is transiting during ToO. "

                if obs_rng[0] < obs_start:
                    free = [[f"FREE PRE-TOO, REPLACE WITH AUX", obs_rng[0], obs_start]]
                    free = pd.DataFrame(free, columns=["Target", "Observation Start", "Observation Stop"])
                    sched_df = pd.concat([sched_df, free], axis=0).reset_index(drop=True)

                sched = [[ToO, obs_start, obs_stop, tf_warning.strip()]]
                sched = pd.DataFrame(sched, columns=["Target", "Observation Start", "Observation Stop", "Comments"])
                
                if sched_df.empty:
                    sched_df = sched.copy()
                else:
                    sched_df = pd.concat([sched_df, sched], axis=0).reset_index(drop=True)
                
                logger.info("Scheduled Target of Opportunity: %s", ToO)
                if tf_warning:
                    logger.warning(tf_warning)

            # Update the observation window
            start = obs_stop
            stop = start + obs_window
            continue
        # ## Next look at each planet and determine if transit occurs in window
        temp_df = helper_codes.check_if_transits_in_obs_window(tracker, temp_df, target_list, start, \
            pandora_start, pandora_stop, sched_start, sched_stop, obs_rng, obs_window, sched_wts, transit_coverage_min)

        ### Check if there's no transit occurring during the observing window
        ### Schedule auxiliary observation if possible
        if len(temp_df) == 0:
            aux_df, log_info, non_primary_obs_time, last_std_obs = Schedule_aux(start, stop, aux_key, \
                non_primary_obs_time=non_primary_obs_time, min_visibility = min_visibility, \
                    deprioritization_limit = deprioritization_limit, last_std_obs = last_std_obs)

            if sched_df.empty:
                sched_df = aux_df.copy()
            else:        
                sched_df = pd.concat([sched_df, aux_df], axis=0)

            # Update observation time for auxiliary targets
            if (len(aux_df) > 0) and (aux_df["Target"].values.all() != "Free Time"):
                for _, row in aux_df.iterrows():
                    if row['Target'] != "Free Time":
                        aux_target = row['Target']
                        aux_duration = row['Observation Stop'] - row['Observation Start']
                        all_target_obs_time[aux_target] = all_target_obs_time.get(aux_target, timedelta()) + aux_duration

            logger.info(log_info, start, stop)
            start = stop
            stop = start + obs_window
            continue                

        else: # THERE ARE TRANSITS DURING THE OBSERVING WINDOW
            if (temp_df["Transit Factor"] <= 2).any():
                temp_df = temp_df.sort_values(by=["Transit Factor"]).reset_index(
                    drop=True
                )
            else:
                # Otherwise, sort by "Quality Factor" (descending) and then by "Transit Factor" (ascending)
                temp_df = temp_df.sort_values(
                    by=["Quality Factor", "Transit Factor"],
                    ascending=[False, True]
                ).reset_index(drop=True)

            planet_name = temp_df["Planet Name"][0]
            ra_tmp, dec_tmp = temp_df["RA"][0], temp_df["DEC"][0]
            star_name = target_list["Star Name"][np.where(target_list["Planet Name"] == planet_name)[0][0]]
            obs_start = temp_df["Obs Start"][0]
            obs_stop = obs_start + timedelta(hours=24)
            trans_cover = temp_df["Transit Coverage"][0]
            saa_cover = temp_df["SAA Overlap"][0]
            s_factor = temp_df["Schedule Factor"][0]
            q_factor = temp_df["Quality Factor"][0]

            # VK BEGIN
            print_every_day = False
            if print_every_day:
                num_sig_digits_ = 2
                if obs_start == pandora_start:
                    logger.debug(
                        "First-day schedule candidate:\n%s",
                        temp_df.head(1).round(num_sig_digits_).to_string(index=False),
                    )
                else:
                    df_string = temp_df.head(1).round(num_sig_digits_).to_string(index=False)# Get the string representation of the DataFrame with headers
                    lines = df_string.split('\n')# Split the string into lines
                    if len(lines) > 1:# Print only the data line (second line, index 1)
                        logger.debug("Daily schedule summary: %s", lines[1])
            # VK END

            if obs_rng[0] < obs_start: # primary target not visible for some time at the start of the visit --> find auxiliary target for that time
                aux_df, log_info, non_primary_obs_time, last_std_obs = Schedule_aux(start, obs_start, aux_key, \
                    non_primary_obs_time=non_primary_obs_time, \
                        min_visibility = min_visibility, deprioritization_limit = deprioritization_limit, last_std_obs = last_std_obs)
                
                if sched_df.empty:
                    sched_df = aux_df.copy()
                else:        
                    sched_df = pd.concat([sched_df, aux_df], axis=0)

                # Update observation time for auxiliary targets
                if len(aux_df) > 0:
                    for _, row in aux_df.iterrows():
                        if row['Target'] != "Free Time":
                            aux_target = row['Target']
                            aux_duration = row['Observation Stop'] - row['Observation Start']
                            all_target_obs_time[aux_target] = all_target_obs_time.get(aux_target, timedelta()) + aux_duration
                logger.info(log_info, start, obs_start)

            sched = [
                [
                    planet_name,
                    ra_tmp, 
                    dec_tmp,
                    obs_start,
                    obs_stop,
                    trans_cover,
                    saa_cover,
                    s_factor,
                    q_factor,
                    np.nan,
                ]
            ]
            sched = pd.DataFrame(
                sched,
                columns=[
                    "Target",
                    "RA",
                    "DEC",
                    "Observation Start",
                    "Observation Stop",
                    "Transit Coverage",
                    "SAA Overlap",
                    "Schedule Factor",
                    "Quality Factor",
                    "Comments",
                ],
            )

            if sched_df.empty:
                sched_df = sched.copy()
            else:        
                sched_df = pd.concat([sched_df, sched], axis=0)

            all_target_obs_time[planet_name] = all_target_obs_time.get(planet_name, timedelta()) + (obs_stop - obs_start)

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

            logger.info(
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
    save_fname = f"Pandora_Schedule_{sched_wts[0]}_{sched_wts[1]}_{sched_wts[2]}_{pandora_start.strftime('%Y-%m-%d')}_to_{pandora_stop.strftime('%Y-%m-%d')}.csv"
    sched_df.to_csv(os.path.join(output_dir, save_fname), sep=",", index=False)

    tracker.to_csv(os.path.join(output_dir, "tracker.csv"), sep=",", index=False)

    # Save tracker
    tracker_path = fname_tracker
    if not os.path.isabs(tracker_path):
        tracker_path = os.path.join(output_dir, tracker_path)
    with open(tracker_path, 'wb') as file:
        pickle.dump(tracker, file)

    # Save all_target_obs_time
    report_name = f"Observation_Time_Report_{pandora_start}.csv"
    helper_codes.save_observation_time_report(
        all_target_obs_time,
        target_list,
        os.path.join(output_dir, report_name),
    )

    return tracker

def Schedule_aux(start, stop, aux_key, non_primary_obs_time, min_visibility, deprioritization_limit, last_std_obs, **kwargs):

    obs_rng = pd.date_range(start, stop, freq = "min")

    obs_std_dur = timedelta(hours = 0.5) 

    # Add standard stars!
    if (start - last_std_obs > timedelta(days = 3.)) and ((start + obs_std_dur) < stop):

        std_fn = f"{PACKAGEDIR}/data/monitoring-standard_targets.csv"
        std_targs = pd.read_csv(std_fn).reset_index(drop=True)
        std_targs = std_targs.sort_values('Priority', ascending=False).reset_index(drop=True)
        std_names = std_targs['Star Name']
        std_ras = std_targs['RA']
        std_decs = std_targs['DEC']
        std_priority = std_targs['Priority']

        for nn in range(len(std_names)):
            vis_file = f"{PACKAGEDIR}/data/aux_targets/{std_names[nn]}/Visibility for {std_names[nn]}.csv"
            vis = pd.read_csv(vis_file, usecols=["Time(MJD_UTC)", "Visible"])

            vis_times = Time(
                vis["Time(MJD_UTC)"].to_numpy(),
                format="mjd",
                scale="utc",
            ).to_datetime()
            vis_times = pd.to_datetime(vis_times)
            time_mask = (vis_times >= start) & (vis_times <= start + obs_std_dur)
            vis_filtered = vis.loc[time_mask]

            if not vis_filtered.empty and vis_filtered['Visible'].all():
                std_name = std_names[nn]
                std_ra, std_dec = std_ras[nn], std_decs[nn]
                priority_std = std_priority[nn]
                print(f"--------> {std_name} scheduled for STD observations with full visibility.")
                logger.info(
                    "%s scheduled for STD observations with full visibility",
                    std_name,
                )
                break
        
        try:
            std_name
        except NameError:
            std_name = 'WARNING: no visible standard star'
            std_ra, std_dec = "NaN", "NaN"
            priority_std = 1.
            print(f"    ----> WARNING: no visible standard star for {start} to {start + obs_std_dur}")
            logger.warning(
                "No visible standard star between %s and %s",
                start,
                start + obs_std_dur,
            )

        # Schedule observations of a standard star 
        STD_obs = [[f"{std_name} STD", start, start + obs_std_dur, std_ra, std_dec]]
        start += obs_std_dur
        last_std_obs = start
        std_priority = priority_std

        from pandas import Timedelta
        existing_STD_time, existing_STD_priority = non_primary_obs_time.get("STD", (timedelta(), 0))
        non_primary_obs_time["STD"] = (Timedelta(existing_STD_time) + obs_std_dur, std_priority)

    if aux_key == None:
        # No aux target scheduling, free time
        free = [["Free Time", start, stop]]
        aux_df = pd.DataFrame(
            free, columns=["Target", "Observation Start", "Observation Stop"]
        )
        log_info = 'Free time, no observation scheduled.'
        return aux_df, log_info, non_primary_obs_time, last_std_obs
        
    for tdf_idx in range(1,len(target_definition_files)):
        aux_fn = f"{PACKAGEDIR}/data/{target_definition_files[tdf_idx]}_targets.csv"
        aux_targs = pd.read_csv(aux_fn).reset_index(drop=True)

        # Create a set of names from aux_targs for efficient lookup
        aux_names = set(aux_targs['Star Name'])

        # Create a dictionary of name-priority pairs from non_primary_obs_time
        non_primary_priorities = {name: priority[-1] for name, (_, priority) in non_primary_obs_time.items() if (not name.endswith(('STD')))}
        mask = (aux_targs['Star Name'].isin(non_primary_priorities.keys())) & \
            (aux_targs['Priority'] != aux_targs['Star Name'].map(non_primary_priorities))

        if aux_key == 'sort_by_tdf_priority':
            if mask.any():
                aux_targs.loc[mask, 'Priority'] = aux_targs.loc[mask, 'Star Name'].map(non_primary_priorities)

            aux_targs = aux_targs.sort_values('Priority', ascending=False).reset_index(drop=True)
            names = aux_targs['Star Name']
            ras = aux_targs['RA']
            decs = aux_targs['DEC']
            aux_priority = aux_targs['Priority']

        vis_all_targs = []
        vis_any_targs = []
        targ_vis = []

        for n in tqdm(range(len(names)), desc=f"Finding visible non-primary target for {str(start)} to {str(stop)} from '{target_definition_files[tdf_idx]}'"):

            try:
                vis_file = f"{PACKAGEDIR}/data/aux_targets/{names[n]}/Visibility for {names[n]}.csv"
                vis = pd.read_csv(vis_file, usecols=["Time(MJD_UTC)", "Visible"])
                vis_times = Time(
                    vis["Time(MJD_UTC)"].to_numpy(),
                    format="mjd",
                    scale="utc",
                ).to_datetime()
                vis_times = pd.to_datetime(vis_times)
                time_mask = (vis_times >= start) & (vis_times <= stop)
                vis_filtered = vis.loc[time_mask]

                if not vis_filtered.empty and vis_filtered['Visible'].all():
                    vis_all_targs.append(n)
                    break
                elif not vis_filtered.empty and vis_filtered['Visible'].any():
                    vis_any_targs.append(n)
                    targ_vis.append(100*(np.sum(vis_filtered['Visible'])/len(vis_filtered)))
            except FileNotFoundError:
                pass

        # Select target
        # If at least one target is visible the entire time, use vis_all_targs
        if len(vis_all_targs) > 0:
            if aux_key == ('sort_by_tdf_priority') or (aux_key == 'closest'):
                idx = 0 
            else:
                idx = np.random.randint(0,len(vis_all_targs))
            name = names[vis_all_targs[idx]]
            ra_tmp, dec_tmp = ras[vis_all_targs[idx]], decs[vis_all_targs[idx]]
            priority_tmp = aux_priority[vis_all_targs[idx]]
            print(f"--------> {name} scheduled for non-primary observations with full visibility from {target_definition_files[tdf_idx]}.")
            logger.info(
                "%s scheduled for non-primary observations with full visibility from %s",
                name,
                target_definition_files[tdf_idx],
            )
            log_info=f"{name} scheduled for non-primary observation with full visibility."
            break
        
        # If at least one target is partially visible, use vis_any_targs
        elif len(vis_any_targs) > 0:
            idx = np.asarray(targ_vis).argmax()
            vis_perc = targ_vis[idx]
            if vis_perc >= 100*min_visibility:
                name = names[vis_any_targs[idx]]
                ra_tmp, dec_tmp = ras[vis_any_targs[idx]], decs[vis_any_targs[idx]]
                priority_tmp = aux_priority[vis_any_targs[idx]]
                print(f"--------> No non-primary target with full visibility; {name} scheduled for non-primary observations with {vis_perc:.2f}% visibility from {target_definition_files[tdf_idx]}.")
                logger.info(
                    "No non-primary target with full visibility; %s scheduled at %.2f%% visibility from %s",
                    name,
                    vis_perc,
                    target_definition_files[tdf_idx],
                )
                log_info=f"No non-primary target with full visibility; {name} scheduled for non-primary observations with {vis_perc:.2f}% visibility from {target_definition_files[tdf_idx]}."
                break
            else:
                print(f"No non-primary target with visibility greater than {100*min_visibility}% from {target_definition_files[tdf_idx]}, try next target list...")
                logger.warning(
                    "No non-primary target with visibility greater than %.2f%% from %s",
                    100 * min_visibility,
                    target_definition_files[tdf_idx],
                )
                continue
        
    try:
        name
    except NameError:
        name = "Free Time"
        priority_tmp = 0.
        log_info = "No fuly or partially visible non-primary targets, Free Time..."
        
    if (not name.endswith(('STD'))):
        existing_time, existing_priority = non_primary_obs_time.get(name, (np.array([]), np.array([])))
        new_time = existing_time + (stop - start) if isinstance(existing_time, timedelta) else stop - start
        if isinstance(existing_time, np.ndarray):
            non_primary_obs_time[name] = (
                np.append(existing_time, new_time),
                np.append(existing_priority, priority_tmp)
            )
        else:
            non_primary_obs_time[name] = (np.array([new_time]), np.array([priority_tmp]))

        total_time = np.sum(non_primary_obs_time[name][0]) if isinstance(non_primary_obs_time[name][0], np.ndarray) else non_primary_obs_time[name][0]

        if total_time > timedelta(hours=deprioritization_limit):
            print(f"----------------------------> Deprioritize {name} <----------------------------")
            logger.warning("Deprioritizing %s due to accumulated auxiliary time", name)
            new_priority = 0.95 * aux_priority[len(names) - 1]
            if isinstance(non_primary_obs_time[name][0], np.ndarray):
                non_primary_obs_time[name] = (
                    np.append(non_primary_obs_time[name][0], new_time),
                    np.append(non_primary_obs_time[name][1], new_priority)
                )
            else:
                non_primary_obs_time[name] = (np.array([new_time]), np.array([new_priority]))

    sched = [[name, start, stop, ra_tmp, dec_tmp]]

    try:
        sched = STD_obs + sched
    except NameError:
        no_std = True

    aux_df = pd.DataFrame(sched, columns=["Target", "Observation Start", "Observation Stop", "RA", "DEC"])

    return aux_df, log_info, non_primary_obs_time, last_std_obs

def Schedule_all_scratch(
    blocks:list,
    pandora_start:str,
    pandora_stop:str,
    primary_targ_list:str, 
    aux_targ_list:str, 
    target_definition_files:list,
    obs_window: timedelta,
    transit_coverage_min: float,
    sched_wts: list,
    aux_key: str,
    fname_tracker: str,
    commissioning_time: int,
    sched_start=None,
    sched_stop=None,
    output_dir: Optional[str] = None,
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
    
    output_dir = output_dir or OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    #Create visibity information for target host stars given avoidance angles 
    #for Sun, Moon, and Earth limb during Pandora's science observation lifetime.

    tdf = np.asarray(target_definition_files)

    for ii1 in range(len(target_definition_files)):
        fn_ = tdf[tdf == tdf[ii1]][0]
        tmp_csv = pd.read_csv(f'{PACKAGEDIR}/data/{fn_}_targets.csv', sep=',')
        tmp_planet_name_lst = tmp_csv['Planet Name']
        tmp_star_name_lst = tmp_csv['Star Name']

        if "exoplanet" in fn_:
            sub_dir = "targets"
        else:
            sub_dir = "aux_targets"

        vis_done = []
        for s in range(len(tmp_star_name_lst)):
            vis_done.append(os.path.exists(f'{PACKAGEDIR}/data/{sub_dir}/{tmp_star_name_lst[s]}/Visibility for {tmp_star_name_lst[s]}.csv'))

        if np.all(vis_done) == False:
            transits.star_vis(blocks[0], blocks[1], blocks[2], pandora_start, pandora_stop, gmat_file, obs_name, \
                save_pth = f'{PACKAGEDIR}/data/{sub_dir}/', targ_list = f'{PACKAGEDIR}/data/{fn_}_targets.csv')

        # Determine primary transits for targets and partners during Pandora's science observation lifetime.
        if "exoplanet" in fn_:
            for s in range(len(tmp_star_name_lst)):
                if not os.path.exists(f'{PACKAGEDIR}/data/{sub_dir}/{tmp_star_name_lst[s]}/{tmp_planet_name_lst[s]}'):
                    transits.transit_timing(f'{fn_}_targets.csv', tmp_planet_name_lst[s], tmp_star_name_lst[s])
    
    targets_prim_csv = pd.read_csv(f'{PACKAGEDIR}/data/{target_definition_files[0]}_targets.csv', sep=',')
    targets_aux_csv = pd.read_csv(f'{PACKAGEDIR}/data/{target_definition_files[1]}_targets.csv', sep=',')
    t_list = targets_prim_csv['Planet Name']
    ts_list = targets_prim_csv['Star Name']
    tp_list = targets_aux_csv['Planet Name']
    ps_list = targets_aux_csv['Star Name']
    for t in tqdm(range(len(t_list))):
        try:
            vis = pd.read_csv(f'{PACKAGEDIR}/data/targets/{ts_list[t]}/{t_list[t]}/Visibility for {t_list[t]}.csv')
            t_over = vis['Transit_Overlap']
        except KeyError:
            transits.Transit_overlap(os.path.basename(primary_targ_list), os.path.basename(primary_targ_list), ts_list[t])

        try:
            vis = pd.read_csv(f'{PACKAGEDIR}/data/targets/{ts_list[t]}/{t_list[t]}/Visibility for {t_list[t]}.csv')
            saa = vis['SAA_Overlap']
        except KeyError:
            transits.SAA_overlap(t_list[t], ts_list[t])
    
    
    #Schedule observations for the scheduling period
    tracker = Schedule(pandora_start, pandora_stop, primary_targ_list, obs_window, transit_coverage_min, sched_wts, min_visibility, deprioritization_limit, \
        aux_key = aux_key, aux_list=aux_targ_list, fname_tracker = fname_tracker, commissioning_time = commissioning_time_, \
            sched_start = sched_start, sched_stop = sched_stop, output_dir=output_dir)
    
# Default values for parameters used
if __name__ == "__main__":

    obs_window = timedelta(hours=24.0)
    pandora_start = "2026-02-05 00:00:00"
    pandora_stop = "2027-02-05 00:00:00"
    sched_start= "2026-02-05 00:00:00"
    sched_stop= "2027-02-05 00:00:00"

    commissioning_time_ = 0  # days

    sched_wts = [0.8, 0.0, 0.2]
    transit_coverage_min = 0.4
    deprioritization_limit = 48  # hours
    min_visibility = 0.5  # for non-primary targets
    
    blocks = [91., 25., 86.]
    gmat_file = 'Pandora-600km-withoutdrag-20251018.txt'
    obs_name = 'Pandora'

    output_dir = OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    output_dir = OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    update_target_list_as_per_json_files = True
    if update_target_list_as_per_json_files:
        target_definition_files = ['exoplanet', 'auxiliary-standard', 'monitoring-standard', 'occultation-standard']

        for keyword_ in target_definition_files:
            fn_tmp = f"{PACKAGEDIR}/data/{keyword_}_targets.csv"
            if os.path.exists(fn_tmp):
                    continue
            else:
                updated_targ_list_df = helper_codes.process_target_files(keyword_)
                updated_targ_list_df.to_csv(fn_tmp, index=False)

        primary_targ_list = f"{PACKAGEDIR}/data/{target_definition_files[0]}_targets.csv"
        aux_targ_list = f"{PACKAGEDIR}/data/occultation-standard_targets.csv"
    fname_tracker = os.path.join(output_dir, f"Tracker_{pandora_start[0:10]}_to_{pandora_stop[0:10]}.pkl")

    if not os.path.exists(f"{PACKAGEDIR}/data/all_targets.csv"):
        create_aux_list = helper_codes.create_aux_list(target_definition_files, PACKAGEDIR)

    aux_key = 'sort_by_tdf_priority'
    # aux_targets = None

    run_ = 'vis_and_schedule'

    if run_ == 'schedule_only':
        Schedule(pandora_start, pandora_stop, primary_targ_list, obs_window, transit_coverage_min, sched_wts, min_visibility, deprioritization_limit, \
            aux_key = aux_key, aux_list=aux_targ_list, fname_tracker = fname_tracker, commissioning_time = commissioning_time_, \
                sched_start = sched_start, sched_stop = sched_stop, output_dir=output_dir)
    elif run_ == 'target_visibility':
        for tt in target_definition_files[0:5]:
            if 'exoplanet' in tt:
                save_path_ = f'{PACKAGEDIR}/data/targets/'
            else:
                save_path_ = f'{PACKAGEDIR}/data/aux_targets/'
            transits.star_vis(blocks[0], blocks[1], blocks[2], pandora_start, pandora_stop, gmat_file, obs_name, \
                save_pth = save_path_, targ_list = f'{PACKAGEDIR}/data/{tt}_targets.csv')
    elif run_ == 'vis_and_schedule':
            Schedule_all_scratch(blocks, pandora_start, pandora_stop, primary_targ_list, aux_targ_list, target_definition_files, \
                obs_window, transit_coverage_min, sched_wts = sched_wts, aux_key=aux_key, \
                    fname_tracker = fname_tracker, commissioning_time=commissioning_time_, \
                        sched_start = sched_start, sched_stop = sched_stop, output_dir=output_dir)
