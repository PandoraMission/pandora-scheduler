import os
import numpy as np
import pandas as pd
from astropy.time import Time
from datetime import datetime, timedelta

# from . import PACKAGEDIR
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))


def Schedule(
    pandora_start: str,
    pandora_stop: str,
    obs_window: timedelta,
    transit_coverage_min: float,
    sched_wts: list,
    sched_start=None,
    sched_stop=None,
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
    sched_start:            string
                            Date and time to start scheduling
                            ex. '2025-04-25 00:00:00'
    sched_stop:             string
                            Date and time to stop scheduling
                            ex. '2026-04-25 00:00:00'

    Returns
    -------
    csv file
        file containing schedule for Pandora
    """
    assert sum(sched_wts) == 1., "Sum of weights should equal 1"

    if sched_start == None:
        sched_start = pandora_start
    if sched_stop == None:
        sched_stop = pandora_stop

    # Convert times to datetime
    pandora_start = datetime.strptime(pandora_start, "%Y-%m-%d %H:%M:%S")
    pandora_stop = datetime.strptime(pandora_stop, "%Y-%m-%d %H:%M:%S")
    sched_start = datetime.strptime(sched_start, "%Y-%m-%d %H:%M:%S")
    sched_stop = datetime.strptime(sched_stop, "%Y-%m-%d %H:%M:%S")

    # Import target list
    target_list = pd.read_csv(f"{PACKAGEDIR}/data/target_list.csv", sep=",")

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
        np.ones(len(planet_names)) * 10, columns=["Transits Needed"]
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

        print(
            planet_name,
            "has",
            len(planet_data),
            "transits with greater transit coverage than",
            transit_coverage_min,
        )

        start_transits = Time(
            planet_data["Transit_Start"], format="mjd", scale="utc"
        ).to_value("datetime")
        end_transits = Time(
            planet_data["Transit_Stop"], format="mjd", scale="utc"
        ).to_value("datetime")
        p_trans = planet_data.index[
            (pandora_start <= start_transits) & (end_transits <= pandora_stop)
        ]
        s_trans = planet_data.index[
            (sched_start <= start_transits) & (end_transits <= sched_stop)
        ]
        pandora_transits_left.append(len(p_trans))
        schedule_transits_left.append(len(s_trans))
    # breakpoint()
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
            "Transit Factor" "Quality Factor",
        ],
    )

    while stop <= sched_stop:
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

            print("Scheduled no phase event", target)
            start = obs_stop
            stop = start + obs_window
            continue

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
                start_transits = planet_data["Transit_Start"]
                end_transits = planet_data["Transit_Stop"]
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
                    start_transits[j] = start_transits[j] - timedelta(
                        seconds=start_transits[j].second,
                        microseconds=start_transits[j].microsecond,
                    )
                    end_transits[j] = end_transits[j] - timedelta(
                        seconds=end_transits[j].second,
                        microseconds=end_transits[j].microsecond,
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
                        gap_time = obs_start - obs_rng[0]
                        s_factor = 1 - (gap_time / obs_window)  # maximize

                        # Calc a quality factor (currently based on transit coverage, SAA crossing, scheduling efficiency)
                        trans_cover = planet_data["Transit_Coverage"][j]  # maximize
                        if trans_cover < transit_coverage_min:
                            print("how", planet_name)
                            print(planet_data)
                            breakpoint()
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
        if len(temp_df) == 0:
            free = [["Free Time", start, stop]]
            free = pd.DataFrame(
                free, columns=["Target", "Observation Start", "Observation Stop"]
            )
            sched_df = pd.concat([sched_df, free], axis=0)

            print("No transits in window!", start, stop)
            start = stop
            stop = start + obs_window
            continue

        else:
            # Check Transit Factor first for planets running out of available transits
            # if (temp_df["Transit Factor"] < 1.5).any():
            if (temp_df["Transit Factor"] < 2).any():
                temp_df = temp_df.sort_values(by=["Transit Factor"]).reset_index(
                    drop=True
                )
                print(temp_df["Planet Name"][0], "Transit Factor Warning")

            # Otherwise sort by Quality Factor and schedule best target
            else:
                temp_df = temp_df.sort_values(
                    by=["Quality Factor"], ascending=False
                ).reset_index(drop=True)

            planet_name = temp_df["Planet Name"][0]
            obs_start = temp_df["Obs Start"][0]
            obs_stop = obs_start + timedelta(hours=24)
            trans_cover = temp_df["Transit Coverage"][0]
            saa_cover = temp_df["SAA Overlap"][0]
            s_factor = temp_df["Schedule Factor"][0]
            q_factor = temp_df["Quality Factor"][0]

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

            print(
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
    # save_fname = "Pandora_Schedule.csv"
    save_fname = f"Pandora_Schedule_{sched_wts[0]}_{sched_wts[1]}_{sched_wts[2]}.csv"
    sched_df.to_csv((f"{PACKAGEDIR}/data/" + save_fname), sep=",", index=False)


if __name__ == "__main__":

    # Specify observing parameters
    obs_window = timedelta(hours=24.0)
    pandora_start = "2025-04-25 00:00:00"
    pandora_stop = "2026-04-25 00:00:00"

    # sched_wts[transit coverage, saa overlap, schedule factor]
    # sched_wts = [0.5, 0.25, 0.25]
    sched_wts = [0.0, 0.0, 1.0]
    # sched_wts = [0., 0., 0.]
    transit_coverage_min = 0.25
    Schedule(pandora_start, pandora_stop, obs_window, transit_coverage_min, sched_wts)
    # sched_start, sched_stop)
