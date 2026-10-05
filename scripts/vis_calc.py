import os

import pandas as pd

from pandorascheduler import transits

fdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "data"))
target_list = "target_list.csv"
partner_list = "target_partner_list.csv"

sun_block = 90.0
moon_block = 40.0
earth_block = 63.0
obs_start = "2025-04-25 00:00:00"
obs_stop = "2026-04-25 00:00:00"

transits.star_vis(sun_block, moon_block, earth_block, obs_start, obs_stop)

targets = pd.read_csv(fdir + "/" + target_list, sep=",")
planet_name = targets["Planet Name"]
star_name = targets["Star Name"]
for i in range(len(targets)):
    print(planet_name[i])
    transits.transit_timing(target_list, planet_name[i], star_name[i])
    transits.SAA_overlap(planet_name[i], star_name[i])

partners = pd.read_csv(fdir + "/" + partner_list, sep=",")
partner_name = partners["Planet Name"]
p_star_name = partners["Star Name"]
for i in range(len(partners)):
    print(partner_name[i])
    transits.transit_timing(partner_list, partner_name[i], p_star_name[i])

stars = star_name.drop_duplicates(
    keep="first",
).reset_index()
for i in range(len(stars)):
    transits.Transit_overlap(target_list, partner_list, stars["Star Name"][i])
