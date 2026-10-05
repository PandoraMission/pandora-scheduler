#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 08:58:02 2023

@author: paul
"""
import os
from datetime import timedelta

import numpy as np
import pandas as pd
import scheduler

starts=['2025-04-25 00:00:00','2025-05-25 00:00:00','2025-06-25 00:00:00','2025-07-25 00:00:00','2025-08-25 00:00:00','2025-09-25 00:00:00','2025-10-25 00:00:00','2025-11-25 00:00:00','2025-12-25 00:00:00','2026-01-25 00:00:00','2026-02-25 00:00:00','2026-03-25 00:00:00','2026-04-25 00:00:00']
stops=['2026-04-25 00:00:00','2026-05-25 00:00:00','2026-06-25 00:00:00','2026-07-25 00:00:00','2026-08-25 00:00:00','2026-09-25 00:00:00','2026-10-25 00:00:00','2026-11-25 00:00:00','2026-12-25 00:00:00','2027-01-25 00:00:00','2027-02-25 00:00:00','2027-03-25 00:00:00','2027-04-25 00:00:00']
obs_window = timedelta(hours=24.0)
sched_wts = [0.0, 0.0, 1.0]
transit_coverage_min = 0.25
blocks=[91.,40.,63.]
target_list='target_list.csv'
target_partner_list='target_partner_list.csv'
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))



for s in range(len(starts)):
    pandora_start=starts[s]
    pandora_stop=stops[s]
    
    
    tr=scheduler.Schedule(pandora_start, pandora_stop, obs_window, transit_coverage_min, sched_wts, commissioning_time=30, aux_key=None)
    tr=tr.sort_values('Planet Name', ignore_index=True)
    transits_needed=pd.DataFrame(np.array([tr['Planet Name'].tolist(),tr['Transits Needed'].tolist()]).T, columns=['Planet Name', 0])
    
    for i in range(99):
        tr=scheduler.Schedule(pandora_start, pandora_stop, obs_window, transit_coverage_min, sched_wts, commissioning_time=30, aux_key=None)
        tr=tr.sort_values('Planet Name', ignore_index=True)
        transits_needed=pd.concat([transits_needed,pd.DataFrame(np.array(tr['Transits Needed'].tolist()).T, columns=[i+1])], axis=1)
    
    transits_needed.to_csv((f"{PACKAGEDIR}/data/lost_time_{pandora_start.split(' ')[0]}.csv"), sep=",", index=False)