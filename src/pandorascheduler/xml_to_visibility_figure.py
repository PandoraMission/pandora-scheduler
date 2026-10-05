import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.time import Time
from lxml import etree

# Use package-relative data path instead of a hardcoded user path
PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(PACKAGEDIR, "data")
fname = os.path.join(DATA_DIR, "Pandora_science_calendar.xml")

parser = etree.XMLParser(load_dtd=True, no_network=False)
tree = etree.parse(fname, parser)
root = tree.getroot()
df = pd.DataFrame(columns=['visitid', 'seqid', 'target', 'start', 'stop'])
idx = 0
for visit in root[1:]:
    visitid = visit[0].text
    for seq in visit[1:]:
        seqid = seq[0].text
        target = seq[1][0].text
        df.loc[idx] = [visitid, seqid, target, Time(seq[1][2][0].text).jd, Time(seq[1][2][1].text).jd]
        idx += 1
df = df.set_index('target')


for ii, row in df[0:].iterrows():
    obs_len = float(24*(row['stop'] - row['start']))
    if obs_len > 1.5:
        print(row['visitid'], row['seqid'], "%0.2f" % float(24*(row['stop'] - row['start'])))


fig, ax = plt.subplots(figsize=(10, 5))
targs = np.unique(np.asarray(df.index))
for idx, targ in enumerate(targs):
    z = np.asarray([df.loc[targ].start, df.loc[targ].stop]).T
    ax.plot(z.T, z.T**0 * idx, c='k')
ax.set(xlabel='Time [JD]', yticks=np.arange(0, len(targs)), yticklabels=targs, title='Pandora Calendar')