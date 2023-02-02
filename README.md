# pandora-scheduler

This tool is used to develop an observation schedule for Pandora SmallSat. The codes primary focus is on scheduling the required number of transiting events for each of the planets on a user provided target list csv file. It also has the functionality to schedule non-phase constrained events that users are interested in including.

For it's scheduling of transits, this tool currently has three factors whose weights can be changed by the user to optimize the schedule:

- Observing efficiency - calculated by the amount of gap time in between observing events
- SAA Crossing Overlap - how much of an observed transit will be captured during an South Atlantic Anomaly crossing of Pandora
- Transit Coverage - how much of a transit will be captured by Pandora while maintaining keep out regions with the Sun, Moon, and Earth


Requires
~~~~~~~~~~
- Internet to query external website for time conversion: https://astroutils.astronomy.osu.edu/time/bjd2utc.html
- General Mission Analysis Tool (GMAT), a free open source tool for space mission design. Required to create orbital position file. More info on GMAT can be found here:
https://software.nasa.gov/software/GSC-17177-1