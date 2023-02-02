# pandora-scheduler

Within transits.py current conversion of time from BJD_tdb to MJD_utc includes use of barycorr.py which requires querying the website: https://astroutils.astronomy.osu.edu/time/bjd2utc.html 
Future work: Current BJD_tdb to JD_utc converter exists as part of EXOFASTv2 package but it is in IDL and needs to be updated to python before implementation.