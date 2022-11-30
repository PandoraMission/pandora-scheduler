"""Functions to calculate transits times from a target list"""
import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from datetime import timedelta
from progressbar import ProgressBar
import barycorr 

from . import PACKAGEDIR


def target_vis(sun_block, moon_block, earth_block, obs_start, obs_stop):
    """ Determine visibility for target(s) host star with Pandora given avoidance angles
    for Sun, Moon, and Earth limb.
        
    Parameters
    ----------
    sun_block:      float
                    Avoidance angle for the Sun
    moon_block:     float
                    Avoidance angle for the Moon
    earth_block:    float
                    Avoidance angle for the Earth limb
    obs_start:      string
                    Date and time of start of Pandora science observing 
                    ex. '2025-04-25 00:00:00'
    obs_stop:       string
                    Date and time of end of Pandora science observing 
                    ex. '2026-04-25 00:00:00'
                                    
    Returns
    -------
    csv file
        file containing target's host star visibility by Pandora 
    """

### Create an array of Pandora's science observation year with exactly 1 min spacings
    # datetime 
    dt_iso_utc = pd.date_range(obs_start, obs_stop, freq='min')
    t_jd_utc   = Time(dt_iso_utc.to_julian_date(), format='jd', scale='utc').value
    t_mjd_utc  = Time(t_jd_utc-2400000.5, format='mjd', scale='utc').value


### Read in GMAT results
    print('Importing GMAT data')
    gmat_data = pd.read_csv(f"{PACKAGEDIR}/data/GMAT_Pandora.txt", sep='\t')

    # Trim dataframe to slightly larger than date range of 
    # Pandora science lifetime defined as obs_start and obs_stop
    gmat_data = gmat_data[(gmat_data['Earth.UTCModJulian']>=(t_jd_utc[0]-2430000.0)-0.0007) & 
            (gmat_data['Earth.UTCModJulian']<=(t_jd_utc[-1]-2430000.0)+0.0007)]


### Convert GMAT times into standard MJD_UTC
    # Note: GMAT uses different offset for it's MJD (uses 2,430,000.0 rather than 2,400,000.5)
    gmat_mjd_utc =np.array(gmat_data['Earth.UTCModJulian']) + 2430000.0 - 2400000.5


### Extract GMAT positions in MJ2000 Earth fixed cartesian coordinates
    # Earth
    gmat_ex = np.array(gmat_data['Earth.EarthMJ2000Eq.X'])
    gmat_ey = np.array(gmat_data['Earth.EarthMJ2000Eq.Y'])
    gmat_ez = np.array(gmat_data['Earth.EarthMJ2000Eq.Z'])

    # Pandora
    gmat_px = np.array(gmat_data['Pandora.EarthMJ2000Eq.X'])
    gmat_py = np.array(gmat_data['Pandora.EarthMJ2000Eq.Y'])
    gmat_pz = np.array(gmat_data['Pandora.EarthMJ2000Eq.Z'])

    # Sun
    gmat_sx = np.array(gmat_data['Sun.EarthMJ2000Eq.X'])
    gmat_sy = np.array(gmat_data['Sun.EarthMJ2000Eq.Y'])
    gmat_sz = np.array(gmat_data['Sun.EarthMJ2000Eq.Z'])

    # Moon
    gmat_mx = np.array(gmat_data['Luna.EarthMJ2000Eq.X'])
    gmat_my = np.array(gmat_data['Luna.EarthMJ2000Eq.Y'])
    gmat_mz = np.array(gmat_data['Luna.EarthMJ2000Eq.Z'])

    # Pandora's Lat and Lon
    gmat_lat = np.array(gmat_data['Pandora.Earth.Latitude'])
    gmat_lon = np.array(gmat_data['Pandora.Earth.Longitude'])


### Interpolate all positions from GMAT to map to 1 minute intervals
    # Earth
    ex = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_ex)
    ey = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_ey)
    ez = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_ez)

    # Pandora
    px = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_px)
    py = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_py)
    pz = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_pz)

    # Moon
    mx = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_mx)
    my = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_my)
    mz = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_mz)

    # Sun
    sx = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_sx)
    sy = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_sy)
    sz = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_sz)

    # Pandora Latitude & Longitude
    gmat_lon_cont = np.copy(gmat_lon)
    for i in range(len(gmat_lon_cont)-1):
        if np.abs(gmat_lon_cont[i+1]-gmat_lon_cont[i]) > 100:
            gmat_lon_cont[i+1:] = gmat_lon_cont[i+1:] - 360
    p_lon_cont = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_lon_cont)
    p_lon = np.copy(p_lon_cont)
    for i in range(len(p_lon)):
        while p_lon[0] < -180:
            p_lon = p_lon + 360        
        if p_lon[i] < -180:
            p_lon[i:] = p_lon[i:] + 360

    p_lon = p_lon * u.deg
    p_lat = np.interp(t_mjd_utc, gmat_mjd_utc, gmat_lat) * u.deg


### Coordinate shift to Pandora reference frame
    #Earth 
    exx = ex-px
    eyy = ey-py
    ezz = ez-pz

    #Pandora
    pxx = px-px
    pyy = py-py
    pzz = pz-pz

    #Sun
    sxx = sx-px
    syy = sy-py
    szz = sz-pz

    #Moon
    mxx = mx-px
    myy = my-py
    mzz = mz-pz


### Create SkyCoord for angular seperation calculations
    ss = SkyCoord(x=sxx, y=syy, z=szz, unit='km', representation_type='cartesian')
    ee = SkyCoord(x=exx, y=eyy, z=ezz, unit='km', representation_type='cartesian')
    mm = SkyCoord(x=mxx, y=myy, z=mzz, unit='km', representation_type='cartesian')
    pp = SkyCoord(x=pxx, y=pyy, z=pzz, unit='km', representation_type='cartesian')


### Define Constraints for each Solar system body
    Sun_constraint   = sun_block   * u.deg
    Moon_constraint  = moon_block  * u.deg
    Earth_constraint = earth_block * u.deg #based on worst case orbital alt of 450km should be 63 deg

    #Could be more robust to pull altitude from each time step but might be overkill
    # Pandora_alt = 450
    # Earth_constraint = np.arctan((1.*u.earthRad)/(1.*u.earthRad+Pandora_alt)).to(u.deg)

### Evaluate at each time step whether Pandora is crossing SAA
    # SAA coordinates at altitude of ~500km
    saa_lat_max =   0. * u.deg
    saa_lat_min = -40. * u.deg
    saa_lon_max =  30. * u.deg
    saa_lon_min = -90. * u.deg
    saa_cross   = np.zeros(len(p_lat))
    for i in range(len(p_lat)):
        if (saa_lat_min <= p_lat[i]) and (p_lat[i] <= saa_lat_max) and \
        (saa_lon_min <= p_lon[i]) and (p_lon[i] <= saa_lon_max):
            saa_cross[i] = 1.


### Import Target list
    target_data = pd.read_csv(f"{PACKAGEDIR}/data/target_list.csv", sep=',')
    
    #Cycle through host star targets
    for i in range(len(target_data['Star Simbad Name'])):
        star_name    = target_data['Star Name'][i]
        star_name_sc = target_data['Star Simbad Name'][i]
        star_sc      = SkyCoord.from_name(star_name_sc)
        print('Analyzing constraints for:', star_name)

        #Evaluate at each time step whether target is blocked by each contraints
        Sun_sep   = np.zeros(len(pp))
        Moon_sep  = np.zeros(len(pp))
        Earth_sep = np.zeros(len(pp))
        Sun_req   = np.zeros(len(pp))
        Moon_req  = np.zeros(len(pp))
        Earth_req = np.zeros(len(pp))

        print('Calculating angular seperation requirements')
        pbar = ProgressBar()
        for i in pbar(range(len(pp))):
            Sun_sep[i]   = ss[i].separation(star_sc).deg
            Sun_req[i]   = Sun_sep[i] * u.deg > Sun_constraint
            Moon_sep[i]  = mm[i].separation(star_sc).deg
            Moon_req[i]  = Moon_sep[i] * u.deg > Moon_constraint
            Earth_sep[i] = ee[i].separation(star_sc).deg
            Earth_req[i] = Earth_sep[i] * u.deg > Earth_constraint
        all_req = Sun_req * Moon_req * Earth_req

        #Check if folder exists for planet and if not create new folder for 
        #output products
        save_dir = f"{PACKAGEDIR}/data/targets/" + star_name + '/'
        if os.path.exists(save_dir) != True:
            os.makedirs(save_dir)
        
        #Save results for each star to csv file
        dt_iso_utc = pd.DataFrame(dt_iso_utc, columns=['Time(datetime)'])
        t_mjd_utc  = pd.DataFrame(t_mjd_utc,  columns=['Time(MJD_UTC)'])

        data = np.vstack((saa_cross, Earth_req, Moon_req, Sun_req, \
            all_req, Earth_sep, Moon_sep, Sun_sep))
        data = data.T.reshape(-1,8)
        data = pd.DataFrame(data, columns = ['SAA_Crossing','Earth_Clear','Moon_Clear','Sun_Clear',\
            'Visible','Earth_Sep','Moon_Sep','Sun_Sep'])
        df = pd.concat([dt_iso_utc, t_mjd_utc, data], axis=1)

        output_file_name = 'Visibility for %s.csv' %star_name
        df.to_csv((save_dir + output_file_name), sep=',', index=False)



def transit_timing(planet_name, star_name, output_dir):
    """ Determine primary transits for target(s) during Pandora's science 
    observation lifetime.
        
    Parameters
    ----------
    planet_name:    string
                    name of target planet
    star_name:      string
                    name of target planet's host star
    output_dir:     string
                    directory where to save output csv and plots
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    """

### Read in Visibility Data
    vis_data = pd.read_csv(f"{PACKAGEDIR}/data/targets/" + star_name + '/' + \
                    'Visibility for ' + star_name + '.csv', sep=',')
    t_mjd_utc = vis_data['Time(MJD_UTC)']
    Visible = np.array(vis_data['Visible'])

    # Convert time to datetime
    T_mjd_utc = Time(t_mjd_utc, format='mjd', scale='utc')
    T_iso_utc = Time(T_mjd_utc.iso, format='iso', scale='utc')
    dt_iso_utc = T_iso_utc.to_value('datetime')


### Extract planet specific parameters from target list
    target_data = pd.read_csv(f"{PACKAGEDIR}/data/target_list.csv", sep=',')
    planet_name_sc = target_data.loc[target_data['Planet Name'] == planet_name,
                                'Planet Simbad Name'].iloc[0]
    planet_sc = SkyCoord.from_name(planet_name_sc)

    transit_dur = target_data.loc[target_data['Planet Name'] == planet_name, 
                            'Transit Duration (hrs)'].iloc[0] * u.hour
    period = target_data.loc[target_data['Planet Name'] == planet_name, 
                            'Period (day)'].iloc[0] *u.day
                            
    epoch_BJD_TDB = target_data.loc[target_data['Planet Name'] == planet_name, 
                'Transit Epoch (BJD_TDB-2400000.5)'].iloc[0]+2400000.5
    epoch_JD_UTC  = barycorr.bjd2utc(epoch_BJD_TDB, planet_sc.ra.degree, planet_sc.dec.degree)
    epoch_JD_UTC  = Time(epoch_JD_UTC, format='jd', scale= 'utc')
    epoch_MJD_UTC = Time(epoch_JD_UTC.mjd, format='mjd', scale='utc')


### Calculate transit times
    # Calculate Pre mid-transit time on target
    half_obs_width = 0.75*u.hour + \
        np.maximum((1.*u.hour+transit_dur/2), transit_dur)

    # Determine minimum number of periods between Epoch and Pandora start of observing run
    min_start_epoch = epoch_MJD_UTC-half_obs_width
    min_pers_start = np.ceil((T_mjd_utc[0]-min_start_epoch)/period)

    # Calculate first transit within Pandora lifetime
    first_transit = epoch_MJD_UTC+(min_pers_start*period)

    # Calc transit times
    Mid_transits = np.arange(first_transit, T_mjd_utc[-1], period)
    for i in range(len(Mid_transits)):
        Mid_transits[i] = Mid_transits[i].mjd
    Mid_transits = (np.array(list(Mid_transits[:])))
    Mid_transits = Time(Mid_transits, format='mjd', scale= 'utc')

    Start_transits = Mid_transits-transit_dur/2
    End_transits   = Mid_transits+transit_dur/2

    start_transits = Start_transits.to_value('datetime')
    end_transits   = End_transits.to_value('datetime')

    # Truncate everything after the minutes place
    for i in range(len(start_transits)):
        start_transits[i] = start_transits[i] - timedelta(seconds=start_transits[i].second,
                                        microseconds=start_transits[i].microsecond)
        end_transits[i]   = end_transits[i] - timedelta(seconds=end_transits[i].second,
                                        microseconds=end_transits[i].microsecond)
    all_transits = np.arange(len(start_transits))


### Calculate which transits are visible to Pandora
    dt_vis_times = [] 
    for i in range(len(dt_iso_utc)):
        if Visible[i] == 1.0:
            dt_vis_times.append(dt_iso_utc[i])

    transit_coverage = np.zeros(len(start_transits))
    for i in range(len(start_transits)):
        tran_len = pd.date_range(start_transits[i], end_transits[i], freq='min')
        tran_len = tran_len.to_pydatetime()
        
        tset = set(tran_len)
        tran_vis = tset.intersection(dt_vis_times)
        
        if len(tran_vis)>0:
            transit_coverage[i] = len(tran_vis)/len(tran_len)
        else:
            continue
    
    #Check if folder exists for planet and if not create new folder for 
    #output products
    save_dir = f"{PACKAGEDIR}/data/targets/" + star_name + '/' + planet_name + '/'
    if os.path.exists(save_dir) != True:
        os.makedirs(save_dir)

### Save transit data to Visibility file
    transit_data = np.vstack((all_transits, Start_transits.value, End_transits.value, transit_coverage))
    transit_data = transit_data.T.reshape(-1, 5)
    df2 = pd.DataFrame(transit_data, columns = ['Transits','Transit_Start','Transit_Stop','Transit_Coverage'])

    output_file_name = 'Visibility for ' + planet_name + '.csv'
    df2.to_csv((save_dir + output_file_name), sep=',', index=False)