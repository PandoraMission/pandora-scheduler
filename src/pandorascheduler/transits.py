"""Functions to calculate transits times from a target list"""
import os
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.time import Time
from datetime import timedelta
# from astropy.time import TimeDelta
# from progressbar import ProgressBar
import logging
# import barycorr
from tqdm import tqdm
# from . import barycorr
import json

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# print(PACKAGEDIR)
# from . import PACKAGEDIR

def star_vis(sun_block:float, moon_block:float, earth_block:float, 
               obs_start:str, obs_stop:str, 
               gmat_file:str,# = 'GMAT_pandora_450_20230713.csv',
               obs_name:str,# = 'Pandora_450km_20230713',
               save_pth:str,# = f'{PACKAGEDIR}/data/targets/',
               targ_list:str):# = f'{PACKAGEDIR}/data/target_list.csv'):
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
    save_pth:       String
                    Path to the directory to save the visibility results
                    default: f'{PACKAGEDIR}/data/targets/'
    targ_list:      String
                    Path to the target list to consider for visibility calculations
                    default: f'{PACKAGEDIR}/data/target_list.csv'
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
    logging.info('Importing GMAT data')
    # gmat_data = pd.read_csv(f'{PACKAGEDIR}/data/{gmat_file}', sep='\t')
    if gmat_file == 'GMAT_pandora_450_20230713.csv':
        gmat_data = pd.read_csv(f'{PACKAGEDIR}/data/{gmat_file}', sep='\t+', engine = 'python')
        keyword_gmat = 'Earth.UTCModJulian'
        obs_name = 'Pandora'
    elif gmat_file == 'GMAT_pandora_600_20240512.txt':
        gmat_data = pd.read_fwf(f'{PACKAGEDIR}/data/{gmat_file}')
        keyword_gmat = 'Earth.UTCModJulian'
    elif gmat_file == 'GMAT_pandora_600_20250706_withdrag.txt':
        gmat_data = pd.read_csv(f'{PACKAGEDIR}/data/{gmat_file}', sep=',', engine = 'python')
        keyword_gmat = 'Earth.A1ModJulian'
        obs_name = 'Pandora'
    elif gmat_file == 'GMAT_pandora_600_20250708_withdrag.txt':
        gmat_data = pd.read_csv(f'{PACKAGEDIR}/data/{gmat_file}', sep=',', engine = 'python')
        keyword_gmat = 'Earth.UTCModJulian'
        obs_name = 'Pandora'
    elif gmat_file == 'GMAT_pandora_600_20250820_withdrag.txt':
        gmat_data = pd.read_csv(f'{PACKAGEDIR}/data/{gmat_file}', sep=',', engine = 'python')
        keyword_gmat = 'Earth.UTCModJulian'
        obs_name = 'Pandora'


    # Trim dataframe to slightly larger than date range of 
    # Pandora science lifetime defined as obs_start and obs_stop
    gmat_data = gmat_data[(gmat_data[keyword_gmat]>=(t_jd_utc[0]-2430000.0)-0.0007) & 
            (gmat_data[keyword_gmat]<=(t_jd_utc[-1]-2430000.0)+0.0007)]

    ### Convert GMAT times into standard MJD_UTC
    # Note: GMAT uses different offset for it's MJD (uses 2,430,000.0 rather than 2,400,000.5)
    gmat_mjd_utc =np.array(gmat_data[keyword_gmat]) + 2430000.0 - 2400000.5

    ### Extract GMAT positions in MJ2000 Earth Centric (EC) cartesian coordinates
    # Earth Centric (EC) coordinates from GMAT
    earth_vectors_gmat = np.asarray(gmat_data[['Earth.EarthMJ2000Eq.X', 'Earth.EarthMJ2000Eq.Y', 'Earth.EarthMJ2000Eq.Z']])
    pandora_vectors_gmat = np.asarray(gmat_data[[f'{obs_name}.EarthMJ2000Eq.X', f'{obs_name}.EarthMJ2000Eq.Y', f'{obs_name}.EarthMJ2000Eq.Z']])
    sun_vectors_gmat = np.asarray(gmat_data[['Sun.EarthMJ2000Eq.X', 'Sun.EarthMJ2000Eq.Y', 'Sun.EarthMJ2000Eq.Z']])
    moon_vectors_gmat = np.asarray(gmat_data[['Luna.EarthMJ2000Eq.X', 'Luna.EarthMJ2000Eq.Y', 'Luna.EarthMJ2000Eq.Z']])


    ### Interpolate all positions from GMAT to map to 1 minute intervals
    earth_vectors_ec = np.empty((len(t_mjd_utc), 3))
    pandora_vectors_ec = np.empty((len(t_mjd_utc), 3))
    sun_vectors_ec = np.empty((len(t_mjd_utc), 3))
    moon_vectors_ec = np.empty((len(t_mjd_utc), 3))
    for i in range(3):
        earth_vectors_ec[:, i] = np.interp(t_mjd_utc, gmat_mjd_utc, earth_vectors_gmat[:, i])
        pandora_vectors_ec[:, i] = np.interp(t_mjd_utc, gmat_mjd_utc, pandora_vectors_gmat[:, i])
        sun_vectors_ec[:, i] = np.interp(t_mjd_utc, gmat_mjd_utc, sun_vectors_gmat[:, i])
        moon_vectors_ec[:, i] = np.interp(t_mjd_utc, gmat_mjd_utc, moon_vectors_gmat[:, i])

    ### Coordinate shift to Pandora Centric (PC) reference frame
    earth_vectors_pc = earth_vectors_ec-pandora_vectors_ec
    sun_vectors_pc = sun_vectors_ec-pandora_vectors_ec
    moon_vectors_pc = moon_vectors_ec-pandora_vectors_ec

    ### Create SkyCoord for angular seperation calculations
    earth_vectors_pc = SkyCoord(earth_vectors_pc, unit='km', representation_type='cartesian') 
    sun_vectors_pc = SkyCoord(sun_vectors_pc, unit='km', representation_type='cartesian')
    moon_vectors_pc = SkyCoord(moon_vectors_pc, unit='km', representation_type='cartesian')  

    ### Define Constraints for each Solar system body
    Sun_constraint   = sun_block   * u.deg
    Moon_constraint  = moon_block  * u.deg
    Earth_constraint = earth_block * u.deg #based on worst case orbital alt of 450km should be 63 deg

    #Could be more robust to pull altitude from each time step but might be overkill
    # Pandora_alt = 450
    # Earth_constraint = np.arctan((1.*u.earthRad)/(1.*u.earthRad+Pandora_alt)).to(u.deg)

    ### Pandora Latitude & Longitude
    gmat_lat = np.array(gmat_data[f'{obs_name}.Earth.Latitude'])
    gmat_lon = np.array(gmat_data[f'{obs_name}.Earth.Longitude'])
    
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
    target_data = pd.read_csv(targ_list, sep=',')
    
    #Cycle through host star targets
    for i in tqdm(range(len(target_data['Star Simbad Name']))):
        
        star_name = target_data['Star Name'][i]
        star_name_sc = target_data['Star Simbad Name'][i]
        # if star_name_sc.startswith("DR3"):
        if "occultation-standard" in targ_list:
            star_name_sc = star_name_sc.replace('G', 'Gaia DR3 ')

        star_sc = SkyCoord.from_name(star_name_sc)
        logging.info('Analyzing constraints for:', star_name)
        
        #Check if folder exists for target and if not create new folder for 
        #output products
        save_dir = save_pth + star_name + '/'
        output_file_name = 'Visibility for %s.csv' %star_name
        save_name = save_dir + output_file_name
        if os.path.exists(save_dir) != True:
            os.makedirs(save_dir)
            
        #Also check if vis calc have already been done for the target
        else:
            if os.path.exists(save_name) == True:
                #Skip this round of the loop if vis calcs have been done
                continue

        #Evaluate at each time step whether target is blocked by each contraints
        Sun_sep   = sun_vectors_pc.separation(star_sc).deg
        Sun_req   = Sun_sep*u.deg > Sun_constraint
        Moon_sep  = moon_vectors_pc.separation(star_sc).deg
        Moon_req  = Moon_sep*u.deg > Moon_constraint
        Earth_sep = earth_vectors_pc.separation(star_sc).deg       
        Earth_req = Earth_sep*u.deg > Earth_constraint
        all_req = Sun_req * Moon_req * Earth_req
        
        
        #Save results for each star to csv file
        data = np.vstack((t_mjd_utc, saa_cross, all_req, Earth_sep, Moon_sep, Sun_sep))
        data = data.T.reshape(-1,6)
        vis_df = pd.DataFrame(data, columns = ['Time(MJD_UTC)', 'SAA_Crossing', \
            'Visible','Earth_Sep','Moon_Sep','Sun_Sep'])
        
        # def custom_float_format(df):
        #     formatters = {}
        #     for col in df.columns:
        #         if col == 'Time(MJD_UTC)':
        #             df[col] = lambda x: f'{x:.6f}'
        #         elif col == 'SAA_Crossing' or col == 'Visible':
        #             df[col] = lambda x: f'{x:.1f}'
        #         else:
        #             df[col] = lambda x: f'{x:.3f}'
        #     return df

        if "exoplanet" not in targ_list:
            vis_df['Time(MJD_UTC)'] = np.round(vis_df['Time(MJD_UTC)'], 6)
        vis_df['SAA_Crossing'] = np.round(vis_df['SAA_Crossing'], 1)
        vis_df['Visible'] = np.round(vis_df['Visible'], 1)
        vis_df['Earth_Sep'] = np.round(vis_df['Earth_Sep'], 3)
        vis_df['Moon_Sep'] = np.round(vis_df['Moon_Sep'], 3)
        vis_df['Sun_Sep'] = np.round(vis_df['Sun_Sep'], 3)

        # formatters = custom_float_format(vis_df)
        # vis_df.to_csv(save_name, sep=',', index=False, float_format=custom_float_format)
        ## vis_df.to_csv(save_name, sep=',', index=False, float_format='%.nf')
        vis_df.to_csv((save_name), sep=',', index=False)


def transit_timing(target_list:str, planet_name:str, star_name:str):
    """ Determine primary transits for target(s) during Pandora's science 
    observation lifetime.
        
    Parameters
    ----------
    target_list:    string
                    name of csv file with list of targets
    planet_name:    string
                    name of target planet
    star_name:      string
                    name of target planet's host star
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    """

    ### Read in Visibility Data
    vis_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/' + star_name + '/' + \
                    'Visibility for ' + star_name + '.csv', sep=',')
    t_mjd_utc = vis_data['Time(MJD_UTC)']
    Visible = np.array(vis_data['Visible'])

    from astropy.time import Time
    # Convert time to datetime
    T_mjd_utc = Time(t_mjd_utc, format='mjd', scale='utc')
    T_iso_utc = Time(T_mjd_utc.iso, format='iso', scale='utc')
    dt_iso_utc = T_iso_utc.to_value('datetime')

    ### Extract planet specific parameters from target list
    target_data = pd.read_csv(f'{PACKAGEDIR}/data/' + target_list, sep=',')
    planet_name_sc = target_data.loc[target_data['Planet Name'] == planet_name,
                                'Planet Simbad Name'].iloc[0]
    planet_sc = SkyCoord.from_name(star_name)

    transit_dur = target_data.loc[target_data['Planet Name'] == planet_name, 
                            'Transit Duration (hrs)'].iloc[0] * u.hour
    period = target_data.loc[target_data['Planet Name'] == planet_name, 
                            'Period (days)'].iloc[0] *u.day
                            
    epoch_BJD_TDB = target_data.loc[target_data['Planet Name'] == planet_name, 
                'Transit Epoch (BJD_TDB-2400000.5)'].iloc[0]+2400000.5

    #epoch_JD_UTC  = barycorr.bjd2utc(epoch_BJD_TDB, planet_sc.ra.degree, planet_sc.dec.degree) 
    # VK note: Jason Eastman's website (barycorr.bjd2utc) is broken as of 2024-01-05. 
    # Instead, use the implementation from https://github.com/shbhuk/barycorrpy/issues/54, i.e.: 

    observer_location = EarthLocation(lat = 0.*u.deg, lon = 0.*u.deg, height=600*u.km)
    bjd_tdb = Time(epoch_BJD_TDB, format='jd', scale='tdb', location=observer_location)

    # From BJD_TDB to JD_TDB
    ltt_bary = bjd_tdb.light_travel_time(planet_sc, kind='barycentric', location=observer_location)
    jd_tdb = bjd_tdb - ltt_bary    
    
    # From JD_TDB to JD_UTC
    epoch_JD_UTC  = Time(jd_tdb, format='jd', scale= 'utc')
    epoch_MJD_UTC = Time(jd_tdb.mjd, format='mjd', scale='utc')

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
    Mid_transits = np.arange(first_transit, T_mjd_utc[-1], period)# --> THIS GIVES A WARNING, TRY WITH LINES BELOW
    # from astropy.time import Time, TimeDelta
    # period_delta = TimeDelta(period)
    # num_periods = int((T_mjd_utc[-1] - first_transit) / period)
    # Mid_transits = first_transit + period_delta * np.arange(num_periods)

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
    save_dir = f'{PACKAGEDIR}/data/targets/' + star_name + '/' + planet_name + '/'
    if os.path.exists(save_dir) != True:
        os.makedirs(save_dir)

    ### Save transit data to Visibility file
    transit_data = np.vstack((all_transits, Start_transits.value, End_transits.value, transit_coverage))
    transit_data = transit_data.T.reshape(-1, 4)
    transit_df = pd.DataFrame(transit_data, columns = ['Transits','Transit_Start','Transit_Stop','Transit_Coverage'])

    output_file_name = 'Visibility for ' + planet_name + '.csv'
    transit_df.to_csv((save_dir + output_file_name), sep=',', index=False)


def Transit_overlap(target_list:str, partner_list:str, star_name:str):
    """ Determine if there is overlap between target planet's transit
    and a companion planet
        
    Parameters
    ----------
    target_list:    string
                    name of target list file
    partner_list:   string
                    name of partner list file
    star_name:      string
                    name of target planet's host star
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    """
    target_data   = pd.read_csv(f'{PACKAGEDIR}/data/' + target_list, sep=',')
    partners_data = pd.read_csv(f'{PACKAGEDIR}/data/' + partner_list, sep=',')

    ### Compile all planets in system to be considered
    total_planets = []
    for i in range(len(target_data)):
        if target_data['Star Name'][i] == star_name:
            total_planets.append(target_data['Planet Name'][i])
    for i in range(len(partners_data)):
        if partners_data['Star Name'][i] == star_name:
            total_planets.append(partners_data['Planet Name'][i])
    total_planets.sort()

    if len(total_planets) > 1:
        total_planets = pd.Series(total_planets).unique().tolist()

    if len(total_planets) < 2:
        logging.info('Only one planet in list around ', star_name)
    else:
        logging.info('Collecting Transit Times for all planets in list around ', star_name)
        for j in range(len(total_planets)):
            planet_name = total_planets[j]
            planet_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/' + star_name + '/' + planet_name + '/' + 
                                        'Visibility for ' + planet_name + '.csv')
            Start_transits = Time(planet_data['Transit_Start'].values, format='mjd', scale='utc')
            End_transits   = Time(planet_data['Transit_Stop'].values, format='mjd', scale='utc')

            start_transits = Start_transits.to_value('datetime')
            end_transits   = End_transits.to_value('datetime')

            # Truncate everything after the minutes place
            for k in range(len(Start_transits)):
                start_transits[k] = start_transits[k] - timedelta(seconds=start_transits[k].second,
                                                microseconds=start_transits[k].microsecond)
                end_transits[k]   = end_transits[k] - timedelta(seconds=end_transits[k].second,
                                                microseconds=end_transits[k].microsecond)

            if j == 0:
                All_start_transits = pd.DataFrame(start_transits, columns=[planet_name])
                All_end_transits   = pd.DataFrame(end_transits, columns=[planet_name])

            else:
                start_transits = pd.DataFrame(start_transits, columns=[planet_name])
                end_transits   = pd.DataFrame(end_transits, columns=[planet_name])
                
                All_start_transits = pd.concat([All_start_transits, start_transits], axis=1) 
                All_end_transits   = pd.concat([All_end_transits, end_transits], axis=1)  

    ###     Find overlaps with other planets in the system
        for j in range(len(total_planets)):
            planet_name = total_planets[j]
            planet_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/' + star_name + '/' + planet_name + '/' + 
                                        'Visibility for ' + planet_name + '.csv')

            overlap = pd.DataFrame(0., index=np.arange(len(All_start_transits[planet_name].dropna())), columns=['Transit_Overlap'])

            for m in range(len(All_start_transits.columns)):
                if All_start_transits.columns[m] == planet_name:
                    logging.info('Analyzing:', planet_name)
                    continue
                else:
                    planet_partner = All_start_transits.columns[m]
                    logging.info('Checking ', planet_name, ' against: ', planet_partner)
                
                for n in range(len(All_start_transits[planet_name].dropna())):
                    
                    for p in range(len(All_start_transits[planet_partner].dropna())):
                        if (All_start_transits[planet_partner][p] < All_start_transits[planet_name][n] and
                            All_end_transits[planet_partner][p] < All_start_transits[planet_name][n]) or \
                            (All_start_transits[planet_partner][p] > All_end_transits[planet_name][n] and
                            All_end_transits[planet_partner][p] > All_end_transits[planet_name][n]):
                            continue
                        else:
                            partner_rng = pd.date_range(All_start_transits[planet_partner][p],
                                                All_end_transits[planet_partner][p], freq='min')
                            partner_rng = partner_rng.to_pydatetime()

                            transit_rng = pd.date_range(All_start_transits[planet_name][n], 
                                                All_end_transits[planet_name][n], freq='min')
                            transit_rng = transit_rng.to_pydatetime()

                            pset = set(partner_rng)
                            tset = set(transit_rng)
                            overlap_times = pset.intersection(tset)
                            transit_overlap = len(overlap_times)/len(transit_rng)
                            # overlap['Transit_Overlap'][n] = transit_overlap
                            overlap.loc[n, 'Transit_Overlap'] = transit_overlap
                       
    ###         Update pandas dataframe and save csv
            if 'Transit_Overlap' in planet_data:
                planet_data.update(overlap)
            else:
                planet_data = pd.concat([planet_data, overlap], axis=1)
            save_dir   = f'{PACKAGEDIR}/data/targets/' + star_name + '/' + planet_name + '/'
            save_fname = 'Visibility for ' + planet_name + '.csv'
            planet_data.to_csv((save_dir + save_fname), sep=',', index=False)


def SAA_overlap(planet_name:str, star_name:str):
    """ Determine if there is overlap between target planet's transit
    and South Atlantic Anomaly crossing
        
    Parameters
    ----------
    planet_name:    string
                    name of target planet
    star_name:      string
                    name of target planet's host star
                                    
    Returns
    -------
    csv file
        file containing target's transits during Pandora's lifetime
    """
    ### Read in SAA crossing data
    star_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/' + star_name + '/' + 
                            'Visibility for ' + star_name + '.csv')
    SAA_time_mjd = Time(star_data['Time(MJD_UTC)'].values, format='mjd', scale='utc')
    SAA_time_iso = Time(SAA_time_mjd.iso, format='iso', scale='utc')
    SAA_time_dt  = SAA_time_iso.to_value('datetime')
    SAA_data = star_data['SAA_Crossing']


    ### Read in planet transit data
    planet_data = pd.read_csv(f'{PACKAGEDIR}/data/targets/' + star_name + '/' + planet_name + '/' + 
                            'Visibility for ' + planet_name + '.csv')
    start_transits = Time(planet_data['Transit_Start'].values, format='mjd', scale='utc')
    end_transits   = Time(planet_data['Transit_Stop'].values, format='mjd', scale='utc')

    start_transits = start_transits.to_value('datetime')
    end_transits   = end_transits.to_value('datetime')


    ### Truncate everything after the minutes place
    for i in range(len(start_transits)):
        start_transits[i] = start_transits[i] - timedelta(seconds=start_transits[i].second,
                                        microseconds=start_transits[i].microsecond)
        end_transits[i]   = end_transits[i] - timedelta(seconds=end_transits[i].second,
                                        microseconds=end_transits[i].microsecond)


    ### Determine SAA overlap for each transit
    saa_overlap = pd.DataFrame(0., index=np.arange(len(planet_data)), columns=['SAA_Overlap'])
    
    for i in range(len(start_transits)):
        transit_rng = pd.date_range(start_transits[i], 
                            end_transits[i], freq='min')
        transit_rng = transit_rng.to_pydatetime()
        saa_rng = []
        for j in range(len(SAA_time_dt)):
            if (transit_rng[0] <= SAA_time_dt[j]) and \
            (SAA_time_dt[j] <= transit_rng[-1]):
                if SAA_data[j] == 1:
                    saa_rng.append(SAA_time_dt[j])
        sset = set(saa_rng)
        tset = set(transit_rng)
        overlap_times = sset.intersection(tset)
        transit_overlap = len(overlap_times)/len(transit_rng)
        # saa_overlap['SAA_Overlap'][i] = transit_overlap
        # saa_overlap.loc[i, "SAA Overlap"] = transit_overlap
        saa_overlap.loc[i] = transit_overlap

    ### Update pandas dataframe and save csv
    if 'SAA_Overlap' in planet_data:
        planet_data.update(saa_overlap)
    else:
        planet_data = pd.concat([planet_data, saa_overlap], axis=1)
    save_dir   = f'{PACKAGEDIR}/data/targets/' + star_name + '/' + planet_name + '/'
    save_fname = 'Visibility for ' + planet_name + '.csv'
    planet_data.to_csv((save_dir + save_fname), sep=',', index=False)