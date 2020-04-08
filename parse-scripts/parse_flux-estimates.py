#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs March 12 17:22:21 2020

@author: Heather Guy
"""

import sys,os
sys.path.append(os.getcwd())
import numpy as np      
import datetime as dt
import pandas as pd
import os
import sys
#sys.path.append('/Users/heather/Desktop/NC_parse_master')
#sys.path.append('/Users/heather/ICECAPS-ACE/DataParse')
from NC_functions_v1 import *
from flux_functions import *
from netCDF4 import Dataset, date2num
import warnings
warnings.filterwarnings("ignore")

####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/processed/'
#out_loc = '/Users/heather/Desktop/temp_out/'

#start='201909010000'
#stop='201909020000'
#avp=30
#level=2

#Example usage: 
# python parse_flux-estimates.py $in_loc $out_loc $start $stop $avp $level
#############################################################

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    # check length of arguments:
    #if len(args_in) != 5:
    #    # get name of the executable:
    #    self_name = os.path.basename(args_in[0])
    #    # print error and exit:
    #    sys.stderr.write('usage: {0} NC_FILE VAR_NAME\n'.format(self_name))
    #    sys.exit()
    # get values:
    
    in_loc = args_in[1]
    out_loc = args_in[2]
    start = dt.datetime.strptime(args_in[3],'%Y%m%d%H%M')
    stop = dt.datetime.strptime(args_in[4],'%Y%m%d%H%M')
    avp = int(args_in[5]) # Averaging time in minutes.
    level=int(args_in[6]) # level 1 for 2m sonic, level 2 for 15 m 

    # return values:
    return in_loc,out_loc,start,stop,avp,level


def main():
    """
    main function generates netCDF and stores in out_loc
    """
    
    # check / get args:
    try:
        in_loc,out_loc,start,stop,avp,level = get_args(sys.argv)
    except:
        print('Input error')
        sys.exit()
        
    # Global attributes
    meta_f = '../metadata/flux_metadata_level%s_%smin.xlsx'%(level,avp)
    meta = pd.read_excel(meta_f)

    var_f_estimates = '../specific_variables/flux-estimates-level%s.xlsx'%level
    var_estimates = pd.read_excel(var_f_estimates)

    var_f_components = '../specific_variables/flux-components-level%s.xlsx'%level
    var_components = pd.read_excel(var_f_components)

    sf = 10   # sampling frequency (10Hz)


    # Days to loop through
    days = pd.date_range(start,stop,freq='1D')
    m = avp * 60 * 10 # Sample length of interval
    # Make sure m is even
    m=np.fix(m/2)*2
    df = sf/m                    # Frequency intervals
    f = np.arange(df,sf/2+df,df) # Frequency timeseries

    # Loop through each day: 

    for day in days:
        day_str = str(day.date()) 
        print(day_str)
    
        # Set up netcdf files.
    
        f1 = 'ace-flux-%s'%level #instrument name
        f2 = 'summit' #platform name  
        f3 = dt.datetime.strftime(day,'%Y%m%d')
        f5 = "v1" #version number
        f6 = ".nc"
        fn_components = out_loc + 'flux-components/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'flux-components' + chr(95) + '%smin'%avp + chr(95) + f5 + f6
        fn_estimates = out_loc + 'flux-estimates/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'flux-estimates' + chr(95) + '%smin'%avp + chr(95) + f5 + f6
        nc_comp = Dataset(fn_components, "w",  format = "NETCDF4_CLASSIC") 
        nc_est = Dataset(fn_estimates, "w",  format = "NETCDF4_CLASSIC")  
   
        len_index = avp * 60 * sf
        len_time = (24 * 60 ) / avp
        NC_Global_Attributes(nc_comp, meta, day,(day + pd.Timedelta(hours=24) - pd.Timedelta(seconds=0.1)))
        NC_Global_Attributes(nc_est, meta, day,(day + pd.Timedelta(hours=24) - pd.Timedelta(seconds=0.1)))

        NC_Dimensions(nc_comp, len_time,len_index)
        NC_Dimensions(nc_est, len_time,len_index)  
   
        time_list = pd.date_range(day,day+pd.Timedelta(days=1),freq='%smin'%avp)[:-1]

        NC_CommonVariables(nc_comp, time_list, np)
        NC_CommonVariables(nc_est, time_list, np)    
    
        NC_SpecificVariables(nc_comp, var_components, np)
        NC_SpecificVariables(nc_est, var_estimates, np)    
      
        # Get the 3D sonic data

        if os.path.isfile(in_loc+'metek/metek%s_%s'%(level,day_str)):
            m1_orig = pd.read_csv(in_loc+'metek/metek%s_%s'%(level,day_str), index_col=0, parse_dates=[0])
            if m1_orig.empty:
                print('Error: Metek File empty, '+day_str)
                continue
        else:
            print('Error: Metek File empty, '+day_str)
            continue
        
        # Get licor data
    
        if os.path.isfile(in_loc+'LiCOR/licor_%s'%day_str):
            licor = pd.read_csv(in_loc+'LiCOR/licor_%s'%day_str, index_col=0, parse_dates=[0])
        else:
            print('Error: Licor File empty, '+day_str)
            continue
        
        # Get HMP 2m T data

        if os.path.isfile(in_loc+'HMP/HMP1_%s'%day_str):
            HMP1 = pd.read_csv(in_loc+'HMP/HMP1_%s'%day_str, index_col=0, parse_dates=[0])
            # Crop to date, time
            HMP1 = HMP1[day:day+pd.Timedelta(hours=24)] 
        else:
            print('Error: HMP File empty, '+day_str)
            continue
    
        # Get SnD data
        
        if os.path.isfile(in_loc+'SnD/snd_%s'%day_str):
            snd = pd.read_csv(in_loc+'SnD/snd_%s'%day_str, index_col=0, parse_dates=[0])
            # Crop to date, time
            snd = snd[day:day+pd.Timedelta(hours=24)] 
        else:
            print('Error: Snd File empty, '+day_str)
            continue        

        
        # Clean metek data 

        m1 = clean_metek(m1_orig)

        # Implement cross-wind temperature correction

        m1['T_corrected'] = Ts_sidewind_correction(m1['T'],m1['x'],m1['y'],m1['z'])

        # Rotate to average streamline for each averaging period. 

        m_rot = rotate_to_run(m1,avp)
    
        # Process licor data
    
        Ta = HMP1['Ta'] + 273.15   # 2m air temperature,  K
        P = licor['P']             # 2m air pressure form licor, Pa
        m_rot['P']=P
        m_rot['height']=snd['depth_Tcorrected'].reindex(m_rot.index,method='nearest')
        Nconc = licor['H2OD']      # H2O number concentration from licor, mol/m3
        m_rot['Nconc']=Nconc
        m_rot['q'],m_rot['PPw'],m_rot['PPd'],m_rot['mmr'] = licor_mmr(Ta,P,Nconc)  # H2O mass mixing ratio, kg/kg
    
        # Implement Ts humidty correction 
        #m_rot['theta'] =  m_rot['T_corrected'] / (1 + (0.51 * m_rot['q']))  
        #Theta = (Ts + 273.15) / (1 + (0.51 * Q))
     
        # Estimate air density, absolute humidity, Cp and Lv  
    
        try:
            m_rot['rho'] = rho_m(Ta,P,m_rot['q'])        # Air density estimation, kg/m3
            m_rot['rho'] = m_rot['rho'].fillna(np.mean(m_rot['rho'][~np.isnan(m_rot['rho'])])) 
        except:
            m_rot['rho']=np.nan
        

        # Calculate absolute humidity (A, kg K / J)
        # T - water and side wind corrected sonic temperature (K)
        # ppwet: water vapor partial pressure (pa
        # C = 2.16679; % g K J-1
        # A = (C * ppwet / T)/1000 # kg K / j
        # Fill empty values of A with mean value of A. 

        try:
            m_rot['A'] = ((2.16679 * m_rot['PPw']) / m_rot['T_corrected'])/1000.0 # kg K /j
            m_rot['A'] = m_rot['A'].fillna(np.mean(m_rot['A'][~np.isnan(m_rot['A'])]))
        except:
            m_rot['A']=np.nan

        # Calculate heat capacity CP
        # cp = 1000 * (1.005 + (1.82*A))  # Cp adjusted for absolute humidity in J/kg/K

        m_rot['Cp'] = 1000.0 * (1.005 + (1.82*m_rot['A']))  # Cp adjusted for absolute humidity in J/kg/K
        m_rot['Cp'][m_rot['Cp'].isnull()]=1003 # for 250K

        # Calculate latent heat of vaporisation
        # Lv = 3147.5 - (2.372 * Ta) # Lv adjusted for temperature (j/g)
        # This equation comes from Table 2.1, A short course in cloud physics, Rogers & Yau.

        m_rot['Lv'] = (3147.5 -( 2.372 * m_rot['T_corrected'])) * 1000
        m_rot['Lv'][m_rot['Lv'].isnull()]= 2500000


        # Split into runs based on averaging time

        m_g = m_rot.groupby(pd.Grouper(freq='%sMin'%avp))    
        keys = list(m_g.groups.keys())
    
        for i in range(0,len(keys)):       
            k=keys[i]
            try:
                m = m_g.get_group(k)
            except:
                # If this part of the file is missing, skip all. 
                continue

            # Interpolate over data gaps that are smaller than around 6 minutes (60% data for 15 min period)

            m = m.interpolate(limit=3600,limit_direction='both')
        
            # Store 10Hz info for each run
        
            nc_comp.variables['sonic_air_temperature'][i,:] = m['T_corrected'].to_numpy()
            #nc_comp.variables['sonic_temperature_theta'][i,:] = m['theta'].to_numpy()
            nc_comp.variables['eastward_wind_rotated_to_run'][i,:] = m['u'].to_numpy()
            nc_comp.variables['northward_wind_rotated_to_run'][i,:] = m['v'].to_numpy()
            nc_comp.variables['upward_air_velocity_rotated_to_run'][i,:] = m['w'].to_numpy()
            if level==1:
                nc_comp.variables['mole_concentration_of_water_vapor_in_air'][i,:] = m['Nconc'].to_numpy()
                nc_comp.variables['specific_humidity'][i,:] = m['q'].to_numpy()
                nc_comp.variables['humidity_mixing_ratio'][i,:] = m['mmr'].to_numpy()
                nc_comp.variables['water_vapour_partial_pressure_in_air'][i,:] = m['PPw'].to_numpy()
        
                if len(m['Nconc'][m['Nconc'].isnull()]) == 0:
                    h2oprime = detrend(m['Nconc'])
                    nc_comp.variables['h2oprime'][i,:] = h2oprime
                else:
                    h2oprime = np.nan
            
                if len(m['q'][m['q'].isnull()]) == 0:
                    qprime = detrend(m['q'])
                    nc_comp.variables['qprime'][i,:] = qprime            
                else:
                    qprime = np.nan      
            
                #thetaprime = detrend(m['theta'])
                #nc_comp.variables['thetaprime'][i,:] = thetaprime 
        
            tsprime = detrend(m['T_corrected'])
            nc_comp.variables['tsprime'][i,:] = tsprime 
        
            uprime = detrend(m['u'])
            nc_comp.variables['uprime'][i,:] = uprime 
        
            uprimeuprime = uprime * uprime
            nc_comp.variables['uprimeuprime'][i,:] = uprimeuprime         
        
            vprime = detrend(m['v'])
            nc_comp.variables['vprime'][i,:] = vprime         
        
            vprimevprime = vprime * vprime
            nc_comp.variables['vprimevprime'][i,:] = vprimevprime          
        
            wprime = detrend(m['w'])
            nc_comp.variables['wprime'][i,:] = wprime          
        
            if level ==1:
                wprimeh2oprime = wprime * h2oprime
                nc_comp.variables['wprimeh2oprime'][i,:] = wprimeh2oprime  
        
                wprimeqprime = wprime * qprime
                nc_comp.variables['wprimeqprime'][i,:] = wprimeqprime 
        
            #  wprimethetaprime = wprime * thetaprime
            #nc_comp.variables['wprimethetaprime'][i,:] = wprimethetaprime
        
            wprimetsprime = wprime * tsprime
            nc_comp.variables['wprimetsprime'][i,:] = wprimetsprime
        
            wprimeuprime = wprime * uprime
            nc_comp.variables['wprimeuprime'][i,:] = wprimeuprime
        
            wprimevprime = wprime * vprime
            nc_comp.variables['wprimevprime'][i,:] = wprimevprime
        
            wprimewprime = wprime * wprime
            nc_comp.variables['wprimewprime'][i,:] = wprimewprime
        
        
            # Store single parameters for each run
 
            air_pressure = (m['P'].mean())/100 # hPa
            nc_comp.variables['air_pressure'][i,:] = (m['P']/100).to_numpy()
        
            nc_comp.variables['number_of_samples_in_run'][i] = len(m)
            nc_est.variables['number_of_samples_in_run'][i] = len(m)        
        
            height_above_surface = m['height'].mean()
            nc_comp.variables['height_above_surface'][i] = height_above_surface
            nc_est.variables['height_above_surface'][i] = height_above_surface
        
            nc_comp.variables['run_length'][i] = (m.index[-1] - m.index[0]).seconds
            nc_est.variables['run_length'][i] = (m.index[-1] - m.index[0]).seconds
        
            nc_comp.variables['start_of_run'][i] = np.float64(date2num(m.index[0],units='seconds since 1970-01-01 00:00:00 UTC'))
            nc_est.variables['start_of_run'][i] = np.float64(date2num(m.index[0],units='seconds since 1970-01-01 00:00:00 UTC'))

        
            sigma_w = np.std(wprime)
            nc_comp.variables['standard_deviation_upward_air_velocity'][i] = sigma_w
        
            if level==1:
                h2obar = m['Nconc'].mean()
                nc_comp.variables['h2obar'][i] = h2obar
        
                qbar = m['q'].mean()
                nc_comp.variables['qbar'][i] = qbar       
        
            #thetabar = m['theta'].mean()
            #nc_comp.variables['thetabar'][i] = thetabar        
        
            tsbar = m['T_corrected'].mean()
            nc_comp.variables['tsbar'][i] = tsbar      
        
            ubar = m['u'].mean()
            nc_comp.variables['ubar'][i] = ubar       
        
            vbar = m['v'].mean()
            nc_comp.variables['vbar'][i] = vbar      
        
            wbar = m['w'].mean()
            nc_comp.variables['wbar'][i] = wbar       
        
            uprimeuprimebar = np.mean(uprimeuprime)
            nc_comp.variables['uprimeuprimebar'][i] = uprimeuprimebar 
        
            vprimevprimebar = np.mean(vprimevprime)
            nc_comp.variables['vprimevprimebar'][i] = vprimevprimebar
        
            wprimewprimebar = np.mean(wprimewprime)
            nc_comp.variables['wprimewprimebar'][i] = wprimewprimebar
        
            if level==1:
                wprimeh2oprimebar = np.mean(wprimeh2oprime)
                nc_comp.variables['wprimeh2oprimebar'][i] = wprimeh2oprimebar 
        
                wprimeqprimebar = np.mean(wprimeqprime)
                nc_comp.variables['wprimeqprimebar'][i] = wprimeqprimebar 
        
            #wprimethetaprimebar = np.mean(wprimethetaprime)
            #nc_comp.variables['wprimethetaprimebar'][i] = wprimethetaprimebar
        
            wprimetsprimebar = np.mean(wprimetsprime)
            nc_comp.variables['wprimetsprimebar'][i] = wprimetsprimebar 
        
            wprimeuprimebar = np.mean(wprimeuprime)
            nc_comp.variables['wprimeuprimebar'][i] = wprimeuprimebar
        
            wprimevprimebar = np.mean(wprimevprime)
            nc_comp.variables['wprimevprimebar'][i] = wprimevprimebar 
        
        
            # Skew
        
            skew_ts = skew(m['T'])
            nc_comp.variables['skew_sonic_air_temperature'][i] = skew_ts
        
            skew_u = skew(m['u'])
            nc_comp.variables['skew_eastward_wind'][i] = skew_u
        
            skew_v = skew(m['v'])
            nc_comp.variables['skew_northward_wind'][i] = skew_u        
        
            skew_w = skew(m['w'])
            nc_comp.variables['skew_upward_air_velocity'][i] = skew_u        
        
        
            # Kurtosis
        
            kurtosis_ts = kurtosis(m['T'])
            nc_comp.variables['kurtosis_sonic_air_temperature'][i] = kurtosis_ts      
        
            kurtosis_u = kurtosis(m['w'])
            nc_comp.variables['kurtosis_eastward_wind'][i] = kurtosis_u
        
            kurtosis_v = kurtosis(m['v'])
            nc_comp.variables['kurtosis_northward_wind'][i] = kurtosis_u 
        
            kurtosis_w = kurtosis(m['w'])
            nc_comp.variables['kurtosis_upward_air_velocity'][i] = kurtosis_u   
        
        
            # Stationarity testing
        
            sst_wts,Cwt, rol_cov_wt = stationarity(m['w'],m['T'])
            nc_comp.variables['sst_wts'][i] = sst_wts
        
            sst_wu, Cwu, rol_cov_wu= stationarity(m['w'],m['u'])
            nc_comp.variables['sst_wu'][i] = sst_wu
        
            sst_wv, Cwv, rol_cov_wv= stationarity(m['w'],m['v'])
            nc_comp.variables['sst_wv'][i] = sst_wv

        
            friction_velocity = (Cwu**2 + Cwv**2)**(1/4)
            nc_comp.variables['friction_velocity'][i] = friction_velocity
        
            obukhov_length = (-np.abs(friction_velocity**3) * np.mean(m['T_corrected'])) / (0.4*9.81*Cwt)
            nc_comp.variables['obukhov_length'][i] = obukhov_length
        
            stability_parameter = height_above_surface / obukhov_length
            nc_comp.variables['stability_parameter'][i] = stability_parameter
        
            # Integral scale test (for turbulence development)
            # theoretical value of sigma_w/ustar - parametrisation after Foken/CarboEurope
        
            if np.abs(stability_parameter) > 1:
                sigma_uw_theory = 2
            elif np.abs(stability_parameter) < 0.032:
                sigma_uw_theory = 1.3
            else:
                sigma_uw_theory = 2 * np.abs(stability_parameter)**0.125  

            itc_w = 100 * ((sigma_uw_theory - (sigma_w/friction_velocity))/sigma_uw_theory) 
            nc_comp.variables['integral_turbulent_characteristic_upward_air_velocity'][i] = itc_w
        
        
            # QC flags
            #qc_flag_itc_class
        
            nc_comp.variables['qc_flag_kurtosis_ts'][i] = kurt_flag(kurtosis_ts)
            nc_est.variables['qc_flag_kurtosis_ts'][i] = kurt_flag(kurtosis_ts)
        
            nc_comp.variables['qc_flag_kurtosis_u'][i] = kurt_flag(kurtosis_u)
            nc_est.variables['qc_flag_kurtosis_u'][i] = kurt_flag(kurtosis_u)
        
            nc_comp.variables['qc_flag_kurtosis_v'][i] = kurt_flag(kurtosis_v)
            nc_est.variables['qc_flag_kurtosis_v'][i] = kurt_flag(kurtosis_v)
        
            nc_comp.variables['qc_flag_kurtosis_w'][i] = kurt_flag(kurtosis_w)
            nc_est.variables['qc_flag_kurtosis_w'][i] = kurt_flag(kurtosis_w)
        
            nc_comp.variables['qc_flag_quality_wts'][i] = flux_devel_test(itc_w, sst_wts)
            nc_est.variables['qc_flag_quality_wts'][i] = flux_devel_test(itc_w, sst_wts)
        
            nc_comp.variables['qc_flag_quality_wu'][i]  = flux_devel_test(itc_w, sst_wu)
            nc_est.variables['qc_flag_quality_wu'][i]  = flux_devel_test(itc_w, sst_wu)
        
            nc_comp.variables['qc_flag_quality_wv'][i] = flux_devel_test(itc_w, sst_wv)
            nc_est.variables['qc_flag_quality_wv'][i] = flux_devel_test(itc_w, sst_wv)
        
            nc_comp.variables['qc_flag_skew_ts'][i] = skew_flag(skew_ts)
            nc_est.variables['qc_flag_skew_ts'][i] = skew_flag(skew_ts)
        
            nc_comp.variables['qc_flag_skew_u'][i] = skew_flag(skew_u)
            nc_est.variables['qc_flag_skew_u'][i] = skew_flag(skew_u)
        
            nc_comp.variables['qc_flag_skew_v'][i] = skew_flag(skew_v)
            nc_est.variables['qc_flag_skew_v'][i] = skew_flag(skew_v)
        
            nc_comp.variables['qc_flag_skew_w'][i] = skew_flag(skew_w)
            nc_est.variables['qc_flag_skew_w'][i] = skew_flag(skew_w)
        
            #qc_flag_sstclass_wts
            #qc_flag_sstclass_wu
            #qc_flag_sstclass_wv       
        
            # Calculate flux estimates
        
            rho = m['rho'].mean()
            if np.isnan(rho):
                rho = 1.4224	 # Density of dry air at -25 C. kg m-3
            cp = m['Cp'].mean()
            lv = m['Lv'].mean()
        
            # Note using ts not theta. Too much missing licor data to use theta. 
            # Also under such cold conditions like at summit, the humidity contribtuion
            # to the sonic temeprature flux is negligible so that sonic temperature
            # flux is a good estimate of sensible heat flux (Andreas et al 2005)
         
            upward_sensible_heat_flux = rho * cp * wprimetsprimebar 
            nc_est.variables['upward_sensible_heat_flux_in_air'][i] = upward_sensible_heat_flux
        
            if level==1:
                upward_latent_heat_flux = rho * lv * wprimeqprimebar
                nc_est.variables['upward_latent_heat_flux_in_air'][i] = upward_latent_heat_flux        
                nc_est.variables['bowen_ratio'][i]  = upward_sensible_heat_flux / upward_latent_heat_flux
                nc_est.variables['kinematic_humidity_flux'][i]  = wprimeqprimebar
        
            nc_est.variables['buoyancy_flux'][i]  = rho * cp * wprimetsprimebar
            nc_est.variables['kinematic_sonic_temperature_flux'][i]  = wprimetsprimebar
            #nc_est.variables['kinematic_heat_flux'][i]  = wprimethetaprimebar
            nc_est.variables['momentum_flux_u'][i]  = - rho * wprimeuprimebar
            nc_est.variables['momentum_flux_v'][i]  = - rho * wprimevprimebar
         
        # Close netcdf files
        
        nc_est.close()
        nc_comp.close()
    
            
if __name__ == '__main__':
    main()          












