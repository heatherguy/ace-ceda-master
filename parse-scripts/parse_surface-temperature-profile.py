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
import sys
#sys.path.append('/Users/heather/Desktop/NC_parse_master')
#sys.path.append('/Users/heather/ICECAPS-ACE/DataParse')
from NC_functions_v1 import *
from fluxtower_parse import *
from netCDF4 import Dataset, num2date


####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/processed/'
#out_loc = '/Users/heather/Desktop/temp_out/'

# Months
#months=[12]
#years = [2019]

# sample interval
#avp = 1 # minutes

#Example usage: 
# python parse_surface-temperature-profile.py in_loc out_loc months years avp
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
    months = map(int, args_in[3].strip('[]').split(','))
    years = map(int, args_in[4].strip('[]').split(','))
    avp = int(args_in[5])

    # return values:
    return in_loc,out_loc,months,years,avp

def main():
    """
    main function generates netCDF and stores in out_loc
    """

    # check / get args:
    try:
        in_loc,out_loc,months,years,avp = get_args(sys.argv)
    except:
        print('Input error')
        sys.exit()
    
    # Global attributes
    meta_f = '../metadata/trh_profile_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '../specific_variables/surface-temperature-profile.xlsx'
    var = pd.read_excel(var_f)


    # Loop through each month: 
    for year in years: 
        for month in months:
            start = dt.datetime(year, month, 1, 0, 0)
            if month==12: 
                stop = dt.datetime(year+1,1,1,0,0)
            else:
                stop = dt.datetime(year,month+1,1,0,0)

            print(start.strftime(format='%Y%m'))
    
            # Set up netcdf files.
            
            f1 = 'ace-tower'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc + 'surface-temperature-profile/' +f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'surface-temperature-profile' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / avp

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='%smin'%avp)[:]
            NC_Dimensions(nc, len(time_list), index=4)  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get HMP155 data
        
            print('Extracting HMP1')
            HMP1 = get_hmp(start,stop,in_loc+'HMP/','HMP1')
            print('Extracting HMP2')
            HMP2 = get_hmp(start,stop,in_loc+'HMP/','HMP2')
            print('Extracting HMP3')
            HMP3 = get_hmp(start,stop,in_loc+'HMP/','HMP3')
            print('Extracting HMP4')
            HMP4 = get_hmp(start,stop,in_loc+'HMP/','HMP4')
            print('Post-processing')
            
            # Minutely averages
        
            hmp1 = HMP1.resample(rule = '1min', how='mean')
            hmp1 = HMP1.reindex(time_list, method='nearest',limit=2)

            hmp2 = HMP2.resample(rule = '1min', how='mean')
            hmp2 = HMP2.reindex(time_list, method='nearest',limit=2)

            hmp3 = HMP3.resample(rule = '1min', how='mean')
            hmp3 = HMP3.reindex(time_list, method='nearest',limit=2)

            hmp4 = HMP4.resample(rule = '1min', how='mean')
            hmp4 = HMP4.reindex(time_list, method='nearest',limit=2)        
 
            # Convert T to Kelvin
        
            t1 = hmp1['Ta'] + 273.15
            t2 = hmp2['Ta'] + 273.15
            t3 = hmp3['Ta'] + 273.15
            t4 = hmp4['Ta'] + 273.15
                
            # Get snowdepth data
        
            snd_nc = Dataset(out_loc + 'snow-height/ACE-tower_Summit_%s_snow-height_v1.nc'%dt.datetime.strftime(start,'%Y%m'),'r')
            snow_height = pd.DataFrame({'time':num2date(snd_nc.variables['time'][:],'seconds since 1970-01-01 00:00:00 0:00'),'snd':snd_nc.variables['distance_to_surface'][:]})
            snow_height.index = pd.DatetimeIndex(snow_height['time'])
            snow_height = snow_height.reindex(time_list, method='nearest',limit=20)

            # Calculate altitude above snow surface
        
            alt_HMP1 = snow_height['snd']
            alt_HMP2 = snow_height['snd'] + 2.5
            alt_HMP3 = snow_height['snd'] + 2.5 + 5.3
            alt_HMP4 = snow_height['snd'] + 2.5 + 5.3 + 3.5
        
            # check QC's
            hmp1['QC'][hmp1['Ta'].isnull()]=0
            hmp1['QC'][hmp1['Ta'].isnull()]=0
            hmp1['QC'][hmp1['Ta'].isnull()]=0
            hmp1['QC'][hmp1['Ta'].isnull()]=0
        
            # Write in data

            nc.variables['air_temperature'][:,0]=t1.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,0]=hmp1['QC'].to_numpy()
            nc.variables['altitude'][:,0]=alt_HMP1.to_numpy()
        
            nc.variables['air_temperature'][:,1]=t2.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,1]=hmp2['QC'].to_numpy()
            nc.variables['altitude'][:,1]=alt_HMP2.to_numpy()
        
            nc.variables['air_temperature'][:,2]=t3.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,2]=hmp3['QC'].to_numpy()
            nc.variables['altitude'][:,2]=alt_HMP3.to_numpy()
        
            nc.variables['air_temperature'][:,3]=t4.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,3]=hmp4['QC'].to_numpy()
            nc.variables['altitude'][:,3]=alt_HMP4.to_numpy()
    
            # Close netcdf file
    
            nc.close()
    
if __name__ == '__main__':
    main()
