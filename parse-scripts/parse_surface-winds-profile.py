#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs March 12 17:22:21 2020

@author: Heather Guy
"""

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
#months=[11]
#years = [2019]

# sample interval
#avp = 1 # minutes

#Example usage: 
# python parse_surface-winds-profile.py $in_loc $out_loc $months $years $avp
#############################################################

# For 3D sonics, use diff = -45
# For 2D sonics, use diff = -18.07
def deg_rot(orig,diff):
    """
    Adds a fixed number of degrees to a wind direction. 
    To account for sonics not being oriented to north. 
    
    Parameters:
        orig:  list or series of wind directions in degrees
        diff:  number of degrees clockwise to rotate (can be negative).
        
    Returns: 
        new: Rotated wind direction in degrees
    """   
    new = orig
    new = new + diff
    new[new<0] = 360 + (orig[new<0] + diff)
    new[new>360] = diff - (360 - orig[new>360])
    new[new==360] = 0

    return new

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
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/wind_profile_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/surface-winds-profile.xlsx'
    var = pd.read_excel(var_f)
    
    # Loop through each ,pmtj: 
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
            fn = out_loc+ 'surface-winds-profile/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'surface-winds-profile' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / avp

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='%smin'%avp)[:]
            NC_Dimensions(nc, len(time_list), index=4)  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get ventus data
        
            print('Extracting v1')
            v1 = get_ventus(start,stop,in_loc+'ventus/','v1', avp='1min')
            v1['wdir_corrected'] = deg_rot(v1['wdir'],-18.07)
            print('Extracting v2')
            v2 = get_ventus(start,stop,in_loc+'ventus/','v2', avp='1min')
            v2['wdir_corrected'] = deg_rot(v2['wdir'],-18.07)
        
            # Get metek data

            print('Extracting m1')    
            m1 = get_metek(start,stop,in_loc+'metek/','metek1', avp='1min')    
            m1['wdir_corrected'] = deg_rot(m1['wdir'],-45)     
            print('Extracting m2')    
            m2 = get_metek(start,stop,in_loc+'metek/','metek2', avp='1min')    
            m2['wdir_corrected'] = deg_rot(m2['wdir'],-45) 
                
            # Get snowdepth data
            
            snd_nc = Dataset(out_loc + 'snow-height/ace-tower_summit_%s_snow-height_v1.nc'%dt.datetime.strftime(start,'%Y%m'),'r')
            snow_height = pd.DataFrame({'time':pd.to_datetime(snd_nc.variables['time'][:],origin='unix',unit='s')),'snd':snd_nc.variables['distance_to_surface'][:]})
            snow_height.index = pd.DatetimeIndex(snow_height['time'])
            snow_height = snow_height.reindex(time_list, method='nearest',limit=20)

            # Calculate altitude above snow surface
        
            alt_1 = snow_height['snd']
            alt_2 = snow_height['snd'] + 2.5
            alt_3 = snow_height['snd'] + 2.5 + 5.3
            alt_4 = snow_height['snd'] + 2.5 + 5.3 + 3.5
        
            # Check QC's for missing data
            v1['QC'][v1.isnull().any(axis=1)]=0
            v2['QC'][v2.isnull().any(axis=1)]=0
            m1['QC'][m1.isnull().any(axis=1)]=0
            m2['QC'][m2.isnull().any(axis=1)]=0
        
            # Write in data

            nc.variables['wind_speed'][:,0]= m1['wsd'].to_numpy()
            nc.variables['wind_from_direction'][:,0]= m1['wdir_corrected'].to_numpy()
            nc.variables['qc_flag'][:,0]= m1['QC'].to_numpy()
            nc.variables['altitude'][:,0]= alt_1.to_numpy()
        
            nc.variables['wind_speed'][:,1]= v1['wsd'].to_numpy()
            nc.variables['wind_from_direction'][:,1]= v1['wdir_corrected'].to_numpy()
            nc.variables['qc_flag'][:,1]= v1['QC'].to_numpy()
            nc.variables['altitude'][:,1]= alt_2.to_numpy()
        
            nc.variables['wind_speed'][:,2]= v2['wsd'].to_numpy()
            nc.variables['wind_from_direction'][:,2]= v2['wdir_corrected'].to_numpy()
            nc.variables['qc_flag'][:,2]= v2['QC'].to_numpy()
            nc.variables['altitude'][:,2]= alt_3.to_numpy()
        
            nc.variables['wind_speed'][:,3]= m2['wsd'].to_numpy()
            nc.variables['wind_from_direction'][:,3]= m2['wdir_corrected'].to_numpy()
            nc.variables['qc_flag'][:,3]= m2['QC'].to_numpy()
            nc.variables['altitude'][:,3]= alt_4.to_numpy()
    
            # Close netcdf file
    
            nc.close()
            
if __name__ == '__main__':
    main()   