#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs November 16 14:23 2023

@authors: Heather Guy, Michael Gallagher

Adapted for use with SAMS downloaded simba data - 05 Aug 2025
"""

import sys, os, re, glob, fnmatch, argparse # builtin python libs
sys.path.append(os.getcwd())
import numpy as np      
import datetime as dt
import pandas as pd
import netCDF4
#sys.path.append('/Users/heather/Desktop/ace-ceda-master/parse-scripts')
from NC_functions_v1 import *
from netCDF4 import Dataset, num2date

# global variables
num_time_vals = 96   # this is hardcoded, heather and grab serial data at 15 min intervals  
num_temp_sensors = 240  # number of sensors on chain 
temp_sensor_spacing = 0.02 # m
height_vals = np.arange(0, -1*num_temp_sensors*temp_sensor_spacing, -1*temp_sensor_spacing)       
nat = np.datetime64('NaT').astype('<M8[ns]') 

####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/'
#out_loc = '/Users/heather/Desktop/temp_out/'

# Months
#months=[12]
#years = [2019]

#Example usage: 
# python parse_snow-temperature-profile.py in_loc out_loc months years
#############################################################

def get_args(args_in):
    """
    check input arguments and return values if all looks o.k.
    """
    # check length of arguments:
    #if len(args_in) != 4:
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

    # return values:
    return in_loc,out_loc,months,years

def round_15(dat):
    return dt.datetime.min + round((dat - dt.datetime.min) / dt.timedelta(minutes=15)) * dt.timedelta(minutes=15)

def compare_arrays(a, b):
    # compare two numpy arrays of the same length
    # where they are not nan, they should be equal, else throw an error
    # where just one is nan, return the good number
    # where both are equal return that value. 
    if a.shape != b.shape:
        raise ValueError("Shape mismatch")
    both_nan = np.isnan(a) & np.isnan(b)
    mismatch = ~(both_nan | (a == b) | np.isnan(a) | np.isnan(b))
    if np.any(mismatch):
        raise ValueError(f"Arrays differ at positions: {np.where(mismatch)}")
    return np.where(np.isnan(a), b, a)

def main():
    """
    main function generates netCDF and stores in out_loc
    """
    # check / get args:
    try:
        in_loc,out_loc,months,years = get_args(sys.argv)
    except:
        print('Input error')
        sys.exit()
    
    # Global attributes

    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/sams-simba_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/snow-temperature-profile.xlsx'
    var = pd.read_excel(var_f)

    # Get simba data
    #dateparse = lambda x: dt.datetime.strptime(x, '%d/%m/%Y %H:%M')
    dateparse = lambda x: dt.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    
    #fpath = '/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/sams_simba/aatestV10TEST2td_2025-10-07_09-03-31.csv'
    #fpath_error = '/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/sams_simba/aatestV10TEST2st_2025-10-07_09-03-34_a.csv'

    fpath = '/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/sams_simba/251008_NCAS_RAW_aatestV10TEST2td_2025-10-08_10-03-24.csv'
    fpath_error = '/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/sams_simba/aatestV10TEST2st_2025-10-07_09-03-34_a.csv'    

    sams_data = pd.read_csv(fpath,index_col=1)
    # sort columns
    new_cols=sams_data.columns[1:]
    sams_data=sams_data[sams_data.columns[:-1]]
    sams_data.columns = new_cols
    
    sams_error = pd.read_csv(fpath_error,index_col=1)
    sams_data.index = pd.to_datetime(sams_data.index)
    sams_error.index = pd.to_datetime(sams_error.index)
    
    # clean:
    #sams_error = sams_error[~sams_error.index.duplicated()]
    #sams_data = sams_data[~sams_data.index.duplicated()]
    sams_error = sams_error.sort_index()
    sams_data = sams_data.sort_index()
    pattern = re.compile(r'^T\d+$')
    T_keys = [s for s in sams_data.columns if pattern.match(s)]
    
    # Loop through each month: 
    for year in years: 
        for month in months:
            start = dt.datetime(year, month, 1, 0, 0)
            if month==12: 
                stop = dt.datetime(year+1,1,1,0,0)
            else:
                stop = dt.datetime(year,month+1,1,0,0)

            print(start.strftime(format='%Y%m'))

            if month < 12: 
                sams_data_month = sams_data[dt.datetime(year,month,1):dt.datetime(year,month+1,1)]
                sams_error_month = sams_error[dt.datetime(year,month,1):dt.datetime(year,month+1,1)]
            else:
                sams_data_month = sams_data[dt.datetime(year,month,1):dt.datetime(year+1,1,1)]
                sams_error_month = sams_error[dt.datetime(year,month,1):dt.datetime(year+1,1,1)]
   
            # Set up netcdf files.
            
            f1 = 'sams-simba'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc + 'snow-temperature-profile/' +f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'snow-temperature-profile' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=15))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=15),freq='15min')[:]
            NC_Dimensions(nc, len(time_list), index=num_temp_sensors)  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            print(f" ... processing month: %s"%month)
            
            # Write in data

            nc.variables['height'][:]=height_vals
            qc_flag = np.ones(np.shape(nc.variables['temperature']))
            all_dtimes =[]
            
            for d in range(0,len(sams_data_month)):
                # Skip if it's a heating sample
                if sams_data_month.iloc[d]['MsgType']!=10:
                    continue
                
                # round time to nearest 15 for index
                d_time = round_15(sams_data_month.index[d].to_pydatetime())

                # Ignore any data prior to sample start time or not in the correct month
                if d_time < dt.datetime(2023,8,7,0):
                    continue
                elif d_time.month != month:
                    continue
                
                i = np.where(nc.variables['time'][:]==netCDF4.date2num(d_time,units='seconds since 1970-01-01 00:00:00 UTC'))[0][0]
                t_data = np.asarray([float(sams_data_month.iloc[d][T_keys[j]])+273.15 for j in range(0,240)]).astype('float32') # Convert to kelvin
                # fill nans
                t_data[t_data < 175] = -999.
                
                # If the time is repeated, you need to compare all lines with the same time stamp, and keep non-nan values. Any values
                # that exists in multiple lines should be duplicates. If that doesn't work, try joining the arrays.

                if d_time in all_dtimes:
                    print('Processing duplicate time stamp')
                    existing_data = nc.variables['temperature'][i,:].filled(-999.)
                    print('Attempting duplicate data method...')
                    try:
                        new_data = compare_arrays(existing_data, t_data)
                        nc.variables['temperature'][i,:] = new_data
                    except:
                        print('Duplicates failed, trying to append')
                        if len(t_data[t_data!=-999.]) + len(existing_data[existing_data!=-999.]) ==240:
                            new_data=np.concatenate([existing_data[existing_data!=-999.],t_data[t_data!=-999.]])
                            nc.variables['temperature'][i,:] = new_data
                        else:
                            print('Failed to establish a good data series')
                            qc_flag[i]=3
                            continue
                else:
                    nc.variables['temperature'][i,:] = t_data
                    all_dtimes.append(d_time)
   
                nc.variables['battery_voltage'][i] = np.float(sams_error['Battery_volts'][d_time-dt.timedelta(days=1):d_time+dt.timedelta(days=1)].mean())
                #nc.variables['qc_flag_temperature'][i] = d[]
            

            # Generate temperature flag qc
            #0 = not_used 
            #1 = good_data 
            #2 = bad_data_temperature_outside_sensor_operational_range 
            #3 = bad_data_unspecified_instrument_error          

            # Check for reasonable values: 
            qc_flag[np.where(nc.variables['temperature'][:] > 295)] = 2
            qc_flag[np.where(nc.variables['temperature'][:] < 195)] = 2

            # Check for nan's
            qc_flag[np.where(nc.variables['temperature'][:].mask==1)] = 3

            # Check that the temperature profile is continuous
            # i.e. if one sensor is wildly different to it's neighbours, flag as unspecified error. 
            diff = np.diff(nc.variables['temperature'][:])
            qc_flag[np.where(np.abs(diff)>5)] = 3

            nc.variables['qc_flag_temperature'][:]=qc_flag

            # Derive valid max and min
            nc.variables['temperature'].valid_min = np.nanmin(nc.variables['temperature'][:][qc_flag==1])
            nc.variables['temperature'].valid_max = np.nanmax(nc.variables['temperature'][:][qc_flag==1])

            nc.variables['sample_start'].valid_min = np.nanmin(nc.variables['sample_start'][:])
            nc.variables['sample_start'].valid_max = np.nanmax(nc.variables['sample_start'][:])

            nc.variables['sample_end'].valid_min = np.nanmin(nc.variables['sample_end'][:])
            nc.variables['sample_end'].valid_max = np.nanmax(nc.variables['sample_end'][:])

            nc.variables['sample_span'].valid_min = np.nanmin(nc.variables['sample_span'][:])
            nc.variables['sample_span'].valid_max = np.nanmax(nc.variables['sample_span'][:])

            nc.variables['battery_voltage'].valid_min = 10.0
            nc.variables['battery_voltage'].valid_max = 20.0

            nc.variables['height'].valid_min = -5
            nc.variables['height'].valid_max = 0

            # Close netcdf file
            nc.close()
    
if __name__ == '__main__':
    main()
