#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs November 16 14:23 2023

@authors: Heather Guy, Michael Gallagher
"""

import sys, os, re, glob, fnmatch, argparse # builtin python libs
sys.path.append(os.getcwd())
import numpy as np      
import datetime as dt
import pandas as pd
import netCDF4
#sys.path.append('/Users/heather/Desktop/ace-ceda-master/parse-scripts')
from NC_functions_v1 import *
from simba_functions import parse_simba_serialstream, get_file_list_serial, extract_serial_samples
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

    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/simba_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/snow-temperature-profile.xlsx'
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
            
            f1 = 'ace-simba'    #instrument
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

            # Get simba data
        
            print('Extracting Simba')
            all_data=[]
            print(f" ... processing month: %s"%month)
            # if the simba daq was restarted today then there are multiple files
            all_files_today = get_file_list_serial(dt.datetime(year,month,1), in_loc+'extracted/simba/')
            file_samples    = extract_serial_samples(all_files_today, '2')
               
            if len(all_files_today)==0: print(f"!!! no data on for month: %s !!!"%month); continue
                
            all_samples = []
            tot_samples = 0
            for file_data in file_samples.items():
                sample_dict = file_data[1][0]; tstamp_dict = file_data[1][1]; 
                parsed_dict = parse_simba_serialstream(sample_dict, tstamp_dict, num_temp_sensors)
                tot_samples += len(parsed_dict.keys())
                for key,value in parsed_dict.items():
                    all_samples.append(value)

            # merge the different data dicts for each file, then write out data
            if tot_samples == 0: 
                print(f"¿¡¡ There was no data for month: %s"%month+ \
                    f"is there a good reason for that?!!")
                continue
            
            print(f"... writing out {tot_samples} total data samples in timeseries !!!")
            
            # Write in data

            nc.variables['height'][:]=height_vals

            for d in all_samples:
                # round time to nearest 15 for index
                d_time = round_15(d['sample_start'])

                # Ignore any data prior to sample start time
                if d_time < dt.datetime(2022,8,22,0):
                    continue

                i = np.where(nc.variables['time'][:]==netCDF4.date2num(d_time,units='seconds since 1970-01-01 00:00:00 UTC'))[0][0]
                
                if len(d['temperature'])!=num_temp_sensors:
                    print('Not enough tempearutre values for %s'%d_time)
                else:
                    for j in range(0,len(d['temperature'])):
                        try: 
                            nc.variables['temperature'][i,j] = float(d['temperature'][j])+273.15 # Conver to kelvin
                        except:
                            print('Cannot convert temperature to float: %s'%d['temperature'][j])
                            continue
                
                nc.variables['sample_start'][i] = netCDF4.date2num(d['sample_start'],units='seconds since 1970-01-01 00:00:00 UTC')
                nc.variables['sample_end'][i] = netCDF4.date2num(d['sample_end'],units='seconds since 1970-01-01 00:00:00 UTC')
                nc.variables['sample_span'][i] = d['sample_span'].total_seconds()
                nc.variables['battery_voltage'][i] = np.float(d['battery_voltage'])
                #nc.variables['qc_flag_temperature'][i] = d[]
            

            # Generate temperature flag qc
            #0 = not_used 
            #1 = good_data 
            #2 = bad_data_temperature_outside_sensor_operational_range 
            #3 = bad_data_unspecified_instrument_error          
            qc_flag = np.ones(np.shape(nc.variables['temperature']))

            # Check for reasonable values: 
            qc_flag[np.where(nc.variables['temperature'][:] > 295)] = 2
            qc_flag[np.where(nc.variables['temperature'][:] < 195)] = 2

            # Check for nan's
            qc_flag[np.where(nc.variables['temperature'][:].mask==1)] = 3

            # Check that the temperature profile is continuous
            # i.e. if one sensor is wildly different to it's neighbours, flag as unspecified error. 
            diff = np.diff(nc.variables['temperature'][:])
            qc_flag[np.where(np.abs(diff)>2)] = 3
            
            nc.variables['qc_flag_temperature'][:]=qc_flag

            # Derive valid max and min
            nc.variables['temperature'].valid_min = np.nanmin(nc.variables['temperature'][:])
            nc.variables['temperature'].valid_max = np.nanmax(nc.variables['temperature'][:])

            nc.variables['sample_start'].valid_min = np.nanmin(nc.variables['sample_start'][:])
            nc.variables['sample_start'].valid_max = np.nanmax(nc.variables['sample_start'][:])

            nc.variables['sample_end'].valid_min = np.nanmin(nc.variables['sample_end'][:])
            nc.variables['sample_end'].valid_max = np.nanmax(nc.variables['sample_end'][:])

            nc.variables['sample_span'].valid_min = np.nanmin(nc.variables['sample_span'][:])
            nc.variables['sample_span'].valid_max = np.nanmax(nc.variables['sample_span'][:])

            nc.variables['battery_voltage'].valid_min = np.nanmin(nc.variables['battery_voltage'][:])
            nc.variables['battery_voltage'].valid_max = np.nanmax(nc.variables['battery_voltage'][:])

            nc.variables['height'].valid_min = -5
            nc.variables['height'].valid_max = 0


            # Write note to netcdf file indicating date used.
            #base_str = 'Platform altitude (h0) is the top of the Met tower. Instrument altitude: HMP1=h0-10.87m, HMP2=h0-9.8m, HMP3=h0-4.5m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4].'
            #nc.setncattr('comment', base_str)
            # Close netcdf file
    
            nc.close()
    
if __name__ == '__main__':
    main()
