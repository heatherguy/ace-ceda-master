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

def round_1h(dat):
    return dt.datetime.min + round((dat - dt.datetime.min) / dt.timedelta(hours=1)) * dt.timedelta(hours=1)

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
    #meta_f = '/Users/heather/Desktop/ace-ceda-master/metadata/sams-simba_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/snow-temperature-profile.xlsx'
    #var_f = '/Users/heather/Desktop/ace-ceda-master/specific_variables/snow-temperature-profile.xlsx'
    var = pd.read_excel(var_f)

    # Get simba data
    #dateparse = lambda x: dt.datetime.strptime(x, '%d/%m/%Y %H:%M')
    dateparse = lambda x: dt.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
    
    fpath_error = '/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/sams_simba/aatestV10TEST2st_2025-10-07_09-03-34_a.csv'
    fpath = '/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/sams_simba/251010_aatestV10TEST2td_2025-10-10_11-40-35.csv'

    sams_data = pd.read_csv(fpath,index_col=1)

    # sort columns
    #new_cols=sams_data.columns[1:]
    #sams_data=sams_data[sams_data.columns[:-1]]
    #sams_data.columns = new_cols
    
    sams_error = pd.read_csv(fpath_error,index_col=1)
    sams_data.index = pd.to_datetime(sams_data.index,dayfirst=True)
    sams_error.index = pd.to_datetime(sams_error.index,dayfirst=True)
    
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
   
            # Skip if it's a heating sample
            sams_data_month = sams_data_month[sams_data_month['MsgType'] == 10]
            sams_data_times = np.asarray([round_1h(d.to_pydatetime()) for d in sams_data_month.index])

            time_list = list(set(sams_data_times))
            time_list.sort()

            # Set up netcdf files.
            
            f1 = 'sams-simba'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc + 'snow-temperature-profile/' +f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'snow-temperature-profile' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(hours=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(hours=1),freq='1h')[:]
            NC_Dimensions(nc, len(time_list), index=num_temp_sensors)  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            print(f" ... processing month: %s"%month)
            
            # Write in data

            nc.variables['height'][:]=height_vals
            qc_flag = np.ones(np.shape(nc.variables['temperature']))
            
            for d in range(0,len(time_list)):
                d_time=time_list[d]
                inds = np.where(sams_data_times==d_time)[0]
                i = np.where(nc.variables['time'][:]==netCDF4.date2num(d_time,units='seconds since 1970-01-01 00:00:00 UTC'))[0][0]

                if len(inds)==0:
                    print('No data found for %s'%d_time)
                    continue

                elif len(inds)==1:
                    temp_data = sams_data_month.iloc[inds]
                    if temp_data["MsgPart"].to_numpy()==1:
                        t_data = np.asarray([float(sams_data_month.iloc[inds[0]][T_keys[j]])+273.15 for j in range(0,240)]).astype('float32') # Convert to kelvin
                        # fill nans
                        t_data[t_data < 175] = -999.
                    else: 
                        # In this case we only have the second part of the message
                        t_data = np.ones(240)*-999.
                        m2_data = np.asarray([float(sams_data_month.iloc[inds[0]][T_keys[j]])+273.15 for j in range(0,240)]).astype('float32') # Convert to kelvin
                        m2_data[m2_data < 175] = np.nan
                        m2_data = m2_data[~np.isnan(m2_data)]
                        t_data[-len(m2_data):]=m2_data
                        # round time to nearest 15 for index
                        d_time = round_1h(sams_data_month.index[d].to_pydatetime())

                else:
                    # make sure we get the messages in the correct order
                    temp_data = sams_data_month.iloc[inds]
                    temp_data = temp_data.sort_values(by="MsgPart")
        
                    # try dropping duplicate rows
                    temp_data = temp_data.drop_duplicates(subset="MsgPart")
        
                    if len(temp_data) == 2:
            
                        p1 = np.asarray([float(temp_data.iloc[0][T_keys[j]])+273.15 for j in range(0,240)]).astype('float32') # Convert to kelvin
                        p2 = np.asarray([float(temp_data.iloc[1][T_keys[j]])+273.15 for j in range(0,240)]).astype('float32') # Convert to kelvin
            
                        # drop nans and check they add up
                        p1[p1<175] = np.nan
                        p2[p2<175] = np.nan
                        p1 = p1[~np.isnan(p1)]
                        p2 = p2[~np.isnan(p2)]
            
                        if (len(p1) + len(p2)) == 240:
                            t_data = np.concatenate((p1,p2))
                        else:
                            print('Total length of data %s, skipping %s'%(len(p1) + len(p2),d_time))
                            qc_flag[i]=3
                    else:
                        print('More duplicate data than expecting, skipping %s'%d_time)
                        qc_flag[i]=3

                # fill nans
                t_data[t_data < 175] = -999.
                t_data[np.isnan(t_data)] = -999.

                nc.variables['temperature'][i,:] = t_data
                nc.variables['battery_voltage'][i] = float(sams_error['Battery_volts'][d_time-dt.timedelta(days=2):d_time+dt.timedelta(days=2)].mean())
                #nc.variables['qc_flag_temperature'][i] = d[]
    
            # Generate temperature flag qc
            #0 = not_used 
            #1 = good_data 
            #2 = bad_data_temperature_outside_sensor_operational_range 
            #3 = bad_data_unspecified_instrument_error          

            # Check for reasonable values: 
            qc_flag[np.where(nc.variables['temperature'][:] > 295)] = 2
            qc_flag[np.where(nc.variables['temperature'][:] < 195)] = 2
            nc.variables['temperature'][:] = np.ma.masked_where(nc.variables['temperature'][:] > 295, nc.variables['temperature'][:])
            nc.variables['temperature'][:] = np.ma.masked_where(nc.variables['temperature'][:] < 195, nc.variables['temperature'][:])

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

            nc.variables['battery_voltage'].valid_min = 10.0
            nc.variables['battery_voltage'].valid_max = 20.0

            nc.variables['height'].valid_min = -5
            nc.variables['height'].valid_max = 0

            # Close netcdf file
            nc.close()
    
if __name__ == '__main__':
    main()
