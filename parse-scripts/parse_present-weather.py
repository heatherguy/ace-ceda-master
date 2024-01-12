#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sunday 16 April 16:59:21 2023

@author: Heather Guy
"""

import numpy as np      
import datetime as dt
import pandas as pd
import sys
import os
import glob
#sys.path.append('/Users/heather/Desktop/NC_parse_master')
#sys.path.append('/Users/heather/ICECAPS-ACE/DataParse')
from NC_functions_v1 import *
from ace_parse import *
from netCDF4 import Dataset, num2date


####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/'
#out_loc = '/Users/heather/Desktop/temp_out/'

# Months
# months=[11]
# years = [2019]

# Example usage: 
# python parse_present-weather.py $in_loc $out_loc $months $years


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

    # return values:
    return in_loc,out_loc,months,years


def main():
    """
    main function generates netCDF and stores in out_loc
    """
    
    # check / get args:
    try:
        in_loc,out_loc,months,years= get_args(sys.argv)
    except:
        print('Input error')
        sys.exit()
        
    # Global attributes
    #meta_f= '/Users/heather/Desktop/ace-ceda-master/metadata/biral_metadata.xlsx'
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/biral_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    #var_f= '/Users/heather/Desktop/ace-ceda-master/specific_variables/present-weather.xlsx'
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/present-weather.xlsx'
    var = pd.read_excel(var_f)
    
    for year in years: 
        for month in months:
            start = dt.datetime(year, month, 1, 0, 0)
            if month==12: 
                stop = dt.datetime(year+1,1,1,0,0)
            else:
                stop = dt.datetime(year,month+1,1,0,0)

            print(start.strftime(format='%Y%m'))
    
            # Set up netcdf files.
    
            f1 = 'ace-biral'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc + 'present-weather/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'present-weather' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / 1 # 1 minute averaging period

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='1min')[:]
            NC_Dimensions(nc, len(time_list))  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get data
    
            biral = extract_biral(time_list[0],time_list[-1],in_loc+'extracted/biral/')
            
            # Write in data
            nc.variables['present_weather_code'][:]=biral['pw'].fillna(-999).astype(int).to_numpy()
            nc.variables['optical_range'][:] = biral['avg_vis'].to_numpy()*1000
            nc.variables['instantaneous_optical_range'][:]= biral['inst_vis'].to_numpy()*1000
                
            # qc
            # no_data = 0
            # good_data = 1
            # suspect_data_failed_self_check = 2
            qc = np.ones(len(biral))
            qc[np.where(biral['avg_vis'].isnull())[0]]=0 # no data
            qc[np.where(biral['qc']==0)[0]]=2 # failed self-check
            nc.variables['qc_flag'][:]=qc

            # Derive valid max and min
            nc.variables['optical_range'].valid_min = np.nanmin(biral['avg_vis'].to_numpy()*1000)
            nc.variables['optical_range'].valid_max = np.nanmax(biral['avg_vis'].to_numpy()*1000)
            nc.variables['instantaneous_optical_range'].valid_min = np.nanmin(biral['inst_vis'].to_numpy()*1000)
            nc.variables['instantaneous_optical_range'].valid_max = np.nanmax(biral['inst_vis'].to_numpy()*1000)
            
            # Close netcdf file
            nc.close()
            
if __name__ == '__main__':
    main() 
