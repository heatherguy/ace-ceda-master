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
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/'
#out_loc = '/Users/heather/Desktop/temp_out/'

# Months
#months=[11]
#years = [2019]

# sample interval
#avp = 1 # minutes

#Example usage: 
# python parse_snow-height.py $in_loc $out_loc $months $years $avp


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
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/snd_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/snow-height.xlsx'
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
    
            f1 = 'ace-tower'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc + 'snow-height/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'snow-height' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / avp

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=10))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=10),freq='%smin'%avp)[:]
            NC_Dimensions(nc, len(time_list))  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get data
    
            snd = extract_snd_data(time_list[0],time_list[-1],in_loc+'raw_extracted/SnD/',in_loc+'processed/HMP/',save=False)
            dat = snd.reindex(time_list, method='nearest',limit=12)
    
            # Write in data
    
            nc.variables['distance_to_surface'][:]=dat['depth_Tcorrected'].to_numpy()
            qc_dat = np.ones(len(dat))
            nan_inds = np.where(dat['depth_Tcorrected'].isnull())[0]
            qc_dat[nan_inds]=0    
            nc.variables['qc_flag_snow_height'][:]=qc_dat
    
            # Close netcdf file
    
            nc.close()
            
if __name__ == '__main__':
    main() 