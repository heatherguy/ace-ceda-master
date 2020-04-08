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
from ace_parse import *
from netCDF4 import Dataset, num2date

####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/ace/processed/'
#out_loc = '/Users/heather/Desktop/temp_out/'

# Months
#months=[11]
#years = [2019]

# sample interval
#avp = 1 # minutes

#Example usage: 
# python parse_aerosol-concentration.py $in_loc $out_loc $months $years $avp
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
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/CPC_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/aerosol-concentration.xlsx'
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
    
            f1 = 'ace-cpc'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc+'aerosol-concentration/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'aerosol-concentration' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / avp

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='%smin'%avp)[:]
            NC_Dimensions(nc, len(time_list))  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get data
        
            CPC = get_cpc(start,stop,in_loc+'CPC/')
        
            # Correct QC's for nans
        
            CPC['QC'][CPC['c/cm3'].isnull()]=0
              
            # Write in data

            nc.variables['number_concentration_of_ambient_aerosol_particles_in_air'][:]=CPC['c/cm3'].to_numpy()
            nc.variables['qc_flag'][:]=CPC['QC'].to_numpy()

    
            # Close netcdf file
    
            nc.close()
            
if __name__ == '__main__':
    main()  