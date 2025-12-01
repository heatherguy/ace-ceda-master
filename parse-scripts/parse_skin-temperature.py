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
#sys.path.append('/Users/heather/Desktop/ace-ceda-master/parse-scripts')
#sys.path.append('/Users/heather/ICECAPS-ACE/DataParse')
from NC_functions_v1 import *
from fluxtower_parse import *
from netCDF4 import Dataset, num2date


####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/processed/'
#out_loc = '/Users/heather/Desktop/temp_out/'

# Months
#months=[6]
#years = [2019]

# sample interval
#avp = 1 # minutes

#Example usage: 
# python parse_skin-temperature.py $in_loc $out_loc $months $years $avp
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
    #meta_f= '/Users/heather/Desktop/ace-ceda-master/metadata/kt_metadata.xlsx'
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/kt_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    #var_f= '/Users/heather/Desktop/ace-ceda-master/specific_variables/skin-temperature.xlsx'
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/skin-temperature.xlsx'
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
            fn = out_loc+'skin-temperature/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'skin-temperature' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / avp

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='%smin'%avp)[:]
            NC_Dimensions(nc, len(time_list))  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get data
        
            KT = get_kt(start,stop,in_loc+'KT15/')
            dat = KT.reindex(time_list, method='nearest',limit=2)
        
            # Convert SkinT to Kelvin
        
            skin_temp = dat['T'] + 273.15
        
            # Make sure QC flag is 1 (good) or 2 (suspect)
        
            qc=np.ones(len(dat)) # good data
            qc[np.where(dat['QC']==0)[0]]=2 # suspect data in log
            qc[np.where(dat['QC'].isnull())[0]]=0 # no data
            qc[np.where(dat['QC']==24)[0]]=3 # Error code 24, under ref limit

            # Check in instrument range of -100 to 100 C
            qc[np.where(dat['T']>100)[0]]=3
            qc[np.where(dat['T']<-100)[0]]=3
    
            # Write in data

            nc.variables['skin_temperature'][:]=skin_temp.to_numpy()
            nc.variables['skin_temperature'].valid_min = skin_temp[dat['QC']!=0].min()
            nc.variables['skin_temperature'].valid_max= skin_temp[dat['QC']!=0].max()
            
            nc.variables['qc_flag_skin_temperature'][:]=qc
    
            # Write note to netcdf file indicating date used.
            if start < dt.datetime(2020,7,1):
                base_str = 'Platform altitude is the top of the Met tower, Instrument altitude is platform altitude minus 10.8 m until 2020-06-26, after which it is platform altitude minus 10.3 m. '
            elif start < dt.datetime(2022,6,1):   
                base_str = 'Platform altitude is the top of the Met tower. Instrument altitude is platform altitude minus 10.3 m. Technician log used to QC data is located here: https://github.com/heatherguy/ace-ceda-master/blob/master/qc-files/KT_bad_dates. Spectral range of sensor is 9.6 to 11.5 um.'
            elif start < dt.datetime(2022,8,1):
                base_str = 'Platform altitude is the top of the Met tower. Instrument altitude is platform altitude minus 10.3 m until 2022-06-06 1000, after which it was raised by 1.1 m to h0-9.2. Intermittent data throughout June and July 2022 due to a damaged cable. Technician log used to QC data is located here: https://github.com/heatherguy/ace-ceda-master/blob/master/qc-files/KT_bad_dates. Spectral range of sensor is 9.6 to 11.5 um.'
            elif start < dt.datetime(2022,9,1): 
                base_str = 'Platform altitude is the top of the Met tower. Instrument altitude is platform altitude minus 9.2 m. Instrument reinstalled after cable repair on 20 August 2022. Technician log used to QC data is located here: https://github.com/heatherguy/ace-ceda-master/blob/master/qc-files/KT_bad_dates. Spectral range of sensor is 9.6 to 11.5 um.'
            elif start < dt.datetime(2023,7,27): 
                base_str = 'Platform altitude is the top of the Met tower. Instrument altitude is platform altitude minus 9.2 m. Technician log used to QC data is located here: https://github.com/heatherguy/ace-ceda-master/blob/master/qc-files/KT_bad_dates. Spectral range of sensor is 9.6 to 11.5 um.'
            else:
                base_str = 'Technician log used to QC data is located here: https://github.com/heatherguy/ace-ceda-master/blob/master/qc-files/KT_bad_dates. Spectral range of sensor is 9.6 to 11.5 um. Platform altitude is the top of the Met tower. Nominal instrument height above snow surface is 2 m'
            
            nc.setncattr('comment', base_str)

            # Close netcdf file
    
            nc.close()
            
if __name__ == '__main__':
    main()  
