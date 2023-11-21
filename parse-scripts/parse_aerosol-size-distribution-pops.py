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
        
    # Tech qc files
    qcf = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/qc-files/pops_qc.txt'
    #qcf = '/Users/heather/Desktop/ace-ceda-master/qc-files/pops_qc.txt'   
    
        
    # Global attributes
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/pops_metadata.xlsx'
    #meta_f = '/Users/heather/Desktop/ace-ceda-master/metadata/pops_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/aerosol-size-distribution_2.xlsx'
    #var_f = '/Users/heather/Desktop/ace-ceda-master/specific_variables/aerosol-size-distribution_2.xlsx'
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
    
            f1 = 'ace-pops'    #instrument
            f2 = 'summit' #platform name  
            f3 = dt.datetime.strftime(start,'%Y%m')
            f5 = "v1" #version number
            f6 = ".nc"
            fn = out_loc+'aerosol-size-distribution/' + f1 + chr(95) + f2 + chr(95) + f3 + chr(95) + 'aerosol-size-distribution' + chr(95) + f5 + f6
            nc = Dataset(fn, "w",  format = "NETCDF4_CLASSIC") 

            len_time = (24 * 60 ) / avp

            NC_Global_Attributes(nc, meta, start,stop - pd.Timedelta(minutes=1))
            time_list = pd.date_range(start,stop - pd.Timedelta(minutes=1),freq='%smin'%avp)[:]
            NC_Dimensions(nc, len(time_list),16)  
            NC_CommonVariables(nc, time_list, np)
            NC_SpecificVariables(nc, var, np)

            # Get data (qc for flights and bad winds already applied in this step)
            df_1min, data, total_conc = extract_pops(start,stop,in_loc+'extracted/pops/')

            if len(df_1min)==0:
                print('No POPS data for %s%s'%(month,year))
                continue

            # Resample to fill one month file
            df_1min = df_1min.reindex(time_list,method='nearest',tolerance='1min')
            total_conc = total_conc.reindex(time_list,method='nearest',tolerance='1min')
            data = data.reindex(time_list,method='nearest',tolerance='1min')

            # Sort QC's
            qc=np.ones(len(total_conc))
            qc[np.where(data['QC']==0)]=2 # 2b is bad due to station pollution
            qc[np.where(data['QC']==2)]=3 # 3b Can't check for station pollution
            
            # 5b, unspecified instrumnet error
            # QC for params
            # flow qc
            # if pops flow rate < 2 cm3 s-1 mark as suspect
            qc[np.where(df_1min[' POPS_Flow']<2)]=5
            # Data status should be 0
            qc[np.where(df_1min['DataStatus']!=0)]=5

            # qc for log
            bad_times = pd.read_csv(qcf,parse_dates={'start_dates':[0],'stop_dates':[1]},header=None)  
            # See if there are any bad dates in this file
            if (bad_times['start_dates'].dt.month==month).any():
                # If yes, set flags
                subset = bad_times[bad_times['start_dates'].dt.month==month]
                
                for i in range(0,len(subset)):
                    start_date = (pd.to_datetime(subset['start_dates'].iloc[i]))
                    stop_date = (pd.to_datetime(subset['stop_dates'].iloc[i]))
                    
                    qc[(total_conc.index >= start_date) & (total_conc.index <= stop_date)] = 4 # Bad data technician log
            
            # Get dNdlogD
            # The default setting for the POPS is 16 bins, and for a 16-bin configuration 
            # the bin boundaries in um are:
            nbins=16
            bounds = np.asarray([115,125,135,150,165,185,210,250,350,475,575,855,1220,1530,1990,2585,3370])/1000 # in um 

            # Write in data

            for i in range(0,len(data)):
                #try:
                mid_points,dNdlogd = get_dist(data[data.columns[0:16]].iloc[i],nbins,bounds)

                nc.variables['ambient_aerosol_particle_diameter'][i,:]=np.asarray(mid_points,dtype='float32')
                if i==0:
                    nc.variables['ambient_aerosol_particle_diameter'].valid_min = np.nanmin(np.asarray(mid_points,dtype='float32'))
                    nc.variables['ambient_aerosol_particle_diameter'].valid_max = np.nanmax(np.asarray(mid_points,dtype='float32'))
        
                nc.variables['ambient_aerosol_size_distribution'][i,:]=dNdlogd.astype('float32')
                
                if isinstance(nc.variables['ambient_aerosol_size_distribution'].valid_min,str):
                    #print(nc.variables['ambient_aerosol_size_distribution'].valid_min)
                    #print(np.nanmin(dNdlogd.astype('float32')))
                    nc.variables['ambient_aerosol_size_distribution'].valid_min=np.nanmin(dNdlogd.astype('float32'))
                    nc.variables['ambient_aerosol_size_distribution'].valid_max=np.nanmax(dNdlogd.astype('float32'))
                if np.isnan(nc.variables['ambient_aerosol_size_distribution'].valid_min):
                    #print(nc.variables['ambient_aerosol_size_distribution'].valid_min)
                    #print(np.nanmin(dNdlogd.astype('float32')))
                    nc.variables['ambient_aerosol_size_distribution'].valid_min=np.nanmin(dNdlogd.astype('float32'))
                    nc.variables['ambient_aerosol_size_distribution'].valid_max=np.nanmax(dNdlogd.astype('float32'))
                if np.nanmin(dNdlogd.astype('float32'))< nc.variables['ambient_aerosol_size_distribution'].valid_min:
                    nc.variables['ambient_aerosol_size_distribution'].valid_min=np.nanmin(dNdlogd.astype('float32'))
                if np.nanmax(dNdlogd.astype('float32'))> nc.variables['ambient_aerosol_size_distribution'].valid_max:
                    nc.variables['ambient_aerosol_size_distribution'].valid_max=np.nanmax(dNdlogd.astype('float32')) 
                #except:
                #    print('Skipping index=%s'%i)
                #    nc.variables['ambient_aerosol_particle_diameter'][i,:]=np.nan
                #    nc.variables['ambient_aerosol_size_distribution'][i,:]=np.nan
            
            nc.variables['qc_flag'][:]=qc
            nc.variables['measurement_temperature'][:]=np.asarray(df_1min['TofP'])+273.15 # in K
            nc.variables['measurement_temperature'].valid_max = np.nanmax(np.asarray(df_1min[' POPS_Flow'],dtype='float32')+273.15)
            nc.variables['measurement_temperature'].valid_mix = np.nanmin(np.asarray(df_1min[' POPS_Flow'],dtype='float32')+273.15)

            nc.variables['sample_flow_rate'][:]=np.asarray(df_1min[' POPS_Flow'])
            nc.variables['sample_flow_rate'].valid_max = np.nanmax(np.asarray(df_1min[' POPS_Flow'],dtype='float32'))
            nc.variables['sample_flow_rate'].valid_mix = np.nanmin(np.asarray(df_1min[' POPS_Flow'],dtype='float32'))
            
            # Close netcdf file
    
            nc.close()
            
if __name__ == '__main__':
    main()  
