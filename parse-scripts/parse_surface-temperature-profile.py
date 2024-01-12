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
from NC_functions_v1 import *
from fluxtower_parse import *
from netCDF4 import Dataset, num2date


####### INPUTS #######
# Data location:
#in_loc = '/Volumes/Data/ICECAPSarchive/fluxtower/'
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
    # meta_f='/Users/heather/Desktop/ace-ceda-master/metadata/trh_profile_metadata.xlsx'
    meta_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/metadata/trh_profile_metadata.xlsx'
    meta = pd.read_excel(meta_f)

    # Specific variables
    # var_f = '/Users/heather/Desktop/ace-ceda-master/specific_variables/surface-temperature-profile.xlsx'
    var_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/specific_variables/surface-temperature-profile.xlsx'
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
            HMP1 = get_hmp(start,stop,in_loc+'processed/HMP/','HMP1')
            print('Extracting HMP2')
            HMP2 = get_hmp(start,stop,in_loc+'processed/HMP/','HMP2')
            print('Extracting HMP3')
            HMP3 = get_hmp(start,stop,in_loc+'processed/HMP/','HMP3')
            print('Extracting HMP4')
            HMP4 = get_hmp(start,stop,in_loc+'processed/HMP/','HMP4')
            print('Post-processing')
            
            # Minutely averages
            if len(HMP1)!=0:
                hmp1 = HMP1.resample(rule = '1min').mean()
                hmp1 = HMP1.reindex(time_list, method='nearest',limit=2)
            else:
                hmp1 = pd.DataFrame({'RH':np.nan,'Ta':np.nan,'QC':np.nan},index = time_list)

            if len(HMP2)!=0:
                hmp2 = HMP2.resample(rule = '1min').mean()
                hmp2 = HMP2.reindex(time_list, method='nearest',limit=2)
            else:
                hmp2 = pd.DataFrame({'RH':np.nan,'Ta':np.nan,'QC':np.nan},index = time_list)


            if len(HMP3)!=0:
                hmp3 = HMP3.resample(rule = '1min').mean()
                hmp3 = HMP3.reindex(time_list, method='nearest',limit=2)
            else:
                hmp3 = pd.DataFrame({'RH':np.nan,'Ta':np.nan,'QC':np.nan},index = time_list)


            if len(HMP4)!=0:
                hmp4 = HMP4.resample(rule = '1min').mean()
                hmp4 = HMP4.reindex(time_list, method='nearest',limit=2)        
            else:
                hmp4 = pd.DataFrame({'RH':np.nan,'Ta':np.nan,'QC':np.nan},index = time_list)

 
            # Convert T to Kelvin
        
            t1 = hmp1['Ta'] + 273.15
            t2 = hmp2['Ta'] + 273.15
            t3 = hmp3['Ta'] + 273.15
            t4 = hmp4['Ta'] + 273.15
                
            # Get snowdepth data
        
            snd_nc = Dataset(out_loc + 'snow-height/ace-tower_summit_%s_snow-height_v1.nc'%dt.datetime.strftime(start,'%Y%m'),'r')
            snow_height = pd.DataFrame({'time':pd.to_datetime(snd_nc.variables['time'][:],origin='unix',unit='s'),'snd':snd_nc.variables['distance_to_surface'][:],'qc':snd_nc.variables['qc_flag_distance_to_surface'][:]})
            # QC
            snow_height['snd'][(snow_height['qc'].astype(int)==0) | (snow_height['qc'].astype(int)==2)| (snow_height['qc'].astype(int)==3)]=np.nan           
            snow_height.index = pd.DatetimeIndex(snow_height['time'])
            snow_height = snow_height.reindex(time_list, method='nearest',limit=20)
            # if qc flag == 4, add note
            #if (snow_height['qc'].astype(int) ==4).any():
            #    base_str = 'Platform height (h0) is the top of the Met tower. Instrument height: HMP1=h0-12.3m, HMP2=h0-9.8m, HMP3=h0-4.5m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4]'
            #    new_str = 'Note: Height above snow surface corresponds to last available data on 2020-04-21.'
            #    nc.setncattr('comment', new_str+base_str)


            # Calculate altitude above snow surface
            # Initialise
            alt_HMP1 = snow_height['snd'] * np.nan
            alt_HMP2 = snow_height['snd'] * np.nan
            alt_HMP3 = snow_height['snd'] * np.nan
            alt_HMP4 = snow_height['snd'] * np.nan
            
            alt_HMP1[dt.datetime(2019,5,1):dt.datetime(2021,7,1)] = snow_height['snd'][dt.datetime(2019,5,1):dt.datetime(2021,7,1)] - 0.4
            alt_HMP2[dt.datetime(2019,5,1):dt.datetime(2021,7,1)] = snow_height['snd'][dt.datetime(2019,5,1):dt.datetime(2021,7,1)] + 1.03
            alt_HMP3[dt.datetime(2019,5,1):dt.datetime(2021,7,1)] = snow_height['snd'][dt.datetime(2019,5,1):dt.datetime(2021,7,1)] + 1.03 + 5.3
            alt_HMP4[dt.datetime(2019,5,1):dt.datetime(2021,7,1)] = snow_height['snd'][dt.datetime(2019,5,1):dt.datetime(2021,7,1)] + 1.03 + 5.3 + 3.5

            # After 1 July 2021 tower raise
            alt_HMP2[dt.datetime(2021,7,1):dt.datetime(2022,6,6,10)] = snow_height['snd'][dt.datetime(2021,7,1):dt.datetime(2022,6,6,10)] + 1.03
            alt_HMP3[dt.datetime(2021,7,1):dt.datetime(2022,6,6,10)] = snow_height['snd'][dt.datetime(2021,7,1):dt.datetime(2022,6,6,10)] + 1.03 + 5.3
            alt_HMP4[dt.datetime(2021,7,1):dt.datetime(2022,6,6,10)] = snow_height['snd'][dt.datetime(2021,7,1):dt.datetime(2022,6,6,10)] + 1.03 + 5.3 + 3.5

            # After 6 June 2022 raise of lower level booms
            alt_HMP2[dt.datetime(2022,6,6,10):dt.datetime(2022,7,22,14)] = snow_height['snd'][dt.datetime(2022,6,6,10):dt.datetime(2022,7,22,14)] + 0.8
            alt_HMP3[dt.datetime(2022,6,6,10):dt.datetime(2022,7,22,14)] = snow_height['snd'][dt.datetime(2022,6,6,10):dt.datetime(2022,7,22,14)] + 0.8 + 4.7
            alt_HMP4[dt.datetime(2022,6,6,10):dt.datetime(2022,7,22,14)] = snow_height['snd'][dt.datetime(2022,6,6,10):dt.datetime(2022,7,22,14)] + 0.8 + 4.7 + 3.5 

            # After 23 July 2022 raise of mid level booms and some lower adjustments
            alt_HMP2[dt.datetime(2022,7,22,14):dt.datetime(2022,8,20,20,30)] = snow_height['snd'][dt.datetime(2022,7,22,14):dt.datetime(2022,8,20,20,30)] + 0.7
            alt_HMP3[dt.datetime(2022,7,22,14):dt.datetime(2022,8,20,20,30)] = snow_height['snd'][dt.datetime(2022,7,22,14):dt.datetime(2022,8,20,20,30)] + 0.7 + 5.3
            alt_HMP4[dt.datetime(2022,7,22,14):dt.datetime(2022,8,20,20,30)] = snow_height['snd'][dt.datetime(2022,7,22,14):dt.datetime(2022,8,20,20,30)] + 0.7 + 5.3 + 2.9

            # After 20 August 2022 HMP1 reinstalled 
            alt_HMP1[dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] = snow_height['snd'][dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] + 0.7
            alt_HMP2[dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] = snow_height['snd'][dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] + 0.7 + 1.8
            alt_HMP3[dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] = snow_height['snd'][dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] + 0.7 + 1.8 + 3.5
            alt_HMP4[dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] = snow_height['snd'][dt.datetime(2022,8,20,20,30):dt.datetime(2023,7,27)] + 0.7 + 1.8 + 3.5 + 2.9
          
            # After 27 July 2023 tower raise and segmed added
            alt_HMP1[dt.datetime(2023,7,27):] = snow_height['snd'][dt.datetime(2023,7,27):] + 1
            alt_HMP2[dt.datetime(2023,7,27):] = snow_height['snd'][dt.datetime(2023,7,27):] + 1 + 2.4
            alt_HMP3[dt.datetime(2023,7,27):] = snow_height['snd'][dt.datetime(2023,7,27):] + 1 + 2.4 + 3.8
            alt_HMP4[dt.datetime(2023,7,27):] = snow_height['snd'][dt.datetime(2023,7,27):] + 1 + 2.4 + 3.8 + 7.2
            
            # do QC's
            # 1b is good data
            qc1 = np.ones(len(t1))
            qc2 = np.ones(len(t2))
            qc3 = np.ones(len(t3))
            qc4 = np.ones(len(t4))
            
            # 2b is outside of operational range
            # temp range: -80C to 60C
            qc1[np.where(t1-273.15< -80)[0]]=2
            qc1[np.where(t1-273.15> 60)[0]]=2
            qc2[np.where(t2-273.15< -80)[0]]=2
            qc2[np.where(t2-273.15> 60)[0]]=2
            qc3[np.where(t3-273.15< -80)[0]]=2
            qc3[np.where(t3-273.15> 60)[0]]=2
            qc4[np.where(t4-273.15< -80)[0]]=2
            qc4[np.where(t4-273.15> 60)[0]]=2
         
            # Add in one off errors here
            # Betwen 18-Aug-20 and 3-Sept-20, the fan on the hmp4 aspirated shield was down. 
            hmp4['QC'][dt.datetime(2020,8,18):dt.datetime(2020,9,3)]=2     
            hmp2.loc[dt.datetime(2021,11,15,14):dt.datetime(2021,11,15,16),'QC']=2
            # Between 7 June 2022 and 20 August 2022, HMP4 aspirated shield was down
            hmp4.loc[dt.datetime(2022,6,7):dt.datetime(2022,8,20),'QC']=2
            # 20th july 2023, hmp1 aspirated shield is down
            hmp1.loc[dt.datetime(2023,7,20):dt.datetime(2023,7,21,13,5),'QC']=2
 
            # 3b is unspecified instrument error
            qc1[np.where(hmp1['QC']==2)[0]]=3
            qc2[np.where(hmp2['QC']==2)[0]]=3
            qc3[np.where(hmp3['QC']==2)[0]]=3
            qc4[np.where(hmp4['QC']==2)[0]]=3
          
            # 0b for no data
            qc1[np.where(hmp1['Ta'].isnull())[0]]=0
            qc2[np.where(hmp2['Ta'].isnull())[0]]=0
            qc3[np.where(hmp3['Ta'].isnull())[0]]=0
            qc4[np.where(hmp4['Ta'].isnull())[0]]=0
        
            # Write in data

            nc.variables['air_temperature'][:,0]=t1.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,0]=qc1
            nc.variables['height_above_snow_surface'][:,0]=alt_HMP1.to_numpy()
        
            nc.variables['air_temperature'][:,1]=t2.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,1]=qc2
            nc.variables['height_above_snow_surface'][:,1]=alt_HMP2.to_numpy()
        
            nc.variables['air_temperature'][:,2]=t3.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,2]=qc3
            nc.variables['height_above_snow_surface'][:,2]=alt_HMP3.to_numpy()
        
            nc.variables['air_temperature'][:,3]=t4.to_numpy()
            nc.variables['qc_flag_surface_temperature'][:,3]=qc4
            nc.variables['height_above_snow_surface'][:,3]=alt_HMP4.to_numpy()
            
            # Derive valid max and min
            nc.variables['air_temperature'].valid_min = np.nanmin([t1.min(), t2.min(),t3.min(),t4.min()])
            nc.variables['air_temperature'].valid_max = np.nanmax([t1.max(), t2.max(),t3.max(),t4.max()])
            
            nc.variables['height_above_snow_surface'].valid_min = np.nanmin([alt_HMP1.min(), alt_HMP2.min(),alt_HMP3.min(),alt_HMP4.min()])
            nc.variables['height_above_snow_surface'].valid_max= np.nanmax([alt_HMP1.max(), alt_HMP2.max(),alt_HMP3.max(),alt_HMP4.max()])          

            # Write note to netcdf file indicating date used.
            if start < dt.datetime(2021,7,1):
                #base_str = 'Platform altitude is the top of the Met tower, sensor height is platform altitude minus 12.3 m until 2020-07-07 1648Z, after this date sensor height is platform altitude minus 10.83 m'
                base_str = 'Platform altitude (h0) is the top of the Met tower. Instrument altitude: HMP1=h0-10.87m, HMP2=h0-9.8m, HMP3=h0-4.5m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4].'
            elif start < dt.datetime(2022,6,1):                
                base_str = 'Platform height (h0) is the top of the Met tower. HMP1 was offline between 1 July 2021 and 22 August 2022. Instrument height: HMP2=h0-9.8m, HMP3=h0-4.5m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4]'
            elif start< dt.datetime(2022,7,1):
                base_str = 'Platform height (h0) is the top of the Met tower. HMP1 was offline between 1 July 2021 and 22 August 2022. Instrument height: HMP2=h0-9.8m until 2022-06-06 10:00, after which it was raised to h0-9.2, HMP3=h0-4.5m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4]'
            elif start < dt.datetime(2022,8,1):
                base_str = 'Platform height (h0) is the top of the Met tower. HMP1 was offline between 1 July 2021 and 22 August 2022. Instrument height: HMP2=h0-9.2, HMP3=h0-4.5m until 2022-07-22 1400,after which it was raised to h0-3.9m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4]'
            elif start < dt.datetime(2022,9,1):
                base_str = 'Platform height (h0) is the top of the Met tower. HMP1 was offline between 1 July 2021 and 22 August 2022. Instrument height: HMP1=h0-9.2 (installed 22 August 2022), HMP2=h0-9.2 until 22 August 2022, after which it was raised to h0-7.4m, HMP3=h0-3.9m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4]'
            elif start<dt.datetime(2023,7,1):
                base_str = 'Platform height (h0) is the top of the Met tower. Instrument height: HMP1=h0-9.2, HMP2=h0-7.4m, HMP3=h0-3.9m, HMP4=h0-1m , Index: [HMP1, HMP2, HMP3, HMP4]'
            else:
                base_str='Index: [HMP1, HMP2, HMP3, HMP4]'
            nc.setncattr('comment', base_str)
                

            nc.setncattr('comment', base_str)
 
    
            # Close netcdf file
    
            nc.close()
    
if __name__ == '__main__':
    main()
