#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 21:39:27 2019

@author: heather

Some useful functions for ICECAPS-ACE data parsing. 
"""


import numpy as np
import pandas as pd
import os
import glob
import tarfile
import datetime as dt
# Supress warnings for sake of log file
import warnings
warnings.filterwarnings("ignore")


def extract_tar(start,stop, dloc,outdir,instr):
    """
    Extract files containting 'instr' from tar.gz files and save in outdir.
    
    Parameters:
        start: Start date (string in format YYMMDD)
        stop: Stop date (string in format YYMMDD)        
        dloc:   Directory to look in
        outdir: Directory to extract to
        instr:  Identifier string
    
    """
    day_range = pd.date_range(dt.datetime.strptime(start,'%y%m%d'),dt.datetime.strptime(stop,'%y%m%d'),freq='1D')
    fnames=[]
    for i in range(0,len(day_range)):
        fnames.append(glob.glob(dloc + r'%s*.tar.gz'%dt.datetime.strftime(day_range[i],'%y%m%d')))
        fnames.append(glob.glob(dloc + r'%s*.tar.gz'%dt.datetime.strftime(day_range[i],'%Y-%m-%d')))

    # flatten list
    fnames_flat = [item for sublist in fnames for item in sublist]
        
    for f in fnames_flat: 
        t = tarfile.open(f,'r')
        mems = t.getmembers()
        for m in mems:
            if instr in m.name:
                if m.isreg():  # skip if the TarInfo is not files
                    m.name = os.path.basename(m.name) # remove the path by reset it
                    t.extract(m,outdir)
        t.close()
        
    return
        

    
def get_NOAA_met(w_dloc):
    """
    Function to extract NOAA weather data from file. 
    
    Parameters:
        w_dloc: input data filepath
        
    Returns: 
        pdf: Dataframe of NOAA weather data
    
    """
    #Fields: Wind D, Wind s (m/s), Wind steadiness, pressure (hPa), 2 m T, 10 m T, tower T, RH, P (mm/hr)
    # Missing values: -999, -999.9, -9, -999.90, -999, -999, -999, -99, -99
    all_dates = []
    all_data =  []
    f = open(w_dloc,mode='r')
    data = f.read().split('\n')
    for i in range(0,len(data)):
        if data[i][4:20]=='':
            continue       
        else:
            try:
                all_dates.append(pd.to_datetime(data[i][4:20],format='%Y %m %d %H %M'))
                freq='Min'
            except:
                all_dates.append(pd.to_datetime(data[i][4:20],format='%Y %m %d %H'))
                freq='H'

            all_data.append(list(map(float, data[i][20:].split())))

    wd = [int(x[0]) for x in all_data]
    ws = [x[1] for x in all_data]
    p = [x[3] for x in all_data]
    T = [x[4] for x in all_data]
    RH = [int(x[7]) for x in all_data]

    # Mask missing values

    wd = np.ma.masked_where(np.asarray(wd)==-999, np.asarray(wd))
    ws = np.ma.masked_where(np.asarray(ws)==-999.9, np.asarray(ws))
    p = np.ma.masked_where(np.asarray(p)==-999.9, np.asarray(p))
    T = np.ma.masked_where(np.asarray(T)==-999.9, np.asarray(T))
    RH = np.ma.masked_where(np.asarray(RH)==-99, np.asarray(RH))

    # Create pandas dataframe
    d = {'date' : all_dates,'wd' : wd, 'ws' : ws, 'pres' : p, 'T' : T, 'RH' : RH}
    pdf = pd.DataFrame(d)
    pdf['date'] = pd.to_datetime(pdf.date,utc=True)
    pdf = pdf.sort_values(by='date')
    pdf = pdf.set_index('date')

    # Delete duplicates and fill any missing data with blank lines   
    #sd = min(all_dates)
    #ed = max(all_dates)
    sd = pdf.index[0]
    ed = pdf.index[-1]  
    
    d_list = pd.date_range(sd, ed,freq=freq)
    pdf = pdf[~pdf.index.duplicated(keep='first')]
    pdf = pdf.reindex(d_list)
    pdf.index = pdf.index.tz_localize(None)
    
    return pdf


def wind_to_uv(wspd,wdir):
    """
    calculated the u and v wind components from wind speed and direction
    Input:
        wspd: wind speed
        wdir: wind direction
    Output:
        u: u wind component
        v: v wind component
    """    
   
    rad = 4.0*np.arctan(1)/180.
    u = -wspd*np.sin(rad*wdir)
    v = -wspd*np.cos(rad*wdir)

    return u,v
        
def wind_uv_to_dir(U,V):
    """
    Calculates the wind direction from the u and v component of wind.
    Takes into account the wind direction coordinates is different than the 
    trig unit circle coordinate. If the wind directin is 360 then returns zero
    (by %360)
    Inputs:
      U = west/east direction (wind from the west is positive, from the east is negative)
      V = south/noth direction (wind from the south is positive, from the north is negative)
    """
    WDIR= (270-np.rad2deg(np.arctan2(V,U)))%360
    return WDIR
    
def wind_uv_to_spd(U,V):
    """
    Calculates the wind speed from the u and v wind components
    Inputs:
      U = west/east direction (wind from the west is positive, from the east is negative)
      V = south/noth direction (wind from the south is positive, from the north is negative)
    """
    WSPD = np.sqrt(np.square(U)+np.square(V))
    return WSPD
<<<<<<< HEAD
=======

def get_nc(dloc,name,var_list,lev_list,years,months):
    from netCDF4 import Dataset, MFDataset,num2date,date2num
    file_list=[]
    for year in years: 
        for month in months: 
            if name=='aerosol-concentration':
                file_list.append(dloc + '%s/ace-cpc_summit_%s%s_%s_v1.nc'%(name,year,str(month).zfill(2),name))
            else:
                file_list.append(dloc + '%s/ace-tower_summit_%s%s_%s_v1.nc'%(name,year,str(month).zfill(2),name))

    if len(file_list)>1:
        nc = MFDataset(file_list,'r',aggdim='time')
    elif len(file_list)==1:
        nc = Dataset(file_list[0],'r')
    else:
        print('No data')

    times = pd.to_datetime(nc.variables['time'][:],origin='unix',unit='s')

    all_vars=[]
    for i in range(0,len(var_list)):  
        if lev_list[i]==-1:
            all_vars.append(nc.variables[var_list[i]][:])
        else:
            all_vars.append(nc.variables[var_list[i]][:,lev_list[i]])
                             
    return times, all_vars
>>>>>>> 46686d08c360a41b32b4cc85ce9016dc9edc09e3
