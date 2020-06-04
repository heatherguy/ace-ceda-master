#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs March 12 17:22:21 2020

@author: Heather Guy

Various functions for parsing data from the ICECAPS-ACE fluxtower suite.
"""

# Import functions
import warnings
warnings.filterwarnings("ignore")
import numpy as np       
import glob
import os
import io
import datetime
import pandas as pd
from scipy.signal import medfilt


def replace_outliers(var,sd):
    """
    Filter for clear outliers - replace with median filtered values
    Limit is 3*standard deviation
    
    Parameters:
        var: input variable
        sd:  standard deviation
        
    Returns:
        var with median filter applied.
    
    """
    var=var.astype(float)
    jj = ~np.isnan(var) # Ignore nans
    temp = var[jj]
    mf = medfilt(temp,11)             # Get median filter
    ii = np.abs(temp - mf) > 3*sd     # Get outliers, where greated than 3 Sd's from median
    temp[ii] = mf[ii]                 # Replace these outliers with the median
    var_clean = var
    var_clean[jj] = temp      # Put back into orginal array/
    
    return var_clean



def extract_KT_data(start,stop,dpath,qcf,save=False):
    """
    Extracts KT15 data from raw output. 
    Resamples to 1 minute averaging time. 
    QC's for bad data from log.
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        qcf:   Filepath for qc log file
        save:  Output directory (optional)
        
    Returns:
        Clean KT15 dataframe
    
    """
    # Raw data format: YYYY MM DD HH MM.mmm TT.tt C
    # TT.tt = temperature, C = celcius

    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*.KT15')  # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]

    # Initialise empty data frames
    KT = pd.DataFrame()

    # Extract the data
    for f in dfs: 
        # Ignore file if it's empty or contains non-ascii characters
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue

        # Filter and report files with non-ascii characters
        content = open(f).read()
        try:
            content.encode('ascii')
        except UnicodeDecodeError:
            print("Error with: %s contains non-ascii characters.\n"%f)  
            continue
        
        # Store good data  
        KT = KT.append(pd.read_csv(f, header=None, delim_whitespace=True,error_bad_lines=False))

    # Sort out the date referencing and columns
    if KT.empty==False:
        try:
            KT[5] = KT[5].astype(int)
            KT['Date'] = pd.to_datetime(KT[0]*10000000000+KT[1]*100000000+KT[2]*1000000+KT[3]*10000+KT[4]*100+KT[5],format='%Y%m%d%H%M%S')
            KT = KT.set_index('Date')
            del KT[0],KT[1],KT[2],KT[3],KT[4],KT[5]
            if np.shape(KT)[1]==2:
                KT.columns = ['T', 'Units']
                KT = KT.sort_values('Date')
                KT.index = pd.DatetimeIndex(KT.index)
                KT = KT[~KT.index.duplicated()]
            
                # Resample to 1 minute averages. 
                new_index = pd.date_range(KT.index[0].round('min'),KT.index[-1].round('min') , freq='min') 
                    
                KT_1min = KT.resample('1min').mean()
                KT_1min = KT_1min.reindex(new_index)
    
                # Crop to datetime
                KT=KT_1min[start:stop]
    
                # QC 
                KT['QC']=np.ones(len(KT))
                bad_times = pd.read_csv(qcf,parse_dates={'start_dates':[0],'stop_dates':[1]},header=None)
                sdate = pd.Timestamp(KT.index[0].date())
        
                # See if there are any bad dates in this file
                if (bad_times['start_dates'].dt.date==sdate).any():
                    # If yes, set flags during flight times to zero.
                    subset = bad_times[bad_times['start_dates'].dt.date==sdate]
                    for i in range(0,len(subset)):
                        start_date = (pd.to_datetime(subset['start_dates'].iloc[i])).date()
                        stop_date = (pd.to_datetime(subset['stop_dates'].iloc[i])).date()
                        start_time = (pd.to_datetime(subset['start_dates'].iloc[i])).time()
                        stop_time = (pd.to_datetime(subset['stop_dates'].iloc[i])).time()
                
                        if start_date == stop_date:
                            KT['QC'][KT.between_time(start_time,stop_time).index]=0
                        else:
                            KT['QC'][KT.between_time(start_time,datetime.time(23,59,59)).index]=0
     
                if save: 
                    KT.to_csv(save+'KT_%s'%(str(start.date())))
            else:
                print('Bad data file')
        except:
            print('Bad data file')
    else:
        print('No KT data found for this period')

    return KT




def extract_snd_data(start,stop,dpath,hmp_dpath,save=False):
    """
    Extracts SnD data from raw output.
    Applies temperature correction. 
    QC's for bad data.
    Saves as .csv if requested
    
    Parameters:
        start:       Start datetime for processing
        stop:        Stop datetime for processing
        dpath:       Raw data filepath
        hmp_dpath:   Filepath for clean HMP1 data for temperature
        save:  Output directory (optional)
        
    Returns:
        Clean SnD dataframe
    
    """
    
    # Raw data format
    #aa;D.DDD;QQQ; VVVVV;CC
    #aa = 33 (serial address of sensor)
    #D.DDD = Distance to target in m (will need temperature adjustment
    #QQQ = Data quality, varies beteen 152-600, 600 is the poorest quality
    #VVVVV = diagnostic tests (only first two are actually something), 1 = pass. 
    #CC = two-character checksum of data packet.
    
    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*.SnD')   # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]

    # Initialise empty data frames
    snd = pd.DataFrame()

    # Extract the data
    for f in dfs: 
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue
        
        try:
            dat_lines =open(f).readlines()
        except:
            print('Data error with: '+f+' \n')
            continue
        good_lines = [y for y in dat_lines if len(y)==48]
        date =[]
        x=[]
        Q=[]
        V=[]
        C=[]
        for i in range(0,len(good_lines)):   
            date.append((pd.to_datetime(good_lines[i][0:23],format='%Y %m %d %H %M %S.%f')).round('min'))
            Q.append(good_lines[i][35:38])
            V.append(good_lines[i][39:41])
            C.append(good_lines[i][45:47])
            try: 
                x.append(float(good_lines[i][29:34]))
            except ValueError:
                x.append(np.nan)

        snd = snd.append(pd.DataFrame({'Date':date,'depth':x,'Q':Q,'V':V,'C':C}))
        
    
    if snd.empty==False:
        snd = snd.set_index('Date')
        snd = snd.sort_values('Date')
        snd.index = pd.DatetimeIndex(snd.index)
        snd = snd[~snd.index.duplicated()]

        # Check diagnostic tests pass.
        snd['depth'][snd['V'].astype(int)!=11]=np.nan

        # Check data quality
        snd['Q'] = snd['Q'].astype(int)
        snd['depth'][snd['Q']<151]=np.nan
        snd['depth']=replace_outliers(snd['depth'],snd['depth'].std())
        
        # Correct for temperature
        HMP1 = get_hmp(start,stop,hmp_dpath,'HMP1')
        tc = HMP1['Ta']
        tk = tc + 273.15 # convert to celcius    
        snd['HMP1_T'] = tk
        snd['depth_Tcorrected'] = snd['depth'] * np.sqrt(snd['HMP1_T']/273.15)

        # Check for crazy values
        snd['depth_Tcorrected']=replace_outliers(snd['depth_Tcorrected'],snd['depth_Tcorrected'].std())
       
        if save: 
            snd.to_csv(save+'snd_%s'%(str(start.date())))
            
    else: 
        print('No data from snd\n')
       

    return snd




def HMP_pdf_sort(df,start,stop):
    """
    Sorts date referencing and columns for raw HMP data.
    Called by extract_HMP_data
           
    Returns:
        Clean df
    
    """
    # Sort out the date referencing and columns
    if df.empty==False:
        df = df.dropna()
        df['Second'] = df['Second'].astype(float)
        df['Second'] = df['Second'].astype(int)
        df['Minute'] = df['Minute'].astype(int)
        df['Hour'] = df['Hour'].astype(int)
        df['Day'] = df['Day'].astype(int)
        df['Month'] = df['Month'].astype(int)
        df['Year'] = df['Year'].astype(int)
        df['Date'] = pd.to_datetime(df['Year']*10000000000+df['Month']*100000000+df['Day']*1000000+df['Hour']*10000+df['Minute']*100+df['Second'],format='%Y%m%d%H%M%S')
        df = df.set_index('Date')
        del df['Year'],df['Month'],df['Day'],df['Hour'],df['Minute'],df['Second'],df['junk']
        df.columns = ['RH', 'Ta', 'Tw', 'Err', 'h']
        df['RH']=df['RH'].astype(float)
        df['Tw']=df['Tw'].astype(float)
        df['Ta']=df['Ta'].astype(float)
        df = df.sort_values('Date')
        df.index = pd.DatetimeIndex(df.index)
        df = df[~df.index.duplicated()]
    else:
        df = pd.DataFrame(columns=['RH', 'Ta', 'Tw', 'Err', 'h'])
    return df

def extract_HMP_data(name,start,stop,dpath,save=False):
    """
    Extracts HMP155 data from raw output. 
    Saves as .csv if requested
    
    Parameters:
        name:  Sensor name ('HMP1', 'HMP2', 'HMP3' or 'HMP4')
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
        
    Returns:
        Clean HMP dataframe
    
    """
    # Raw data format: YYYY MM DD HH MM.mmm TT:TT:TT RH Ta Tw Err hs
    # Ta = seperate probe T, Tw = wetbulb t, hs = heating status

    os.chdir(dpath)                     # Change directory to where the data is
    all_files = glob.glob('*.%s'%name)  # List all data files
    
    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]

    # Initialise empty data frames
    HMP = pd.DataFrame()

    # Extract the data
    for f in dfs:
        # Ignore file if it's empty 
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.')
            continue
        
        fd = io.open(f,"r",errors='replace')
        f_dat = fd.readlines()
        clean_dat = [i for i in f_dat if len(i)>=60 and len(i)<=63]
        pdf = pd.DataFrame(clean_dat)
        pdf[1] = pdf[0].str.split()
        final_df = pd.DataFrame(pdf[1].values.tolist(), columns=['Year','Month','Day','Hour','Minute','Second','junk','RH','Ta', 'Tw', 'Err', 'h'])
        # Store good data
        HMP = HMP.append(final_df)
     
    # Sort out the date referencing and columns
    HMP = HMP_pdf_sort(HMP,start,stop)

    if save: 
        HMP.to_csv(save+'%s_%s'%(name,str(start.date())))

    return HMP



def ventus_pdf_sort(df,start,stop):
    """
    Sorts date referencing and columns for raw 2D sonic data.
    Called by extract_ventus_data
           
    Returns:
        Clean df
    
    """
    # Sort out the date referencing and columns
    if df.empty==False:
        try:
            df[5] = df[5].astype(int)
            df['Date'] = pd.to_datetime(df[0]*10000000000+df[1]*100000000+df[2]*1000000+df[3]*10000+df[4]*100+df[5],format='%Y%m%d%H%M%S')
            df = df.set_index('Date')
            del df[0],df[1],df[2],df[3],df[4],df[5]
            df.columns = ['wsd', 'wdir', 'T', 'Checksum']
            df = df.sort_values('Date')
            new_idx = pd.date_range(start,stop,freq='1s' )
            df.index = pd.DatetimeIndex(df.index)
            df = df[~df.index.duplicated()]
            df= df.reindex(new_idx, fill_value=np.NaN)
        except:
            df = pd.DataFrame()
            print('Bad data')
    return df


def extract_ventus_data(start,stop,dpath,save=False):
    """
    Extracts 2D sonic data from raw output. 
    QC's for bad data. 
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
        
    Returns:
        v1, v2: Clean dataframes of lower 2D sonic (v1) and upper (v2)
    
    """
    # Raw Data format:
    #<STX>SS.S DDD +TT.T xx*XX<CR><ETX>
    #        SS.S = wind speed (m/s)
    #        DDD = wind direction
    #        +TT.T = signed virtual temperature
    #        xx = status
    #        XX = checksum

    os.chdir(dpath)                     # Change directory to where the data is
    all_files = glob.glob('*.ventus*')  # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]

    # Initialise empty data frames
    v1 = pd.DataFrame()
    v2 = pd.DataFrame()

    # Extract the data
    for f in dfs: 
        # Ignore file if it's empty or contains non-ascii characters
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.')
            continue
        # Filter and report files with non-ascii characters
        try:
            content = open(f).read()
        except UnicodeDecodeError:
            print("Error with: %s contains non-ascii characters."%f)  
            continue
        try:
            content.encode('ascii')
        except UnicodeDecodeError:
            print("Error with: %s contains non-ascii characters."%f)  
            continue
            
        if f[-1]=='1':
            v1 = v1.append(pd.read_csv(f, header=None, delim_whitespace=True, error_bad_lines=False))
        if f[-1]=='2':
            v2 = v2.append(pd.read_csv(f, header=None, delim_whitespace=True, error_bad_lines=False))
        
    # Sort out the date referencing and columns
    v1 = ventus_pdf_sort(v1,start,stop)
    v2 = ventus_pdf_sort(v2,start,stop)
    
    # Convert to numeric
    if v1.empty==False:
        v1['wsd'] = pd.to_numeric(v1['wsd'].str.split(pat='\x02',expand=True)[1],errors='coerce')
        v1['wdir'] = pd.to_numeric(v1['wdir'],errors='coerce')
        v1['T'] = pd.to_numeric(v1['T'],errors='coerce')

        # Add qc column
        v1['QC'] = np.ones(len(v1))

        # Check for bad data
        # Check checksum, status should be 08 (heating on) or 00
        v1['status']= v1['Checksum'].str.split(pat='*',expand=True)[0]
        v1['QC'][(v1['status']!='08') & (v1['status']!='00')]=2
        v1['QC'][v1.isnull().any(axis=1)]=0
    else:
        print('No data for v1')
    
        # Convert to numeric
    if v2.empty==False:
        v2['wsd'] = pd.to_numeric(v2['wsd'].str.split(pat='\x02',expand=True)[1],errors='coerce')
        v2['wdir'] = pd.to_numeric(v2['wdir'],errors='coerce')
        v2['T'] = pd.to_numeric(v2['T'],errors='coerce')

        # Add qc column
        v2['QC'] = np.ones(len(v2))

        # Check for bad data
        # Check checksum, status should be 08 (heating on) or 00
        v2['status']= v2['Checksum'].str.split(pat='*',expand=True)[0]
        v2['QC'][(v2['status']!='08') & (v2['status']!='00')]=2
        v2['QC'][v2.isnull().any(axis=1)]=0
    else:
        print('No data for v2')

        
    #print('Data parse finished\n')
    if save: 
        v1.to_csv(save+'v1_%s'%(str(start.date())))
        v2.to_csv(save+'v2_%s'%(str(start.date())))
        #print('Saved csv')

    return v1,v2




def extract_metek_data(start,stop,dpath,save=False):
    """
    Extracts 3D sonic data from raw output. 
    Cleans bad data.
    Converts to SI units.
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
        
    Returns:
        m1, m2: Clean dataframe for lower 3D sonic (m1) and upper (m2)
    
    """
    # raw data format:
    #2019 04 02 16 23 41.734 M:x =    14 y =    -1 z =    12 t =  2357
    # M = measured data heater off
    # H = measured data heater on
    # D = measured data heater defect
    # x,y,z componenets of wind in cm/s
    # t = acoustic temperature 2357 = 23.57C
    
    os.chdir(dpath)                    # Change directory to where the data is
    all_files = glob.glob('*.metek*')  # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]

    # Initialise empty data frames
    m1 = pd.DataFrame()
    m2 = pd.DataFrame()

    # Get data.
    for f in dfs: 
        
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue
        
        # Filter out bad data files
        try:
            fd = io.open(f,"r",errors='replace')
            f_dat = fd.readlines()
        except:
            print('Data error with %s'%f)
            continue   

        # Iterate through all data lines. 
        # Delete duplicated error lines (E)
        # Replace bad data lines with nans. 
        # Iterate backwards so you can delete in place.
        # Replace any lines of non standard length with nan line.

        nan_line = 'x =   nan y =   nan z =   nan t = nan  \n'
        for i in range(len(f_dat) -1 , -1, -1):
            if len(f_dat[i])!= 66:
                f_dat[i] = f_dat[i][0:26] + nan_line
            if f_dat[i][24]=='D':
                f_dat[i] = f_dat[i][0:26] + nan_line
            elif f_dat[i][24]=='E':
                f_dat[i-1] = f_dat[i-1][0:26] + nan_line
                del f_dat[i]

        # Filter out incomplete data files
        if len(f_dat) < 36000:
            print('Incomplete file %s'%f)
            continue
        
        # Ignore extra data points
        f_dat = f_dat[0:36000]
        
        # Now split up strings and put in pandas df.
        pdf = pd.DataFrame(f_dat)
        pdf[1] = pdf[0].str.split()
        df = pd.DataFrame(pdf[1].values.tolist(), columns=['Year','Month','Day','Hour','Minute','Second','Status','junk1','x','junk2','junk3','y','junk4','junk5','z','junk6','junk7','T'])
        df['Logger_Date'] = pd.to_datetime(df[['Year','Month','Day','Hour','Minute','Second']])
        del df['junk1'],df['junk2'],df['junk3'],df['junk4'],df['junk5'],df['junk6'],df['junk7'],df['Year'],df['Month'],df['Day'],df['Hour'],df['Minute'],df['Second']

        # Make a standard 10Hz time series. 
        # Start time
        st = pd.to_datetime(f[0:16],format='%y%m%d_%H%M%S.%f')
        # End time
        et = st + pd.Timedelta(hours=1)
        et = et - pd.Timedelta(seconds=0.1)
        # 10 Hz time series (36000 data points in one hour)
        df['Date'] = pd.date_range(st,et,periods=36000)
        
        # Tidy and sort units.
        df = df.set_index('Date')
        df['T'] = pd.to_numeric(df['T'], errors='coerce')
        df['T']=df['T'].astype(float)
        df['x']=df['x'].astype(float)
        df['y']=df['y'].astype(float)
        df['z']=df['z'].astype(float)
        df['T']=df['T']/100.0

        df = df.sort_values('Date')
        df.index = pd.DatetimeIndex(df.index)
        #df = df[~df.index.duplicated()]
               
        # Change units from cm/s to m/s, (* 0.01), and T to kelvin
        df['x']=df['x']*0.01
        df['y']=df['y']*0.01
        df['z']=df['z']*0.01
        df['T']= df['T']+ 273.15

        # Append to relavent dataframe
        if f[-1]=='1':
            try:              
                m1 = m1.append(df)
            except:
                print('Data error with %s'%f)
                continue
 
        if f[-1]=='2':
            try:              
                m2 = m2.append(df)
            except:
                print('Data error with %s'%f)
                continue   
            
    # crop data for date/time
    m1=m1[start:stop]
    m2=m2[start:stop]
        
    if save: 
        m1.to_csv(save+'metek1_%s'%str(start.date()))
        m2.to_csv(save+'metek2_%s'%str(start.date()))
   
    return m1,m2  




def extract_licor_data(start,stop,dpath,save=False):
    """
    Extracts licor data from raw output. 
    Corrects for data transmission delay. 
    QC's licor data. 
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
        
    Returns:
        Clean licor dataframe
    
    """
    # Raw data format: 
    #2019 04 03 11 11 56.453 89	189	0.16469	35.4518	0.04404	297.105	20.74	99.0	1.5224
    # Date, Ndx, DiagVal, CO2R, CO2D, H2OR, H2OD, T, P, cooler

    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*.licor')  # List all data files

    # Get start and stop filenames
    start_f = int(start.strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))

    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]

    # Initialise empty data frames
    licor = pd.DataFrame()

    # Function to convert to float otherwise write nan
    def convert_float(x):
        try:
            return np.float(x)
        except:
            return np.nan

    # Extract the data
    for f in dfs: 

        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue 
        
        try:
            li = pd.read_csv(f, header=None, sep='\t',error_bad_lines=False)
        except:
            print('Data Error with: '+f+'\n')
            continue        

        # Filter out incomplete data files
        if len(li) < 36000:
            print('Incomplete file %s'%f)
            continue
        
        # Ignore extra data points
        li = li[0:36000]
 
        # Sort out the date referencing
        li['Logger_Date']=li[0].str[0:23]
        li = li[li['Logger_Date'].map(len) ==23]
        li['Logger_Date'] = pd.to_datetime(li['Logger_Date'],format='%Y %m %d %H %M %S.%f')

        # Make a standard 10Hz time series. 

        # DONT FORGET, LICOR DATA TRANSMISSION IS DELAYED BY 200MS (TWO DATA POINTS)
        # EACH DATA POINT SHOULD BE SHIFTED BACK IN TIME BY THIS AMOUNT. ••

        # Start time
        st = pd.to_datetime(f[0:16],format='%y%m%d_%H%M%S.%f') - pd.Timedelta(seconds=0.2) # •• Implemented here. 
        # End time
        et = st + pd.Timedelta(hours=1)
        et = et - pd.Timedelta(seconds=0.1)
        # 10 Hz time series (36000 data points in one hour)
        li['Date'] = pd.date_range(st,et,periods=36000)

        # Sor the rest of the columns
        li['Ndx']=li[0].str[24:].astype('int')
        li['DiagV']=li[1].astype('int')
        li['CO2R'] = li[2]
        li['CO2R']=li['CO2R'].apply(convert_float)
        li['CO2D'] = li[3]
        li['CO2D']=li['CO2D'].apply(convert_float)/1000 # mol/m3
        li['H2OR'] = li[4]
        li['H2OR']=li['H2OR'].apply(convert_float)
        li['H2OD'] = li[5]
        li['H2OD']=li['H2OD'].apply(convert_float)/1000 # mol/m3
        li['T'] = li[6].astype('float')+273.15      # K
        li['P'] = li[7].astype('float')*1000        # Pa
        li['cooler'] = li[8].astype('float')
        del li[0],li[1],li[2],li[3],li[4],li[5],li[6],li[7],li[8],li[9]                 
        li = li.sort_values('Date')
        li = li.set_index('Date')
         
        licor = licor.append(li)
    
 
    # OC licor data
    # Interp diagnostic bitmap
    #The cell diagnostic value (diag) is a 1 byte unsigned integer (value between 0 and 255) with the
    #following bit map:
    #bit 7 bit 6 bit 5 bit 4 bit 3 bit 2 bit 1 bit 0
    #Chopper Detector PLL Sync <------------------------ ----AGC / 6.25 ------------------------>
    #1=ok 1=ok 1=ok 1=ok
    #Example: a value is 125 (01111101) indicates Chopper not ok, and AGC = 81% (1101 is 13,
    #times 6.25)
    try: 
        #   For now 1 is good, zero is bad, 2=no met data
        licor['QC']=np.ones(len(licor))
        diag_list = licor['DiagV'].to_list()
        chopper = []
        detector = []
        pll = []
        sync = []
        agc = []
        for i in range(0,len(diag_list)):
            binar = format(diag_list[i], "008b")
            chopper.append(int(binar[0]))
            detector.append(int(binar[1]))
            pll.append(int(binar[2]))
            sync.append(int(binar[3]))
            agc_temp = binar[4:]
            agc.append(int(agc_temp, 3)* 6.25)
    
        licor['chopper']=chopper
        licor['detector']=detector
        licor['pll']=pll
        licor['sync']=sync
        licor['agc']=agc

        # Filter out bad values
        #licor['QC'][licor['chopper']==0] = 0 
        #licor['QC'][licor['detector']==0] = 0
        licor['QC'][licor['pll']==0] = 0
        licor['QC'][licor['sync']==0] = 0
        licor['QC'][licor['agc']>90] = 0
        
        # Delete now useless params
        #del licor['chopper'],licor['detector'],licor['pll'],licor['sync']

        # Clean pressure for values over 80000 or under 50000
        # Clean temperature for values under 200 or over 300
        licor['QC'][licor['P']<50000]=0
        licor['QC'][licor['P']>80000]=0
        licor['QC'][licor['T']<200]=0
        licor['QC'][licor['T']>300]=0

        # Filter clear outliers.
        jj = ~np.isnan(licor['H2OD']) # Not nan indices
        sd = np.std(licor['H2OD'][jj]) # standard deviation 
        licor['H2OD']=replace_outliers(licor['H2OD'],sd)
        
        # Remove negative concentrations
        licor['QC'][licor['H2OD']<0]=0
        licor['QC'][licor['H2OR']<0]=0

        if save: 
            licor.to_csv(save+'licor_%s'%str(start.date()))
    
    except:
        print('QC error. No data. \n')

    return licor    

def get_kt(start,stop,d_loc):
    """
    Retrieves KT15 from .csv output of extract_kt_data
    
    Parameters:
        start:      Start datetime for processing
        stop:       Stop datetime for processing
        d_loc:      data filepath
        
    Returns:
        KT dataframe
    
    """
    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    KT_out = pd.DataFrame(columns=['T','Units','QC'])
    for date in f_date_list:
        f = d_loc + r'KT_%s'%(str(date.date()))
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0])
        except:
            print('No data for %s'%str(date.date()))
            continue

        KT_out = KT_out.append(data,sort=True)    
    
    # Get rid of any duplicates
    KT_out = KT_out[~KT_out.index.duplicated()]   
    # Crop to datetime
    KT_out=KT_out[start:stop]
    
    # Fill any missing minutes with nans
    if len(KT_out)!=0:
        new_index = pd.date_range(KT_out.index[0],KT_out.index[-1] , freq='min')
        KT_out = KT_out.reindex(new_index)
    else:
        new_index = pd.date_range(start,stop, freq='min')
        KT_out = pd.DataFrame(columns=['T','Units','QC'],index=new_index)
    
    return KT_out


def get_hmp(start,stop,d_loc,hmp_name):
    """
    Retrieves HMP155 from .csv output of extract_hmp_data
    QC's hmp data
    resamples to one minutely averages.
    
    Parameters:
        start:      Start datetime for processing
        stop:       Stop datetime for processing
        d_loc:      data filepath
        hmp_name:   sensor name (HMP1, HMP2, HMP3, or HMP4)

        
    Returns:
        HMP dataframe
    
    """
    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    out = pd.DataFrame(columns=['Date','RH','Ta','Tw','Err','h'])
    out2 = pd.DataFrame(columns=['RH','Ta','QC'])
    for date in f_date_list:
        f = d_loc + r'%s_%s'%(hmp_name,str(date.date()))
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0])
        except:
            print('No data for %s'%str(date.date()))
            continue

        out = out.append(data,sort=True)    
    
    # Get rid of any duplicates
    out = out[~out.index.duplicated()]   
    # Crop to datetime
    out=out[start:stop]
    qc = np.ones(len(out))
    qc[out['Err']!=0]= 2
    qc[out['Ta'].isnull()] = 0
    out['QC']=qc
    
    # Resample to one minutely
    if len(out)!=0:
        out2['Ta'] = out['Ta'].resample(rule = '1min').mean()
        out2['RH'] = out['RH'].resample(rule = '1min').mean()
        out2['QC'] = out['QC'].resample(rule = '1min').max()
        # Fill any missing minutes with nans
        new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq='1min')
        out2 = out2.reindex(new_index)
    
        return out2
    else:
        return out


def get_ventus(start,stop,d_loc,name, avp='1min'):
    """
    Retrieves 2D sonic data from .csv output of extract_ventus_data
    Resamples to requested averaging time (default=1 minutely)
    
    Parameters:
        start:  Start datetime for processing
        stop:   Stop datetime for processing
        d_loc:  data filepath
        name:   sensor name, lower='v1', upper='v2'
        avp:    Resamples to avp averaging period in minutes (default = '1min')
        
    Returns:
        2D sonic dataframe
    
    """
    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    out = pd.DataFrame(columns=['Date','wsd','wdir','T','Checksum','QC','status'])
    out2 = pd.DataFrame(columns=['wsd','wdir','T','QC'])    
    for date in f_date_list:
        f = d_loc + r'%s_%s'%(name,str(date.date()))
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0])
        except:
            print('No data for %s'%str(date.date()))
            continue

        out = out.append(data,sort=True)    
    
    # Get rid of any duplicates
    out = out[~out.index.duplicated()]   
    # Crop to datetime
    out=out[start:stop]
    out.index = pd.DatetimeIndex(out.index)    

    # Resample to one minutely
    if len(out)!=0:
        out2['T'] = out['T'].resample(rule = '1min').mean()   
        out2['QC'] = out['QC'].resample(rule = '1min').max()
    
        # Convert to vectors
        u,v = wind_to_uv(out['wsd'],out['wdir'])
    
        # average vectors
        u_mean = u.resample(rule=avp).mean()
        v_mean = v.resample(rule=avp).mean()
    
        # convert back to wsd and wdir
        out2['wsd'] = wind_uv_to_spd(u_mean,v_mean)
        out2['wdir'] = wind_uv_to_dir(u_mean,v_mean)


        # Fill any missing minutes with nans
        new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq=avp)
        out2 = out2.reindex(new_index)
    
        return out2
    else:
        new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq=avp)
        out2 = pd.DataFrame(columns=['wsd','wdir','T','QC'],index=new_index)   
        
        return out2



def get_metek(start,stop,d_loc,name, avp='1min'):
    """
    Retrieves 3D sonic data from .csv output of extract_metek_data
    Resamples to requested averaging time (default=1 minutely)
    
    Parameters:
        start:  Start datetime for processing
        stop:   Stop datetime for processing
        d_loc:  data filepath
        name:   sensor name, lower='vm', upper='m2'
        avp:    Resamples to avp averaging period in minutes (default = '1min')
        
    Returns:
        3D sonic dataframe
    
    """

    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    out = pd.DataFrame(columns=['Date','Status','x','y','z','T','Logger_Date'])
    out2 = pd.DataFrame(columns=['wsd','wdir','T','QC'])    
    for date in f_date_list:
        f = d_loc + r'%s_%s'%(name,str(date.date()))
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0])
        except:
            print('No data for %s'%str(date.date()))
            continue

        out = out.append(data,sort=True)    
    
    # Get rid of any duplicates
    out = out[~out.index.duplicated()]   
    # Crop to datetime
    out=out[start:stop]
    del out['Date']
    # QC
    out['QC']=np.ones(len(out))
    out['QC'][out.isnull().any(axis=1)]=0
    
    if len(out)!=0:
        # Resample to one minutely
        out2['T'] = out['T'].resample(rule = avp).mean()  
        out2['QC'] = out['QC'].resample(rule = avp).max() 
       
        # average vectors
        u_mean = out['x'].resample(rule=avp).mean() 
        v_mean = out['y'].resample(rule=avp).mean() 
    
        # convert to wsd and wdir
        out2['wsd'] = wind_uv_to_spd(u_mean,v_mean)
        out2['wdir'] = wind_uv_to_dir(u_mean,v_mean)

        # Fill any missing minutes with nans
        new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq=avp)
        out2 = out2.reindex(new_index)
    
        return out2
    else:
        new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq=avp)       
        out2 = pd.DataFrame(columns=['wsd','wdir','T','QC'],index=new_index)   
        return out2


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

