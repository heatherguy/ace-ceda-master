#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 22:32:23 2019

@author: heather
"""

import sys,os
sys.path.append(os.getcwd())
import numpy as np
import datetime as dt
import pandas as pd
import os
import glob
from scipy import io
from utils import * 


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# QC aerosol data: Remove low winds, MSF winds, TAWO winds and flight days. 
def qc_aerosol(qc_in):
    """
    Function to QC aerosol data. 
    Flags low winds
    Flags winds coming from across station
    Flags flight days.
   
    Parameters:
        qc_in:      Data to qc
        
    Returns:
        dataframe with additional qc flag
        0 = bad data
        1 = good data
        2 = no met data available for wind flagging.
    
    """

    sdate = pd.Timestamp(qc_in.index[0].date())

    #   For now 1 is good, zero is bad, 2=no met data
    qc_in['QC']=np.ones(len(qc_in))

    # Get flight dates
    #flight_dates_f = '/Users/heather/Desktop/ace-ceda-master/qc-files/flight_days.csv'    
    flight_dates_f = '/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/qc-files/flight_days.csv'
    flight_times = pd.read_csv(flight_dates_f,parse_dates={'Dates':[0],'ondeckUTC':[0,3],'offdeckUTC':[0,4]})

    # See if there are any flight dates in this file
    if (flight_times['Dates']==sdate).any():
        # If yes, set flags during flight times to zero.
        subset = flight_times[flight_times['Dates']==pd.Timestamp(qc_in.index[0].date())]
        for i in range(0,len(subset)):
            try:
                start_time = (pd.to_datetime(subset['ondeckUTC'].iloc[i]) - pd.Timedelta(hours=1)).time()
            except:
                start_time = dt.time(0,0)
            try:
                stop_time = (pd.to_datetime(subset['offdeckUTC'].iloc[i]) + pd.Timedelta(hours=1)).time()
            except:
                stop_time = dt.time(23,59)
                
            qc_in['QC'][qc_in.between_time(start_time,stop_time).index]=0
           
    # Try to get Noaa Met data     
    w_dloc = '/gws/nopw/j04/ncas_radar_vol1/heather/Summit_Met/met_sum_insitu_1_obop_minute_%s_%s.txt'%(sdate.year,str(sdate.month))
    try:
        met = get_NOAA_met(w_dloc)
        # try the qc
        # Change qc direction when aerosol instruments moved to the AWO
        if sdate < dt.datetime(2022,7,1):
            qc_in['QC'][ace_winds['ws']<1]=0
            qc_in['QC'][ace_winds['wdir']>270]=0
        else: 
            # updated north winds criteria between 345 and 55 degrees
            qc_in['QC'][ace_winds['ws']<1]=0
            qc_in['QC'][ace_winds['wdir']>345]=0  
            qc_in['QC'][ace_winds['wdir']<55]=0
    except:
        try:
            # Try to get met data from ace netcdfs
            #nc_loc = '/Volumes/Data/ICECAPSarchive/ACE_netcdfs/'
            nc_loc = '/gws/nopw/j04/ncas_radar_vol1/heather/final_nc/'
        
            # Get 2m ws & direction
            var_list=['wind_speed','wind_from_direction','qc_flag']
            lev_list=[0,0,0]
            var_alts = ['height_above_snow_surface']
            wind_times,[ws,wdir,qc] = get_nc(nc_loc,'surface-winds-profile',var_list,lev_list,[sdate.year],[str(sdate.month)])
            # Put in dataframe and qc
            ace_winds = pd.DataFrame({'ws':ws,'wdir':wdir,'qc':qc},index=wind_times)
            ace_winds[ace_winds['qc']!=1]=np.nan

            # try the qc
            # Change qc direction when aerosol instruments moved to the AWO
            if sdate < dt.datetime(2022,7,1):
                qc_in['QC'][ace_winds['ws']<1]=0
                qc_in['QC'][ace_winds['wdir']>270]=0
            else: 
                # updated north winds criteria between 345 and 55 degrees
                qc_in['QC'][ace_winds['ws']<1]=0
                qc_in['QC'][ace_winds['wdir']>345]=0  
                qc_in['QC'][ace_winds['wdir']<55]=0
        
        except:
        
            # If no met data, QC=2
            print('No good met data')
            qc_in['QC']=2

    return qc_in



def extract_cpc(start,stop,dpath,save=False):
    """
    Extracts cpc data from raw output. 
    QC's for bad data. 
    Resamples to 1 minutely medians
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
    
    """
    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*CPC*')
    file_dates = np.asarray([(dt.datetime.strptime(f[-14:-4], '%Y-%m-%d')).date() for f in all_files]) 
    idxs = np.where(np.logical_and(file_dates>=start.date(), file_dates<=stop.date()))[0]
    dfs = [all_files[i] for i in idxs]
    for f in dfs: 
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue 
        cpc = pd.read_csv(f,sep=',',error_bad_lines=False,header=None,parse_dates={'Dates' : [0,1,2,3,4,5]})
        cpc.Dates = pd.to_datetime(cpc.Dates,format='%Y %m %d %H %M %S')
        cpc = cpc.sort_values('Dates')
        cpc = cpc.set_index(cpc['Dates'])
        cpc.index = pd.DatetimeIndex(cpc.index)
        del cpc['Dates']
        cpc = cpc.rename(columns={6:'c/cm3'})

        #  to minutly average
        new_index = pd.date_range(start,start + pd.Timedelta(minutes=(24*60)-1), freq='min')      
        cpc_1min = cpc.resample('1min').median()
        cpc_1min = cpc_1min.reindex(new_index)

        # QC
        cpc_qcd = qc_aerosol(cpc_1min)
        cpc_qcd['QC'][cpc_qcd['c/cm3']==0]=0
        cpc_qcd['QC'][cpc_qcd['c/cm3']==np.nan]=0        
        # Save if neccesary
        if save:
            cpc_qcd.to_csv(save+'CPC_%s'%(str(start.date())))

    return


def extract_skyopc(start,stop,dpath,save=False):
    """
    Extracts skyopc data from raw output. 
    QC's for bad data. 
    resamples to 1 minutely mean
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
    
    """
    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*skyopc')
    file_dates = np.asarray([(dt.datetime.strptime(f[0:9], '%y%m%d_%H')) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start, file_dates<=stop))[0]
    dfs = [all_files[i] for i in idxs]
    c=np.nan
    skyopc = pd.DataFrame()
    # Extract the data
    for f in dfs: 
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue 
        f_data = open(f)
        d = f_data.readlines()
        f_data.close()
        for i in range(0,len(d)):
            line=d[i].split()
            if line[0+6] =='P':
                if len(line)!=17+6:
                    c=0
                    datetime=np.nan
                    continue
                #Year Mon Day Hr Min Loc 4Tmp Err pA/p pR/p UeL Ue4 Ue3 Ue2 Ue1 Iv 
                datetime = dt.datetime(int(line[1+6])+2000,int(line[2+6]),int(line[3+6]),int(line[4+6]),int(line[5+6]))
                #datetime = dt.datetime.strptime('20'+line[1]+line[2]+line[3]+line[4]+line[5],'%Y%m%d%H%M')
                quad_Tmp = int(line[7+6])
                Err = int(line[8+6])
                pAp = int(line[9+6])
                pRp = int(line[10+6])
                Int = int(line[16+6])
                c=0
            
            elif len(line)!=9+6:
                continue
                
            elif c==0: 
                ch1=int(line[1+6])
                ch2=int(line[2+6])
                ch3=int(line[3+6])
                ch4=int(line[4+6])
                ch5=int(line[5+6])
                ch6=int(line[6+6])
                ch7=int(line[7+6])
                ch8=int(line[8+6])
                c = c+1    
            elif c ==1:
                ch9=int(line[1+6])
                ch10=int(line[2+6])
                ch11=int(line[3+6])
                ch12=int(line[4+6])
                ch13=int(line[5+6])
                ch14=int(line[6+6])
                ch15=int(line[7+6])
                ch16=int(line[8+6])
                c = c+1
            elif c == 2:
                ch17=int(line[1+6])
                ch18=int(line[2+6])
                ch19=int(line[3+6])
                ch20=int(line[4+6])
                ch21=int(line[5+6])
                ch22=int(line[6+6])
                ch23=int(line[7+6])
                ch24 =int(line[8+6])
                c= c+1
            elif c==3:
                ch25=int(line[1+6])
                ch26=int(line[2+6])
                ch27=int(line[3+6])
                ch28=int(line[4+6])
                ch29=int(line[5+6])
                ch30=int(line[6+6])
                ch31=int(line[7+6])
                ch32=int(line[8+6])
                c = 0
                n = int(line[0+6][-2])
                if isinstance(datetime,dt.datetime):
                    skyopc = skyopc.append(pd.Series([datetime+dt.timedelta(seconds=n*6), ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9, ch10, ch11, ch12, ch13, ch14, ch15, ch16, ch17, ch18, ch19, ch20, ch21, ch22, ch23, ch24, ch25, ch26, ch27, ch28, ch29, ch30, ch31, ch32, quad_Tmp,Err,pAp,pRp,Int]),ignore_index=True)

    #try: 
        # remove repeated channel 16
    del skyopc[16]
        # Correct counts for size bins 'all counts above lower threshold.'        
    for i in range(2,16):
        skyopc[i-1]=skyopc[i-1]-skyopc[i]
    for i in range(18,33):
        skyopc[i-1]=skyopc[i-1]-skyopc[i]
    
    skyopc=skyopc.rename(columns={0: 'Date',1:'ch1' ,2: 'ch2', 3: 'ch3',4: 'ch4',5: 'ch5',6: 'ch6',7: 'ch7',8: 'ch8',9: 'ch9',10: 'ch10',11: 'ch11',12: 'ch12',13: 'ch13',14: 'ch14',15: 'ch15',16: 'ch16',17: 'ch17',18: 'ch18',19: 'ch19', 20:'ch20',21: 'ch21',22: 'ch22',23: 'ch23',24: 'ch24',25: 'ch25',26: 'ch26',27: 'ch27',28: 'ch28',29: 'ch29',30: 'ch30',31: 'ch31',32: 'ch32',33: 'quad_Tmp',34:'Err',35:'pAp',36:'pRp',37:'Int'})
    skyopc.dropna(inplace=True)
    skyopc = skyopc.set_index('Date')
    skyopc = skyopc.sort_values('Date')
    skyopc.index = pd.DatetimeIndex(skyopc.index)
    skyopc = skyopc[~skyopc.index.duplicated()]
        
    # resample to 1 minute before calculating concentrations
    skyopc_counts = skyopc[skyopc.columns[0:31]]
    skyopc_counts = skyopc_counts.apply(pd.to_numeric, errors='coerce') # Counts in counts/ 6 seconds
    skyopc_params = skyopc[skyopc.columns[31:]]
    skyopc_counts =skyopc_counts / 100.0 # convert from counts/100ml to counts/cm3    

    # Resample to minutly average
    new_index = pd.date_range(pd.to_datetime(start).round('min'),pd.to_datetime(stop).round('min'), freq='min') 

    skyopc_1min = skyopc_counts.resample('1min').mean()
    skyopc_1min = skyopc_1min.reindex(new_index)
    params_1min = skyopc_params.resample('1min').mean()
    params_1min = params_1min.reindex(new_index)
    
    # QC for north winds/ flight days 
    #skyopc_qcd = qc_aerosol(skyopc_1min)
    skyopc_qcd=skyopc_1min.copy()
        
    # QC for params
    skyopc_qcd['qc_pap']=np.ones(len(skyopc_qcd))
    skyopc_qcd['qc_pap'][params_1min['pAp']>120]=0
        
    skyopc_qcd['qc_prp']=np.ones(len(skyopc_qcd))
    skyopc_qcd['qc_prp'][params_1min['pRp']>120]=0
        
    skyopc_qcd['qc_err']=np.ones(len(skyopc_qcd))
    skyopc_qcd['qc_err'][params_1min['Err']!=0]=0
        
    skyopc_qcd['qc_int']=np.ones(len(skyopc_qcd))
    skyopc_qcd['qc_int'][params_1min['Int']!=6]=0
        
    if save: 
        skyopc_qcd.to_csv(save+'SKYOPC_%s'%(str(start.date())))
        skyopc_params.to_csv(save+'SKYOPC_params_%s'%(str(start.date())))

    #except:
        #print('no data')
    try:
        return skyopc_qcd, skyopc_params
    except:
        print('Cant find any suitable data')
        return pd.DataFrame(columns=['ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8', 'ch9', 'ch10',
       'ch11', 'ch12', 'ch13', 'ch14', 'ch15', 'ch17', 'ch18', 'ch19', 'ch20',
       'ch21', 'ch22', 'ch23', 'ch24', 'ch25', 'ch26', 'ch27', 'ch28', 'ch29',
       'ch30', 'ch31', 'ch32', 'qc_pap', 'qc_prp', 'qc_err', 'qc_int']),pd.DataFrame(columns=['quad_Tmp', 'Err', 'pAp', 'pRp', 'Int'])



    
def extract_opc(opc_n,start,stop,dpath,save=False):
    """
    Extracts alphasense-N3 opc data from raw output. 
    QC's for bad data. 
    resamples to 1 minutely median
    Saves as .csv if requested
    
    Parameters:
        ocp_n: name, 'MSF' or 'TAWO'
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
    
    """
    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*%s*OPC*'%opc_n)
    file_dates=[]
    for f in all_files: 
        try:
            file_dates.append((dt.datetime.strptime(f[-12:-4], '%Y%m%d')).date())
        except:
            file_dates.append((dt.datetime.strptime(f[-14:-4], '%Y-%m-%d')).date())

    file_dates = np.asarray(file_dates)
           
    idxs = np.where(np.logical_and(file_dates>=start.date(), file_dates<=stop.date()))[0]
    dfs = [all_files[i] for i in idxs]
    opc = pd.DataFrame()
    
    if len(dfs)!=0: 
        # Extract the data
        for f in dfs: 
            # Ignore file if it's empty
            if os.path.getsize(f)==0:
                print('Error with: '+f+' this file is empty.\n')
                continue 
            opc = opc.append(pd.read_csv(f, skiprows=4,sep=',',error_bad_lines=False))  
            opc['Dates'] = pd.to_datetime(opc['time'],format='%Y-%m-%d %H:%M:%S')
            opc = opc.sort_values('Dates')
            opc = opc.set_index(opc['Dates'])
            opc.index = pd.DatetimeIndex(opc.index)
            #opc = opc[~opc.index.duplicated()]
            del opc['time'], opc['Dates']

        # Convert flow rate from L/min to cm3/s
        # 1 L/min = 16.66667 cm3/s
        opc.FlowRate = opc.FlowRate/100 * 16.66667

        opc_counts = opc[opc.columns[0:24]]
        opc_counts = opc_counts.apply(pd.to_numeric, errors='coerce')
        opc_params = opc[opc.columns[24:]]
    
        # Convert counts/interval to total counts/s
        opc.period = opc.period/100.0 # period in s
        opc_counts = opc_counts.divide(opc.period, axis=0)
        # Convert total counts/second to counts/cm3
        opc_counts = opc_counts.divide(opc.FlowRate, axis=0)
    
        # QC for OPC params. 
        opc_counts[opc_params['period']>500] = np.nan
        opc_counts[opc_params['period']<400] = np.nan

        # Resample to minutly average
        new_index = pd.date_range(opc_counts.index[0].round('min'),opc_counts.index[-1].round('min') , freq='min')      
        opc_1min = opc_counts.resample('1min').mean()
        opc_1min = opc_1min.reindex(new_index)
    
        # QC for winds and flight days
        opc_qc = qc_aerosol(opc_1min)
    
        # Save if neccessary
        if save: 
            opc_qc.to_csv(save+'%s_OPC_%s'%(opc_n,str(start.date())))
            opc_params.to_csv(save+'%s_OPC_params_%s'%(opc_n,str(start.date())))
            
        return
    else:
        print('No data for %s'%str(start))
        return



def get_clasp(d_loc,d1,d2,claspn,calfile,save=False):
    """
    Extracts clasp data from raw output. 
    QC's for bad data. 
    QC-ing for CLASP parameters is currently not implemented. 
    resamples to 1 minutely median
    Saves as .csv if requested
    
    Parameters:
        d1:     Start datetime for processing
        d2:     Stop datetime for processing
        d_loc:  Raw data filepath
        claspn: CLASP name, i.e. 'CLASP_F'
        save:   Output directory (optional)
    
    """
    # Function to convery interger to binary.
    get_bin = lambda x, n: format(x, 'b').zfill(n)
    os.chdir(d_loc)                  # Change directory to where the data is
    all_files = glob.glob('*%s*'%claspn)
    file_dates = np.asarray([(dt.datetime.strptime(f[-14:-4], '%Y-%m-%d')).date() for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=d1.date(), file_dates<=d2.date()))[0]
    dfs = [all_files[i] for i in idxs]

    # Loop through, extract and sort data into the dataframes initialised above
    #for i in range(0,np.shape(data_block)[0]):
    for f in dfs:
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue 
        
        CLASP = pd.read_csv(f,header=None,delim_whitespace=True, parse_dates={'Dates' : [0,1,2,3,4,5]},error_bad_lines=False)
        if np.shape(CLASP)[1]!=20:
            print('Issues with CLASP file %s'%f)
            continue

        CLASP['Dates'] = pd.to_datetime(CLASP['Dates'],format='%Y %m %d %H %M %S')
        CLASP = CLASP.sort_values('Dates')
        CLASP = CLASP.set_index(CLASP['Dates'])
        CLASP.index = pd.DatetimeIndex(CLASP.index)
        CLASP = CLASP[~CLASP.index.duplicated()]
        
        del CLASP['Dates']
        CLASP.columns = ['status', 'value','overflow','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16']
    
        # Extract, and convert staus addresses
        #statusbyte=float(split[6])
        binary=[]
        statusaddr=[]
        heater=[]
        for i in range(0,len(CLASP)):
            try:
                bi = get_bin(int(CLASP['status'].iloc[i]),8)

                # Check overflow flags and correct histogram
                # Overflow is for channels 1 to 8 only
                obin=get_bin(int(CLASP['overflow'].iloc[i]),8)
                for n in range(0,8):
                    if int(obin[n])>0:
                        CLASP[str(n)].iloc[i] = CLASP[str(n)].iloc[i] + int(obin[n]) * 256
                
                statusaddr.append(int(bi[4:8],2))
                heater.append(int(bi[2]))
                
            except:
                CLASP.iloc[i]=np.nan
                heater.append(np.nan)
                binary.append(np.nan)
                statusaddr.append(np.nan)
                    
        statusaddr = np.asarray(statusaddr)
        CLASP['heater']=heater
    
        # Arrange parameters into a neat dataframe  
        param_df = pd.DataFrame()
        param_len = int(len(statusaddr)/10)
        parameter = CLASP['value']
        param_df['rejects'] = parameter[np.where(statusaddr==0)[0][0:param_len]].resample('1min').mean()
        param_df['threshold'] = parameter[np.where(statusaddr==1)[0][0:param_len]].resample('1min').mean()
        param_df['ThisFlow']= parameter[np.where(statusaddr==2)[0][0:param_len]].resample('1min').mean()
        param_df['FlowPWM'] = parameter[np.where(statusaddr==3)[0][0:param_len]].resample('1min').mean()
        param_df['PumpCurrent'] = parameter[np.where(statusaddr==4)[0][0:param_len]].resample('1min').mean()
        param_df['SensorT'] = parameter[np.where(statusaddr==5)[0][0:param_len]].resample('1min').mean()
        param_df['HousingT'] = parameter[np.where(statusaddr==6)[0][0:param_len]].resample('1min').mean()
        param_df['PumpT'] = parameter[np.where(statusaddr==7)[0][0:param_len]].resample('1min').mean()
        param_df['SupplyV'] = parameter[np.where(statusaddr==8)[0][0:param_len]].resample('1min').mean()
        param_df['LaserR']  = parameter[np.where(statusaddr==9)[0][0:param_len]].resample('1min').mean()
        param_df['heater'] = CLASP['heater'].resample('1min').max()

        del CLASP['status'],CLASP['value'],CLASP['overflow'],CLASP['heater']
   
        # Resample to minutly average
        new_index = pd.date_range(CLASP.index[0].round('min'),CLASP.index[-1].round('min'), freq='min')      
        CLASP_1min = CLASP.resample('1min').median()
        CLASP_1min = CLASP_1min.reindex(new_index)

        # Apply flow corrections and quality flags, 
        # convert raw counts to concentrations in particles per ml.
        # Get calibration data
        cal_dict=io.loadmat(calfile,matlab_compatible=True)
        TSIflow = cal_dict['calibr'][0][0][8][0]         # array of calibration flow rates from TSI flow meter
        realflow = cal_dict['calibr'][0][0][9][0]        # array of measured A2D flow rates matching TSflow

        # Do flow correction and convert to concentations
        # TSI flow is from the TSI flowmeter, realflow is the flow the CLASP records internally
        P = np.polyfit(realflow,TSIflow,2) # These are from the flow calibration - fit a polynomial
        flow = np.polyval(P,param_df['ThisFlow']) # flow in L/min
        param_df['Sample volume (ml/s)'] = ((flow/60)*1000)/1

        # Now to plot concentrations in counts/cm3, just need to divide the counts/s by the sample volume
        CLASP_1min = CLASP_1min.apply(pd.to_numeric, errors='coerce')
        clasp_conc = CLASP_1min.divide(param_df['Sample volume (ml/s)'],axis=0)
   
        # QC for params
        # Currently not implemented. 
        
        # QC for winds/ flight days
        clasp_qcd = qc_aerosol(clasp_conc)
    
        # Save if neccessary
        if save: 
            clasp_qcd.to_csv(save+'%s_%s'%(claspn,str(d1.date())))
            param_df.to_csv(save+'%s_params_%s'%(claspn,str(d1.date())))
        
    
    return

def extract_pops(start,stop,dpath,save=False):
    """
    Extracts POPS data from Housekeeping csv files. 
    Resamples to 1 minutely mean
    Saves as .csv if requested
    
    Parameters:
        start:     Start datetime for processing
        stop:     Stop datetime for processing
        dpath:  Filepath for housekeeping files
        save:   Output directory (optional)
    
    """
    os.chdir(dpath)
    all_files = glob.glob('HK_*.csv')
    all_files.sort()
    file_dates = np.asarray([(dt.datetime.strptime(f[3:11], '%Y%m%d')).date() for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start.date(), file_dates<=stop.date()))[0]
    fils = [all_files[i] for i in idxs]
    
    for fil in fils:
        print(str(fil))
        if fil == fils[0]:
            df = pd.read_csv(fil)
            df['time'] = pd.to_datetime(df['DateTime'],unit='s')
            df.set_index('time',inplace=True)
            #df.index = df['time'].dt.round('1s')
            df = df[~df.index.duplicated()]
            df.sort_index(inplace=True)
    
            # resample to 1-minutely means 
            new_index = pd.date_range(df.index[0].round('1min'),df.index[-1].round('1min'), freq='min')      
            df_1min = df.resample('1min').mean()
            df_1min = df_1min.reindex(new_index,fill_value=np.NaN)
        else:
            df_temp = pd.read_csv(fil)
            df_temp['time'] = pd.to_datetime(df_temp['DateTime'],unit='s')
            df_temp.set_index('time',inplace=True)
            #df.index = df['time'].dt.round('1s')
            df_temp = df_temp[~df_temp.index.duplicated()]
            df_temp.sort_index(inplace=True)
    
            # resample to 1-minutely means 
            new_index = pd.date_range(df_temp.index[0].round('1min'),df_temp.index[-1].round('1min'), freq='min')      
            df_1min_temp = df_temp.resample('1min').mean()
            df_1min_temp = df_1min_temp.reindex(new_index,fill_value=np.NaN)   
            
            df_1min = df_1min.append(df_1min_temp)
    
    
    # PartCt = particle counds / s, Note this is the same as 'HistSum'
    # PartCon = particle concentration cm-3
    # POPS_Flow = current sample flow rate (cm3 /s)
    # LD_Mon = Laser diode output power mointor
    # BatV = Battery DV power voltage
    # Bins: b0 .... nbins
    
    df_1min  = df_1min[~df_1min.index.duplicated()]
    df_1min.sort_index(inplace=True)

    nbins = list(set(df_1min['nbins'].dropna()))
    if len(nbins)!=1: 
        print('Caution, bin length changed.')
    
    data_keys = ['b%s'%int(s) for s in np.arange(0,nbins[0])]
    data = df_1min[data_keys]
    total_conc = df_1min['PartCon']
  
    # QC for winds/ flight days
    df_1min = qc_aerosol(df_1min)
    data = qc_aerosol(data)

    # Save if neccessary
    if save: 
        df_1min.to_csv(save+'POPS_%s'%(str(d1.date())))

    return df_1min, data, total_conc



def get_dist(df,nbins,bounds,sum_hist=False):
    """
    Calculateds dNdlogd from aerosol size distribution
    
    Parameters:
        df:     Size distribution
        nbins:  total number of bins
        bounds: Bin boundaries (list of length nbins+1) 
        
    Returns: 
        mid_points: Bin mid-points
        dNdlogd
    
    """
    if len(bounds)!=nbins+1:
        print('Error bounds')
        return

    mid_points = [(bounds[i+1]+bounds[i])/2 for i in range(0,nbins)]
    logd = np.log10(bounds)   
    dlogd = [logd[i+1]-logd[i] for i in range(0,len(mid_points))]
    # Sum columns
    if sum_hist:
        df = df.sum(axis=0)
    
    dNdlogd = df/dlogd
    return mid_points,dNdlogd


def plot_dist(dists,labels,xlims):
    """
    Plots aerosol size distributions
    
    Parameters:
        dists:  list of size distributions to plot (dNdlogd)
        labels: Names for legend
        xlims:  x-axis (particle diameter) limits 
        
    Returns: 
        fig: Size distribution plot
    
    """
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    ax.grid(True)
    for i in range(0,len(dists)):
        ax.loglog(dists[i][0],dists[i][1],label=labels[i])

    ax.set_xlim(xlims[0],xlims[1])
    ax.set_xlabel('Diameter (d) / $\mu$m')
    ax.set_ylabel('dN/dlogd (cc$^{-3}$)')
    ax.legend(loc='best',fontsize=10)
    fig.tight_layout()
    
    return fig


def get_cpc(start,stop,d_loc):
    """
    Retrieves CPC data from .csv output of extract_cpc_data
    
    Parameters:
        start:      Start datetime for processing
        stop:       Stop datetime for processing
        d_loc:      data filepath
        
    Returns:
        CPC dataframe
    
    """
    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    CPC_out = pd.DataFrame(columns=['c/cm3', 'QC'])
    for date in f_date_list:
        f = d_loc + r'CPC_%s'%(str(date.date()))
        
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0])
        except:
            print('No data for %s'%str(date.date()))
            continue

        CPC_out = CPC_out.append(data,sort=True)    
    
    # Get rid of any duplicates
    CPC_out = CPC_out[~CPC_out.index.duplicated()]
    
    # Fill any missing minutes with nans
    new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq='min')
    CPC_out = CPC_out.reindex(new_index)
    
    # Crop to datetime
    CPC_out=CPC_out[start:stop]
    
    return CPC_out


def get_skyopc(start,stop,d_loc):
    """
    Retrieves skyopc data from .csv output of extract_skyopc_data
    
    Parameters:
        start:      Start datetime for processing
        stop:       Stop datetime for processing
        d_loc:      data filepath
        
    Returns:
        skyopc dataframe
    
    """
    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    # changed column number from 20 to 30 after inlet change 2022,8,21
    SKY_out = pd.DataFrame(columns=list(np.arange(0,31,1))+['qc_pap','qc_prp','qc_err','qc_int'])
    for date in f_date_list:
        f = d_loc + r'SKYOPC_%s'%(str(date.date()))
        
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0],skiprows=1,names=list(np.arange(0,31,1))+['qc_pap','qc_prp','qc_err','qc_int'])
        except:
            print('No data for %s'%str(date.date()))
            continue
        
        SKY_out = SKY_out.append(data,sort=True)    
    
    # Get rid of any duplicates
    SKY_out = SKY_out[~SKY_out.index.duplicated()]
    
    # Fill any missing minutes with nans
    new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq='min')
    SKY_out = SKY_out.reindex(new_index)
    
    # Crop to datetime
    SKY_out=SKY_out[start:stop]
    
    return SKY_out

def get_opc(start,stop,d_loc):
    """
    Retrieves alphasense opc data from .csv output of extract_opc_data
    
    Parameters:
        start:      Start datetime for processing
        stop:       Stop datetime for processing
        d_loc:      data filepath
        
    Returns:
        opc dataframe
    
    """
    f_date_list = pd.date_range(start.date(),stop.date(),freq='1D')
    out = pd.DataFrame(columns=list(np.arange(0,24,1))+['QC'])
    for date in f_date_list:
        f = d_loc + r'MSF_OPC_%s'%(str(date.date()))
        
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0],skiprows=1,names=list(np.arange(0,24,1))+['QC'])
        except:
            try: 
                f = d_loc + r'TAWO_OPC_%s'%(str(date.date()))
                data = pd.read_csv(f,parse_dates=[0],index_col=[0],skiprows=1,names=list(np.arange(0,24,1))+['QC'])
            except:            
                print('No data for %s'%str(date.date()))
                continue
        
        out = out.append(data,sort=True)    
    
    # Get rid of any duplicates
    out = out[~out.index.duplicated()]
    
    # Fill any missing minutes with nans
    new_index = pd.date_range(start,stop-pd.Timedelta(minutes=1), freq='min')
    out = out.reindex(new_index)
    
    # Crop to datetime
    out=out[start:stop]
    
    return out

def extract_biral(start,stop,dpath,save=False):
    """
    Extracts biral present weather sensor data from raw output. 
    QC's for bad data. 
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
    
    """
    
    headers = ['time','adr','int','avg_vis','xx_a','pw','xx_b','inst_vis','self-test']
    custom_date_parse = lambda x: dt.datetime.strptime(x,'%Y %m %d %H %M %S.%f SWS100')
    weather_codes={10:'No weather',20:'Haze',30:'Fog',40:'Indeterminate precip',50:'Drizzle',60:'Rain',70:'Snow'}
    
    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*.biral')
    
    # Get start and stop filenames
    start_f = int((start-dt.timedelta(days=1)).strftime("%y%m%d"))
    stop_f = int(stop.strftime("%y%m%d"))
    
    # Extract daterange
    file_dates = np.asarray([int(f[0:6]) for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start_f, file_dates<=stop_f))[0]
    dfs = [all_files[i] for i in idxs]
    dfs.sort()
    biral_df = pd.DataFrame(columns=headers)
    
    for f in dfs: 
        # Ignore file if it's empty
        if os.path.getsize(f)==0:
            print('Error with: '+f+' this file is empty.\n')
            continue 
        try:
            temp = pd.read_csv(f,delimiter=',',parse_dates=[0],date_parser=custom_date_parse,header=None,index_col=[0],names=headers,error_bad_lines=False,warn_bad_lines=True)
            # convert strings to float
            avg_vis = [np.float(np.asarray(temp['avg_vis'])[ind][0:5]) for ind in range(0,len(temp['avg_vis']))]
            inst_vis = [np.float(np.asarray(temp['inst_vis'])[ind][0:5]) for ind in range(0,len(temp['inst_vis']))]
            temp['avg_vis']=avg_vis
            temp['inst_vis']=inst_vis
            biral_df = pd.concat([biral_df,temp],sort=True)
        except:
            print('Skipping: %s'%files[i])
    
    biral_df.sort_index(inplace=True)
    
    # sort dates
    if len(biral_df)>0:
        resample_dates = pd.date_range(biral_df.index[0].replace(second=0,microsecond=0),biral_df.index[-1].replace(second=0,microsecond=0),freq='1min')
        biral_df = biral_df[~biral_df.index.duplicated()]
        biral_df = biral_df.reindex(resample_dates,method='nearest',tolerance='30s')
    else:
        print('No data')
        return biral_df
    
    # Get bad data times
    fail_self_check_indices = biral_df[((biral_df['self-test']!='OOO') & (biral_df['self-test']!='XOO'))].index
    # round to nearest minute
    fail_self_check = [fsc.round('min') for fsc in fail_self_check_indices]
    
    # Get rid of 'xx' values from present weather:
    biral_df['pw'] = pd.to_numeric(biral_df['pw'],errors='coerce')
    
    # change pw=4 (haze) to 20 for the sake of the plot and pw=0 (no weather) to 10
    biral_df['pw_adjusted'] = biral_df['pw']
    biral_df.loc[biral_df['pw']==4,'pw_adjusted'] = 20
    biral_df.loc[biral_df['pw']==0,'pw_adjusted'] = 10
    
    biral_df['qc']=1
    biral_df.loc[fail_self_check,'qc']=0
    
    # reindex to daterange
    biral_df = biral_df.reindex(pd.date_range(start,stop,freq='1min'))
    
    
    # Save if neccesary
    if save:
        biral_df.to_csv(save+'biral_%s'%(str(start.date())))
    
    return biral_df

