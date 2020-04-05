#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 22:32:23 2019

@author: heather
"""

import numpy as np
import datetime as dt
import pandas as pd
import os
import glob
from scipy import io
from utils import * 

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
    flight_dates_f = '/Volumes/Data/ICECAPSarchive/qc_files/flight_days.csv'
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
           
    # Get Met data     
    w_dloc = '/Volumes/Data/ICECAPSarchive/Summit_Met/met_sum_insitu_1_obop_minute_%s_%s.txt'%(sdate.year,str(sdate.month).zfill(2))
    try:
        met = get_NOAA_met(w_dloc)
        qc_in['QC'][met['ws']<1]=0
        qc_in['QC'][met['wd']>270]=0
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

        # Resample to minutly average
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
    resamples to 1 minutely median
    Saves as .csv if requested
    
    Parameters:
        start: Start datetime for processing
        stop:  Stop datetime for processing
        dpath: Raw data filepath
        save:  Output directory (optional)
    
    """
    os.chdir(dpath)                  # Change directory to where the data is
    all_files = glob.glob('*SKYOPC*')
    file_dates = np.asarray([(dt.datetime.strptime(f[-14:-4], '%Y-%m-%d')).date() for f in all_files])
    idxs = np.where(np.logical_and(file_dates>=start.date(), file_dates<=stop.date()))[0]
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
            if line[0] =='P':
                if len(line)!=17:
                    c=0
                    datetime=np.nan
                    continue
                #Year Mon Day Hr Min Loc 4Tmp Err pA/p pR/p UeL Ue4 Ue3 Ue2 Ue1 Iv 
                datetime = dt.datetime(int(line[1])+2000,int(line[2]),int(line[3]),int(line[4]),int(line[5]))
                #datetime = dt.datetime.strptime('20'+line[1]+line[2]+line[3]+line[4]+line[5],'%Y%m%d%H%M')
                quad_Tmp = int(line[7])
                Err = int(line[8])
                pAp = int(line[9])
                pRp = int(line[10])
                Int = int(line[16])
                c=0
            
            elif len(line)!=9:
                continue

            elif c==0: 
                ch1=int(line[1])
                ch2=int(line[2])
                ch3=int(line[3])
                ch4=int(line[4])
                ch5=int(line[5])
                ch6=int(line[6])
                ch7=int(line[7])
                ch8=int(line[8])
                c = c+1    
            elif c ==1:
                ch9=int(line[1])
                ch10=int(line[2])
                ch11=int(line[3])
                ch12=int(line[4])
                ch13=int(line[5])
                ch14=int(line[6])
                ch15=int(line[7])
                ch16=int(line[8])
                c = c+1
            elif c == 2:
                ch17=int(line[1])
                ch18=int(line[2])
                ch19=int(line[3])
                ch20=int(line[4])
                ch21=int(line[5])
                ch22=int(line[6])
                ch23=int(line[7])
                ch24 =int(line[8])
                c= c+1
            elif c==3:
                ch25=int(line[1])
                ch26=int(line[2])
                ch27=int(line[3])
                ch28=int(line[4])
                ch29=int(line[5])
                ch30=int(line[6])
                ch31=int(line[7])
                ch32=int(line[8])
                c = 0
                n = int(line[0][-2])
                if isinstance(datetime,dt.datetime):
                    skyopc = skyopc.append(pd.Series([datetime+dt.timedelta(seconds=n*6), ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8, ch9, ch10, ch11, ch12, ch13, ch14, ch15, ch16, ch17, ch18, ch19, ch20, ch21, ch22, ch23, ch24, ch25, ch26, ch27, ch28, ch29, ch30, ch31, ch32, quad_Tmp,Err,pAp,pRp,Int]),ignore_index=True)

    try: 
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
    
        # QC for params
        skyopc[skyopc['pAp']>120]=np.nan
        skyopc[skyopc['pRp']>120]=np.nan  
        skyopc[skyopc['Err']!=0]=np.nan
        skyopc[skyopc['Int']!=6]=np.nan
    
        skyopc_counts = skyopc[skyopc.columns[0:31]]
        skyopc_counts = skyopc_counts.apply(pd.to_numeric, errors='coerce') # Counts in counts/ 6 seconds
        skyopc_counts =skyopc_counts / 100.0 # convert from counts/100ml to counts/cm3    
        skyopc_params = skyopc[skyopc.columns[31:]]
    
        # Resample to minutly average
        new_index = pd.date_range(skyopc.index[0].round('min'),skyopc.index[-1].round('min'), freq='min')      
        skyopc_1min = skyopc_counts.resample('1min').median()
        skyopc_1min = skyopc_1min.reindex(new_index)
    
        # QC
        skyopc_qcd = qc_aerosol(skyopc_1min)
        
        if save: 
            skyopc_qcd.to_csv(save+'SKYOPC_%s'%(str(start.date())))
            skyopc_params.to_csv(save+'SKYOPC_params_%s'%(str(start.date())))

    except:
        print('no data')
    
    return 



    
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
    if opc_n=='TAWO':
        file_dates = np.asarray([(dt.datetime.strptime(f[-12:-4], '%Y%m%d')).date() for f in all_files])
    else:
        file_dates = np.asarray([(dt.datetime.strptime(f[-14:-4], '%Y-%m-%d')).date() for f in all_files])
           
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
        opc_1min = opc_counts.resample('1min').median()
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




def get_dist(df,nbins,bounds):
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
    logd = np.log(bounds)   
    dlogd = [logd[i+1]-logd[i] for i in range(0,len(mid_points))]
    # Sum columns
    hist = df.sum(axis=0)
    dNdlogd = hist/dlogd
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
    SKY_out = pd.DataFrame(columns=list(np.arange(0,31,1))+['QC'])
    for date in f_date_list:
        f = d_loc + r'SKYOPC_%s'%(str(date.date()))
        
        try:
            data = pd.read_csv(f,parse_dates=[0],index_col=[0],skiprows=1,names=list(np.arange(0,31,1))+['QC'])
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

