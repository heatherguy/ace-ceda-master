#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:08:16 2020

@author: heather

Function to pre-process raw data. 

Inputs: 
    in_loc: location of extracted raw data
    out_loc: location to save pre-processed data
    instrument: Name of instrument to process
    start: Start date to process
    stop: stop date to process
    option_qc: Optional path to qc file if neccessary 
            (required for CLASP, KT and SnD)
    
Example usage: 
   python $in_loc $out_loc $instrument $all_start $all_stop $option_qc
"""

import sys,os
sys.path.append(os.getcwd())
from fluxtower_parse import *
from ace_parse import *
import numpy as np      
import datetime as dt
import pandas as pd

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
    instrument=args_in[3]
    all_start = dt.datetime.strptime(args_in[4],'%Y%m%d%H%M')
    all_stop = dt.datetime.strptime(args_in[5],'%Y%m%d%H%M')
    if instrument=='CLASP_F':
        calfile = args_in[6]
        return in_loc,out_loc,instrument,all_start,all_stop, calfile
    elif instrument=='CLASP_G':
        calfile = args_in[6]
        return in_loc,out_loc,instrument,all_start,all_stop, calfile
    elif instrument=='KT':
        KT_qcf = args_in[6]
        return in_loc,out_loc,instrument,all_start,all_stop, KT_qcf
    elif instrument=='SnD':
        hmp_dpath= args_in[6]
        return in_loc,out_loc,instrument,all_start,all_stop, hmp_dpath
    elif instrument=='SKYOPC':
        inlet_cal = args_in[6]
        return in_loc,out_loc,instrument,all_start,all_stop, inlet_cal
    else:
        # return values:
        return in_loc,out_loc,instrument,all_start,all_stop

def main():

    # check / get args:
    try:
        in_loc,out_loc,instrument,all_start,all_stop, calfile = get_args(sys.argv)
    except:
        try:
            in_loc,out_loc,instrument,all_start,all_stop = get_args(sys.argv)
        except:
            print('Input error')
            sys.exit()
        
    # Extract and process all data, and save to 'processed'

    all_days = pd.date_range(all_start,all_stop,freq='1D')

    # Loop through and split into daily files
    for i in range(0,len(all_days)-1):
        start = all_days[i]
        stop = all_days[i+1] - pd.Timedelta(seconds=0.1)
        print(str(start) + ' to ' + str(stop))
        
        if instrument=='KT':
            extract_KT_data(start,stop,in_loc,calfile,save=out_loc)
            #except:
            #    print('Unable to parse')
            #    continue
            
        elif instrument=='ventus':
            try:
                v1,v2s = extract_ventus_data(start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse')
                continue
            
        elif instrument=='metek':
            #try:
            m1,m2 = extract_metek_data(start,stop,in_loc,save=out_loc)
            #except:
            #    print('Unable to parse')
            #    continue
            
        elif instrument=='licor':
            try:
                licor = extract_licor_data(start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse')
                continue            
            
        elif instrument=='SnD':
           # try:
            snd = extract_snd_data(start,stop,in_loc,calfile,save=out_loc)
           # except:
             #   print('Unable to parse')
             #   continue

        elif instrument=='CPC':
            try:
                extract_cpc(start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse')
                continue
            
        elif instrument=='SKYOPC':
            #try:
            extract_skyopc(start,stop,in_loc,calfile,save=out_loc)
            #except:
            #    print('Unable to parse')
            #    continue
            
        elif instrument=='TAWO_OPC':
            try:
                extract_opc('TAWO',start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse')
                continue
            
        elif instrument=='MSF_OPC':
            try:
                extract_opc('MSF',start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse')
                continue
            
        elif instrument=='HMP':
            try:
                HMP1 = extract_HMP_data('HMP1',start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse HMP1')
            try:
                HMP2 = extract_HMP_data('HMP2',start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse HMP2')            
            try:
                HMP3 = extract_HMP_data('HMP3',start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse HMP3')            
            try:
                HMP4 = extract_HMP_data('HMP4',start,stop,in_loc,save=out_loc)
            except:
                print('Unable to parse HMP4')            
            
        elif instrument=='CLASP_F':
            try:
                get_clasp(in_loc,start,stop,'CLASP_F',calfile,save=out_loc)
            except:
                print('Unable to parse')  
                
        elif instrument=='CLASP_G':
            try:
                get_clasp(in_loc,start,stop,'CLASP_G',calfile,save=out_loc) 
            except:
                print('Unable to parse')  
            
        else:
            print('Instrument identifier failed')
            
if __name__ == '__main__':
    main()              
