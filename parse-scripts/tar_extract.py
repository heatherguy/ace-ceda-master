#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:06:02 2020

@author: heather

Function to extract files from tar.gz based on an identifier string

Inputs: 
    in_loc: location of tar.gaz files
    out_loc: Directory to save extracted files
    id_str: String identifier
    
Extracted files are saved in out_loc.

Example usage:
    
python tar_extract.py $in_loc $out_loc $id_str

"""

import sys,os
sys.path.append(os.getcwd())
from utils import *
import numpy as np       
import datetime as dt  


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
    
    start=args_in[1]
    stop=args_in[2]
    in_loc = args_in[3]
    out_loc = args_in[4]
    id_str = args_in[5]

    return start,stop,in_loc,out_loc,id_str

def main():

    # Some examples:
    # Data names: 
    # HMP1: 2m T/RH
    # HMP2: 4m T/RH
    # HMP3: 8m T/RH
    # HMP4: 15m T/RH
    # metek1: 2m 3D sonic
    # metek2: 15m 3D sonic
    # licor: H2O and CO2 analyzer
    # KT15: Snow surface temp
    # ventus1: 4m 2D sonic
    # ventus2: 8m 2D sonic
    # SnD: Snow depth sensor.

    # What data do you want to get?
    #dname = 'SKYOPC'
    #dname = 'Summit_MSF_ICECAPSACE_OPCN3'
    #dname = 'CPC'

    # check / get args:
    try:
        start,stop,in_loc,out_loc,id_str = get_args(sys.argv)
    except:
        print('Input error')
        sys.exit()

    # Extract all tar files. 
    print('running extract tar')
    extract_tar(start,stop,in_loc,out_loc,id_str)
    
    return

if __name__ == '__main__':
    main()  
