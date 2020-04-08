#!/bin/bash

# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Extract data files from .raw
ext_start='190121'
ext_stop='200401'
raw_dir='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/fluxtower/raw/'
#'/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/ace/raw/'
extract_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/'

#python tar_extract.py $ext_start $ext_stop $raw_dir $extract_dir 'HMP'



# Pre-process raw data

proc_dir='/gws/nopw/j04/ncas_radar_vol1/heather/processed/hmp/'
all_start='202001210000'
all_stop='202003012359'
instrument='HMP'

#python process_raw.py $extract_dir $proc_dir $instrument $all_start $all_stop



# Generate netcdf file
meta_dir='/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master'
netcdf_out='/gws/nopw/j04/ncas_radar_vol1/heather/final_nc/'
months=[6,7,8,9,10,11,12]
years=[2019]
avp=1

python parse_surface-temperature-profile.py $proc_dir $netcdf_out $months $years $avp


