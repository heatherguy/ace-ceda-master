#!/bin/bash

# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Extract data files from .raw
ext_start='190601'
ext_stop='200401'
raw_dir='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/fluxtower/raw/'
#'/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/ace/raw/'
extract_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/SnD/'

#python tar_extract.py $ext_start $ext_stop $raw_dir $extract_dir 'SnD'



# Pre-process raw data

proc_dir='/gws/nopw/j04/ncas_radar_vol1/heather/processed/SnD'
all_start='201906010000'
all_stop='201912312359'
instrument='SnD'
hmp_path='/gws/nopw/j04/ncas_radar_vol1/heather/processed/HMP/'

#python process_raw.py $extract_dir $proc_dir $instrument $all_start $all_stop $hmp_path



# Generate netcdf file
meta_dir='/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master'
netcdf_out='/gws/nopw/j04/ncas_radar_vol1/heather/final_nc/'
months=[6,7,8,9,10,11,12]
years=[2019]
avp=1


in_loc='/gws/nopw/j04/ncas_radar_vol1/heather/'
#python parse_snow-height.py $in_loc $netcdf_out $months $years $avp
python parse_surface-temperature-profile.py $in_loc $netcdf_out $months $years $avp


