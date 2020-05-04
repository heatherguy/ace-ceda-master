#!/bin/bash

# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Extract data files from .raw
ext_start='190201'
ext_stop='190601'
raw_dir_flux='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/fluxtower/raw/'
raw_dir_ace='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/ace/raw/'
extract_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/'


#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'SnD'
##python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'HMP'

#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'KT'

#echo 'trying cpc'
python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'CPC'

#echo 'trying skyopc'
#python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'SKYOPC'

#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'licor'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'metek'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'ventus'


