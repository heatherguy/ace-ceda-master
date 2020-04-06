#!/bin/bash

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Extract data files from .raw

raw_dir='/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/fluxtower/raw/'
#'/gws/nopw/j04/ncas_radar_vol2/data/ICECAPSarchive/ace/raw/'

extract_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/'

python tar_extract.py $raw_dir $extract_dir 'HMP'



# Pre-process raw data

proc_dir='/gws/nopw/j04/ncas_radar_vol1/heather/processed/'
all_start='201906010000'
all_stop='201912312359'
instrument='HMP'

python process_raw.py $extract_dir $proc_dir $instrument $all_start $all_stop



# Generate netcdf file

#netcdf_out='/gws/nopw/j04/ncas_radar_vol1/heather/processed/final_nc/surface-temperature-profile/'
#months=[6,7,8,9,10,11,12]
#years=[2019]
#avp=1

#python parse_surface-temperature-profile.py $proc_dir $netcdf_out $months $years $avp


