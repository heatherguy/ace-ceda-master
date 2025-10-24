#!/bin/bash
#SBATCH --job-name="ICECAPS-extract_main"
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=icecaps
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH --time 24:00:00



# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Extract data files from .raw
ext_start='230701'
ext_stop='250101'

raw_dir_flux='/gws/ssde/j25b/icecaps/ICECAPSarchive/fluxtower/raw/'
#raw_dir_flux='/gws/nopw/j04/ncas_radar_vol1/heather/fluxtower_temp/'

raw_dir_ace='/gws/ssde/j25b/icecaps/ICECAPSarchive/ace/raw/'
extract_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/'


#echo 'Extracting SnD'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'SnD'

#echo 'Extracting biral'
#python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'biral'

#echo 'Extracting HMP'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'HMP'

#echo 'Extracting KT'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'KT'

#echo 'Extracting CPC'
#python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'CPC'

#echo 'Extracting SKYOPC'
#python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'skyopc'

echo 'Extracting POPS'
python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'F2024'

#echo 'Extracting msf-opc'
#python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'MSF'

#echo 'Extracting tawo-opc'
#python tar_extract.py $ext_start $ext_stop $raw_dir_ace $extract_dir 'TAWO'

#echo 'Extracting licor'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'licor'

#echo 'Extracting metek'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'metek'

#echo 'Extracting ventus'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'ventus'

#echo 'Extracting simba'
#python tar_extract.py $ext_start $ext_stop $raw_dir_flux $extract_dir 'simba'
