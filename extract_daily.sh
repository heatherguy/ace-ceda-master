#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 24:00:00
# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Extract data files from .raw
yesterday=`date -d "1 day ago" '+%y%m%d'`
today=`date '+%y%m%d'`

# try the last 10 days incase some got missed. 
# Clean up data over a week old
old_date=`date -d "20 days ago" '+%y%m%d'`
flux_extract_date=$old_date
ace_extract_date=$old_date

raw_dir_flux='/gws/nopw/j04/ncas_radar_vol1/heather/fluxtower_temp/'
raw_dir_ace='/gws/nopw/j04/ncas_radar_vol1/heather/fluxtower_temp/'
extract_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/'

while [ "$flux_extract_date" != "$today" ]; do
    echo "$flux_extract_date"
    # If date is not in log file, then run extract
	if grep -Fx "$flux_extract_date" /gws/nopw/j04/ncas_radar_vol1/heather/logs/fluxtower_extract_list.log
	then
    	echo "File already processed"
	else
    	echo "Trying to run fluxtower extract..."

    	if test -f "/gws/nopw/j04/ncas_radar_vol1/heather/fluxtower_temp/"$extract_date"_fluxtower.tar.gz"; then
  			echo 'Extracting SnD'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'SnD'

			echo 'Extracting HMP'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'HMP'

			echo 'Extracting KT'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'KT'

			echo 'Extracting licor'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'licor'

			echo 'Extracting metek'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'metek'

			echo 'Extracting ventus'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'ventus'

			echo 'Extracting simba'
			python tar_extract.py $extract_date $extract_date $raw_dir_flux $extract_dir 'simba'

			echo 'Moving data into appropriate directories...'
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.SnD /gws/nopw/j04/ncas_radar_vol1/heather/extracted/SnD/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.HMP* /gws/nopw/j04/ncas_radar_vol1/heather/extracted/HMP/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.KT15 /gws/nopw/j04/ncas_radar_vol1/heather/extracted/KT15/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.licor /gws/nopw/j04/ncas_radar_vol1/heather/extracted/licor/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.metek* /gws/nopw/j04/ncas_radar_vol1/heather/extracted/metek/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.ventus* /gws/nopw/j04/ncas_radar_vol1/heather/extracted/ventus/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.simba /gws/nopw/j04/ncas_radar_vol1/heather/extracted/simba/
			echo  "$flux_extract_date" >> /gws/nopw/j04/ncas_radar_vol1/heather/logs/fluxtower_extract_list.log
    		flux_extract_date=$(date -d "$flux_extract_date + 1 day" '+%y%m%d')
    	else
    		echo "can't find fluxtower file"
    	fi
    	flux_extract_date=$(date -d "$flux_extract_date + 1 day" '+%y%m%d')
    fi
done

while [ "$ace_extract_date" != "$today" ]; do
    echo "$ace_extract_date"
    # If date is not in log file, then run extract
	if grep -Fx "$ace_extract_date" /gws/nopw/j04/ncas_radar_vol1/heather/logs/ace_extract_list.log
	then
    	echo "File already processed"
	else
    	echo "Trying to run ace extract..."
    	if test -f "/gws/nopw/j04/ncas_radar_vol1/heather/fluxtower_temp/"$extract_date"_ACE.tgz"; then
			echo 'Extracting SKYOPC'
			python tar_extract.py $extract_date $extract_date $raw_dir_ace $extract_dir 'skyopc'

			echo 'Extracting tawo-opc'
			python tar_extract.py $extract_date $extract_date $raw_dir_ace $extract_dir 'TAWO_AQRPI8_OPCN3'

			echo 'Extracting biral'
			python tar_extract.py $extract_date $extract_date $raw_dir_ace $extract_dir 'biral'

			echo 'Extracting pops'
			python tar_extract.py $extract_date $extract_date $raw_dir_ace $extract_dir 'F20'

			echo 'Moving data into appropriate directories...'
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.skyopc /gws/nopw/j04/ncas_radar_vol1/heather/extracted/SKYOPC/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*TAWO*OPC*.csv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/tawo-opc/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*.biral /gws/nopw/j04/ncas_radar_vol1/heather/extracted/biral/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/*POPS* /gws/nopw/j04/ncas_radar_vol1/heather/extracted/pops/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/Peak* /gws/nopw/j04/ncas_radar_vol1/heather/extracted/pops/
			mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/HK* /gws/nopw/j04/ncas_radar_vol1/heather/extracted/pops/
            mv /gws/nopw/j04/ncas_radar_vol1/heather/extracted/Log_*.txt /gws/nopw/j04/ncas_radar_vol1/heather/extracted/pops/
			echo  "$ace_extract_date" >> /gws/nopw/j04/ncas_radar_vol1/heather/logs/ace_extract_list.log
    	else
    		echo "can't find ace file"
    	fi
    	ace_extract_date=$(date -d "$ace_extract_date + 1 day" '+%y%m%d')
	fi
done
