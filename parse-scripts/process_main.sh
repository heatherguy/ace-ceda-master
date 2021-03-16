#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 05:00:00

# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Pre-process raw data

snd_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/SnD/'
kt_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/KT15/'
v_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/ventus/'
m_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/metek/'
proc_dir='/gws/nopw/j04/ncas_radar_vol1/heather/processed/'
li_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/licor/'
cpc_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/CPC/'
sky_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/SKYOPC/'
hmp_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/HMP/' 
opc_dir='/gws/nopw/j04/ncas_radar_vol1/heather/extracted/msf-opc/'

all_start='202010150000'
all_stop='202101010000'
hmp_path='/gws/nopw/j04/ncas_radar_vol1/heather/processed/HMP/'
kt_qcf='/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/qc-files/KT_bad_dates'
skyopc_qc='/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/qc-files/skyopc_percent_loss_data.csv'

python process_raw.py $hmp_dir $proc_dir 'HMP' $all_start $all_stop
python process_raw.py $snd_dir $proc_dir 'SnD' $all_start $all_stop $hmp_path
#python process_raw.py $kt_dir $proc_dir 'KT' $all_start $all_stop $kt_qcf
python process_raw.py $v_dir $proc_dir 'ventus' $all_start $all_stop
python process_raw.py $m_dir $proc_dir 'metek' $all_start $all_stop
python process_raw.py $li_dir $proc_dir 'licor' $all_start $all_stop
#python process_raw.py $cpc_dir $proc_dir 'CPC' $all_start $all_stop
python process_raw.py $sky_dir $proc_dir 'SKYOPC' $all_start $all_stop $skyopc_qc
python process_raw.py $opc_dir $proc_dir 'MSF_OPC' $all_start $all_stop





