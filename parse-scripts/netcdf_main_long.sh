#!/bin/bash
#SBATCH --partition=long-serial
#SBATCH -o /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.out 
#SBATCH -e /gws/nopw/j04/ncas_radar_vol1/heather/logs/%j.err
#SBATCH -t 72:00:00



# Check we're in the right directory
cd /gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master/parse-scripts

# activate python environment
source /home/users/guyh/miniconda3/etc/profile.d/conda.sh
conda activate guyh

# Generate netcdf file
meta_dir='/gws/nopw/j04/ncas_radar_vol1/heather/ace-ceda-master'
netcdf_out='/gws/nopw/j04/ncas_radar_vol1/heather/final_nc/'
months=[1,2,3,4,5]
years=[2021]
avp=1


in_loc='/gws/nopw/j04/ncas_radar_vol1/heather/'
in_loc_proc='/gws/nopw/j04/ncas_radar_vol1/heather/processed/'

#python parse_snow-height.py $in_loc $netcdf_out $months $years $avp
#python parse_surface-temperature-profile.py $in_loc $netcdf_out $months $years $avp
#python parse_skin-temperature.py $in_loc_proc $netcdf_out $months $years $avp
#python parse_surface-moisture-profile.py $in_loc $netcdf_out $months $years $avp
#python parse_surface-winds-profile.py $in_loc $netcdf_out $months $years $avp
#python parse_aerosol-concentration.py $in_loc_proc $netcdf_out $months $years $avp
#python parse_aerosol-size-distribution.py $in_loc_proc $netcdf_out $months $years $avp
#python parse_aerosol-opc.py $in_loc_proc $netcdf_out $months $years $avp

start_dat='202107200000'
stop_dat='202201010000'

#python parse_flux-estimates.py $in_loc_proc $netcdf_out $start_dat $stop_dat 30 1
#python parse_flux-estimates.py $in_loc_proc $netcdf_out $start_dat $stop_dat 30 2
python parse_flux-estimates.py $in_loc_proc $netcdf_out $start_dat $stop_dat 15 1
python parse_flux-estimates.py $in_loc_proc $netcdf_out $start_dat $stop_dat 15 2

