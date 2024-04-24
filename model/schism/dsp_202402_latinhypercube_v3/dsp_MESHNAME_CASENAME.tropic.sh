#!/bin/sh
set -e

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)

# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
gen_elev2d --outfile {meshname}.{cname}.elev2D.th.nc --hgrid=hgrid.gr3 --stime={year_start}-{month_start}-{day_start} --etime={year_end}-{month_end}-{day_end} --slr 0.0 ../pt_reyes{tidal_pert}.csv ../monterey{tidal_pert}.csv
ln -sf {meshname}.{cname}.elev2D.th.nc elev2D.th.nc

# CREATE OTHER SYMBOLIC LINKS ----------------------------
# shared inputs
ln -sf station.in station.in

# shared spatial inputs
ln -sf bctides.in.2d bctides.in
ln -sf vgrid.in.2d vgrid.in

# modified TH inputs
{linked_th_file_strings}

# inputs specific to this setup
ln -sf param.nml.tropic param.nml
# ln -sf ../launch.tropic.pbs launch.pbs # only necessary with HPC4

# CHECK OUTPUT DIRECTORY IS PRESENT ----------------------

mkdir -p outputs # only creates outputs if it's not preset

# RUN MODEL ----------------------------------------------

# qsub launch.pbs # HPC4
sbatch slurm_{meshname}_{cname}.tropic.sh # HPC5