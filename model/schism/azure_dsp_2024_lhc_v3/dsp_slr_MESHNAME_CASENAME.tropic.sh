#!/bin/sh
set -e

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)

# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
gen_elev2d --outfile {meshname}.{cname}.elev2D.th.nc --hgrid=hgrid.gr3 --stime={year_start}-{month_start}-{day_start} --etime={year_end}-{month_end}-{day_end} --slr 1.07 pt_reyes{tidal_pert}.csv monterey{tidal_pert}.csv
ln -sf {meshname}.{cname}.elev2D.th.nc elev2D.th.nc

# CREATE OTHER SYMBOLIC LINKS ----------------------------
# shared spatial inputs
ln -sf bctides.in.2d bctides.in
ln -sf vgrid.in.2d vgrid.in

# modified TH inputs
{linked_th_file_strings}

# inputs specific to this setup
ln -sf SAL_nu_roms.nc SAL_nu.nc
ln -sf salinity_nudge_roms.gr3 SAL_nudge.gr3
ln -sf TEM_nu_roms.nc TEM_nu.nc
ln -sf temperature_nudge_roms.gr3 TEM_nudge.gr3
ln -sf param.nml.tropic param.nml
# ln -sf ../launch.tropic.pbs launch.pbs # only necessary with HPC4

# CHECK OUTPUT DIRECTORY IS PRESENT ----------------------

mkdir -p outputs # only creates outputs if it's not preset

# model runs with yml file for Azure setup