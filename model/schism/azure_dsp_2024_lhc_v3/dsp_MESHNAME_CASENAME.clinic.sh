#!/bin/sh
set -e

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)
# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
gen_elev2d --outfile {meshname}.{cname}.elev2D.th.nc --hgrid=hgrid.gr3 --stime={year_start}-{month_start}-{day_start} --etime={year_end}-{month_end}-{day_end} --slr 0.0 pt_reyes{tidal_pert}.csv monterey{tidal_pert}.csv
# link .th.nc file
ln -sf {meshname}.{cname}.elev2D.th.nc elev2D.th.nc

# run script to create uv3d.th.nc
ulimit -s unlimited

# START !no-interp commands
[ "$1" != "no-interp" ] && 
# CHANGE OUTPUTS DIR -------------------------------------
mkdir -p outputs_tropic && 
mv outputs/* outputs_tropic && 
mkdir -p outputs &&
	

# CREATE OCEAN BOUNDARY ----------------------------------

cd outputs_tropic &&
rsync -avz ../interpolate_variables.in . &&

# link necessary files
ln -sf ../hgrid.gr3 bg.gr3 &&
ln -sf ../hgrid.gr3 fg.gr3 &&
ln -sf ../vgrid.in.2d vgrid.bg &&
ln -sf ../vgrid.in.3d vgrid.fg &&

interpolate_variables8 && # this takes quite a while
cp uv3D.th.nc ../uv3D.th.nc &&
cd ../
# END !no-interp commands

# Re link barotropic files -------------------------------

# modified TH inputs
{linked_th_file_strings}

# CREATE CLIINIC SYMBOLIC LINKS ----------------------------

# add new links
ln -sf SAL_nu_roms.nc SAL_nu.nc
ln -sf salinity_nudge_roms.gr3 SAL_nudge.gr3
ln -sf TEM_nu_roms.nc TEM_nu.nc
ln -sf temperature_nudge_roms.gr3 TEM_nudge.gr3
ln -sf hotstart_{case_year}.nc hotstart.nc

# change existing links
ln -sf bctides.in.3d bctides.in
ln -sf vgrid.in.3d vgrid.in
ln -sf param.nml.clinic param.nml

# model runs with Azure yml file