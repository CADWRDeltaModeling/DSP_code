#!/bin/sh
set -e

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)

# CHANGE OUTPUTS DIR -------------------------------------
mkdir -p outputs_tropic
mv outputs/* outputs_tropic
mkdir -p outputs

# CREATE OCEAN BOUNDARY ----------------------------------

cd outputs_tropic
rsync -avz ../interpolate_variables.in .

# link necessary files
ln -sf ../hgrid.gr3 bg.gr3
ln -sf ../hgrid.gr3 fg.gr3
ln -sf ../vgrid.in.2d vgrid.bg
ln -sf ../vgrid.in.3d vgrid.fg

# run script to create uv3d.th.nc
ulimit -s unlimited
interpolate_variables8 # this takes quite a while
cd ../
cp ./outputs_tropic/uv3D.th.nc uv3D.th.nc

# CREATE CLIINIC SYMBOLIC LINKS ----------------------------	

# add new links
ln -sf /SAL_nu_roms.nc SAL_nu.nc
ln -sf /salinity_nudge_roms.gr3 SAL_nudge.gr3
ln -sf /TEM_nu_roms.nc TEM_nu.nc
ln -sf /temperature_nudge_roms.gr3 TEM_nudge.gr3
ln -sf hotstart_{case_year}.nc hotstart.nc

# change existing links
ln -sf bctides.in.3d bctides.in
ln -sf vgrid.in.3d vgrid.in
ln -sf param.nml.clinic param.nml

# model runs with Azure yml file