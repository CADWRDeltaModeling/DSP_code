#!/bin/sh
set -e

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)

# CHANGE OUTPUTS DIR -------------------------------------
mv outputs outputs_tropic
mkdir -p outputs

# CREATE OCEAN BOUNDARY ----------------------------------

cd outputs_tropic
rsync -avz ../interpolate_variables.in .

# link necessary files
ln -sf ../hgrid.gr3 bg.gr3
ln -sf ../hgrid.gr3 fg.gr3
ln -sf ../../vgrid.in.2d vgrid.bg
ln -sf ../vgrid.in.3d vgrid.fg

# run script to create uv3d.th.nc
module purge
# module load schism/5.11.0 # HPC4
module load intel/2024.0 openmpi/5.0.2 hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1 schism/5.11.0 # HPC5
ulimit -s unlimited
interpolate_variables8 # this takes quite a while
cd ../
ln -sf ./outputs_tropic/uv3D.th.nc uv3D.th.nc

# CREATE CLIINIC SYMBOLIC LINKS ----------------------------	

# add new links
ln -sf salinity_nudge.gr3 SAL_nudge.gr3
ln -sf temperature_nudge.gr3 TEM_nudge.gr3
ln -sf TEM_1.th temp.th
ln -sf SAL_1.th salt.th
ln -sf SAL_nu_roms.nc SAL_nu.nc
ln -sf TEM_nu_roms.nc TEM_nu.nc
ln -sf hotstart/hotstart_2008.nc hotstart.nc

# change existing links
ln -sf bctides.in.3d bctides.in
ln -sf vgrid.in.3d vgrid.in
ln -sf param.nml.clinic param.nml
ln -sf launch.clinic.pbs launch.pbs

# RUN MODEL ----------------------------------------------

# qsub launch.pbs # HPC4
sbatch slurm_batch_file.sh # HPC5