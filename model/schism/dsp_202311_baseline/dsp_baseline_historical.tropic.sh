#!/bin/sh

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)

# MAKE SFLUX LINKS --------------------------------------

mkdir -p sflux

# remove existing links
find ./sflux/ -name '*.nc' -delete

# make new links
cd sflux
python ../make_links_historical.py
cd ..

# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
gen_elev2d --outfile base.hist.elev2D.th.nc --hgrid=hgrid_dsp_baseline.gr3 --stime=2009-1-1 --etime=2009-7-1 --slr 0.0 noaa_download/pryc_filled_data.csv noaa_download/noaa_mtyc1_9413450_water_level_2006_2017.csv
ln -sf base.hist.elev2D.th.nc elev2d.th.nc

# CREATE OTHER SYMBOLIC LINKS ----------------------------
ln -sf bctides.in.2d bctides.in
ln -sf hgrid_dsp_baseline.gr3 hgrid.gr3
ln -sf vgrid.in.2d vgrid.in
ln -sf param.nml.tropic param.nml
ln -sf launch.tropic.pbs launch.pbs # only necessary with HPC4

# CHECK OUTPUT DIRECTORY IS PRESENT ----------------------

mkdir -p outputs

# RUN MODEL ----------------------------------------------

qsub launch.pbs # HPC4
# sbatch slurm_batch_file.sh # HPC5