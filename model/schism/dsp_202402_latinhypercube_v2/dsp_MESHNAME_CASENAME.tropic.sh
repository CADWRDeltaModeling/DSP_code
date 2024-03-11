#!/bin/sh
set -e

# LOAD CONDA ENVIRONMENT PRIOR TO RUNNING THIS BASH SCRIPT
# conda activate schism-dms (for LAT HPC4)

# MAKE SFLUX LINKS --------------------------------------

mkdir -p sflux

# remove existing links
find ./sflux/ -name '*.nc' -delete

# make new links
cd sflux
python ../../make_links_{casename}.py
cd ..

# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
gen_elev2d --outfile {meshname}.{casename}.elev2D.th.nc --hgrid=hgrid_dsp_{meshname}.gr3 --stime=2006-11-14 --etime=2008-12-31 --slr 0.0 ../../dsp_202311_baseline/noaa_download/pryc_filled_data.csv ../../dsp_202311_baseline/noaa_download/noaa_mtyc1_9413450_water_level_2006_2017.csv
ln -sf {meshname}.{casename}.elev2D.th.nc elev2D.th.nc

# CREATE OTHER SYMBOLIC LINKS ----------------------------
# shared inputs
ln -sf ../../sflux/sflux_inputs.txt ./sflux/sflux_inputs.txt
ln -sf ../schism.sh schism.sh
ln -sf ../station.in station.in

# shared spatial inputs
ln -sf ../bctides.in.2d bctides.in
ln -sf ../vgrid.in.2d vgrid.in

# shared TH inputs
ln -sf ../ccfb_gate.th ccfb_gate.th
ln -sf ../delta_cross_channel.th delta_cross_channel.th
ln -sf ../flux.th flux.th
ln -sf ../grantline_barrier.th grantline_barrier.th
ln -sf ../grantline_culvert.th grantline_culvert.th
ln -sf ../grantline_weir.th grantline_weir.th
ln -sf ../midr_culvert_l.th midr_culvert_l.th
ln -sf ../midr_culvert_r.th midr_culvert_r.th
ln -sf ../midr_weir.th midr_weir.th
ln -sf ../montezuma_boat_lock.th montezuma_boat_lock.th
ln -sf ../montezuma_flash.th montezuma_flash.th
ln -sf ../montezuma_radial.th montezuma_radial.th
ln -sf ../msource.th msource.th
ln -sf ../oldr_head_barrier.th oldr_head_barrier.th
ln -sf ../oldr_tracy_culvert.th oldr_tracy_culvert.th
ln -sf ../oldr_tracy_weir.th oldr_tracy_weir.th
ln -sf ../SAL_1.th SAL_1.th
ln -sf ../TEM_1.th TEM_1.th
ln -sf ../tom_paine_sl_culvert.th tom_paine_sl_culvert.th
ln -sf ../vsink.th vsink.th
ln -sf ../vsource.th vsource.th
ln -sf ../west_false_river_barrier_leakage.th west_false_river_barrier_leakage.th

# inputs specific to this geometry
ln -sf hgrid_dsp_{meshname}.gr3 hgrid.gr3
ln -sf hgrid_{meshname}.ll hgrid.ll
ln -sf hydraulics_{meshname}.in hydraulics.in
ln -sf source_sink_{meshname}.in source_sink.in
ln -sf fluxflag_{meshname}.prop fluxflag.prop

# inputs specific to this setup
ln -sf param.nml.tropic param.nml
ln -sf launch.tropic.pbs launch.pbs # only necessary with HPC4

# CHECK OUTPUT DIRECTORY IS PRESENT ----------------------

mkdir -p outputs

# RUN MODEL ----------------------------------------------

qsub launch.pbs # HPC4
# sbatch slurm_batch_file.sh # HPC5