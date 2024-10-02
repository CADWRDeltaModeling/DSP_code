#!/bin/bash
DOCKER_IMAGE=suxarray:v2024.09.0
docker run --rm -it -v .:/workdir -v /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_4:/data $DOCKER_IMAGE python calculate_depth_average.py --path_study /data

cp -f out2d_1.nc depth_averaged_monthly_salinity_franks.nc;
