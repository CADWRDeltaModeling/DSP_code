#!/bin/bash
DOCKER_IMAGE=suxarray:v2024.09.0
docker run --rm -it -v .:/workdir -v /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/suisun_lhc_4:/data $DOCKER_IMAGE python calculate_depth_average.py --path_study /data

cp -f /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/suisun_lhc_4/outputs/zCoordinates_1.nc ./;
