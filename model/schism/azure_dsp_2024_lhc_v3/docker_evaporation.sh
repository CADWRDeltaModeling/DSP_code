#!/bin/env bash
extract_simulation_day()
{
   # Extract the simulation day from the filename
    # The filename is of the form prefix_salinity_100.nc
    # The simulation day is 100
    fname=$1
    # Remove the directory part of the filename
    fname=$(basename $fname)
    # Remove the prefix and suffix, leaving only the number
    fname=${fname#*_}
    fname=${fname%.*}
    echo $fname
}

# Function ---------------------------------------------
get_start_date()
{
    param="$(readlink -f $1)"
    start_year=$(cat $param | grep 'start_year = ' | grep -Eo '[0-9]*')
    start_month=$(cat $param | grep 'start_month = ' | grep -Eo '[0-9]*')
    start_day=$(cat $param | grep 'start_day = ' | grep -Eo '[0-9]*')

    echo "${start_year}-${start_month}-${start_day}"
}

# Run all ---------------------------------------------------

simulation_start_date=`get_start_date param.nml` # of the form 2005-01-24

# Get highest salinity file in outputs
highest=-1
for file in outputs/salinity_*.nc
do
  if [[ $file =~ salinity_([0-9]+).* ]]
  then
    [[ "${BASH_REMATCH[1]}" -gt "$highest" ]] && highest=${BASH_REMATCH[1]}
  fi
done

echo "Simulation start date from param.nml is $simulation_start_date"

for NUMBER in $(seq 1 $highest);
do
  #!/bin/bash
  DOCKER_IMAGE=suxarray:v2024.09.0
  docker run --rm -it -v .:/workdir -v ./:/data $DOCKER_IMAGE python evaporation_integrate.py --path_out2d /data/outputs/out2d_${NUMBER}.nc --path_polygons /data/evap/evap_poly.shp --polygon_domain baseline

done

mkdir -p evap
mv evap*.csv evap/
    