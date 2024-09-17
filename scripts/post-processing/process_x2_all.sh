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

add_days_to_start_date()
{
    # Add days to the start date
    # The start date is of the form 2005-01-24
    # The number of days to add is 100
    start_date=$1
    days=$2
    days=$((days-1))
    new_date=$(date -d "$start_date + $days days" +%Y-%m-%d)
    echo $new_date
}

# Function ------------------------------------------
process_x2()
{
    # file_num is the first arg to this script
    # This script assumes it is the outputs directory where the files are located.
    file_num=$1
    date_to_process=$2
    station_bp=$3 # the station.bp filename
    echo "Processing X2 for file number $file_num for date $date_to_process"
    # run for the file_num
    ulimit -s unlimited
    ln -sf ../$station_bp.bp station.bp
    ln -sf ../vgrid.in.3d vgrid.in
    read_output10_xyz <<EOF
1
1
salinity
1
1
$file_num $file_num
1
EOF
    python $BAY_DELTA_SCHISM_HOME/bdschism/bdschism/x2_time_series.py --salt_data_file fort.18 --start $date_to_process --x2route station.bp --output x2_$file_num.csv 
    # we expect 1 line of output per file and we want to append to x2.out
    tail -1 x2_$file_num.csv >> $station_bp.out
    # cleanup
    rm fort.18 x2_$file_num.csv fort.20
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

cd outputs;
echo "Simulation start date from param.nml is $simulation_start_date"

for NUMBER in $(seq 1 $highest); # all salinity
do
    fname="salinity_${NUMBER}.nc" # of the form salinity_100.nc
    
    file_num=`extract_simulation_day $fname`
    date_to_process=`add_days_to_start_date $simulation_start_date $file_num`
    echo "Model date to be processed is $date_to_process"

    echo "Processing Bay-SJR X2 route............................................"
    process_x2 "$file_num" "$date_to_process" x2_bay_sjr

    echo "Processing Bay-NY-SJR X2 route............................................"
    process_x2 "$file_num" "$date_to_process" x2_bay_nysjr

    echo "Processing Bay-Suisun X2 route............................................"
    process_x2 "$file_num" "$date_to_process" x2_bay_mzm

    echo "Processing Bay-Sacramento X2 route............................................"
    process_x2 "$file_num" "$date_to_process" x2_bay_sac
done

    