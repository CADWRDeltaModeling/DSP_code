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
    new_date=$(date -d "$start_date + $days days" +%Y-%m-%d)
    echo $new_date
}



# Function ------------------------------------------
process_x2()
{
    # file_num is the first arg to this script
    # This script assumes it is the outputs directory where the files are located.
    file_num=$1
    start_date=$2
    station_bp=$3 # the station.bp filename
    x2_csv_out=$4 # the x2 csv filename that get used as the .out file
    rid=$5 # sjr, sac, or mzm for extract_x2_station_xyz>x2_route2_bp
    echo "Processing X2 for file number $file_num for date $start_date"
    # run for the file_num
    ulimit -s unlimited
    ln -sf $station_bp.bp station.bp
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
    python $BAY_DELTA_SCHISM_HOME/bdschism/bdschism/x2_time_series.py --salt_data_file fort.18 --start $start_date --x2route $x2_csv_out.csv --east_route_id $rid --output x2_temp.csv 
    # we expect 1 line of output per file and we want to append to x2.out
    tail -1 x2_temp.csv >> $x2_csv_out.out
    # cleanup
    rm fort.18 x2_temp.csv fort.20
}
simulation_start_date="$1" # of the form 2005-01-24
fname="$2" # of the form salinity_100.nc
#
file_num=`extract_simulation_day $fname`
date_to_process=`add_days_to_start_date $simulation_start_date $file_num`

echo "Processing Bay-SJR X2 route............................................"
process_x2 "$file_num" "$date_to_process" x2_bay_sjr_station x2route_bay_sjr sjr

echo "Processing Bay-NY-SJR X2 route............................................"
process_x2 "$file_num" "$date_to_process" x2_bay_ny_sjr_station x2route_bay_ny_sjr sjr

echo "Processing Bay-Suisun X2 route............................................"
process_x2 "$file_num" "$date_to_process" x2_bay_suisun_station x2route_bay_suisun mzm

echo "Processing Bay-Sacramento X2 route............................................"
process_x2 "$file_num" "$date_to_process" x2_bay_sac_station x2route_bay_sac sac