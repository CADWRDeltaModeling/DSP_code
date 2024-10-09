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
    poly_file=$3 # the polygon yaml file
    evap_csv_out=$4 # the x2 csv filename that get used as the .out file
    echo "Processing Evaporation for file number $file_num for date $start_date"
    # run for the file_num
    ulimit -s unlimited
    ln -sf $poly_file evap_polygons.yaml
    read_output10_xyz <<EOF
1
1
salinity
1
1
$file_num $file_num
1
EOF
    python #TODO: Add calls here
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

# This should call the script which subsets and summarizes evap polygons into csv file(s)