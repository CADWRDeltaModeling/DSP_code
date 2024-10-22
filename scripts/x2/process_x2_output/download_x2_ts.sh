study_dir="process_x2_output"
routes=("bay_sjr" "bay_nysjr" "bay_mzm" "bay_sac")

for NUMBER in $(seq 1 7); # all cases
do
    for ROUTE in "${routes[@]}"; # all cases
    do
        # echo "${AZLINK}${study_dir}/baseline_lhc_${NUMBER}/evap_ts_${POLY}.csv" # for troubleshooting

        # First attempt to copy the file
        azcopy copy "${AZLINK}${study_dir}/baseline_lhc_${NUMBER}/x2_ts_${ROUTE}.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/x2/process_x2_output/baseline_lhc_${NUMBER}_x2_ts_${ROUTE}.csv" #--dry-run

        # Check if the previous command failed
        if [ $? -ne 0 ]; then
            # If it failed, try copying the file with the alternate name
            azcopy copy "${AZLINK}${study_dir}/baseline_lhc_${NUMBER}/x2_ts_${ROUTE}_missing_data.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/x2/process_x2_output/baseline_lhc_${NUMBER}_x2_ts_${ROUTE}_missing_data.csv" #--dry-run
        fi

    done
done

# Troubleshooting one file:
# azcopy copy "${AZLINK}${study_dir}/baseline_lhc_1/x2_ts_bay_sjr.out?${sas}" "/scratch/tomkovic/DSP_code/scripts/x2/process_x2_output" --dry-run
    