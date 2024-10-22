study_dir="process_evap_output"
polygons=("legal_delta" "full_domain" "suisun_bay")
# 'all': ['legal_delta', 'suisun_bay', 'grizzly_restoration', 'prospect_restoration', 'lookout_restoration', 'egbert_restoration', 'full_domain'],
# 'baseline': ['legal_delta', 'suisun_bay', 'full_domain'],
# 'suisun': ['legal_delta', 'suisun_bay', 'grizzly_restoration', 'full_domain'],
# 'cache': ['legal_delta', 'suisun_bay', 'prospect_restoration', 'lookout_restoration', 'egbert_restoration', 'full_domain']

for NUMBER in $(seq 1 7); # all cases
do
    for POLY in "${polygons[@]}"; # all cases
    do
        # echo "${AZLINK}${study_dir}/baseline_lhc_${NUMBER}/evap_ts_${POLY}.csv" # for troubleshooting

        azcopy copy "${AZLINK}${study_dir}/baseline_lhc_${NUMBER}/evap_ts_${POLY}.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/process_evap_output/baseline_lhc_${NUMBER}_evap_ts_${POLY}.csv"

    done
done

# Troubleshooting one file:
# azcopy copy "${AZLINK}${study_dir}/baseline_lhc_1/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/process_evap_output" --dry-run
    