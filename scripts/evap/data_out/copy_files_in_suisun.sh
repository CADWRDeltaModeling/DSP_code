# Run line by line

# suisun 1
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/suisun_lhc_1/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_suisun_lhc_1.csv --dry-run

# suisun 2
azcopy copy "${AZLINK}process_evap_output/suisun_lhc_2/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_suisun_lhc_2.csv" --recursive --preserve-symlinks --dry-run

# suisun 3
azcopy copy "${AZLINK}process_evap_output/suisun_lhc_3/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_suisun_lhc_3.csv" --recursive --preserve-symlinks --dry-run

# suisun 4
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/suisun_lhc_4/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_suisun_lhc_4.csv

# suisun 5
azcopy copy "${AZLINK}process_evap_output/suisun_lhc_5/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_suisun_lhc_5.csv" --recursive --preserve-symlinks --dry-run

# suisun 6
azcopy copy "${AZLINK}process_evap_output/suisun_lhc_6/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_suisun_lhc_6.csv" --recursive --preserve-symlinks --dry-run

# suisun 7
azcopy copy "${AZLINK}process_evap_output/suisun_lhc_7/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_suisun_lhc_7.csv" --recursive --preserve-symlinks --dry-run