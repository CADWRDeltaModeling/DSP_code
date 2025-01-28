# Run line by line

# baseline 1
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_1/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_lhc_1.csv

# baseline 2
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_2/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_lhc_2.csv --dry-run

# baseline 3
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_3/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_lhc_3.csv --dry-run

# baseline 4
azcopy copy "${AZLINK}process_evap_output/baseline_lhc_4/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_baseline_lhc_4.csv" --recursive --preserve-symlinks --dry-run

# baseline 5
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_5/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_lhc_5.csv

# baseline 6
azcopy copy "${AZLINK}process_evap_output/baseline_lhc_6/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_baseline_lhc_6.csv" --recursive --preserve-symlinks --dry-run

# baseline 7
azcopy copy "${AZLINK}process_evap_output/baseline_lhc_7/evap_ts_legal_delta.csv?${sas}" "/scratch/tomkovic/DSP_code/scripts/evap/data_out/evap_ts_legal_delta_baseline_lhc_7.csv" --recursive --preserve-symlinks --dry-run

# mss
cp /scratch/tomkovic/DSP_code/model/schism/mss_th_dsp_geom/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_mss_2021.csv