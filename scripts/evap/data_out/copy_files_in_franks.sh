# Run line by line

# franks 1
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_1/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_1.csv

# franks 2
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_2/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_2.csv

# franks 3
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_3/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_3.csv

# franks 4
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_4/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_4.csv --dry-run

# franks 5
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_5/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_5.csv --dry-run

# franks 6
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_6/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_6.csv --dry-run

# franks 7
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_7/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_franks_lhc_7.csv --dry-run
