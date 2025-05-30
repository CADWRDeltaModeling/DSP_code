# Run line by line

# cache 1
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_1/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_1.csv

# cache 2
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_2/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_2.csv

# cache 3
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_3/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_3.csv --dry-run

# cache 4
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_4/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_4.csv --dry-run

# cache 5
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_5/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_5.csv --dry-run

# cache 6
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_6/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_6.csv --dry-run

# cache 7
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_7/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_cache_lhc_7.csv --dry-run