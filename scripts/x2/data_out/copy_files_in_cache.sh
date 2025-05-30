# Run line by line

# cache 1
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_1/x2/x2_ts_bay_sac_missing_data.csv ./x2_ts_bay_sac_cache_lhc_1.csv

# cache 2
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_2/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_cache_lhc_2.csv

# cache 3
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_3/x2/x2_ts_bay_sac_missing_data.csv ./x2_ts_bay_sac_cache_lhc_3.csv --dry-run

# cache 4
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_4/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_cache_lhc_4.csv --dry-run

# cache 5
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_5/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_cache_lhc_5.csv --dry-run

# cache 6
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_6/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_cache_lhc_6.csv --dry-run

# cache 7
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/cache_lhc_7/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_cache_lhc_7.csv --dry-run