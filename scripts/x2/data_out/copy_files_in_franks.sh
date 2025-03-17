# Run line by line

# franks 1
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_1/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_1.csv

# franks 2
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_2/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_2.csv

# franks 3
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_3/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_3.csv

# franks 4
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_4/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_4.csv --dry-run

# franks 5
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_5/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_5.csv --dry-run

# franks 6
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_6/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_6.csv --dry-run

# franks 7
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/franks_lhc_7/x2/x2_ts_bay_sac.csv ./x2_ts_bay_sac_franks_lhc_7.csv --dry-run