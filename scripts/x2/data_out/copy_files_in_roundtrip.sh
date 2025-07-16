# Run line by line

# suisun-base
cp /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base/x2/x2_ts_bay_sac_missing_data.csv ./x2_ts_bay_sac_baseline_suisun-base.csv

# suisun-suisun
cp /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun/x2/x2_ts_bay_sac_missing_data.csv ./x2_ts_bay_sac_suisun_suisun-suisun.csv

# slr-base
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/slr-base/x2/x2_ts_bay_sac_missing_data.csv ./x2_ts_bay_sac_baseline_slr-base.csv --dry-run

# slr-slr
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/slr-slr/x2/x2_ts_bay_sac_missing_data.csv ./x2_ts_bay_sac_baseline_slr-slr.csv --dry-run