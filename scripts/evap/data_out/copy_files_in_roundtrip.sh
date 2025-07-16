
# suisun-base
cp /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_suisun-base.csv

# suisun-suisun
cp /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_suisun-suisun.csv

# slr-base
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/slr-base/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_slr-base.csv --dry-run

# slr-slr
rsync -avz tomkovic@10.3.80.41:/scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/slr-slr/evap/evap_ts_legal_delta.csv ./evap_ts_legal_delta_baseline_slr-slr.csv --dry-run
