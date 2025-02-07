import os
import pandas as pd
from dms_datastore.read_ts import read_ts

os.chdir(os.path.dirname(os.path.abspath(__file__)))

obs_dir = "W:/repo/continuous/screened"
out_dir = "./data_out"

obs_fns = {"des_bdl_49_ec_*.csv":"des_bdl_obs_ec_1990_2025.csv",
           "des_cse@upper_22_ec_*.csv":"des_cse_obs_ec_1990_2025.csv",
           "ncro_emm2_b91120_ec_*.csv":"ncro_emm2_obs_ec_1990_2025.csv",
           "usbr_jer_jer_ec_*.csv":"usbr_jer_obs_ec_1990_2025.csv",
           "dwr_rsl_rsl_ec_*.csv":"dwr_rsl_obs_ec_1990_2025.csv"}

for obs_fn_in, obs_fn_out in obs_fns.items():
    df = read_ts(os.path.join(obs_dir,obs_fn_in),
                start=pd.to_datetime('1990-1-1 00:00'),
                end=pd.to_datetime('2025-02-05 00:00'))

    df.to_csv(os.path.join(out_dir, obs_fn_out), float_format="%.2f", index=True, index_label='datetime')
