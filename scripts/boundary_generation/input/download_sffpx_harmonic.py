import os
import pandas as pd
from dms_datastore.read_ts import read_ts

os.chdir(os.path.dirname(os.path.abspath(__file__)))

df = read_ts("W:/repo/continuous/screened/noaa_sffpx_9414290_predictions_*.csv",
             start=pd.to_datetime('1990-1-1 00:00'),
             end=pd.to_datetime('2025-02-05 00:00'))

df.to_csv('./noaa_sffpx_9414290_prediction_1990_2025.csv', float_format="%.2f", index=True, index_label='datetime')