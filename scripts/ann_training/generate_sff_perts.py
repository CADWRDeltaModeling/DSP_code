# DSM2-output-specific code to generate the SFFPX tidal energy and filter data needed for the ANN

import pandas as pd
import datetime as dt

from vtools.functions.filter import cosine_lanczos
from pydelmod.create_ann_inputs import get_dss_data

from schimpy import schism_yaml

import os
import re

import matplotlib.pyplot as plt

import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# shifted subtide
def shift_subtide(df, shift_val, out_col):
    filt_df = cosine_lanczos(df.copy(),cutoff_period='40H', padtype='odd')
    shift_df = filt_df.copy()
    shift_df.index = shift_df.index + pd.Timedelta(days=shift_val)
    # shift_df = shift_df.to_frame()
    shift_df.columns = ['shift_filter']

    shift_df['inst'] = df
    shift_df['eg_filter'] = filt_df

    shift_df[out_col] = shift_df['inst'] - shift_df['eg_filter'] + shift_df['shift_filter']
    
    return shift_df[[out_col]]

# shifted instantaneous
def shift_ts(df, shift_val, out_col):
    shift_df = df.copy()
    shift_df.index = shift_df.index + pd.Timedelta(days=shift_val)
    # shift_df = shift_df.to_frame()
    shift_df.columns = [out_col]
    
    return shift_df[[out_col]]

def calc_filt_nrg(df):
    filt_df = cosine_lanczos(df.copy(),cutoff_period='40H', padtype='odd') # <z>
    nrg_df = cosine_lanczos((df - filt_df)**2, cutoff_period ='40H', padtype='odd') # = < (z- <z>)^2 >
    filt_df.columns = ['sf_tidal_filter']
    nrg_df.columns = ['sf_tidal_energy']
    filt_df = filt_df.resample('D', closed='right').mean()
    nrg_df = nrg_df.resample('D', closed='right').mean()

    out_df = pd.merge(nrg_df, filt_df, left_index=True, right_index=True, how='outer')

    return out_df
    
sf_fn = "../boundary_generation/input/9414290_gageheight.txt"
sf_raw_ts = pd.read_csv(sf_fn, index_col=0, parse_dates=[0], 
                        usecols=["Date Time", "Water Level"], 
                        skipinitialspace=True,
                        dtype={"Water Level": str},
                        sep=',', comment='#')
sf_raw_ts = sf_raw_ts.dropna()
sf_raw_ts["Water Level"] = pd.to_numeric(sf_raw_ts["Water Level"], errors="coerce")
sf_raw_ts = sf_raw_ts.resample('15min').ffill()
sf_raw_ts.index.freq = pd.infer_freq(sf_raw_ts.index)

# load lhc_v4
# lhc_fn = "../boundary_generation/data_out/lhc_v3.csv"
# case_nums = range(1001,1008)
# case_nums = [1001]

lhc_fn = "../boundary_generation/data_out/lhc_v4.csv"
# case_nums=range(1,108)
case_nums=range(107,108)
# case_nums=[5,12,17,21,27,33,34,41,42,44,48,49,50,55,70,71,72,73,75,78,82,88,89,98,99] # shifted - 50d cases
# case_nums=[55,70,71,72,73,75,78,82,88,89,98,99]

run_wait = False
lhc_df = pd.read_csv(lhc_fn)

# casanntra dir
cas_dir = '../../../casanntra/data'
csv_fmt = "dsm2_base_{case_num}.csv"

# Regular
reg_out_df = calc_filt_nrg(sf_raw_ts)

# Shifted - 50d
shift_minus50_df = shift_ts(sf_raw_ts, -50, 'shift-50')
shift_minus50_df = calc_filt_nrg(shift_minus50_df)

# Shifted + 100d
shift_plus100_df = shift_ts(sf_raw_ts, 100, 'shift+100')
shift_plus100_df = calc_filt_nrg(shift_plus100_df)

# Subtidal Pert
pert_fn = "../boundary_generation/data_out/perturb_historical_subtide_v1.csv"
mrz_fn = "../boundary_generation/input/dsm2_mrz_stage_ft_15min_clean.csv "
mrz_filt_fn = "../boundary_generation/data_out/dsm2_mrz_stage_ft_15min_clean_filtered.csv"
mrz_pert_df = pd.read_csv(pert_fn, index_col=[0], parse_dates=[0], header=None)
mrz_df = pd.read_csv(mrz_fn, index_col=[0], parse_dates=[0])
mrz_filt_df = pd.read_csv(mrz_filt_fn, index_col=[0], parse_dates=[0])
mrz_tidal_df = mrz_df - mrz_filt_df
mrz_subtide_pert = mrz_pert_df[1] - mrz_tidal_df['stage (ft)']

sf_tidal_df = sf_raw_ts["Water Level"] - reg_out_df['sf_tidal_filter'].resample('15min').ffill()
sf_tidal_df = sf_tidal_df.dropna()

sf_subtide_pert = sf_tidal_df + mrz_subtide_pert
sf_subtide_pert = sf_subtide_pert.dropna()
sf_subtide_pert = sf_subtide_pert.resample('15min').ffill()
sf_subtide_pert.index.freq = pd.infer_freq(sf_subtide_pert.index)
sf_subtide_pert = sf_subtide_pert.to_frame()
sf_subtide_pert.columns = ['stage']

subtide_pert_df = calc_filt_nrg(sf_subtide_pert)

# Shifted Subtide + 7 days
shift_subtide_plus7_df = shift_subtide(sf_raw_ts, 7, 'stage')
shift_subtide_plus7_df = calc_filt_nrg(shift_subtide_plus7_df)

# Shifted Subtide - 7 days
shift_subtide_minus7_df = shift_subtide(sf_raw_ts, -7, 'stage')
shift_subtide_minus7_df = calc_filt_nrg(shift_subtide_minus7_df)

# Subtidal Pert + Shifted Subtide + 7 days
subtide_pert_shift_plus7_df = shift_subtide(sf_subtide_pert, 7, 'stage')
subtide_pert_shift_plus7_df = calc_filt_nrg(subtide_pert_shift_plus7_df)

# Subtidal Pert + Shifted Subtide - 7 days
subtide_pert_shift_minus7_df = shift_subtide(sf_subtide_pert, -7, 'stage')
subtide_pert_shift_minus7_df = calc_filt_nrg(subtide_pert_shift_minus7_df)

col_order = ['model','scene','case',
             'northern_flow','sac_flow','sjr_flow','exports','cu_flow','ndo',
             'dcc','smscg',
             'vern_ec',
             'mrz_tidal_energy','mrz_tidal_filter',
             'sf_tidal_energy','sf_tidal_filter', 
             'anc','anh','bac','bdl','bdt','bet','cll','cse',
             'dsj','emm2','frk','god','gys','gzl','hll','hol2',
             'ibs','jer','mal','mtz','nsl2','obi','oh4','old',
             'pct','ppt','rri2','rsl','sal','snc','srv','sss',
             'tms','trp','tss','uni','vcu','vol','wci',
             'x2']

# now write those into dsm2 dataframes
for index, row in lhc_df.iterrows():
    case_num = re.search(r'(\d+)$', row['case']).group(1)
    if int(case_num) in case_nums:
        print(row['case'])
        casanntra_casefile = os.path.join(cas_dir, csv_fmt.format(case_num=case_num))

        # to run this in parallel with ongoing/overnight check if the model is finished running
        if run_wait:
            next_ann_fn = casanntra_casefile.replace(f'_{case_num}', f'_{int(case_num)+1}')
            while not os.path.exists(next_ann_fn):
                print(f"Waiting for file {next_ann_fn} to appear so that case {case_num} is done...")
                time.sleep(120)  # Wait for 30 seconds before checking again

            print(f"File {next_ann_fn} is now available! Post-processing model results for case {case_num}")
        
        cdf = pd.read_csv(casanntra_casefile,
                        parse_dates=[0],
                        index_col=0)
        
        
        if 'sf_tidal_filter' in cdf.columns:
            cdf.drop(columns=['sf_tidal_filter','sf_tidal_energy'], inplace=True)
        
        if row['tide'] == 'Shifted + 100d':
            merge_df = pd.merge(cdf, shift_plus100_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')

        elif row['tide'] == 'Subtidal Pert':
            merge_df = pd.merge(cdf, subtide_pert_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')
            
        elif row['tide'] == 'Regular':
            merge_df = pd.merge(cdf, reg_out_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')
            
        elif row['tide'] == 'Shifted - 50d':
            merge_df = pd.merge(cdf, shift_minus50_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')
            
        elif row['tide'] == 'Shifted Subtide + 7 days':
            merge_df = pd.merge(cdf, shift_subtide_minus7_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')
            
        elif row['tide'] == 'Shifted Subtide - 7 days':
            merge_df = pd.merge(cdf, shift_subtide_plus7_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')
            
        elif row['tide'] == 'Subtidal Pert + Shifted Subtide + 7 days':
            merge_df = pd.merge(cdf, subtide_pert_shift_plus7_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')
            
        elif row['tide'] == 'Subtidal Pert + Shifted Subtide - 7 days':
            merge_df = pd.merge(cdf, subtide_pert_shift_minus7_df.loc[cdf.index[0]:cdf.index[-1]], left_index=True, right_index=True, how='outer')

        else:
            print("NO METHOD")
            
        merge_df = merge_df[[col for col in col_order if col in merge_df.columns]]
        if isinstance(merge_df.index, pd.DatetimeIndex):
            merge_df.index = merge_df.index.to_period('D')
        merge_df.to_csv(casanntra_casefile, float_format="%.2f", index=True, index_label='datetime')