import pandas as pd
import datetime as dt

from vtools.functions.filter import cosine_lanczos
from pydelmod.create_ann_inputs import get_dss_data

from schimpy import schism_yaml

import os
import re

os.chdir(os.path.dirname(os.path.abspath(__file__)))

casanntra_dir = "../../../casanntra/data"
cases = range(1,101)
csv_fmt = "dsm2_base_{case_num}.csv"
in_fname = "./input/ann_config_lathypcub_v4_dsm2.yaml"
with open(in_fname, 'r') as f:
    # loader = RawLoader(stream)
    inputs = schism_yaml.load(f)
experiment = inputs.get('experiment').format(**locals())
model_folder = inputs.get('model_folder').format(**locals())
case_setup = pd.read_csv(inputs.get('case_setup'))

col_order = ['model','scene','case',
             'sac_flow','sjr_flow','exports','cu_flow','ndo',
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

for index, row in case_setup.iterrows():
    print(row['case'])
    case_num = re.search(r'(\d+)$', row['case']).group(1)

    casanntra_casefile = os.path.join(casanntra_dir, csv_fmt.format(case_num=case_num))
    cdf = pd.read_csv(casanntra_casefile,
                      parse_dates=[0],
                      index_col=0)
    
    hist_dss_file = inputs.get('hist_dss_file').format(**locals())
    
    b_part_dss_filename_dict={'RSAC054': hist_dss_file}
    b_part_c_part_dict={'RSAC054': 'STAGE'}
    df_mtz_stage = get_dss_data(b_part_dss_filename_dict, 'b_part', \
        primary_part_c_part_dict=b_part_c_part_dict, daily_avg=False)
    
    df_sff_stage = df_mtz_stage.copy()
    df_sff_stage.index = df_sff_stage.index - pd.Timedelta(hours=1.5)

    # df_sff_stage = df_sff_stage.loc[pd.to_datetime(row['start'])-pd.Timedelta(days=30):pd.to_datetime(row['end']),:] # allow 30 days to lose track of initial sf (t+2hours)   
    # cutoff = df_sff_stage.index[-1] - pd.Timedelta(hours=2)

    # # Filter rows starting from `start_time` and iterate backwards
    # for idx, rw in df_sff_stage.loc[:cutoff].iterrows():
    #     # sf(t) = 2 * mrz(t + 2hr) - sf(t + 1hr)
    #     # df_sff_stage.loc[idx:idx+pd.Timedelta(hours=2)]
    #     # df_mtz_stage.loc[idx:idx+pd.Timedelta(hours=2)]
    #     df_sff_stage.loc[idx] = 2*df_mtz_stage.loc[idx + pd.Timedelta(hours=2)] - df_sff_stage.loc[idx + pd.Timedelta(hours=1)]
    #     # df_sff_stage.loc[idx:idx+pd.Timedelta(hours=2)]

    df_filter = cosine_lanczos(df_sff_stage.copy(), cutoff_period ='40H', padtype='odd')
    df_nrg = cosine_lanczos((df_sff_stage-df_filter)**2, cutoff_period ='40H', padtype='odd') # = < (z- <z>)^2 >
    df_sff_tidal_filter = df_filter.resample('D', closed='right').mean()
    df_sff_tidal_filter.columns=['tidal_filter']
    df_sff_tidal_filter.index = df_sff_tidal_filter.index.to_timestamp()
    df_sff_tidal_energy = df_nrg.resample('D', closed='right').mean()
    df_sff_tidal_energy.columns=['tidal_energy']
    df_sff_tidal_energy.index = df_sff_tidal_energy.index.to_timestamp()

    cdf['sf_tidal_energy'] = df_sff_tidal_energy['tidal_energy']
    cdf['sf_tidal_filter'] = df_sff_tidal_filter['tidal_filter']

    cdf = cdf[[col for col in col_order if col in cdf.columns]]

    cdf.to_csv(casanntra_casefile, float_format="%.2f", index=True, index_label='datetime')
