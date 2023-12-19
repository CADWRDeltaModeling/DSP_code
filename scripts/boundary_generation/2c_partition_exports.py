import os

import numpy as np
import pandas as pd
from pydelmod.create_ann_inputs import get_dss_data
import pyhecdss


output_dir = r"D:\projects\delta_salinity\scripts\DSP_code\scripts\boundary_generation\data_out"

# the modified input file generated from step 2a
mod_in_file = os.path.join(output_dir, "Exports_2006-2016_perturb_historical.csv")

dsp_home = "D:/projects/delta_salinity" # LAT Computer: D:/projects/delta_salinity, EA computer: f:/projects/ann_dsp

# Define DSS data/parameters
in_dss_dir = f"{dsp_home}/model/dsm2/2021DSM2FP_202301/timeseries"
hist_dss = "hist.dss"
hist_dss_file = os.path.join(in_dss_dir, hist_dss)

# Get pathnames function
def get_pathname(dss_filename, b_part, c_part, e_part=None, f_part=None, filter_b_part_numeric=None):
    with pyhecdss.DSSFile(dss_filename) as d:
        catdf = d.read_catalog()
        dss_file_parts = dss_filename.split('/')
        dfilename = dss_file_parts[len(dss_file_parts)-1]
        filtered_df = None
        if b_part is not None:
            filtered_df = filtered_df[(catdf.B==b_part)] if filtered_df is not None else catdf[(catdf.B==b_part)]
        if c_part is not None:
            filtered_df = filtered_df[(catdf.C==c_part)] if filtered_df is not None else catdf[(catdf.C==c_part)]
        if e_part is not None:
            filtered_df = filtered_df[(catdf.E==e_part)] if filtered_df is not None else catdf[(catdf.E==e_part)]
        if f_part is not None:
            filtered_df = filtered_df[(catdf.F==f_part)] if filtered_df is not None else catdf[(catdf.F==f_part)]
        if filter_b_part_numeric:
            filtered_df = filtered_df[(catdf.B.str.isnumeric())]
        path_list = d.get_pathnames(filtered_df)

    return path_list

# access dss file for export parts
primary_pathname_part_dss_filename_dict = {'CHSWP003': hist_dss_file, 
                                           'CHDMC004': hist_dss_file, 
                                           'CHCCC006': hist_dss_file, 
                                           'ROLD034': hist_dss_file, 
                                           'CHVCT001': hist_dss_file}
primary_part_c_part_dict = {'CHSWP003': 'FLOW-EXPORT', 
                            'CHDMC004': 'FLOW-EXPORT', 
                            'CHCCC006': 'FLOW-DIVERSION', 
                            'ROLD034': 'FLOW-EXPORT', 
                            'CHVCT001': 'FLOW-EXPORT'}
unit_part_dict = {'CHSWP003': 'CFS',
                  'CHDMC004': 'CFS',
                  'CHCCC006': 'CFS',
                  'ROLD034': 'CFS',
                  'CHVCT001': 'CFS'}
primary_pathname_part = 'b_part'

df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part, 
                        primary_part_c_part_dict=primary_part_c_part_dict,
                        primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)

exports_list = ['CHSWP003', 'CHDMC004', 'CHCCC006', 'ROLD034', 'CHVCT001']

# read in modified timeseries
mod_in_ts = pd.read_csv(mod_in_file, header=None, index_col=0, names=['Modified'])
mod_in_ts.set_index(pd.DatetimeIndex(mod_in_ts.index).to_period(freq='D'), inplace=True)
mod_in_ts.Modified = mod_in_ts.Modified.astype(float)

# scale flows per day
in_ts = df_input[exports_list]
scale_in_ts = in_ts.copy()
scale_in_ts = scale_in_ts.apply(lambda x: x.div(x.sum()), axis=1)

# Join the scaled dataframe with the modified input timeseries 
merge = pd.merge(mod_in_ts,scale_in_ts, how='inner', left_index=True, right_index=True)
merge.iloc[:, 1:] =  merge.iloc[:, 1:].multiply(merge.iloc[:,0], axis='index') # distribute the modified timeseries across columns
merge.index = merge.index.to_timestamp()

# Loop through components and write out to DSS
for b_part in exports_list:
    in_dss = primary_pathname_part_dss_filename_dict[b_part]
    pathname = get_pathname(in_dss, b_part, primary_part_c_part_dict[b_part])[0]
    
    d_part = pathname.split('/')[1:7][3]
    e_part = pathname.split('/')[1:7][4]
    f_part = pathname.split('/')[1:7][5]
    
    with pyhecdss.DSSFile(in_dss) as d_in:
        df = None
        units = None
        ptype = None
        df, units, ptype = d_in.read_rts(pathname)
        
        out_ts = merge[[b_part]].copy()
        out_ts.columns = ['Modified']
        # add historical timeseries to beginning and end using margin_ts (copy of df)
        if 'Timestamp' not in str(type(out_ts.index[0])):
            out_ts.index = out_ts.index.to_timestamp()
        margin_ts = df.copy()
        if 'Timestamp' not in str(type(margin_ts.index[0])):
            margin_ts.index = margin_ts.index.to_timestamp()
        margin_ts = margin_ts[(margin_ts.index < out_ts.index[0]) | (margin_ts.index > out_ts.index[-1])]
        margin_ts.columns = ['Modified']
        out_ts = pd.concat([out_ts, margin_ts])
        out_ts.sort_index(inplace=True)
        # out_ts.index = out_ts.index + pd.Timedelta(days=1) # need to shift by one day because of DSS writing timestamp issues

    print(f'Writing out {b_part}')
    out_ts.to_csv(os.path.join(output_dir,f"{b_part}_perturb_historical.csv"), 
                    index=True,
                    header=False)
