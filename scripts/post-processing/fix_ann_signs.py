# Code to fix signs on already generated ANN data. This is paired with a change in the code that makes this step no longer necessary.

import pandas as pd
import datetime as dt

from vtools.functions.filter import cosine_lanczos
from pydelmod.create_ann_inputs import get_dss_data

from schimpy import schism_yaml

import os
import re

import matplotlib.pyplot as plt

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# casanntra dir
cas_dir = '../../../casanntra/data'
csv_fmt = "dsm2_base_{case_num}.csv"

# load lhc_v4
lhc_fn = "../boundary_generation/data_out/lhc_v4.csv"
lhc_df = pd.read_csv(lhc_fn)

# columns_to_flip = ['exports']
# cases_to_fix = range(1,108)
columns_to_flip = []
cases_to_fix = range(1001,1008)

rename = {'scene':{'baseline':'base'}} #{COLUMN_NAME: {STRINGTOREPLACE: REPLACEMENT}}

# now write those into dsm2 dataframes
for index, row in lhc_df.iterrows():
    print(row['case'])
    case_num = re.search(r'(\d+)$', row['case']).group(1)
    if int(case_num) in cases_to_fix:

        casanntra_casefile = os.path.join(cas_dir, csv_fmt.format(case_num=case_num))
        cdf = pd.read_csv(casanntra_casefile,
                        parse_dates=[0],
                        index_col=0)
        
        for flip_col in columns_to_flip:
            cdf[flip_col] = -cdf[flip_col]

        for col, dict in rename.items():
            cdf[col] = cdf[col].replace(dict)
            
        cdf.to_csv(casanntra_casefile, float_format="%.2f", index=True, index_label='datetime')