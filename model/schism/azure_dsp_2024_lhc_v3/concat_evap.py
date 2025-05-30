from pathlib import Path
import re
import numpy as np
import xarray as xr
import geopandas as gpd
import suxarray as sx
import pandas as pd

import argparse
from dateutil.parser import parse

import re
import os


def create_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_dir', default=None, type=Path,
                        help='')
    parser.add_argument('--polygon_domain', default=None,  type=str,
                        help='')
    return parser


def concat_daily_evap(csv_dir, poly_ids):

    # loop through relevant polygons (shouldn't have a polygon without cells in the hgrid)
    for poly in poly_ids:

        # get all csvs for this domain
        matching_files = [f for f in os.listdir(
            csv_dir) if re.match(f'evap_{poly}_\d+\.csv', f)]
        highest_file = max([int(re.search(r'_(\d+)\.csv', filename).group(1))
                           for filename in matching_files])
        print(f'Going to read file list for {poly} until day {highest_file}.')

        out_df = pd.DataFrame()
        out_df.index = pd.DatetimeIndex([])
        out_fn = os.path.join(csv_dir, f'evap_ts_{poly}.csv')

        # print out any missing files
        missing_inds = []
        for f in range(1, highest_file+1):
            day_file = os.path.join(csv_dir, f'evap_{poly}_{f}.csv')
            if not os.path.exists(day_file):
                missing_inds.append(f)
        if len(missing_inds) > 0:
            raise ValueError(f"There are missing files: {missing_inds}")

        for f in range(1, highest_file+1):
            day_in = pd.read_csv(os.path.join(
                csv_dir, f'evap_{poly}_{f}.csv'), parse_dates=[0], index_col='time')
            out_df = pd.concat([out_df, day_in])

        out_df.to_csv(out_fn, float_format="%.1f")


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    csv_dir = args.csv_dir
    polygon_domain = args.polygon_domain

    # Polygon Name options
    polygon_dict = {'all': ['legal_delta', 'suisun_bay', 'grizzly_restoration',
                            'prospect_restoration', 'lookout_restoration', 'egbert_restoration', 'full_domain'],
                    'baseline': ['legal_delta', 'suisun_bay', 'full_domain'],
                    'suisun': ['legal_delta', 'suisun_bay',
                               'grizzly_restoration', 'full_domain'],
                    'cache': ['legal_delta', 'suisun_bay', 'prospect_restoration',
                              'lookout_restoration', 'egbert_restoration', 'full_domain']}

    poly_ids = polygon_dict[polygon_domain]

    concat_daily_evap(csv_dir, poly_ids)


# def main_hardwire():

#     os.chdir(os.path.dirname(__file__))
#     csv_dir = './evap'

#     # Baseline (but input will be fed as "baseline" for polygon_domain)
#     base_polys = ['legal_delta', 'suisun_bay', 'full_domain']

#     concat_daily_evap(csv_dir, base_polys)


if __name__ == "__main__":
    # main_hardwire()
    main()
