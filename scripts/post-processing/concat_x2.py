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
                        help='Directory where the csv files from daily x2')
    return parser


def concat_daily_x2(csv_dir, route_ids):

    # loop through relevant polygons (shouldn't have a polygon without cells in the hgrid)
    for route in route_ids:

        # get all csvs for this route
        matching_files = [f for f in os.listdir(
            csv_dir) if re.match(f'x2_{route}_\d+\.out', f)]
        highest_file = max([int(re.search(r'_(\d+)\.out', filename).group(1))
                           for filename in matching_files])
        print(f'Going to read file list for {route} until day {highest_file}.')

        out_df = pd.DataFrame()
        out_df.index = pd.DatetimeIndex([])
        out_fn = os.path.join(csv_dir, f'x2_ts_{route}.csv')

        # print out any missing files
        missing_inds = []
        for f in range(1, highest_file+1):
            day_file = os.path.join(csv_dir, f'x2_{route}_{f}.out')
            if not os.path.exists(day_file):
                missing_inds.append(f)
        if len(missing_inds) > 0:
            raise ValueError(f"There are missing files: {missing_inds}")

        print('Reading the .out files to concatenate to csv output')
        empty_ids = []
        for f in range(1, highest_file+1):
            day_file = os.path.join(csv_dir, f'x2_{route}_{f}.out')
            if os.path.getsize(day_file) == 0:
                empty_ids.append(f)
            else:
                day_in = pd.read_csv(day_file, parse_dates=[
                    0], index_col=[0], header=None)
                out_df = pd.concat([out_df, day_in])
        if len(empty_ids) > 0:
            out_fn = out_fn.replace('.csv', '_missing_data.csv')
            print(f"ERROR: Missing files for route {route} are {empty_ids}")

        out_df.to_csv(out_fn, float_format="%.1f")


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    csv_dir = args.csv_dir

    # Polygon Name options
    route_ids = ['bay_sjr', 'bay_nysjr', 'bay_mzm', 'bay_sac']

    concat_daily_x2(csv_dir, route_ids)


# def main_hardwire():

#     os.chdir(os.path.dirname(__file__))
#     csv_dir = './evap'

#     # Baseline (but input will be fed as "baseline" for polygon_domain)
#     base_polys = ['legal_delta', 'suisun_bay', 'full_domain']

#     concat_daily_evap(csv_dir, base_polys)


if __name__ == "__main__":
    # main_hardwire()
    main()
