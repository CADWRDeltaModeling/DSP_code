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
    parser.add_argument("--path_out2d", default=None, type=Path, help="")
    parser.add_argument("--path_polygons", default=None, type=Path, help="")
    parser.add_argument("--polygon_domain", default=None, type=str, help="")
    return parser


def extract_evap_daily(
    path_out2d, path_polygons, poly_ids, varname="evaporationRate", chunks={"time": 12}
):

    # Read the out2d files and catalog number
    ds_out2d = xr.open_mfdataset(
        path_out2d,
        chunks=chunks,
        parallel=True,
        engine="h5netcdf",
        data_vars="minimal",
        coords="minimal",
        compat="override",
    ).astype(np.float64)
    match = re.search(r"out2d_(\d+)\.nc", str(path_out2d))

    if match:
        num_file = match.group(1)

    # Create a grid object
    grid = sx.open_grid(path_out2d, chunks=chunks)
    # Create a dataset for the variable
    sxds = sx.read_schism_nc(
        grid,
        ds_out2d[
            [
                varname,
            ]
        ],
    )

    # Read polygons from a shapefile
    gdf = gpd.read_file(path_polygons)

    # Convert the coordinate system to what our SCHISM uses.
    crs_target = 26910
    gdf = gdf.to_crs(crs_target)

    # Show what we have read.
    # print(gdf[["comments"]])

    # loop through relevant polygons (shouldn't have a polygon without cells in the hgrid)
    for poly in poly_ids:

        polygon = gdf[gdf.Name == poly].geometry
        # comment = gdf.loc[poly].comments

        # Slice the data using the polygon
        da_subset = sxds[varname].subset.bounding_polygon(polygon)
        # Integrate node values to each element
        da_integrated = da_subset.integrate()

        # If we want to convert the unit from kg/s to cfs.
        kg_to_cfs = 0.035314667
        # cfs per timestep within polygon
        da_out = da_integrated.sum(dim="n_face") * kg_to_cfs
        # add 8 hours to timestamps there's an 8-hr shift issue from xarray converting to UTC.
        # SCHISM has time information attached so xarray shifts by 8 hours
        da_out["time"] = da_out.time + pd.to_timedelta(8, "h")

        # calculate daily average water evaporation cfs per day by averageing da_out over time (use xarray)
        da_out_mean = da_out.resample(time="1D", closed="right").mean()
        # print(da_out_mean)
        out_df = da_out_mean.to_dataframe()
        out_df = out_df.head(1)
        # print(out_df)

        # write da_out
        out_fn = os.path.join(
            os.path.dirname(path_polygons), f"evap_{poly}_{num_file}.csv"
        )
        print(out_fn)
        out_df.to_csv(out_fn, float_format="%.1f")


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    path_out2d = args.path_out2d
    path_polygons = args.path_polygons
    polygon_domain = args.polygon_domain

    # Polygon Name options
    polygon_dict = {
        "all": [
            "legal_delta",
            "suisun_bay",
            "grizzly_restoration",
            "prospect_restoration",
            "lookout_restoration",
            "egbert_restoration",
            "full_domain",
        ],
        "baseline": ["legal_delta", "suisun_bay", "full_domain"],
        "franks": ["legal_delta", "suisun_bay", "full_domain"],
        "suisun": ["legal_delta", "suisun_bay", "grizzly_restoration", "full_domain"],
        "cache": [
            "legal_delta",
            "suisun_bay",
            "prospect_restoration",
            "lookout_restoration",
            "egbert_restoration",
            "full_domain",
        ],
    }

    poly_ids = polygon_dict[polygon_domain]

    extract_evap_daily(path_out2d, path_polygons, poly_ids)


# def main_hardwire():
#     model_out_dir = "/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_7/outputs"
#     os.chdir(model_out_dir)
#     path_out2d = os.path.join(model_out_dir, 'out2d_700.nc')
#     path_polygons = '/scratch/tomkovic/DSP_code/scripts/evap/evap_poly.shp'

#     # Baseline
#     base_polys = ['legal_delta', 'suisun_bay', 'full_domain']

#     extract_evap_daily(path_out2d, path_polygons, base_polys)


if __name__ == "__main__":
    # main_hardwire()
    main()
