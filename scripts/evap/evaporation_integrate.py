from pathlib import Path
import re
import numpy as np
import xarray as xr
import geopandas as gpd
import suxarray as sx
# import click
import pandas as pd


def natural_sort_key(s, _nsre=re.compile(r"(\d+)")):
    """Sort function to sort file names naturally"""
    return [
        int(text) if text.isdigit() else text.lower() for text in _nsre.split(str(s))
    ]


# @click.command()
# @click.option("--path_study", type=Path)
# @click.option("--path_input", type=Path)
def main():
    # path_study = Path(
    #    "/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_5/"
    # )

    # The paths to the study and the input are hard-coded to /data, assuming
    # that those are mounted via Docker.
    path_study = Path("/data")
    # path_inputs = Path("/home/knam/support/dsp-60")
    path_input = Path("/input")
    path_out2d = sorted(
        (path_study / "outputs").glob("out2d_?.nc"), key=natural_sort_key
    )
    # We are going to process the following variable.
    varname = "evaporationRate"
    chunks = {"time": 12}

    # Read polygons from a shapefile
    path_polygons = path_input / "dsp_evaporation_polygons.shp"
    gdf = gpd.read_file(path_polygons)

    # Convert the coordinate system to what our SCHISM uses.
    crs_target = 26910
    gdf = gdf.to_crs(crs_target)

    # Show what we have read.
    print(gdf[["comments"]])
    # In this case, we will uses the second polygon from the shapefile.
    # loop through relevant polygons (shouldn't have a polygon without cells in the hgrid)
    polygon = gdf.iloc[1].geometry
    comment = gdf.iloc[1].comments

    # Read the out2d files.
    ds_out2d = xr.open_mfdataset(
        path_out2d,
        chunks=chunks,
        parallel=True,
        engine="h5netcdf",
        data_vars="minimal",
        coords="minimal",
        compat="override",
    ).astype(np.float64)

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

    # Slice the data using the polygon
    da_subset = sxds[varname].subset.bounding_polygon(polygon)
    # Integrate node values to each element
    da_integrated = da_subset.integrate()
    # If we want to convert the unit from kg/s to cfs.
    kg_to_cfs = 0.035314667
    # Write the output
    # cfs per timestep within polygon
    da_out = da_integrated.sum(dim="n_face") * kg_to_cfs
    # da_out = da_integrated.sum(dim="n_face") # in kg/s
    attrs = da_out.attrs
    attrs["units"] = "cfs"
    # attrs["units"] = "kg s-1"
    attrs["long_name"] = "water_evaporation_rate_over_area"
    attrs["area"] = f"{comment}"
    da_out = da_out.assign_attrs(attrs)
    path_out = f"{varname}.nc"
    da_out.to_netcdf(path_out)

    # calculate daily average water evaporation cfs per day by averageing da_out over time (use xarray)
    da_out.index = "SHIFT"  # add 8 hours to timestamps there's an 8-hr shift issue from xarray converting to UTC. SCHISM has time information attached so xarray shifts by 8 hours
    # maybe not use closed right?
    da_out = da_out.resample('D', closed='right').mean()
    1    .01235
    2    .13590913i5
    ...

    return


if __name__ == "__main__":
    main()
