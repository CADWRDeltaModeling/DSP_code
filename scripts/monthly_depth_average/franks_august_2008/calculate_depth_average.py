#!/usr/bin/env python
# Use suxarray v2024.09.0
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import suxarray as sx
from suxarray.utils.computing import _depth_average
from shapely.geometry import Point
from dask.distributed import Client
import click


@click.command()
@click.option("--path_study", type=Path)
def main(path_study: Path):

    path_output = path_study / "outputs_protected"

    client = Client(n_workers=20, threads_per_worker=1)

    start_day = 626 # 687 is October 1, 2008; 626 is August 1, 2008
    end_day = 656 # 717 is October 31, 2008; 656 is August 31, 2008

    salinity_files = [f for f in path_output.glob("salinity_*.nc") if start_day <= int(f.stem.split('_')[1]) <= end_day]
    ds_salinity = xr.open_mfdataset(salinity_files, engine="h5netcdf", parallel=True).astype(np.float64)
    
    files_out2d = [f for f in path_output.glob("out2d_*.nc") if start_day <= int(f.stem.split('_')[1]) <= end_day]
    files_zcoords = [f for f in path_output.glob("zCoordinates_*.nc") if start_day <= int(f.stem.split('_')[1]) <= end_day]
    grid = sx.open_grid(files_out2d, files_zcoords)
    print(grid)

    sxds = sx.read_schism_nc(grid, ds_salinity)
    sxda_salinity_depth_average = sxds.salinity.depth_average()
    
    print(sxda_salinity_depth_average)

    sxda_salinity_depth_average_monthly = sxda_salinity_depth_average.resample(
        time=f"{end_day-start_day+1}D", origin="start").mean(dim="time") # 31 days in october (COUlD BE GENERIC)
    sxda_salinity_depth_average_monthly.to_netcdf("salinity_depth_averaged_monthly_franks.nc")
    
    sx.write_schism_nc(sxda_salinity_depth_average_monthly)

    client.close()


if __name__ == "__main__":
    main()
