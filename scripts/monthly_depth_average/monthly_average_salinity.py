# Import packages
from pathlib import Path
import numpy as np
import xarray as xr
import suxarray as sx

# Where the test data is located
dir_data = Path('../../model/schism/azure_dsp_2024_lhc_v3/simulations/suisun_lhc_4')

files_out2d = [str(p) for p in [dir_data / 'out2d_1.nc', dir_data / 'out2d_2.nc']]
files_zcoords = [str(p) for p in [dir_data / 'zCoordinates_1.nc', dir_data / 'zCoordinates_2.nc']]
grid = sx.open_grid(files_out2d, files_zcoords)

grid.plot().opts(width=600, height=400)