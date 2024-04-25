# -*- coding: utf-8 -*-
"""
Nudging example for a multiple scenario schism run with TEMP and SALT as tracers.
"""

import os
from schimpy import nudging
from schimpy import schism_yaml
from schimpy.schism_mesh import write_mesh, read_mesh
import xarray as xr
import datetime

from schimpy.geo_tools import ll2utm
from shapely.ops import cascaded_union
from shapely.geometry import Polygon

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os

def make_nudging_nc(yaml_fn, date_range, interp_source, model_dir):
    rnday = (date_range[1] - date_range[0]).days
    date_start = date_range[0]

    with open(yaml_fn, 'r') as file:
        yaml_in = file.read()
    with open(temp_yml, 'w') as file:
        file.write(yaml_in.format_map(locals())) # uses rnday, date_start, interp_source, model_dir

    # create a nudging file for SCHISM
    ndge = nudging.Nudging(temp_yml,proj4 ='EPSG:26910')
    ndge.read_yaml()
    ndge.create_nudging()
    
    # delete temp_yml
    # os.remove(temp_yml)


if __name__ == '__main__':

    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    lhc_setup_yml = '../../../../../scripts/boundary_generation/input/lathypcub_v3_setup.yaml'
    with open(lhc_setup_yml, 'r') as f:
        inputs = schism_yaml.load(f)

    date_range_dict = {}
    interp_dict = {}

    for cp in inputs['cases']:
        year = cp['model_year']
        if not year in date_range_dict.keys():
            # set up interpolant based on data availability
            if cp['case_start'] > datetime.date(2014,1,1):
                interp_dict[year] = ['cencoos_ca_roms/processed/cencoos_hourly_pst'] # available 2014-01-01 to 2022-09-14
                date_range_dict[year] = [cp['case_start'],cp['case_end']] # set the model run period based on lhc setup file
            elif cp['case_end'] <= datetime.date(2014,5,1):
                interp_dict[year] = ['sfb_forecast/processed/ca_subSFB_fcst_hourly_pst'] # available 2005-01-01 to 2014-04-29
                date_range_dict[year] = [cp['case_start'],cp['case_end']] # set the model run period based on lhc setup file
            elif (cp['case_start'] <= datetime.date(2014,1,1)) & (cp['case_end'] > datetime.date(2014,5,1)):
                interp_dict[year] = ['sfb_forecast/processed/ca_subSFB_fcst_hourly_pst',
                                     'cencoos_ca_roms/processed/cencoos_hourly_pst']
                date_range_dict[year] = [[cp['case_start'], datetime.date(2014,1,1)], 
                                     [datetime.date(2014,1,2), cp['case_end']]] # set the model run period based on lhc setup file

    model_dir = os.path.join(os.getcwd(),'../')
    yaml_fn = os.path.join(os.getcwd(),"./nudging_fmt.yaml")
            
    temp_yml = './temp_yml.yml'     

    for year, date_ranges in date_range_dict.items():

        if not os.path.exists(str(year)):
            os.mkdir(str(year))
        os.chdir(str(year))

        if len(interp_dict[year]) == 1:
            date_range = date_ranges
            interp_source = interp_dict[year][0]
            make_nudging_nc(yaml_fn, date_range, interp_source, model_dir)
            print(f'writing nudging files for {year}')

        else:
            # the simulation period spans two interpolant data sources, need to run two different nudging calls
            print(f'Need to create combined netCDF for {year}')
            interp_dirs = []
            for i, interp_source in enumerate(interp_dict[year]):
                interp_dir = interp_source.split('/')[0]
                interp_dirs.append(interp_dir)
                if not os.path.exists(interp_dir):
                    os.mkdir(interp_dir)
                os.chdir(interp_dir)
                date_range = date_ranges[i]
                make_nudging_nc(yaml_fn, date_range, interp_source, model_dir)
                os.chdir('../')

            # TODO: stitch together two nudging netCDF files
            for var in ['TEM','SAL']:
                nc1 = xr.open_dataset(os.path.join(interp_dirs[0],f'{var}_nu_roms.nc'))
                nc2 = xr.open_dataset(os.path.join(interp_dirs[1],f'{var}_nu_roms.nc'))
                nc2.assign(time=nc2.time + ((date_ranges[1][0]-date_ranges[0][0]).days * 86400)) # shift the second time dataset by the difference between datasteps
                ncout = xr.concat([nc1,nc2], dim='time')
                ncout.to_netcdf(f'{var}_nu_roms.nc')
                print(f"Wrote combined netCDF: {var}_nu_roms.nc")
        
        os.chdir('../')
