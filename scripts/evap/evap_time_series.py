# -*- coding: utf-8 -*-
"""

This script reads evap_polygons.yaml 
It queries the evaporation rates in all the cells contained in each polygon in the yaml file to generate a single 


usage: --start 2004-04-18 --evap_polys evap_polygons.yaml
output: evaporation rate on surface of all cells in polygons time series in csv format,
        evap_out.csv

"""


from osgeo import ogr
from shapely.geometry import shape, Point
import json
import xarray as xr
from math import isnan
import numpy as np
import os
import re
import shutil
from dms_datastore.logging_config import logger
import argparse
import netCDF4

# read in shapefile 
def records(file):
    # generator
    reader = ogr.Open(file)
    layer = reader.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        yield json.loads(feature.ExportToJson())

def process_out2d_var(outputs_dir, time_range, partition_shp, model_start, output_file, variable):
    """Process  evaporation rates into a time series
    """
    # time period
    # TODO:generate the time period and specify which out2d files will be looped through as well as creating the dataframe with the timestamp values on them.


    # load data from first out2d_file
    out2d_src = netCDF4.Dataset(out2d_files[0])
    node_x = out2d_src.variables['SCHISM_hgrid_node_x'][:]
    node_y = out2d_src.variables['SCHISM_hgrid_node_y'][:]
    n_hgrid = len(node_x)
    
    # load shapefile and store regions
    regions = []
    reg_names = []
    print('Reading shapefile and out2d file...')
    try:
        mpoly = records(partition_shp)
        for mp in mpoly:
            regions.append(mp)
            reg_names.append(mp['properties']['name'])
    except BaseException as err:
        print(f"Error reading partition shape file {err}, {type(err)}")

    # Create empty index dictionary
    nodesinROI = {key:[] for key in reg_names}

    # Find which nodes are in which regions of interest (polygons)
    print('Determining nodes in ROI...')
    for ii in range(n_hgrid):
        point = Point(node_x[ii], node_y[ii])
        for mp in regions:
            reg_name = mp['properties']['name']
            if point.within(shape(mp['geometry'])):
                nodesinROI[reg_name].append(ii)

    for out2d_file in out2d_files:
        out2d_src = netCDF4.Dataset(out2d_file)

        # get out2d variables
        var_out = out2d_src.variables[variable][:,:]
        time_dim  = out2d_src.dimensions["time"]
        
        

    print('hi')

if __name__ == "__main__":

    simulation = 'baseline_lhc_5'
    model_dir = f'/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/{simulation}/outputs'
    os.chdir(model_dir)

    evap_polygons = '/scratch/tomkovic/DSP_code/scripts/post-processing/input/evaporation_polygons.shp'

    output_file = f'/scratch/tomkovic/DSP_code/scripts/post-processing/data_out/evaporation/{simulation}_evap.csv'

    out2d_file = os.path.join(model_dir,'outputs/out2d_2.nc')

    model_start = 'none'


    # File names
    out2d_files = [os.path.join(model_dir, f"out2d_{fnum}.nc") for fnum in range(2,5)]

    variable = 'evaporationRate'
    
    process_out2d_var(out2d_files, evap_polygons, model_start, output_file, variable)


