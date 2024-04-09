# -*- coding: utf-8 -*-
"""
Hotstart example for a basic schism run with TEMP and SALT as tracers.
"""

import schimpy.schism_hotstart as sh
import numpy as np
import xarray as xr
import uxarray as ux
import suxarray as sx
import suxarray.helper

import os

global date_start

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def make_hotstart_nc(yaml_fn, year, date_start, modules=['TEM','SAL']):

    hotstart_fn = f"hotstart_{year}.nc" # output hotstart file

    # create a hotstart file for SCHISM
    h = sh.hotstart(temp_yml,modules=modules,
                    crs ='EPSG:26910')
    h.create_hotstart()
    hnc = h.nc_dataset
    hnc.to_netcdf(hotstart_fn)   
        

def visit_ncs(yaml_fn, year, modules=['TEM','SAL']):
    h = sh.hotstart(yaml_fn,modules=modules,
                    crs ='EPSG:26910')
    h.read_yaml()
    hotstart_fn = f"hotstart_{year}.nc" # output hotstart file

    #%% converting hotstart file to schism output format so that it can be viewd by VisIt
    sh.hotstart_to_outputnc(hotstart_fn,str(h.date),hgrid_fn='../hgrid_dsp_baseline.gr3', 
                            vgrid_fn='../vgrid.in.3d',vgrid_version=h.vgrid_version,
                            outname=f"schout_hotstart_{year}.nc")
    

def _depth_average(v, zs, k, dry):
    # Select the values with the last index from the bottom index
    # array, k
    z_bottom = np.take_along_axis(zs, k.reshape(1, -1, 1), axis=-1)
    depth = zs[:, :, -1] - z_bottom.squeeze(axis=-1)
    # Mask nan values
    v = np.ma.masked_invalid(v, copy=False)
    zs = np.ma.masked_invalid(zs, copy=False)
    result = np.trapz(v, x=zs, axis=-1) / depth
    result[dry == 1.] = np.nan
    return result


def hotst_to_xy(hotstart_fn, out2d_fn, out_name):

    hnc = xr.open_dataset(hotstart_fn)
    onc = xr.open_dataset(out2d_fn)

    sal_nd = hnc.tr_nd[:,:,1]
    z_coords = hnc.z

    dryFlagNode = hnc.idry
    bottom_index_node = onc.bottom_index_node

    del hnc, onc
    
    da_da = xr.apply_ufunc(_depth_average,
                            sal_nd,
                            z_coords,
                            bottom_index_node - 1,
                            dryFlagNode,
                            input_core_dims=[["nVert",],
                                            ["nVert",],
                                            [],
                                            []],
                            dask='parallelized',
                            output_dtypes=[float])

    with open(out_name, 'w') as outcsv:
        print('hi')


if __name__ == '__main__':

    yaml_fn = "./hotstart_fmt.yaml"
    date_starts = {'2008': '2006-11-14',
                '2010': '2009-10-29',
                '2012': '2011-10-18',
                '2014': '2013-10-24'
                }
            
    temp_yml = 'temp_yml.yml'     

    for year, date_start in date_starts.items():

        with open(yaml_fn, 'r') as file:
            yaml_in = file.read()
        with open(temp_yml, 'w') as file:
            file.write(yaml_in.format(**locals()))

        # # make_hotstart_nc(temp_yml, year, date_start)
        # os.remove('hgrid.nc')
        # visit_ncs(temp_yml, year, date_start)

        # # delete temp_yml
        # os.remove(temp_yml)

        # create xy spatial depth-averaged output
        hotstart_fn = f"hotstart_{year}.nc" # output hotstart file
        out2d_fn = "../outputs/out2d_1.nc"
        out_name = hotstart_fn.replace(".nc",".csv")
        out_name = out_name.replace("schout_","")
        hotst_to_xy(hotstart_fn, out2d_fn, out_name)

