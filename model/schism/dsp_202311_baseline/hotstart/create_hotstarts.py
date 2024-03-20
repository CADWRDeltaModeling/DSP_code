# -*- coding: utf-8 -*-
"""
Hotstart example for a basic schism run with TEMP and SALT as tracers.
"""

import schimpy.schism_hotstart as sh
import numpy as np
import xarray as xr
# import uxarray as ux
# import suxarray as sx
# import suxarray.helper
# import matplotlib.pyplot as plt
# from dms_datastore.read_ts import *
# from dms_datastore.dstore_config import *
# from dms_datastore.read_multi import *
# from schimpy.unit_conversions import *
# from vtools.data.vtime import days
# import glob
import pandas as pd
import os

global date_start

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def make_hotstart_nc(yaml_fn, year, date_start, modules=['TEM','SAL']):

    hotstart_fn = f"hotstart_{year}.nc" # output hotstart file

    # create a hotstart file for SCHISM
    h = sh.hotstart(yaml_fn,modules=modules,
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
    sh.hotstart_to_outputnc(hotstart_fn,str(h.date),hgrid_fn='../hgrid.gr3', 
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

def create_csv_data(t0):

    # t0 = pd.Timestamp(2006,11,14)
    # nudgelen =days(30)


    stations = ['anh','benbr','hsl','bts','snc','ibs','cyg','hun','bdl',
                'fmb','msl','cll','gzl','ryc','hon',
                'c24','pct','flt','mrz','pct','mal','pts','carqb','benbr',
                'co5','ssi','emm','sdi','blp','jer','sjj',
                'dsj','frp','fal','bet','hol2','hll','orq2','frk','holm','bac','mdm',
                'dbi','ori','oh4','mab','pri','ppt','trn','rindg','sjc','rri','sjg','bdt',
                'dvi','sjr','orx','pdup','tpp','uni','sga','gle','old','twa','orm','oad',
                'trp','glc2','wci','vcu','mab','mtb','rri2','mdmzq','sdc','ges','swe','gss',
                'nmr','sus','sss','sut','snod','gln','rye','ryf','rvb','mir','dws','lib','ucs','has',
                'srh','awb','afo','hst','ist','ssw','von','few','fre','wlk','gys','god','sal']

    # Stations where an "upper" and "lower" sublocation occur and we must distinguish the upper
    add_upper = ["anh","cll","mrz","emm","mal","pts"] 

    ##  


    station_df=station_dbase()

    buf = days(5)
    sdata = t0 - buf
    edata = t0 + nudgelen + buf

    station_df = station_df.loc[stations]

    repo = "//cnrastore-bdo/Modeling_Data/continuous_station_repo_beta/formatted_1yr"
    no_such_file = []
    tndx = pd.date_range(t0,t0+nudgelen,freq='H')
    all_vars = ["temperature","salinity"]
    used_stations = set()
    nudging_dfs = {}
    accepted_loc= [] 
    for label_var in all_vars:
        var = {"temperature":"temp","salinity":"ec"}[label_var] # working variable for data
        print(f"Working on variable: {label_var},{var}")
        vals = []
        accepted = {}

        for ndx,row in station_df.iterrows():
            x = row.x
            y = row.y
            fndx = ndx+"@upper" if ndx in add_upper else ndx
            
            pat = f"*_{fndx}_*_{var}*_20??.csv"
            pat = os.path.join(repo,pat)
            matches = glob.glob(pat)
            if len(matches) == 0:
                no_such_file.append((ndx,var))
                continue
            try:
                subloc='upper' if ndx in add_upper else None
                ts = read_ts_repo(ndx,var,subloc = subloc,
                    repo=repo,
                    src_priority='infer')
                ts = ts.loc[sdata:edata]
                ts = ts.interpolate(limit=4)
                if ts.shape[1] >1:
                    ts=ts.mean(axis=1)
                    ts.name="value"
                else:
                    ts = ts.squeeze()
                if var == "temp":
                    topquant = ts.quantile(q=0.25)
                    if (topquant > 35):            
                        print("Transforming F to C based on 25% qyantuke > 35deg")
                        print("Transforming F to C based on 25% qyantuke > 35deg")
                        ts = fahrenheit_to_celsius(ts)
                    if ndx in ['clc'] and (ts<0.).all():
                        ts = celsius_to_farenheit(ts)
                elif var == "ec":
                    ts = ec_psu_25c(ts)
                else:
                    raise ValueError(f"Haven't worked out transforms needed except for {var}, only salt/temp")
                
            
                val = ts.at[t0]
                if not np.isnan(val): 
                    vals.append((ndx,x,y,val))
                    
                # This is the fraction of missing data
                ts = ts.reindex(tndx)
                gap_frac = ts.isnull().sum()/len(ts)
                print(f"Fraction of mssing data for {ndx} {var} is {gap_frac}")                
                if gap_frac < 0.25:
                    print(f"Accepted {ndx} {var}")
                    ts.columns=[ndx]
                    ts = ts.fillna(-9999.)
                    accepted[ndx]=ts
                    if ndx not in used_stations:
                        accepted_loc.append((ndx,x,y))
                        used_stations.add(ndx)
                
                        
            except Exception as err:
                print("Exception")
                print(str(err))
                print(ndx,var)
                print(ts.iloc[0:5])
                print(err)
        var_df = pd.DataFrame(data = vals,columns=("station","x","y",f"{label_var}"))
        var_df.set_index("station")
        var_df.to_csv(f"hotstart_data_{label_var}.csv",sep=",",float_format="%.2f")
    
        nudging_df = pd.concat(accepted,axis=1)
        nudging_df.index.name='datetime'
        nudging_dfs[label_var] = nudging_df
        print(nudging_df)

    obs_xy = pd.DataFrame(data=accepted_loc,columns=["site","x","y"])
    print("reindexing and printing")

    for label_var in all_vars:
        nudging_dfs[label_var].to_csv(f"nudging_data_{label_var}.csv",sep=",",float_format="%.2f")
        # Deprecated
        #nudging_dfs[label_var].reindex(columns=obs_xy.site).to_csv(f"nudging_data_{label_var}_b.csv",sep=",",float_format="%.2f")


    obs_xy = obs_xy.set_index("site",drop=True)
    obs_xy.to_csv(f"obs_xy.csv",sep=",",float_format="%.2f")

    print("No such file")
    for item in no_such_file:
        print(item)
        



if __name__ == '__main__':

    yaml_fn = "./hotstart_fmt.yaml"
    date_starts = {'2008': '2006-11-14',
                '2010': '2009-10-29',
                '2012': '2011-10-18',
                '2014': '2013-10-24'
                }
            
    temp_yml = 'temp_yml.yml'     

    # TODO: need to modify the download of hotstart_data csv for each year

    for year, date_start in date_starts.items():

        with open(yaml_fn, 'r') as file:
            yaml_in = file.read()
        with open(temp_yml, 'w') as file:
            file.write(yaml_in.format(**locals()))

        make_hotstart_nc(temp_yml, year, date_start)
        # os.remove('hgrid.nc')
        # visit_ncs(temp_yml, year, date_start)

        # delete temp_yml
        os.remove(temp_yml)

        # # create xy spatial depth-averaged output
        # hotstart_fn = f"hotstart_{year}.nc" # output hotstart file
        # out2d_fn = "../outputs/out2d_1.nc"
        # out_name = hotstart_fn.replace(".nc",".csv")
        # out_name = out_name.replace("schout_","")
        # hotst_to_xy(hotstart_fn, out2d_fn, out_name)

