####!/usr/bin/env python
# -*- coding: utf-8 -*-import pandas as pd
import matplotlib.pyplot as plt
from dms_datastore.read_ts import *
from dms_datastore.dstore_config import *
from dms_datastore.read_multi import *
from schimpy.unit_conversions import *
from vtools.data.vtime import days
import glob
import pandas as pd
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

## Items you may want to change

t0 = pd.Timestamp(2006,11,1)
tend = pd.Timestamp(2015,1,1)
nudgelen = days((tend-t0).days)

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

station_df=station_dbase()

buf = days(5)
sdata = t0 - buf
edata = t0 + nudgelen + buf

station_df = station_df.loc[stations]

repo = "/nasbdo/modeling_data/repo/continuous/formatted"
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
            print(f"Fraction of missing data for {ndx} {var} is {gap_frac}")                
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
  
    nudging_df = pd.concat(accepted,axis=1)
    nudging_df.index.name='datetime'
    nudging_dfs[label_var] = nudging_df
    print(nudging_df)

obs_xy = pd.DataFrame(data=accepted_loc,columns=["site","x","y"])
print("reindexing and printing")

for label_var in all_vars:
    nudging_dfs[label_var].to_csv(f"nudging_data_{label_var}.csv",sep=",",float_format="%.2f")

print("No such file")
for item in no_such_file:
    print(item)
    
