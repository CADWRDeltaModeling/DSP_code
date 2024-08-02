#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to gather input/output data from SCHISM for training ANN


"""

import os
from schimpy.station import *
from schimpy.model_time import file_to_timestamp
from vtools.functions.filter import cosine_lanczos
import datetime
import warnings
import string
import re

def get_single_staout(outputs_fpath, station_inpath, time_basis, loc, depth='upper'):

    all_ts = read_staout(outputs_fpath, station_inpath, time_basis)
    station_df = read_station_in(station_inpath)

    if isinstance(loc, list):
        loc = loc[0]

    if depth == 'upper':
        loc_outs = [stc for stc in station_df.index if loc in stc]
        if len(loc_outs) == 0:
            raise ValueError(f"There is no station output for {loc}")
        elif len([stc for stc in loc_outs if 'upper' in stc]) > 0:
            col = f'{loc}_upper'
        else:
            col = f'{loc}_default'
    else:
        raise ValueError('Currently only programmed to get the upper layer by default')

    out_ts = all_ts[col]

    return out_ts

def get_start_date_from_param(param_in):

    with open(param_in, 'r') as param:
        for line in param.readlines():
            if 'start_year' in line:
                sy =  int(re.findall(r'\b\d+\b', line)[0])
            elif 'start_month' in line:
                sm =  int(re.findall(r'\b\d+\b', line)[0])
            elif 'start_day' in line:
                sd =  int(re.findall(r'\b\d+\b', line)[0])
                
    start_date = datetime.datetime(sy, sm, sd)
    
    return start_date

outputs_fpath = '/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_5/outputs_tropic/staout_1'
station_inpath = "/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_5/station.in" 
param_fpath = "/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_5/param.nml.tropic"
loc = "SFO"

time_basis = get_start_date_from_param(param_fpath) 

out_ts = get_single_staout(outputs_fpath, station_inpath, time_basis, loc)

out_fn = ""

out_ts.to_csv(out_fn, index=True)