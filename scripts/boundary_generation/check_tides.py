
import pandas as pd
import numpy as np
import os, sys

from pydelmod.create_ann_inputs import get_dss_data
import pyhecdss
from dms_datastore.read_ts import read_noaa
from schimpy.station import read_staout
from vtools.functions.filter import cosine_lanczos

import matplotlib.pyplot as plt

os.chdir(os.path.dirname(os.path.abspath(__file__)))

def tidal_energy(df):

    df = df * 3.28084  # ft/m
    df = cosine_lanczos((df-cosine_lanczos(df.copy(),cutoff_period='40H', padtype='odd'))**2, cutoff_period='40H', padtype='odd')  # = < (z- <z>)^2 >
    if not isinstance(df.index, pd.DatetimeIndex):
        df.index = df.index.to_timestamp()
    df = df.resample('D', closed='right').mean()
    df.columns = ['tidal_energy']
    df = df.resample('15min').ffill()

    return df

plot_dates = pd.to_datetime(['2008-01-01','2008-06-01'])

rma_sfo_bdry_filename = '../../../../data/_fromRMA/GGtide_shift100d/GGtide_shift100d.dss'
b_part_dss_filename_dict = {'SAN FRANCISCO': rma_sfo_bdry_filename}
b_part_c_part_dict = {'SAN FRANCISCO': 'STAGE-NAVD88'}
b_part_f_part_dict = {'SAN FRANCISCO': 'NOAA SHIFTED + 100D-MERGED'}
rma_sfo_100d_shift = get_dss_data(b_part_dss_filename_dict, 'b_part', primary_part_c_part_dict=b_part_c_part_dict, 
                            primary_part_f_part_dict=b_part_f_part_dict, daily_avg=False)
rma_sfo_100d_shift = rma_sfo_100d_shift.loc[plot_dates[0]:plot_dates[1]]
rma_sfo_100d_shift.index = rma_sfo_100d_shift.index.to_timestamp()

# sfo NOAA station
# download_noaa --syear 1989 --eyear 2020 --param water_level sfo_station.txt # in input

sffpx_elev = './input/noaa_download/noaa_sffpx_9414290_water_level_1989_2021.csv'

sffpx = read_noaa(sffpx_elev, force_regular=True)
sffpx_shift = sffpx.copy()
sffpx_shift.index = sffpx_shift.index + pd.Timedelta(days=100)
sffpx = sffpx.loc[plot_dates[0]:plot_dates[1]]
sffpx_shift = sffpx_shift.loc[plot_dates[0]:plot_dates[1]]

# SCHISM sffpx outputs
schism_sfo4 = read_staout('./input/staout_1.base_4', './input/station.in', pd.to_datetime('2006-11-14'))
schism_sfo4 = schism_sfo4['sffpx_default']
schism_sfo5 = read_staout('./input/staout_1.base_5', './input/station.in', pd.to_datetime('2006-11-14'))
schism_sfo5 = schism_sfo5['sffpx_default']

# # Plot Results
# fig, ax = plt.subplots(1)
# ax.plot(sffpx.index,sffpx.values)
# ax.plot(sffpx_shift.index,sffpx_shift.values)
# # ax.plot(rma_sfo_100d_shift.index,rma_sfo_100d_shift.values)
# ax.plot(schism_sfo4.index,schism_sfo4.values)
# ax.plot(schism_sfo5.index,schism_sfo5.values)
# ax.legend(['noaa sffpx','python shift','dss shift','schism','schism shifted'])    
# plt.show()

# plot energies
schism_sfo5_nrg = tidal_energy(schism_sfo5)
rma_sfo_100d_shift_nrg = tidal_energy(rma_sfo_100d_shift)
schism_sfo5_calc_nrg = pd.read_csv('../../../casanntra/data/schism_base_5.csv', index_col=0, parse_dates=['datetime'])
schism_sfo5_calc_nrg = schism_sfo5_calc_nrg['sf_tidal_energy']
rma_mod_sfo5_calc_nrg = pd.read_csv('../../../casanntra/data/rma_base_5.csv', index_col=0, parse_dates=['datetime'])
rma_mod_sfo5_calc_nrg = rma_mod_sfo5_calc_nrg['sf_tidal_energy']

# Plot Results
fig, ax = plt.subplots(1)
ax.plot(schism_sfo5_nrg.index,schism_sfo5_nrg.values)
ax.plot(rma_sfo_100d_shift_nrg.index,rma_sfo_100d_shift_nrg.values)
ax.plot(schism_sfo5_calc_nrg.index,schism_sfo5_calc_nrg.values)
ax.plot(rma_mod_sfo5_calc_nrg.index,rma_mod_sfo5_calc_nrg.values)
ax.legend(['schism filtered','rma filtered','schism mod out','rma mod out'])    
plt.show()