##!/usr/bin/env python
# -*- coding: utf-8 -*

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
import pyhecdss
from pydsm.functions import tsmath
# from dms_datastore.read_ts import read_ts
# from dms_datastore.read_multi import read_ts_repo

from pydelmod.create_ann_inputs import get_dss_data
import os


def perturb_binary(orig_ts, P00=0.988, P11=0.9625):
    """ Given a series of  ones and zeros, perturb to  the complementary state. "Unperturbed" means original value.
    Parameters 
    ----------
    ts : Series
        Series to perturb

    P00 : float
        Markov transition probability for unoperated staying unoperated. Probability of transition to perturbed
        is 1.0-P00. Typically > 0.98 for daily values unless you want a lot of chatter

    P11 : float
        Transition probability from perturbed to perturbed. Probability of reverting to unperturbed is 1-P11. Typically 
        a little bit smaller than P00 if you want the original series to dominate.

    Returns
    -------
    Frame with columns "perturbed" and "orig"
    """

    ts = orig_ts.to_frame()
    ts.columns = ["orig"]
    ts["perturbed"] = ts["orig"]

    ts['use_orig'] = True

    USEORIGNDX = 2  # ts.columns.index('use_orig')

    # Transitions
    rnum = np.random.rand(len(ts))
    for i in range(1, len(ts)):
        if ts.iloc[i-1, USEORIGNDX]:
            ts.iloc[i, USEORIGNDX] = True if rnum[i-1] < P00 else False
        else:
            ts.iloc[i, USEORIGNDX] = False if rnum[i-1] < P11 else True

    ts["perturbed"] = ts.orig.where(ts.use_orig, (ts.perturbed+1.) % 2)
    ts.drop("use_orig", axis=1, inplace=True)
    ts.drop("orig", axis=1, inplace=True)
    return ts


def sample_binary_ts():
    ndx = pd.date_range(pd.Timestamp(2007, 1, 1), freq='d', periods=1000)
    ts = pd.Series(index=ndx, data=None)
    ts[:] = 1.0

    off_periods = [("2007-6-01", "2007-09-15"), ("2008-03-01", "2008-05-05"),
                   ("2008-10-01", "2008-11-15"), ("2009-03-01", "2009-05-05")]

    for off in off_periods:
        ts.loc[pd.to_datetime(off[0]):pd.to_datetime(off[1])] = 0

    ts.name = "pos"
    return ts

# example:
# ts_orig = sample_binary_ts()
# P00=0.988
# P11=0.9625
# ts = perturb_binary(ts_orig,P00,P11)  # Defaults are fine
# print(ts)
# ts["perturbed"]*=0.75

# ts.plot() #drawstyle="steps-post")
# plt.legend()
# plt.show()


def perturb_suisun_ops(radial_df, flash_df, boat_df, P00=0.990, P11=0.990, column_order=['radial', 'flashboard', 'boat_lock']):
    """ Given a series of  ones and zeros, perturb to  the complementary state. "Unperturbed" means original value.
    Parameters 
    ----------
    ts : Series
        Series to perturb

    P00 : float
        Markov transition probability for unoperated staying unoperated. Probability of transition to perturbed
        is 1.0-P00. Typically > 0.98 for daily values unless you want a lot of chatter

    P11 : float
        Transition probability from perturbed to perturbed. Probability of reverting to unperturbed is 1-P11. Typically 
        a little bit smaller than P00 if you want the original series to dominate.

    Returns
    -------
    Frame with columns "perturbed" and "orig"
    """
    # three separate time series that get read in and modified similar to perturb_binary but with some extra checks.

    ts = radial_df.index.to_list()

    radial_pert_df = radial_df.copy()
    flash_pert_df = flash_df.copy()
    boat_pert_df = boat_df.copy()

    # check that time series all overlap
    if not (radial_df.index.min() == flash_df.index.min() == boat_df.index.min()):
        print(f"""START TIMES NOT ALIGNED! radial: {radial_df.index.min()} 
                         flashboard: {flash_df.index.min()} 
                         boat lock: {boat_df.index.min()}""")
        ts = [t for t in ts if t >= max(
            radial_df.index.min(), flash_df.index.min(), boat_df.index.min())]
    if not (radial_df.index.max() == flash_df.index.max() == boat_df.index.max()):
        print(f"""END TIMES NOT ALIGNED! radial: {radial_df.index.max()} 
                       flashboard: {flash_df.index.max()} 
                       boat lock: {boat_df.index.max()}""")
        ts = [t for t in ts if t <= min(
            radial_df.index.max(), flash_df.index.max(), boat_df.index.max())]
    print(f"Time series to be perturbed: {ts[0]} - {ts[-1]}")

    all_df = pd.DataFrame(index=ts, columns=column_order)
    all_df['radial'] = radial_pert_df['op_up']
    all_df['flashboard'] = flash_pert_df['op_up']
    all_df['boat_lock'] = boat_pert_df['op_up']

    # Transitions
    rnum = np.random.rand(len(ts))
    for i in range(1, len(ts)):
        # print(radial_pert_df.loc[radial_pert_df.index[i-1],'op_up'])
        if radial_pert_df.loc[radial_pert_df.index[i-1], 'op_up'] == 0:
            # radial gate is operational
            if rnum[i-1] >= P00:
                radial_pert_df.loc[radial_pert_df.index[i], ['op_up', 'op_down']] = [
                    1.0, 1.0]  # change to open if above P00 threshold
                all_df.loc[radial_pert_df.index[i], 'radial'] = 1.0
            else:
                radial_pert_df.loc[radial_pert_df.index[i], ['op_up', 'op_down']] = [
                    0.0, 1.0]  # keep operational if below P00 threshold
                all_df.loc[radial_pert_df.index[i], 'radial'] = 0.0
        else:
            # radial gate is open (op_up!=0)
            if rnum[i-1] >= P11:
                radial_pert_df.loc[radial_pert_df.index[i], ['op_up', 'op_down']] = [
                    0.0, 1.0]  # change to operational if above P00 threshold
                all_df.loc[radial_pert_df.index[i], 'radial'] = 0.0
            else:
                radial_pert_df.loc[radial_pert_df.index[i], ['op_up', 'op_down']] = [
                    1.0, 1.0]  # keep operational if below P00 threshold
                all_df.loc[radial_pert_df.index[i], 'radial'] = 1.0
        # Check that boat and flash are operated appropriately
        if radial_pert_df.loc[radial_pert_df.index[i-1], 'op_up'] == 0:
            flash_pert_df.loc[radial_pert_df.index[i], ['op_up', 'op_down']] = [
                0.0, 1.0]  # flashboards need to be IN
            all_df.loc[radial_pert_df.index[i], 'flashboard'] = 0.0
            boat_pert_df.loc[radial_pert_df.index[i], ['op_up', 'op_down']] = [
                1.0, 1.0]  # boatlock needs to be OPEN
            all_df.loc[radial_pert_df.index[i], 'boat_lock'] = 1.0

    # remove any rows that are unecessary
    radial_pert_clean = clean_df(radial_pert_df)
    flash_pert_clean = clean_df(flash_pert_df)
    boat_pert_clean = clean_df(boat_pert_df)
    all_clean = clean_df(all_df)

    return radial_pert_clean, flash_pert_clean, boat_pert_clean, all_clean


def clean_df(in_df):
    keep_rows = [0]
    for r in range(1, len(in_df.index)):
        if not (list(in_df.iloc[r, :]) == list(in_df.iloc[r-1, :])):
            # keep the row where something changes and the row before (for plotting)
            if not r-1 in keep_rows:
                keep_rows.append(r-1)
            keep_rows.append(r)
    out_df = in_df.iloc[keep_rows, :]

    return out_df


def read_th(infile):

    in_df = pd.read_table(infile,
                          delim_whitespace=True,
                          index_col='datetime',
                          comment="#")
    in_df.index = pd.to_datetime(in_df.index, format="%Y-%m-%dT%H:%M")

    # monotonic increase check
    # True for monotonically-increasing data
    mon_inc = all(x < y for x, y in zip(in_df.index, in_df.index[1:]))
    if not mon_inc:
        # prints the row(s) where monotonicity is broken
        print(in_df.loc[in_df.index.to_series().diff()
              < pd.to_timedelta('0 seconds')])

    in_df = in_df.reindex(
        pd.date_range(start=in_df.index.min(),
                      end=in_df.index.max(),
                      freq='D'),
        method='ffill')

    print(
        f"TH FILE ==== {infile} min: {in_df.index.min()} max:{in_df.index.max()}")
    return in_df


def write_th(in_df, outfile):
    in_df.index = in_df.index.strftime("%Y-%m-%dT%H:%M")

    in_df.to_csv(outfile,
                 sep=' ',
                 header=True,
                 index=True,
                 index_label='datetime')


def write_DSM2_to_SCHISM(dsm2_dss, radial_df, flash_df, boat_df, out_dir, case, ts_range):
    # primary_pathname_part_dss_filename_dict = {
    #     'MTZSL': dsm2_dss}
    c_parts = {'RADIAL_OP': radial_df,
               'FLASHBOARD_OP': flash_df,
               'BOATLOCK_OP': boat_df}
    c_part = 'RADIAL_OP'

    ts = radial_df.index.to_list()

    # check that time series all overlap
    if not (radial_df.index.min() == flash_df.index.min() == boat_df.index.min()):
        print(f"""START TIMES NOT ALIGNED! radial: {radial_df.index.min()} )
                         flashboard: {flash_df.index.min()} 
                         boat lock: {boat_df.index.min()}""")
        ts = [t for t in ts if t >= max(
            radial_df.index.min(), flash_df.index.min(), boat_df.index.min())]
    if not (radial_df.index.max() == flash_df.index.max() == boat_df.index.max()):
        print(f"""END TIMES NOT ALIGNED! radial: {radial_df.index.max()} 
                       flashboard: {flash_df.index.max()} 
                       boat lock: {boat_df.index.max()}""")
        ts = [t for t in ts if t <= min(
            radial_df.index.max(), flash_df.index.max(), boat_df.index.max())]
    print(f"Time series to be perturbed: {ts[0]} - {ts[-1]}")

    with pyhecdss.DSSFile(dsm2_dss) as d:
        catdf = d.read_catalog()
        dss_df, units, ptype = d.read_its(
            f'/HIST+GATE/MTZSL/{c_part}/{(ts_range[0]- pd.DateOffset(months=6)).strftime("%d%b%Y").upper()} - {ts_range[1].strftime("%d%b%Y").upper()}/IR-DECADE/DSP_{case.upper()}/')
        # dss_df = tsmath.per_aver(dss_df, '1D')
        dss_df.columns = ['MTZSL']
        dss_df.loc[dss_df['MTZSL'] == 1, 'MTZSL'] = 1
        dss_df.loc[dss_df['MTZSL'] == -10, 'MTZSL'] = 0
        dss_df.index = dss_df.index.normalize()
        # dss_df.index = dss_df.index.to_timestamp()

    # add in data from DSM2 DSS
    closest_start = radial_df.index[radial_df.index.searchsorted(
        ts_range[0], side='right') - 1]
    radial_pert_df = radial_df.loc[closest_start:ts_range[1]]
    start_op = radial_pert_df.loc[closest_start, 'op_up']
    radial_pert_df = tsmath.per_aver(radial_pert_df, '1D')
    radial_pert_df.index = radial_pert_df.index.to_timestamp()
    radial_pert_df.op_up = np.nan
    radial_pert_df.loc[radial_pert_df.index.isin(
        dss_df.index), 'op_up'] = dss_df['MTZSL']

    radial_pert_df.loc[closest_start, 'op_up'] = start_op
    radial_pert_df = radial_pert_df.dropna()
    radial_pert_df = radial_pert_df.loc[(
        radial_pert_df != radial_pert_df.shift()).any(axis=1)]  # gets rid of duplicate rows
    radial_pert_df.loc[:, 'install'] = radial_pert_df['install'].astype(float).map(
        '{:.0f}'.format)
    radial_pert_df.loc[:, 'ndup'] = radial_pert_df['ndup'].astype(
        float).map('{:.0f}'.format)

    # fix op_up/op_down (when op_up=0, op_down=1, else op_down=1)
    radial_pert_df['op_down'] = 1

    write_th(radial_pert_df, os.path.join(
        out_dir, f"montezuma_radial_{case}.th"))

    # set boat_lock and flashboard dfs:
    flash_pert_df = flash_df.head(len(radial_pert_df))
    flash_pert_df.loc[:, 'op_up'] = np.nan
    flash_pert_df.loc[:, 'op_up'] = radial_pert_df['op_up'].values
    flash_pert_df.loc[:, 'op_down'] = flash_pert_df['op_up']
    flash_pert_df.loc[:, 'install'] = flash_pert_df['install'].astype(float).map(
        '{:.0f}'.format)
    flash_pert_df.loc[:, 'ndup'] = flash_pert_df['ndup'].astype(
        float).map('{:.0f}'.format)
    flash_pert_df.index = pd.to_datetime(radial_pert_df.index)

    write_th(flash_pert_df, os.path.join(
        out_dir, f"montezuma_flash_{case}.th"))

    boat_pert_df = boat_df.head(len(radial_pert_df))
    boat_pert_df.loc[:, 'op_up'] = np.nan
    boat_pert_df.loc[:, 'op_up'] = -radial_pert_df['op_up'].values
    # switch so boatlock is always opposite the radial gate
    boat_pert_df.loc[boat_pert_df['op_up'] == 0, 'op_up'] = 1
    boat_pert_df.loc[boat_pert_df['op_up'] == -1, 'op_up'] = 0
    boat_pert_df.loc[:, 'op_down'] = boat_pert_df['op_up']
    boat_pert_df.index = pd.to_datetime(radial_pert_df.index)
    boat_pert_df.loc[:, 'install'] = boat_pert_df['install'].astype(float).map(
        '{:.0f}'.format)
    boat_pert_df.loc[:, 'ndup'] = boat_pert_df['ndup'].astype(
        float).map('{:.0f}'.format)

    write_th(boat_pert_df, os.path.join(
        out_dir, f"montezuma_boat_lock_{case}.th"))


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # DCC Gates ---------------------------------------------------------------------------------------

    if False:  # switches so I can just re-write the necessary component
        bddsm_dir = "D:/dsm2/DSM2v821/timeseries"
        primary_pathname_part_dss_filename_dict = {
            'RSAC128': os.path.join(bddsm_dir, "gates-v8-201712.dss")}
        primary_part_c_part_dict = {'RSAC128': 'POS'}
        unit_part_dict = {'RSAC128': 'UNSPECIF'}
        primary_pathname_part = 'b_part'

        P00 = 0.988
        P11 = 0.9625

        df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part,
                                primary_part_c_part_dict=primary_part_c_part_dict,
                                primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)

        output_dir = r"D:/projects/delta_salinity/scripts/DSP_code/scripts/boundary_generation/data_out"

        print("DCC")

        # Read in the gate operation data in th format
        in_ts = df_input['RSAC128']

        # need to fill NaNs and convert from 0-2 float to 0 or 1 binary
        in_ts = in_ts.ffill()

        in_ts = (in_ts > 0.1).astype(int)
        out_ts = perturb_binary(in_ts, P00=P00, P11=P11)

        # Save to output
        out_ts.to_csv(os.path.join(output_dir, f"RSAC128_markov_pert.csv"),
                      index=True,
                      header=False)

    # Suisun Gates -------------------------------------------------------------------------------------

    if False:  # switches so I can just re-write the necessary component

        bds_dir = '/home/tomkovic/BayDeltaSCHISM/data/time_history'
        out_dir = './data_out/suisun_gates/'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        smscg_files = {'boat': os.path.join(bds_dir, 'montezuma_boat_lock.th'),
                       'flash': os.path.join(bds_dir, 'montezuma_flash.th'),
                       'radial': os.path.join(bds_dir, 'montezuma_radial.th')}

        smscg_dfs = {}

        # header for all_clean df in perturb_suisun_ops ['radial','flashboard','boat_lock']
        dsm2_header = ['/HIST+GATE/MTZSL/RADIAL_OP//IR-DECADE/DWR-ESO/',
                       '/HIST+GATE/MTZSL/FLASHBOARD_OP//IR-DECADE/DWR-ESO/',
                       '/HIST+GATE/MTZSL/BOATLOCK_OP//IR-DECADE/DWR-ESO/']

        for gp in smscg_files:
            smscg_dfs[gp] = read_th(smscg_files[gp])

        for v in [1, 2]:
            radial_pert_clean, flash_pert_clean, boat_pert_clean, all_clean = perturb_suisun_ops(
                smscg_dfs['radial'], smscg_dfs['flash'], smscg_dfs['boat'], P00=0.99, P11=0.99)

            if False:
                # plot results
                plt_dr = [datetime.date(2006, 1, 1), datetime.date(2010, 1, 1)]
                f, axs = plt.subplots(2)
                axs[0] = smscg_dfs['radial'].plot(
                    y='op_up', use_index=True, label='Original')
                axs[1] = radial_pert_clean.plot(
                    y='op_up', use_index=True, label='Perturbed')
                axs[0].set_xlim(plt_dr)
                axs[1].set_xlim(plt_dr)
                plt.show()

            # Export back to *.th files
            write_th(radial_pert_clean, os.path.join(
                out_dir, f'montezuma_radial_lhc_v{v}.th'))
            write_th(flash_pert_clean, os.path.join(
                out_dir, f'montezuma_flash_lhc_v{v}.th'))
            write_th(boat_pert_clean, os.path.join(
                out_dir, f'montezuma_boat_lock_lhc_v{v}.th'))

            # change radial on/off to match DSM2 operations (-10 is operating tidally (0), !=-10 is open (1))
            all_clean.loc[all_clean['radial'] != 0, 'radial'] = 1
            all_clean.loc[all_clean['radial'] == 0, 'radial'] = -10
            all_clean.to_csv(f'./data_out/MTZSL_markov_pert_v{v}.csv',
                             header=dsm2_header,
                             index=True)

    if True:  # just write out SCHISM form DSM2
        bds_dir = '/home/tomkovic/BayDeltaSCHISM/data/time_history'
        out_dir = './data_out/suisun_gates/'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        smscg_files = {'boat': os.path.join(bds_dir, 'montezuma_boat_lock.th'),
                       'flash': os.path.join(bds_dir, 'montezuma_flash.th'),
                       'radial': os.path.join(bds_dir, 'montezuma_radial.th')}

        smscg_dfs = {}

        for gp in smscg_files:
            smscg_dfs[gp] = read_th(smscg_files[gp])

        dss_files = {'lhc_1': r"./data_out/latinhypercube_v3_bundles/lhc_1/lhc_1_gates.dss",
                     'lhc_2': r"./data_out/latinhypercube_v3_bundles/lhc_2/lhc_2_gates.dss",
                     'lhc_5': r"./data_out/latinhypercube_v3_bundles/lhc_5/lhc_5_gates.dss",
                     'lhc_6': r"./data_out/latinhypercube_v3_bundles/lhc_6/lhc_6_gates.dss",
                     'lhc_7': r"./data_out/latinhypercube_v3_bundles/lhc_7/lhc_7_gates.dss"}

        ts_ranges = {'lhc_1': [pd.Timestamp(2013, 10, 24), pd.Timestamp(2015, 12, 31)],
                     'lhc_2': [pd.Timestamp(2011, 10, 18), pd.Timestamp(2012, 12, 31)],
                     'lhc_5': [pd.Timestamp(2006, 11, 14), pd.Timestamp(2008, 12, 31)],
                     'lhc_6': [pd.Timestamp(2011, 10, 18), pd.Timestamp(2012, 12, 31)],
                     'lhc_7': [pd.Timestamp(2009, 10, 29), pd.Timestamp(2011, 12, 31)]}

        for case in dss_files:

            write_DSM2_to_SCHISM(dss_files[case],  smscg_dfs['radial'], smscg_dfs['flash'],
                                 smscg_dfs['boat'], out_dir, case, ts_ranges[case])
