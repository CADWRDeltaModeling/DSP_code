##!/usr/bin/env python
# -*- coding: utf-8 -* 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime
# from dms_datastore.read_ts import read_ts
# from dms_datastore.read_multi import read_ts_repo

from pydelmod.create_ann_inputs import get_dss_data
import os


def perturb_binary(orig_ts,P00=0.988,P11=0.9625):
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
    ts.columns=["orig"]
    ts["perturbed"] = ts["orig"]

    ts['use_orig'] = True
    
    USEORIGNDX = 2 #ts.columns.index('use_orig')

    # Transitions
    rnum = np.random.rand(len(ts))
    for i in range(1,len(ts)):
        if ts.iloc[i-1,USEORIGNDX]:
            ts.iloc[i,USEORIGNDX] = True if rnum[i-1] < P00 else False
        else:
            ts.iloc[i,USEORIGNDX] = False if rnum[i-1] < P11 else True


    ts["perturbed"] = ts.orig.where(ts.use_orig,(ts.perturbed+1.)%2)
    ts.drop("use_orig",axis=1,inplace=True)
    ts.drop("orig",axis=1,inplace=True)
    return ts


def sample_binary_ts():
    ndx = pd.date_range(pd.Timestamp(2007,1,1),freq='d',periods=1000)
    ts = pd.Series(index=ndx,data=None)
    ts[:]=1.0


    off_periods = [ ("2007-6-01","2007-09-15"),("2008-03-01","2008-05-05"),
                ("2008-10-01","2008-11-15"),("2009-03-01","2009-05-05")]

    for off in off_periods:
        ts.loc[pd.to_datetime(off[0]):pd.to_datetime(off[1])] = 0

    ts.name = "pos"
    return ts

#### example:
# ts_orig = sample_binary_ts()
# P00=0.988
# P11=0.9625
# ts = perturb_binary(ts_orig,P00,P11)  # Defaults are fine
# print(ts)
# ts["perturbed"]*=0.75

# ts.plot() #drawstyle="steps-post")
# plt.legend()
# plt.show()

def perturb_suisun_ops(radial_df, flash_df, boat_df, P00=0.990,P11=0.990):
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
    if not (radial_df.index.min()==flash_df.index.min()==boat_df.index.min()): 
        print(f"""START TIMES NOT ALIGNED! radial: {radial_df.index.min()} 
                         flashboard: {flash_df.index.min()} 
                         boat lock: {boat_df.index.min()}""")
        ts = [t for t in ts if t>=max(radial_df.index.min(),flash_df.index.min(),boat_df.index.min())]
    if not (radial_df.index.max()==flash_df.index.max()==boat_df.index.max()): 
        print(f"""END TIMES NOT ALIGNED! radial: {radial_df.index.max()} 
                       flashboard: {flash_df.index.max()} 
                       boat lock: {boat_df.index.max()}""")
        ts = [t for t in ts if t<=min(radial_df.index.max(),flash_df.index.max(),boat_df.index.max())]
    print(f"Time series to be perturbed: {ts[0]} - {ts[-1]}")

    all_df = pd.DataFrame(index=ts, columns=['radial','flashboard','boat_lock'])
    all_df['radial'] = radial_pert_df['op_up']
    all_df['flashboard'] = flash_pert_df['op_up']
    all_df['boat_lock'] = boat_pert_df['op_up']

    # Transitions
    rnum = np.random.rand(len(ts))
    for i in range(1,len(ts)):
        if radial_pert_df.loc[radial_pert_df.index[i-1],'op_up']==0:
            # radial gate is operational
            if rnum[i-1]>=P00:
                radial_pert_df.loc[radial_pert_df.index[i],['op_up','op_down']] = [1.0,1.0] # change to open if above P00 threshold
                all_df.loc[radial_pert_df.index[i],'radial'] = 1.0
        else:
            # radial gate is open (op_up!=0)
            if rnum[i-1] >= P11:
                radial_pert_df.loc[radial_pert_df.index[i],['op_up','op_down']] = [0.0,0.0] # change to operational if above P00 threshold
                all_df.loc[radial_pert_df.index[i],'radial'] = 0.0
        # Check that boat and flash are operated appropriately
        if radial_pert_df.loc[radial_pert_df.index[i-1],'op_up']==0:
            flash_pert_df.loc[radial_pert_df.index[i],['op_up','op_down']]=[0.0,0.0] # flashboards need to be IN
            all_df.loc[radial_pert_df.index[i],'flashboard'] = 0.0
            boat_pert_df.loc[radial_pert_df.index[i],['op_up','op_down']]=[1.0,1.0] # boatlock needs to be OPEN
            all_df.loc[radial_pert_df.index[i],'boat_lock'] = 1.0

    # remove any rows that are unecessary
    radial_pert_clean = clean_df(radial_pert_df)
    flash_pert_clean = clean_df(flash_pert_df)
    boat_pert_clean = clean_df(boat_pert_df)
    all_clean = clean_df(all_df)

    return radial_pert_clean, flash_pert_clean, boat_pert_clean, all_clean

def clean_df(in_df):
    keep_rows = [0]
    for r in range(1,len(in_df.index)):
        if not (list(in_df.iloc[r,:])==list(in_df.iloc[r-1,:])):
            # keep the row where something changes
            keep_rows.append(r)
    out_df = in_df.iloc[keep_rows,:]

    return out_df

def read_th(infile):

    in_df = pd.read_table(infile,
                          delim_whitespace=True,
                          index_col='datetime',
                          comment="#")
    in_df.index = pd.to_datetime(in_df.index,format="%Y-%m-%dT%H:%M")

    # monotonic increase check
    mon_inc = all(x<y for x, y in zip(in_df.index, in_df.index[1:])) # True for monotonically-increasing data
    if not mon_inc: print(in_df.loc[in_df.index.to_series().diff() < pd.to_timedelta('0 seconds')]) # prints the row(s) where monotonicity is broken

    in_df = in_df.reindex(
        pd.date_range(start=in_df.index.min(),
                      end=in_df.index.max(),
                      freq='D'),
                      method='ffill')
    
    print(f"TH FILE ==== {infile} min: {in_df.index.min()} max:{in_df.index.max()}")
    return in_df

def write_th(in_df, outfile):
    in_df.index = in_df.index.strftime("%Y-%m-%dT%H:%M")
    
    in_df.to_csv(outfile,
                 sep=' ',
                 header=True,
                 index=True,
                 index_label='datetime')



if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # DCC Gates ---------------------------------------------------------------------------------------

    if False:
        bddsm_dir = "D:/dsm2/DSM2v821/timeseries"
        primary_pathname_part_dss_filename_dict = {'RSAC128': os.path.join(bddsm_dir, "gates-v8-201712.dss")}
        primary_part_c_part_dict = {'RSAC128': 'POS'}
        unit_part_dict = {'RSAC128': 'UNSPECIF'}
        primary_pathname_part = 'b_part'
        
        P00 = 0.988
        P11 = 0.9625

        df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part, 
                            primary_part_c_part_dict=primary_part_c_part_dict,
                            primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)
        
        output_dir = r"D:\projects\delta_salinity\scripts\DSP_code\scripts\boundary_generation\data_out"

        print("DCC")
        
        ## Read in the gate operation data in th format
        in_ts = df_input['RSAC128']

        # need to fill NaNs and convert from 0-2 float to 0 or 1 binary
        in_ts = in_ts.ffill()

        in_ts = (in_ts>0.1).astype(int)
        out_ts = perturb_binary(in_ts,P00=P00,P11=P11)

        # Save to output
        out_ts.to_csv(os.path.join(output_dir,f"RSAC128_markov_pert.csv"), 
                    index=True,
                    header=False)

    # Suisun Gates -------------------------------------------------------------------------------------
    bds_dir = '/home/tomkovic/BayDeltaSCHISM/data/time_history'
    out_dir = './data_out/suisun_gates/'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    smscg_files = {'boat': os.path.join(bds_dir, 'montezuma_boat_lock.th'),
                   'flash':os.path.join(bds_dir, 'montezuma_flash.th'),
                   'radial':os.path.join(bds_dir, 'montezuma_radial.th')}

    smscg_dfs = {}

    for gp in smscg_files:
        smscg_dfs[gp] = read_th(smscg_files[gp])

    for v in [1,2]:
        radial_pert_clean, flash_pert_clean, boat_pert_clean, all_clean = perturb_suisun_ops(smscg_dfs['radial'], smscg_dfs['flash'], smscg_dfs['boat'], P00=0.99555,P11=0.99333)

        if False:
            # plot results
            plt_dr = [datetime.date(2006,1,1), datetime.date(2010,1,1)]
            f, axs = plt.subplots(2)
            axs[0] = smscg_dfs['radial'].plot(y='op_up', use_index=True, label='Original')
            axs[1] = radial_pert_clean.plot(y='op_up', use_index=True, label='Perturbed')
            axs[0].set_xlim(plt_dr)
            axs[1].set_xlim(plt_dr)
            plt.show()

        # Export back to *.th files
        write_th(radial_pert_clean, os.path.join(out_dir,f'montezuma_radial_lhc_v{v}.th'))
        write_th(flash_pert_clean, os.path.join(out_dir,f'montezuma_flash_lhc_v{v}.th'))
        write_th(boat_pert_clean, os.path.join(out_dir,f'montezuma_boat_lock_lhc_v{v}.th'))

        all_clean.to_csv(f'./data_out/MTZSL_markov_pert_v{v}.csv',
                         header=False, 
                         index=True)
