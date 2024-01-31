##!/usr/bin/env python
# -*- coding: utf-8 -* 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
        Markov transition probability for unperturbed staying unperturbed. Probability of transition to perturbed
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


if __name__ == "__main__":

    bddsm_dir = "D:/dsm2/DSM2v821/timeseries"
    primary_pathname_part_dss_filename_dict = {'RSAC128': os.path.join(bddsm_dir, "gates-v8-201712.dss"),
                                               'MTZSL': os.path.join(bddsm_dir, "gates-v8-201712.dss")}
    primary_part_c_part_dict = {'RSAC128': 'POS',
                                'MTZSL': 'RADIAL_OP'}
    unit_part_dict = {'RSAC128': 'UNSPECIF',
                      'MTZSL':'UNSPECIF'}
    primary_pathname_part = 'b_part'
    
    P00 = 0.988
    P11 = 0.9625

    df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part, 
                        primary_part_c_part_dict=primary_part_c_part_dict,
                        primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)
    
    output_dir = r"D:\projects\delta_salinity\scripts\DSP_code\scripts\boundary_generation\data_out"

    for gate in primary_pathname_part_dss_filename_dict.keys():
        print(gate)
        
        ## Read in the gate operation data in th format
        in_ts = df_input[gate]

        # need to fill NaNs and convert from 0-2 float to 0 or 1 binary
        in_ts = in_ts.ffill()

        # two versions of perturbed binary for Suisun
        if 'MTZSL' == gate:
            in_ts[in_ts==-10] = 1
            in_ts[in_ts!=-10] = 0

            # v1
            out_ts = perturb_binary(in_ts,P00=P00,P11=P11)
            out_ts[out_ts==1] = -10
            # Save to output
            out_ts.to_csv(os.path.join(output_dir,f"{gate}_markov_pert_v1.csv"), 
                        index=True,
                        header=False)
            
            # v2
            out_ts = perturb_binary(in_ts,P00=P00,P11=P11)
            out_ts[out_ts==1] = -10
            # Save to output
            out_ts.to_csv(os.path.join(output_dir,f"{gate}_markov_pert_v2.csv"), 
                        index=True,
                        header=False)
        else:
            in_ts = (in_ts>0.1).astype(int)
            out_ts = perturb_binary(in_ts,P00=P00,P11=P11)

            # Save to output
            out_ts.to_csv(os.path.join(output_dir,f"{gate}_markov_pert.csv"), 
                        index=True,
                        header=False)
