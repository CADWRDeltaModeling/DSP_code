##!/usr/bin/env python
## -*- coding: utf-8 -* 
import pandas as pd
import numpy as np
import os
import datetime as dt
import matplotlib.pyplot as plt
from pydelmod.create_ann_inputs import *



def adjust_src_sink(src,sink,perturb,sinkfrac=0.001):
    """ Distributes a perturbation to dataframe sc and sink
    Parameters
    ----------

    src : DataFrame
        DataFrame of positive sources

    sink : DataFrame
        DataFrame of sinks. Must be all of the same sign. If positive, they are assumed to represent the
    the magnitude of the (latent) sinks. If they are negative, they are the sink values themselves. Return
    value will match the convention of the input.

    perturb : Series
        Series representing a perturbation to be distributed proportionally over src and sink. Sign convention
    agrees with that of src. When perturb is negative, it augments sink. When it is positive 
    it reduces the sink until sinkfrac fraction of the sink
    is exhausted, after which the remainder augments the source.

    sinkfrac : float
        Fraction of sink that can be reduced in order to achieve the perturbation by reducing the sink. Remainder
    is applied by augmenting src. Has no effect on time steps when perturb is negative.
            
    """


    # convert so that sink is a magnitude and switch the sign so that it contributes to sink
    # deal with magnitudes and in terms of consumptive use as positive,
    # will undo at the end to match convention of the inputs
    

    if (sink<0).any(axis=None): 
        negsink = True
        hasneg = sink.loc[(sink<0).any(axis=1),:]
        colswithneg = (hasneg<0).any(axis=0)
        hasneg.loc[:,colswithneg].to_csv("test.csv")
        print("negsink")
        print(sink.loc[(sink<0).any(axis=1),:])
    else: 
        negsink = False
        print("not",(sink>=0).all(axis=None),)
        sink = -sink 

    src_total = src.sum(axis=1)
    sink_total = sink.sum(axis=1)
    dcd_total = src_total + sink_total

    if 'Timestamp' not in str(type(perturb.index[0])):
        perturb.index = perturb.index.to_timestamp()
    pert = perturb.reindex(dcd_total.index) # using only the indices in dcd_total (from the src/sink inputs)
    pert
    neg = pert <= 0.
    pos = ~neg

    # For negative (e.g. more depletion) , achieve change by scaling all sinks
    scaleneg = (sink_total+pert)/sink_total # all quantities negative, so scaleneg > 1
    scaleneg = scaleneg.where(neg,1.)      # Ensure No effect in pos case 

    # For positive adjustments, remove up to FRAC fraction of total sinks, then ... 
    adjustmag = pert.where(pos,0.) # total positive adjustment required
    FRAC = sinkfrac
    sinklimit = -FRAC*sink_total     # max amount done by reducing sink (this is pos) 
    from_sink =  adjustmag.clip(upper=sinklimit) 
    residual = adjustmag - from_sink # This should be nonnegative for pertinent cases where ~pos
    residual = residual.where(pos,0.)  # Just to clarify
    resscale = (residual + src_total)/src_total
    src = src.mul(resscale,axis=0)
    
    # Now put together the adjustments to sinks under both cases, pos and ~pos
    scalepos = (from_sink + sink_total)/sink_total    # for positive (source increase)
    scale = scalepos.where(pos,scaleneg)              # choose case
    sink=sink.mul(scale,axis=0)
    if not negsink:
        sink = -sink
        #sink.loc[sink==0.] = -0.0
    return src,sink

def get_net_srcsnk(dcd_dss_file):
    start=pd.Timestamp(1980,1,1)
    end=pd.Timestamp(2022,10,1)        
    div_seep_dcd_c_part_dss_filename_dict = {'DIV-FLOW': dcd_dss_file, 'SEEP-FLOW': dcd_dss_file}
    drain_dcd_c_part_dss_filename_dict = {'DRAIN-FLOW': dcd_dss_file}
    df_div_seep_dcd = get_dss_data(div_seep_dcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True).loc[start:end]
    df_div_seep_dcd = df_div_seep_dcd.clip(lower=0.0)
    df_drain_dcd = get_dss_data(drain_dcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True).loc[start:end]
    
    df_div_seep_dcd.columns = [strip_dpart(col) for col in df_div_seep_dcd.columns]
    df_drain_dcd.columns = [strip_dpart(col) for col in df_drain_dcd.columns]
    
    use_cols = [col for col in df_div_seep_dcd.columns if "total" not in col.lower()]
    df_div_seep_dcd=df_div_seep_dcd[use_cols]
    use_cols = [col for col in df_drain_dcd.columns if "total" not in col.lower()]
    df_drain_dcd=df_drain_dcd[use_cols]

    src0=df_drain_dcd.sum(axis=1)
    sink0=df_div_seep_dcd.sum(axis=1)
    serieslen = len(src0)

    return src0, sink0

def dcd_from_dsm2_pert(dcd_dss_file, dsm2_dcd_dss_file, schism_in, out_dir, version):
    
    src0, sink0 = get_net_srcsnk(dcd_dss_file) # source and sink for historical data
    psrc0, psink0 = get_net_srcsnk(dsm2_dcd_dss_file) # perturbed source and sink for this version

    net = src0 - sink0 # net flow for source/sink for historical
    pnet = psrc0 - psink0 # net flow for source/sink for perturbed data for this version

    perturb = pnet - net # the difference in net flows to be applied to the schism input/output data
    perturb = perturb.div(35.3147) # convert to cms

    sch_src = pd.read_csv(os.path.join(schism_in,'vsource_dated.th'), sep=' ', header=3, index_col=0, parse_dates=['datetime'], dtype=np.float64)
    sch_sink = pd.read_csv(os.path.join(schism_in,'vsink_dated.th'), sep=' ', header=5, index_col=0, parse_dates=['datetime'], dtype=np.float64)

    src,sink = adjust_src_sink(sch_sink, sch_src, perturb) # create the adjusted source/sink values to be used for this version in SCHISM
    src.index = src.index.strftime('%Y-%m-%dT00:00')
    sink.index = sink.index.strftime('%Y-%m-%dT00:00')

    fn_src = os.path.join(out_dir, f'vsource_{version}_dated.th')
    fn_sink = os.path.join(out_dir, f'vsink_{version}_dated.th')

    sink = -sink
    
    src.to_csv(fn_src, sep=' ', float_format="%.2f")
    sink.to_csv(fn_sink, sep=' ', float_format="%.2f")


def strip_dpart(colname):
    parts = colname.split("/")
    parts[4]=""
    return "/".join(parts)

if __name__ == "__main__":
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    dcd_dss_file = "../../model/dsm2/2021DSM2FP_202301/timeseries/DCD_hist_Lch5.dss"
    dsm2_dcd_dss_file = "../../model/dsm2/DSP_DSM2_202307/latinhypercube_v2/timeseries/lhc_1_dcd.dss"
    version = 'v1'

    schism_in = "/home/tomkovic/BayDeltaSCHISM/data/channel_depletion"
    out_dir = "./data_out/schism_dcd/" 

    dcd_from_dsm2_pert(dcd_dss_file, dsm2_dcd_dss_file, schism_in, out_dir, version)
    
    dsm2_dcd_dss_file = "../../model/dsm2/DSP_DSM2_202307/latinhypercube_v2/timeseries/lhc_2_dcd.dss"
    version = 'v2'
    
    dcd_from_dsm2_pert(dcd_dss_file, dsm2_dcd_dss_file, schism_in, out_dir, version)
