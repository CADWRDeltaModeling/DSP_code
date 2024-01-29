##!/usr/bin/env python
## -*- coding: utf-8 -* 
import pandas as pd
import numpy as np
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
    pert = perturb.reindex(dcd_total.index)
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



def dcd_from_dss(dcd_dss_file, out_dir, version, pvar1, pvar2, ptvar1, ptvar2):
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

    net0 = src0 - sink0
    perturb = pvar1*np.sin(2.*np.pi*np.arange(serieslen)/ptvar1) +  pvar2*np.sin(2.*np.pi*np.arange(serieslen)/ptvar1)

    perturb = net0*0.+perturb # convert to dataframe
    src,sink=adjust_src_sink(df_drain_dcd,df_div_seep_dcd,perturb)
    
    src.to_csv(os.path.join(out_dir,os.path.basename(dcd_dss_file).replace('.dss',f'_{version}_source.csv')))
    sink.to_csv(os.path.join(out_dir,os.path.basename(dcd_dss_file).replace('.dss',f'_{version}_sink.csv')))

    # fig,(ax0,ax1,ax2)=plt.subplots(3,sharex=True)
    
    # print("div seep ==============")
    # print (df_div_seep_dcd)
    # print("drain-----------------")
    # print(df_drain_dcd)
    # src0.plot(ax=ax0)
    # src_tot=src.sum(axis=1)
    # src_tot.plot(ax=ax0)
    # perturb.plot(ax=ax0)   
    # ax0.legend(["Drain (Orig. Src)","Src"])
 
    # sink0.plot(ax=ax1)
    # sink_tot=sink.sum(axis=1)
    # sink_tot.plot(ax=ax1)
    # ax1.legend(["Div/Seep (Orig Snk)","sink"])

    # net0.plot(ax=ax2)
    # newnet = src_tot - sink_tot
    # diff = newnet - net0
    # newnet.plot(ax=ax2)
    # diff.plot(ax=ax2)
    # ax2.legend(["net","net perturbed","diff"])
    # ax0.grid();ax1.grid();ax2.grid()
    # plt.show()


def strip_dpart(colname):
    parts = colname.split("/")
    parts[4]=""
    return "/".join(parts)

    
def test_adjust():
    dr = pd.date_range(pd.Timestamp(1990,1,1),freq='D',periods=4)
    source_arr = np.array([[3.,2.5,5.],[0.,5.,1.],[0.,1.,2.],[4.,5.,6.]])
    src = pd.DataFrame(data=source_arr,index=dr,columns=["N0","N1","N3"])
    sink_arr = np.array([[1.,2.,5.],[0.,12.,3.],[1.,8.,3.],[4.,3.,4.]])
    snk = -pd.DataFrame(data=sink_arr,index=dr,columns=["N0","N1","N3"])    
    dr = pd.date_range(pd.Timestamp(1989,12,31),freq='D',periods=6)
    adjust = -pd.Series(data=np.array([0.,5.,-3.,-4.,-20.,37.]),index=dr)
    
    print("Adjust")
    print(adjust)
    print("Original src")
    print(src)
    print("Original sink")
    print(snk)
    newsrc,newsink = adjust_src_sink(src,snk,adjust)
    print("*****")
    print("New src")
    print(newsrc)
    print("New sink")
    print(newsink)

    total = src.sum(axis=1).to_frame()
    total.columns=["src"]
    total["sink"]=snk.sum(axis=1) 
    total["src_new"]=newsrc.sum(axis=1)
    total["sink_new"]=newsink.sum(axis=1)

    total["adjust"]=adjust
    print("Marginal totals")
    print(total)

if __name__ == "__main__":
    
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    hist_fn = "../../model/dsm2/2021DSM2FP_202301/timeseries/DCD_hist_Lch5.dss"

    out_dir = "./data_out"

    version = 'v1'
    pvar1 = 1200.
    pvar2 = 200.
    ptvar1 = 100.
    ptvar2 = 16.

    dcd_from_dss(hist_fn, out_dir, version, pvar1, pvar2, ptvar1, ptvar2)
    
    version = 'v2'
    pvar1 = 1200.
    pvar2 = 500.
    ptvar1 = 200.
    ptvar2 = 30.

    dcd_from_dss(hist_fn, out_dir, version, pvar1, pvar2, ptvar1, ptvar2)






