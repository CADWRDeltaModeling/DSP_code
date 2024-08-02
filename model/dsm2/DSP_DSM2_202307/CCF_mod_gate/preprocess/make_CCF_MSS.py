# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:06:50 2024

@author: smunger

This script will create the dss input for the Clifton Court adaptive gate operation. 
You'll need to either specify a constant swp export level or a time serie in a dss format. 
You'll also need to provide the predicted tide at San Francisco for the duration of your simulation


"""

import shutil, os, glob
import pyhecdss
import pydsm
from pydsm.functions import tidalhl
import pandas as pd
import matplotlib.pyplot as plt
from vtools.functions.unit_conversions import CMS2CFS, FT2M
from dms_datastore.read_ts import read_noaa
import numpy as np

#%% Functions

def tlmax(arr):
    '''return HH(1) or LH (0)
    '''
    idx = np.argmax(arr) # only the first occurence of the maxima is return
    #print(arr,idx)
    return idx
def tlmin(arr):
    '''return LL(1) or HL(0)'''
    idx = np.argmin(arr)
    return idx

def get_tidal_hh_lh(sh):
    sth=sh.rolling(2).apply(tlmax,raw=True)
    sth.iloc[0]=0 if sth.iloc[1][0] > 0 else 1 # fill in the first value based on next value
    return sth.iloc[:,0].map({np.nan:'', 0:'LH',1:'HH'}).astype(str)

def get_tidal_ll_hl(sl):
    stl=sl.rolling(2).apply(tlmin,raw=True)
    stl.iloc[0]=0 if stl.iloc[1][0] > 0 else 1 # fill in the first value based on next value
    return stl.iloc[:,0].map({np.nan:'', 0:'HL',1:'LL'}).astype(str)

def flow_to_priority(flow, 
                     breaks = [-100, 2000, 4000.,9000.,99999.],
                     labels=[1,2,3,4]):
    """Convert export flows to priorities based on numerical brackets with breaks. 
       Labels must be integers"""
    priority = pd.cut(
      flow,
      breaks,   # These are the boundaries between priorities
      labels=labels
    ).astype(int)
    priority.name="priority"
    return priority

def flow_to_max_gate(flow, 
                     breaks = [-100, 400, 1200, 3000.,4000,99999.],
                     labels=[3,5,8,10,16]):
    """Convert export flows to max gate height on numerical brackets with breaks."""
    gmax = pd.cut(flow,breaks,labels=labels)
    gmax.name="max_gate"
    return gmax

def create_priority_series(p1,p2,p3,p4,priority, stime, etime):
    """Choose priorities day-by-day based on the value of the priority argument"""
    pgate = pd.concat([p1,p2,p3,p4],axis=1)[stime:etime]
    pgate.columns = pgate.columns.map(int) # has to be integer
    priority2 = priority.loc[pgate.index.date]
    pgate = pgate.ffill()
    pgate2 = pgate.stack()
    lookup = pd.MultiIndex.from_arrays([pgate.index,priority2.values])
    pgate2.name="op"
    pgate3 = pgate2.reindex(lookup).dropna()
    pgate3 = pgate3.loc[pgate3 != pgate3.shift(1)]
    pgate4=pgate3.reset_index()
    #pgate4=pgate4.set_index('level_0').rename(columns={'level_1':'priority'})
    pgate4=pgate4.set_index('DATETIME').rename(columns={'level_1':'priority'})
    pgate4.index.names = ['Datetime']
    return pgate4

def make_3_prio(input_tide, stime, etime):
    
    '''
    Function that makes the priorities schedule time serie based on the predicted tide at San Francisco. 
    Input: a time serie of the Predicted tide at SF in LST
    The time serie should be taken from the datastore (water level in m, NAVD88). Headers: datetime,predicted_m,elev_m. 

    Output: 3 irregular time serie that contain the schedule for the priority 1, 2 and 3.
    '''
    
    #input_tide = r'D:\Projects\DSM2\Gate_prio_for_ines\SFO_dms_data.csv' 
    
    # # range of dates for the output files
    # stime = '1997-01-01'
    # etime = '2024-12-31'
    
    if False: # os.path.exists('./prio_ts'): # suppress re-reading for now. Not working
        print('Retreaving the priorities already created from the <prio_ts> directory')
        
        prio1_ts = pd.read_csv(os.path.join('prio_ts','p1.csv'), parse_dates=True, index_col = [0])
        prio2_ts = pd.read_csv(os.path.join('prio_ts','p2.csv'), parse_dates=True, index_col = [0])
        prio3_ts = pd.read_csv(os.path.join('prio_ts','p3.csv'), parse_dates=True, index_col = [0])
        p4 = prio1_ts.rename(columns={1:4}).resample('1D').mean()*0+1
        s = 'dummy_var' # dummy variable for the return. No need to reload the csv with the tide, it's a really big file
        
    else: 
        
        print('Making priorities from tide')
    
        # wl_df = pd.read_csv(input_tide, parse_dates=True, index_col = [0])
        input_tide.columns = ['value']
        sf = input_tide/FT2M # from m to ft.
        lag = 7.25 # in hours to for the best fit with local tide
        sf_in = sf.resample('1min').interpolate(method='cubic') # resample and interpolate to fine resolution so that the lag can be in minutes
        sf_lagged = sf_in.shift(int(lag*60)).resample('5min').mean() # lag and downsample back to 5min.
        s = sf_lagged[stime:etime] # Input tidal signal for priority development
        
        
        # Find minimum and maximums
        sh,sl=tidalhl.get_tidal_hl(s) #Get Tidal highs and lows
        sh=pd.concat([sh,get_tidal_hh_lh(sh)],axis=1)
        sh.columns=['max','max_name']
        sl=pd.concat([sl,get_tidal_ll_hl(sl)],axis=1)
        sl.columns=['min','min_name']
        
       
        # --------  flagg open and close portions - OG Priority 3
        # when it opens
        idx1=sl[sl['min_name']=='LL'].index+pd.Timedelta('1h')
        idx2=sh[sh['max_name']=='HH'].index-pd.Timedelta('1h')
        idx3=sh[sh['max_name']=='LH'].index-pd.Timedelta('1h') #This is so it opens during the HL-LH-HL sequence
        ci = idx1.union(idx2).union(idx3)
        opens=pd.DataFrame(data=np.ones(len(ci)),index=ci)
        
        # when it closes
        idx1=sl[sl['min_name']=='HL'].index+pd.Timedelta('2h')
        idx2=sl[sl['min_name']=='LL'].index-pd.Timedelta('2h')
        closes=pd.DataFrame(data=np.zeros(len(sl)),index=idx1.union(idx2))
        
        # combined the df
        open_close=pd.concat([opens,closes],axis=1)
        
        prio3_ts=open_close.iloc[:,0].combine_first(open_close.iloc[:,1]).to_frame()
        prio3_ts.columns=[3]
        prio3_ts.index.name='DATETIME'
        prio3_ts = prio3_ts[prio3_ts[3] != prio3_ts[3].shift()]
        prio3_ts.head()
        
        
        #------- flagg open and close portions - Priority 2
        
        # when it opens
        idx1=sl[sl['min_name']=='LL'].index+pd.Timedelta('1h')
        idx2=sh[sh['max_name']=='HH'].index-pd.Timedelta('1h')
        idx3=sh[sh['max_name']=='LH'].index-pd.Timedelta('1h') #This is so it opens during the HL-LH-HL sequence
        ci = idx1.union(idx2).union(idx3)
        opens=pd.DataFrame(data=np.ones(len(ci)),index=ci)
        opens.head()
        
        # when it closes
        idx1=sl[sl['min_name']=='HL'].index-pd.Timedelta('1h')
        idx2=sl[sl['min_name']=='LL'].index-pd.Timedelta('2h')
        closes=pd.DataFrame(data=np.zeros(len(sl)),index=idx1.union(idx2))
        
        # combined the df
        open_close=pd.concat([opens,closes],axis=1)
        prio2_ts=open_close.iloc[:,0].combine_first(open_close.iloc[:,1]).to_frame()
        prio2_ts.columns=[2]
        prio2_ts.index.name='DATETIME'
        prio2_ts = prio2_ts[prio2_ts[2] != prio2_ts[2].shift()]
        prio2_ts.head()
    
        #------ flagg open and close portions - Priority 1
        
        # when it opens
        idx1=sh[sh['max_name']=='LH'].index+pd.Timedelta('1h')
        idx2=sh[sh['max_name']=='HH'].index+pd.Timedelta('1h')
        ci = idx1.union(idx2)
        opens=pd.DataFrame(data=np.ones(len(ci)),index=ci)
        
        # when it closes
        idx1=sl[sl['min_name']=='HL'].index-pd.Timedelta('1h')
        idx2=sl[sl['min_name']=='LL'].index-pd.Timedelta('2h')
        closes=pd.DataFrame(data=np.zeros(len(sl)),index=idx1.union(idx2))
        
        # combined the df
        open_close=pd.concat([opens,closes],axis=1)
        prio1_ts=open_close.iloc[:,0].combine_first(open_close.iloc[:,1]).to_frame()
        prio1_ts.columns=[1]
        prio1_ts.index.name='DATETIME'
        prio1_ts = prio1_ts[prio1_ts[1] != prio1_ts[1].shift()]
        prio1_ts.head()
        
        p4 = prio1_ts.rename(columns={1:4}).resample('1D').mean()*0+1 # Priority 4 is when exports are so large that gates are always open. 1 value per day.
        
        save_prio_ts('prio_ts',s, prio1_ts, prio2_ts, prio3_ts, p4)
    
    return s, prio1_ts, prio2_ts, prio3_ts, p4

def save_prio_ts(tsdir, tide_lagged, p1,p2,p3,p4):

    if not os.path.exists(tsdir):
        os.makedirs(tsdir)

    print('saving prio time serie to',tsdir)
    p1.to_csv(os.path.join(tsdir,'p1.csv'))
    p2.to_csv(os.path.join(tsdir,'p2.csv'))
    p3.to_csv(os.path.join(tsdir,'p3.csv'))
    p4.to_csv(os.path.join(tsdir,'p4.csv'))
    tide_lagged.to_csv(os.path.join(tsdir,'tide_lagged.csv'))
    

def write_priority_ts_to_dss(dss_name, prio1, prio2, prio3, prio4):
    #this is use in the case where export at SWP are constant. Then only one priority type is used 
    
    #dss_name = 'ccfb_3prios_planning.dss'
    with pyhecdss.DSSFile(dss_name,create_new=True) as d:
        d.write_its('/HIST+GATE/CCFB/GATE-PRIORITY//IR-YEAR/TYPE-4/',prio4,'POS','INST-VAL')
        d.write_its('/HIST+GATE/CCFB/GATE-PRIORITY//IR-YEAR/TYPE-3/',prio3,'POS','INST-VAL')
        d.write_its('/HIST+GATE/CCFB/GATE-PRIORITY//IR-YEAR/TYPE-2/',prio2,'POS','INST-VAL')
        d.write_its('/HIST+GATE/CCFB/GATE-PRIORITY//IR-YEAR/TYPE-1/',prio1,'POS','INST-VAL')

def write_gate_op_to_dss(dss_name, priority_df, max_gate_height,export_DM, export_15m):
    '''
    this is use in the case where export at SWP is a time serie and priority changes in time
    use function create_priority_series to create the df
    '''
    with pyhecdss.DSSFile(dss_name,create_new=True) as d:
        d.write_its('/HIST+GATE/CCFB/GATE-PRIORITY//IR-YEAR/DWR-MSS/',priority_df.priority.astype('float64'),'POS','INST-VAL')
        d.write_its('/HIST+GATE/CCFB/GATE-OP//IR-YEAR/DWR-MSS/',priority_df.op,'POS','INST-VAL')
        d.write_rts('/HIST+GATE/CCFB/GATE-MAX//1DAY/DWR-MSS/',max_gate_height,'POS','INST-VAL')
        d.write_rts('/HIST+GATE/CHSWP003/FLOW-EXPORT//1DAY/DWR-DMS/', export_DM,'CFS','INST-VAL')
        d.write_rts('/HIST+GATE/CHSWP003/FLOW-EXPORT//15MIN/DWR-DMS/', export_15m,'CFS','INST-VAL')

def get_exports(dss_in):

    data = list(pyhecdss.get_ts(dss_in,'/FILL+CHAN/CHSWP003/FLOW-EXPORT//1DAY//'))[0].data
    new_ind = data.index.astype('datetime64[ns]')
    exports_df = data.set_index(new_ind)
    exports_df.rename(columns={exports_df.columns[0]: "exports" }, inplace = True)
    return exports_df
    

    
def export_lookup(x): 

    xp = [-100, 400, 1200, 2000, 3000, 4000, 9000, 99999]
    p = [1, 1, 1, 2 ,2, 3, 4]
    max_g = [3,5,8,8,10,16,16]
    
    if (x >= xp[0]) & (x < xp[1]):  
        prio = p[0] 
        max_gate = max_g[0]
    if (x >= xp[1]) & (x < xp[2]):  
        prio = p[1] 
        max_gate = max_g[1]
    if (x >= xp[2]) & (x < xp[3]):  
        prio = p[2] 
        max_gate = max_g[2]
    if (x >= xp[3]) & (x < xp[4]):  
        prio = p[3]
        max_gate = max_g[3]
    if (x >= xp[4]) & (x < xp[5]):  
        prio = p[4] 
        max_gate = max_g[4]
    if (x >= xp[5]) & (x < xp[6]):  
        prio = p[5] 
        max_gate = max_g[5]    
    if (x >= xp[6]) & (x < xp[7]):  
        prio = p[6] 
        max_gate = max_g[6]
    
    print('Export is',x,'CFS--> Priority =',prio,' Max GH =',max_gate)
    
    return prio, max_gate
    
def make_prio_for_constant_exports(input_tide, dss_output,export_level):
    wl_df = pd.read_csv(input_tide, parse_dates=True, index_col = [0])
    stime = wl_df.index[0]
    etime = wl_df.index[-1]
    tide_lag, p1, p2, p3, p4 = make_3_prio(input_tide,stime, etime) # this is slow could be written in csv instead of recreating it everytime
   
    prio, max_gate = export_lookup(export_level) # look up corresponding gate height and prio
    
    prio_dict = {1: p1, 2:p2, 3:p3, 4:p4} 
    prio_df = pd.DataFrame(data=prio_dict[prio]) # retreive correct prio for export level
    prio_df.rename(columns={prio_df.columns[0]: "op" }, inplace = True)
    prio_df['priority']=prio
    
    max_gate_height = p4*0+max_gate
    export_DM = p4*0+export_level
    export_15m = export_DM.resample('15T').ffill()

    write_gate_op_to_dss(dss_output, prio_df, max_gate_height,export_DM, export_15m)

def make_prio_for_varying_exports(input_tide, dss_output, dss_in):
    
    export_df= get_exports(dss_in) # read the .th file and write in dss. Output swp exports for later
    stimee=export_df.index[0]
    etimee=export_df.index[-1]
    
    # wl_df = pd.read_csv(input_tide, parse_dates=True, index_col = [0])
    input_tide
    stime = max(input_tide.index[0],stimee)
    etime = min(input_tide.index[-1],etimee)
    
    export_15min = export_df.squeeze()
    export_1day = export_15min.resample('1D').mean()
        
    tide_lag, p1, p2, p3, p4 = make_3_prio(input_tide, stime, etime) 
         
    priority = flow_to_priority(export_1day) # export flows to priorities type
    priority_df = create_priority_series(p1,p2,p3,p4,priority, stime, etime) # make the priority schedule irr ts
    
    max_gate_height = flow_to_max_gate(export_1day).astype('float64') # assign max gate heights base on export level (sipping)
    export_15min =  export_1day.resample('15T').ffill()
    
    write_gate_op_to_dss(dss_output, priority_df, max_gate_height,export_1day, export_15min)
    print('writting to dss')


    
#%%

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    tidal_infile = "./data/noaa_sffpx_shift_forward_100d.csv"
    tidal_df = read_noaa(tidal_infile,force_regular=True)

    # For dsm2 runs with export that are time varying
    # point to the dss that has the swp export time series. 
    # you might have to adjust the path in the <get_exports> function
    # dss_output = 'CCF_gate_op_5.dss'
    # dss_in = '../timeseries/lhc_5_hist.dss' 
    # make_prio_for_varying_exports(tidal_df, dss_output, dss_in)

    # shutil.copyfile(dss_output, f"../CCF_inputs/CCF_gate_op_5.dss")


    dss_output = 'CCF_gate_op_3.dss'
    dss_in = '../timeseries/lhc_3_hist.dss' 
    make_prio_for_varying_exports(tidal_df, dss_output, dss_in)

    shutil.copyfile(dss_output, f"../CCF_inputs/CCF_gate_op_3.dss")