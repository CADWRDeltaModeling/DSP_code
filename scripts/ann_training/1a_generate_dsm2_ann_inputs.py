import pandas as pd
from vtools.functions.filter import cosine_lanczos
import pyhecdss
import os

from pydelmod.create_ann_inputs import get_dss_data


def process_gate_data(dss_filename, b_part, c_part, map_zero_one=['df==0','df==1'], startDateStr=None, endDateStr=None):
    '''
    Read delta cross-channel gate operation data
    Create daily time series indicating fraction of maximum gate opening (100% means both gates open all day).
    '''
    
    with pyhecdss.DSSFile(dss_filename) as d:
        fdname, generated = d._check_condensed_catalog_file_and_recatalog(condensed=True)
        catdf = pyhecdss.DSSFile._read_catalog_dsd(fdname)

        filtered_df = catdf[(catdf.B == b_part) & (catdf.C == c_part)]
        p = d.get_pathnames(filtered_df)[0]

        if d.parse_pathname_epart(p).startswith('IR-'):
            df, units, ptype = d.read_its(p, startDateStr=startDateStr, endDateStr=endDateStr)
        else:
            df,units,ptype=d.read_rts(p)

    df[eval(map_zero_one[0])] = 0
    df[eval(map_zero_one[1])] = 1

    # resample to 1 minute, then fill forward (with last value)
    df_1min = df.resample('T', closed='right').ffill()
    # now find daily averages of one minute data
    df_daily_avg = df_1min.resample('D', closed='right').mean()
    df_daily_avg = df_daily_avg.rename(columns={df_daily_avg.columns[0]:'gate_pos'})

    return df_daily_avg

def generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder,
                        b_part_e_part_dict = None):
    '''
    1. Northern flow = Sum(Sac, Yolo, Moke, CSMR, Calaveras, -NBA)
    2. San Joaquin River flow (the model input time series)
    3. Exports: Sum(Banks, Jones, CCC plants(Rock Sl, Middle R (actually old river), Victoria))
    4. DCC gate operation as daily percentage
    5. Net Delta CU, daily (DIV+SEEP-DRAIN) for DCD and SMCD
    6. Tidal Energy: daily max-daily min
    7. SJR inflow salinity at vernalis, daily
    8. Sacramento River EC
    9. EC Output for various locations
    '''

    if not os.path.exists(output_folder): os.makedirs(output_folder)

    print('hist_dss_file='+hist_dss_file)
    print('gate_dss_file='+gate_dss_file)
    print('dcd_dss_file='+dcd_dss_file)
    print('smcd_dss_file='+smcd_dss_file)
    print('model_ec_file='+model_ec_file)
    print('output_folder='+output_folder)

    #################
    # Northern Flow #
    #################
    print('northern flow: hist_dss_file='+hist_dss_file)
    b_part_dss_filename_dict = {'RSAC155': hist_dss_file, 'BYOLO040': hist_dss_file, 'RMKL070': hist_dss_file, \
        'RCSM075': hist_dss_file, 'RCAL009': hist_dss_file, 'SLBAR002': hist_dss_file}
    df_northern_flow = get_dss_data(b_part_dss_filename_dict, 'b_part')
    print('northern flow columns='+str(df_northern_flow.columns))
    df_northern_flow.fillna(0, inplace=True)
    df_northern_flow['northern_flow'] = df_northern_flow['RSAC155'] + df_northern_flow['BYOLO040']+df_northern_flow['RMKL070'] +\
        df_northern_flow['RCSM075'] + df_northern_flow['RCAL009']-df_northern_flow['SLBAR002']
    
    out_df = df_northern_flow['northern_flow'].to_frame()
    out_df['sac_flow'] = df_northern_flow['RSAC155']

    #############
    # SJR Flow  #
    #############
    b_part_dss_filename_dict = {'RSAN112': hist_dss_file}
    b_part_c_part_dict = {'RSAN112': 'FLOW'}
    df_sjr_flow = get_dss_data(b_part_dss_filename_dict, 'b_part', b_part_c_part_dict)
    
    df_sjr_flow.columns = ['sjr_flow']
    out_df['sjr_flow'] = df_sjr_flow

    ###############################################################################################
    # 3. Exports: Sum(Banks, Jones, CCC plants(Rock Sl, Middle R (actually old river), Victoria)) #
    ###############################################################################################
    b_part_dss_filename_dict = {'CHSWP003': hist_dss_file, 'CHDMC004': hist_dss_file, 'CHCCC006': hist_dss_file,\
        'ROLD034': hist_dss_file, 'CHVCT001': hist_dss_file}
    df_exports_flow = get_dss_data(b_part_dss_filename_dict, 'b_part')
    df_exports_flow['exports'] = df_exports_flow['CHSWP003']+df_exports_flow['CHDMC004']+df_exports_flow['CHCCC006']+\
        df_exports_flow['ROLD034']+df_exports_flow['CHVCT001']
    
    out_df['exports'] = df_exports_flow['exports']

    #############################################
    # 4. DCC gate operation as daily percentage #
    #############################################
    b_part = 'RSAC128'
    c_part = 'POS'
    dcc_op_df = process_gate_data(gate_dss_file, b_part, c_part, map_zero_one=['df<2','df==2'])
    out_df['dcc'] = dcc_op_df
    
    ################################################
    # 5. Suisun gate operation as daily percentage #
    ################################################
    b_part = 'MTZSL'
    c_part = 'RADIAL_OP'
    suisun_op_df = process_gate_data(gate_dss_file, b_part, c_part, map_zero_one = ['df!=-10','df==-10'], startDateStr='01JAN1953', endDateStr='01JAN2030')
    out_df['smscg'] = suisun_op_df

    ############################################################
    # 6. Net Delta CU, daily (DIV+SEEP-DRAIN) for DCD and SMCD #
    ############################################################
    div_seep_dcd_c_part_dss_filename_dict = {'DIV-FLOW': dcd_dss_file, 'SEEP-FLOW': dcd_dss_file}
    div_seep_smcd_c_part_dss_filename_dict = {'DIV-FLOW': smcd_dss_file, 'SEEP-FLOW': smcd_dss_file}
    drain_dcd_c_part_dss_filename_dict = {'DRAIN-FLOW': dcd_dss_file}
    drain_smcd_c_part_dss_filename_dict = {'DRAIN-FLOW': smcd_dss_file}

    df_div_seep_dcd = get_dss_data(div_seep_dcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True)
    df_div_seep_smcd = get_dss_data(div_seep_smcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True)
    df_drain_dcd = get_dss_data(drain_dcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True)
    df_drain_smcd = get_dss_data(drain_smcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True)

    df_div_seep_dcd['dcd_divseep_total']=df_div_seep_dcd[df_div_seep_dcd.columns].sum(axis=1)
    df_div_seep_smcd['smcd_divseep_total']=df_div_seep_smcd[df_div_seep_smcd.columns].sum(axis=1)

    df_drain_dcd['dcd_drain_total']=df_drain_dcd[df_drain_dcd.columns].sum(axis=1)
    df_drain_smcd['smcd_drain_total']=df_drain_smcd[df_drain_smcd.columns].sum(axis=1)

    cu_total_dcd = pd.merge(df_div_seep_dcd, df_drain_dcd, how='left', left_index=True, right_index=True)
    cu_total_smcd = pd.merge(df_div_seep_smcd, df_drain_smcd, how='left', left_index=True, right_index=True)
    cu_total = pd.merge(cu_total_dcd, cu_total_smcd, how='left', left_index=True, right_index=True)

    cu_total['cu_total']=cu_total['dcd_divseep_total']+cu_total['smcd_divseep_total']-cu_total['dcd_drain_total']-cu_total['smcd_drain_total']

    out_df['delta_cu'] = cu_total['cu_total']

    # Net Delta Outflow RMA:  (Sac, SJR, Yolo Bypass, Mokelumne, Cosumnes, Calaveras), exports (SWP, CVP, Contra Costa, NBA) and net DICU
    out_df['ndo'] = out_df['northern_flow'] + out_df['sjr_flow'] - out_df['exports'] - out_df['delta_cu']

    ########################################
    # 7. Tidal Energy: daily max-daily min #
    ########################################
    b_part_dss_filename_dict={'RSAC054': hist_dss_file}
    b_part_c_part_dict={'RSAC054': 'STAGE'}
    df_mtz_stage = get_dss_data(b_part_dss_filename_dict, 'b_part', \
        primary_part_c_part_dict=b_part_c_part_dict, daily_avg=False)
    df_filter = cosine_lanczos(df_mtz_stage.copy(), cutoff_period ='40H', padtype='odd')
    df_nrg = cosine_lanczos((df_mtz_stage-df_filter)**2, cutoff_period ='40H', padtype='odd') # = < (z- <z>)^2 >
    df_mrz_tidal_filter = df_filter.resample('D', closed='right').mean()
    df_mrz_tidal_filter.columns=['tidal_filter']
    df_mrz_tidal_energy = df_nrg.resample('D', closed='right').mean()
    df_mrz_tidal_energy.columns=['tidal_energy']

    out_df['mrz_tidal_energy'] = df_mrz_tidal_energy['tidal_energy']
    out_df['mrz_tidal_filter'] = df_mrz_tidal_filter['tidal_filter']

    ######################################
    # 8. EC Output for various locations #
    ######################################
    b_part_dss_filename_dict = {'CHDMC006': model_ec_file, 'CHSWP003': model_ec_file,\
        'CHVCT000': model_ec_file, 'OLD_MID': model_ec_file, 'ROLD024': model_ec_file,
        'ROLD059': model_ec_file, 'RSAC064': model_ec_file, 'RSAC075': model_ec_file,
        'RSAC081': model_ec_file, 'RSAC092': model_ec_file, 'RSAC101': model_ec_file,
        'RSAN007': model_ec_file, 'RSAN018': model_ec_file, 'RSAN032': model_ec_file,
        'RSAN037': model_ec_file, 'RSAN058': model_ec_file, 'RSAN072': model_ec_file,
        'RSMKL008': model_ec_file, 'SLCBN002': model_ec_file, 'SLDUT007': model_ec_file,
        'SLMZU011': model_ec_file, 'SLMZU025': model_ec_file, 'SLSUS012': model_ec_file,
        'SLTRM004': model_ec_file, 'SSS': model_ec_file, 'RSAC054': hist_dss_file}
    b_part_c_part_dict = {'CHDMC006': 'EC', 'CHSWP003': 'EC',\
        'CHVCT000': 'EC', 'OLD_MID': 'EC', 'ROLD024': 'EC',
        'ROLD059': 'EC', 'RSAC064': 'EC', 'RSAC075': 'EC',
        'RSAC081': 'EC', 'RSAC092': 'EC', 'RSAC101': 'EC',
        'RSAN007': 'EC', 'RSAN018': 'EC', 'RSAN032': 'EC',
        'RSAN037': 'EC', 'RSAN058': 'EC', 'RSAN072': 'EC',
        'RSMKL008': 'EC', 'SLCBN002': 'EC', 'SLDUT007': 'EC',
        'SLMZU011': 'EC', 'SLMZU025': 'EC', 'SLSUS012': 'EC',
        'SLTRM004': 'EC', 'SSS': 'EC', 'RSAC054': 'EC'}
    if b_part_e_part_dict == None:
        b_part_e_part_dict = {'CHDMC006': '15MIN', 'CHSWP003': '15MIN',\
            'CHVCT000': '15MIN', 'OLD_MID': '15MIN', 'ROLD024': '15MIN',
            'ROLD059': '15MIN', 'RSAC064': '15MIN', 'RSAC075': '15MIN',
            'RSAC081': '15MIN', 'RSAC092': '15MIN', 'RSAC101': '15MIN',
            'RSAN007': '15MIN', 'RSAN018': '15MIN', 'RSAN032': '15MIN',
            'RSAN037': '15MIN', 'RSAN058': '15MIN', 'RSAN072': '15MIN',
            'RSMKL008': '15MIN', 'SLCBN002': '15MIN', 'SLDUT007': '15MIN',
            'SLMZU011': '15MIN', 'SLMZU025': '15MIN', 'SLSUS012': '15MIN',
            'SLTRM004': '15MIN', 'SSS': '15MIN', 'RSAC054': '1HOUR'}
    df_model_ec = get_dss_data(b_part_dss_filename_dict, 'b_part', \
        primary_part_c_part_dict=b_part_c_part_dict, primary_part_e_part_dict=b_part_e_part_dict)
    
    # now add model output ec near CCC intakes
    b_part_dss_filename_dict = {'ROLD034': model_ec_file, 'SLRCK005': model_ec_file}
    b_part_c_part_dict = {'ROLD034': 'EC', 'SLRCK005': 'EC'}
    df_model_ec_2 = get_dss_data(b_part_dss_filename_dict, 'b_part', b_part_c_part_dict)
    df_model_ec = pd.merge(df_model_ec, df_model_ec_2, how='outer', left_index=True, right_index=True)
    
    # now add another copy of Mtz ec
    df_model_ec['RSAC054'] = df_model_ec['RSAC054']
        
    # now rename some of the columns
    col_rename_dict = {'CHDMC006': 'trp', 'CHSWP003': 'wci', 'CHVCT000': 'vcu', 
                       'OLD_MID': 'uni', 'ROLD024': 'rsl', 'ROLD059': 'old', 
                       'RSAC064': 'pct', 'RSAC075': 'mal', 'RSAC081': 'cse',
                       'RSAC092': 'emm2', 'RSAC101': 'srv', 'RSAN007': 'anc', 
                       'RSAN018': 'jer', 'RSAN032': 'sal', 'RSAN037': 'ppt', 
                       'RSAN058': 'rri2', 'RSAN072': 'bdt', 'RSMKL008': 'lps',
                       'SLCBN002': 'snc', 'SLDUT007': 'dsj', 'SLMZU011': 'bdl', 
                       'SLMZU025': 'nsl2', 'SLSUS012': 'vol', 'SLTRM004': 'tss', 
                       'SSS': 'sss', 'ROLD034': 'oh4', 'SLRCK005': 'inb', 'RSAC054': 'mrz'}

    df_model_ec.rename(columns=col_rename_dict, inplace=True)

    out_df = pd.merge(out_df, df_model_ec, left_index=True, right_index=True, how='outer')

    return out_df

def run_ann_input(in_fname, e_part_dict=None):

    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)

    experiment = inputs.get('experiment')

    base_study_folder = inputs.get('base_study_folder').format(**locals())
    model_input_folder = inputs.get('model_input_folder').format(**locals())
    model_folder = inputs.get('model_folder').format(**locals())
    model_output_folder = inputs.get('model_output_folder').format(**locals())
    smcd_dss_file = inputs.get('smcd_dss_file').format(**locals())
    
    col_order = ['model','scene','case',
                 'northern_flow','sac_flow','sjr_flow','exports','delta_cu','ndo',
                 'dcc','smscg',
                 'mrz_tidal_energy','mrz_tidal_filter',
                #  'x2',
                 'trp','wci','vcu','uni','rsl','old','pct','mal',
                 'emm2','srv','anc','jer','sal','ppt','rri2','bdt','lps',
                 'snc','dsj','bdl','nsl2','vol','tss','sss','cse',
                 'oh4','rsl','vcu','god',
                #  'frk','bac','hol','mtz','cll','tms','anh','gzl'
                 ]

    cases = list(eval(inputs.get('cases')))
    for case_num in cases:
        hist_dss_file = inputs.get('hist_dss_file').format(**locals())
        gate_dss_file = inputs.get('gate_dss_file').format(**locals())
        model_ec_file = inputs.get('model_ec_file').format(**locals())
        output_folder = inputs.get('output_folder').format(**locals())
        dcd_dss_file = inputs.get('dcd_dss_file').format(**locals())

        # generate aggregated ANN inputs from DSM2 outputs
        case_out_df = generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder,
                            b_part_e_part_dict=e_part_dict)
        case_out_df['case'] = case_num
        case_out_df['model'] = 'dsm2'
        case_out_df['scene'] = 'baseline'
        
        case_out_df = case_out_df[[col for col in col_order if col in case_out_df.columns]]
        
        case_out_df.to_csv(os.path.join(output_folder,f'dsm2_base_{case_num}.csv'), float_format="%.2f", index=True)


if __name__ == '__main__':

    from schimpy import schism_yaml
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "./input/ann_config_lathypcub_v4_dsm2.yaml"
    e_part_dict = None


    run_ann_input(in_fname, e_part_dict=e_part_dict)