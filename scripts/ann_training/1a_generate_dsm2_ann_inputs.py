import pandas as pd
from vtools.functions.filter import cosine_lanczos
import pyhecdss
import os

from pydelmod.create_ann_inputs import get_dss_data

from bdschism.x2_time_series import find_x2

from vtools.functions.unit_conversions import ec_psu_25c

import re
import time


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
    
    out_df['exports'] = -df_exports_flow['exports']

    #############################################
    # 4. DCC gate operation as daily percentage #
    #############################################
    b_part = 'RSAC128'
    c_part = 'POS'
    dcc_op_df = process_gate_data(gate_dss_file, b_part, c_part, map_zero_one=['df<2','df==2'])
    dcc_op_df.index = dcc_op_df.index.to_period()
    out_df['dcc'] = dcc_op_df['gate_pos']
    
    ################################################
    # 5. Suisun gate operation as daily percentage #
    ################################################
    b_part = 'MTZSL'
    c_part = 'RADIAL_OP'
    suisun_op_df = process_gate_data(gate_dss_file, b_part, c_part, 
                                     map_zero_one = ['df!=-10','df==-10'], startDateStr='01JAN1953', endDateStr='01JAN2030')
    suisun_op_df.index = suisun_op_df.index.to_period()
    out_df['smscg'] = suisun_op_df['gate_pos']

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

    out_df['cu_flow'] = cu_total['cu_total']

    # Net Delta Outflow RMA:  (Sac, SJR, Yolo Bypass, Mokelumne, Cosumnes, Calaveras), exports (SWP, CVP, Contra Costa, NBA) and net DICU
    out_df['ndo'] = out_df['northern_flow'] + out_df['sjr_flow'] + out_df['exports'] - out_df['cu_flow']

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
    ec_locs = ['anc','anh','bac','bdl','bdt','bet','cll','cse',
                 'dsj','emm2','frk','god','gys','gzl','hll','hol2',
                 'ibs','jer','mal','mtz','nsl2','obi','oh4','old',
                 'pct','ppt','rri2','rsl','sal','snc','srv','sss',
                 'tms','trp','tss','uni','vcu','vol','wci','ver']
    
    
    b_part_dss_filename_dict = {name.upper():model_ec_file for name in ec_locs}
    b_part_c_part_dict = {name.upper():'EC' for name in ec_locs}
    b_part_e_part_dict = {name.upper():'1HOUR' for name in ec_locs}

    df_model_ec = get_dss_data(b_part_dss_filename_dict, 'b_part', \
        primary_part_c_part_dict=b_part_c_part_dict, primary_part_e_part_dict=b_part_e_part_dict)
    df_model_ec.columns = df_model_ec.columns.str.lower()
    
    col_rename_dict = {'ver': 'vern_ec'}

    df_model_ec.rename(columns=col_rename_dict, inplace=True)
    
    out_df = pd.merge(out_df, df_model_ec, left_index=True, right_index=True, how='outer')

    return out_df

def calc_x2(dss_file, names):

    b_part_dss_filename_dict = {str(name):dss_file for name in names}
    b_part_c_part_dict = {str(name):'EC' for name in names}
    b_part_e_part_dict = {str(name):'1HOUR' for name in names}

    ec_df = get_dss_data(b_part_dss_filename_dict, 'b_part', \
        primary_part_c_part_dict=b_part_c_part_dict, primary_part_e_part_dict=b_part_e_part_dict)
    
    for col_ec in ec_df.columns:
        ec_df[col_ec] = ec_psu_25c(ec_df.loc[:, col_ec])

    x2_df = pd.DataFrame(index=ec_df.index, columns=['x2'])

    for index, row in ec_df.iterrows():
        row.index = row.index.astype(float)
        x2_df.loc[index,'x2'] = find_x2(row)

    return x2_df
    

def run_ann_input(in_fname, case_nums=range(0,9999)):

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
                 'sac_flow','sjr_flow','exports','cu_flow','ndo',
                 'dcc','smscg',
                 'vern_ec','mrz_tidal_energy','mrz_tidal_filter',
                 'anc','anh','bac','bdl','bdt','bet','cll','cse',
                 'dsj','emm2','frk','god','gys','gzl','hll','hol2',
                 'ibs','jer','mal','mtz','nsl2','obi','oh4','old',
                 'pct','ppt','rri2','rsl','sal','snc','srv','sss',
                 'tms','trp','tss','uni','vcu','vol','wci',
                 'x2']
    
    float_columns = ['sac_flow','sjr_flow','exports','cu_flow','ndo',
                     'vern_ec','mrz_tidal_energy','mrz_tidal_filter',
                     'anc','anh','bac','bdl','bdt','bet','cll','cse',
                     'dsj','emm2','frk','god','gys','gzl','hll','hol2',
                     'ibs','jer','mal','mtz','nsl2','obi','oh4','old',
                     'pct','ppt','rri2','rsl','sal','snc','srv','sss',
                     'tms','trp','tss','uni','vcu','vol','wci',
                     'x2']

    x2_csv_infile = inputs.get('x2_csv_infile').format(**locals())
    x2_names = pd.read_csv(x2_csv_infile, comment='#')
    x2_names = x2_names['distance'].tolist()[::10] # only use every 10 columns out of the ~1000

    case_setup = pd.read_csv(inputs.get('case_setup'))

    for index, row in case_setup.iterrows():
        case_num = re.search(r'(\d+)$', row['case']).group(1)
        mod_case_num = False
        if int(case_num) in case_nums:
            if int(case_num) > 1000:
                mod_case_num = True
                case_num = int(case_num) - 1000

            hist_dss_file = inputs.get('hist_dss_file').format(**locals())
            gate_dss_file = inputs.get('gate_dss_file').format(**locals())
            model_ec_file = inputs.get('model_ec_file').format(**locals())
            model_x2_ec_file = inputs.get('model_x2_ec_file').format(**locals())
            output_folder = inputs.get('output_folder').format(**locals())
            casanntra_folder = inputs.get('casanntra_folder').format(**locals())
            dcd_dss_file = inputs.get('dcd_dss_file').format(**locals())

            # # to run this in parallel with ongoing/overnight check if the model is finished running
            # while not os.path.exists(model_x2_ec_file):
            #     print(f"Waiting for file {model_x2_ec_file} to appear...")
            #     time.sleep(120)  # Wait for 30 seconds before checking again

            # print(f"File {model_x2_ec_file} is now available! Post-processing model results")

            # calculate x2
            case_x2 = calc_x2(model_x2_ec_file, x2_names) # takes a while

            # generate aggregated ANN inputs from DSM2 outputs
            case_out_df = generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder)
            
            if mod_case_num:
                case_num = case_num + 1000
            case_out_df['case'] = case_num
            case_out_df['model'] = 'dsm2'
            case_out_df['scene'] = 'base'
            case_out_df = case_out_df.loc[pd.to_datetime(row['start']):(pd.to_datetime(row['end'])-pd.Timedelta(days=1)),:]

            case_out_df['x2'] = case_x2.loc[case_out_df.index, 'x2']
            
            case_out_df = case_out_df[[col for col in col_order if col in case_out_df.columns]]
            case_out_df[float_columns] = case_out_df[float_columns].astype(float)
            
            case_out_df.to_csv(os.path.join(casanntra_folder,f'dsm2_base_{case_num}.csv'), float_format="%.2f", index=True, index_label='datetime')


if __name__ == '__main__':

    from schimpy import schism_yaml
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "./input/ann_config_lathypcub_v3_dsm2.yaml"

    # run_ann_input(in_fname, case_nums=range(1001,1008))
    run_ann_input(in_fname, case_nums=range(1007,1008))