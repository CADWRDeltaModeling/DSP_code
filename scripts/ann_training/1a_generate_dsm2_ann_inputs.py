import pandas as pd
from vtools.functions.filter import cosine_lanczos
import pyhecdss
import os

from pydelmod.create_ann_inputs import get_dss_data


def process_gate_data(dss_filename, output_file, b_part, c_part, map_zero_one=['df==0','df==1'], startDateStr=None, endDateStr=None):
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
    df_daily_avg.to_csv(output_file)

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
    df_northern_flow.to_csv(output_folder + '/df_northern_flow.csv')

    #############
    # SJR Flow  #
    #############
    b_part_dss_filename_dict = {'RSAN112': hist_dss_file}
    b_part_c_part_dict = {'RSAN112': 'FLOW'}
    df_sjr_flow = get_dss_data(b_part_dss_filename_dict, 'b_part', b_part_c_part_dict)
    df_sjr_flow.to_csv(output_folder + '/df_sjr_flow.csv')

    ###############################################################################################
    # 3. Exports: Sum(Banks, Jones, CCC plants(Rock Sl, Middle R (actually old river), Victoria)) #
    ###############################################################################################
    b_part_dss_filename_dict = {'CHSWP003': hist_dss_file, 'CHDMC004': hist_dss_file, 'CHCCC006': hist_dss_file,\
        'ROLD034': hist_dss_file, 'CHVCT001': hist_dss_file}
    df_exports_flow = get_dss_data(b_part_dss_filename_dict, 'b_part')
    # df_exports_flow.fillna(0, inplace=True)
    df_exports_flow['exports'] = df_exports_flow['CHSWP003']+df_exports_flow['CHDMC004']+df_exports_flow['CHCCC006']+\
        df_exports_flow['ROLD034']+df_exports_flow['CHVCT001']
    df_exports_flow.to_csv(output_folder+'/df_exports_flow.csv')

    #############################################
    # 4. DCC gate operation as daily percentage #
    #############################################
    b_part = 'RSAC128'
    c_part = 'POS'
    gate_output_file = output_folder+'/dcc_gate_op.csv'
    process_gate_data(gate_dss_file, gate_output_file, b_part, c_part, map_zero_one=['df<2','df==2'])
    
    ################################################
    # 5. Suisun gate operation as daily percentage #
    ################################################
    b_part = 'MTZSL'
    c_part = 'RADIAL_OP'
    gate_output_file = output_folder+'/suisun_gate_op.csv'
    process_gate_data(gate_dss_file, gate_output_file, b_part, c_part, map_zero_one = ['df!=-10','df==-10'], startDateStr='01JAN1953', endDateStr='01JAN2030')

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
    # now only save the grand total column to csv
    cu_total[['cu_total']].to_csv(output_folder+'/df_cu_total.csv')

    ########################################
    # 7. Tidal Energy: daily max-daily min #
    ########################################
    b_part_dss_filename_dict={'RSAC054': hist_dss_file}
    b_part_c_part_dict={'RSAC054': 'STAGE'}
    df_mtz_stage = get_dss_data(b_part_dss_filename_dict, 'b_part', \
        primary_part_c_part_dict=b_part_c_part_dict, daily_avg=False)
    df_nrg = cosine_lanczos((df_mtz_stage-cosine_lanczos(df_mtz_stage.copy(), 
                                                         cutoff_period ='40H', padtype='odd'))**2, 
                            cutoff_period ='40H', padtype='odd') # = < (z- <z>)^2 >
    df_mtz_tidal_energy = df_nrg.resample('D', closed='right').mean()
    df_mtz_tidal_energy.columns=['tidal_energy']

    df_mtz_stage.to_csv(output_folder+'/df_mtz_stage.csv')    
    df_mtz_tidal_energy.to_csv(output_folder+'/df_mtz_tidal_energy.csv')

    #############################################
    # 8. SJR inflow salinity at vernalis, daily #
    #############################################
    b_part_dss_filename_dict = {'RSAN112': hist_dss_file}
    b_part_c_part_dict = {'RSAN112': 'EC'}
    df_sjr_ec = get_dss_data(b_part_dss_filename_dict, 'b_part', primary_part_c_part_dict=b_part_c_part_dict)
    df_sjr_ec.to_csv(output_folder + '/df_sjr_ec.csv')

    ##########################
    # 9. Sacramento River EC #
    ##########################
    b_part_dss_filename_dict = {'RSAC139': hist_dss_file}
    b_part_c_part_dict = {'RSAC139': 'EC'}
    df_sac_ec = get_dss_data(b_part_dss_filename_dict, 'b_part', primary_part_c_part_dict=b_part_c_part_dict)
    df_sac_ec.to_csv(output_folder + '/df_sac_ec.csv')

    ######################################
    # 10. EC Output for various locations #
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
    # df_model_ec = df_model_ec.resample('D').mean()

    # now add duplicate columns
    duplication_dict = {'RSAN007': 'Antioch_dup', 'CHSWP003': 'CCFB_Intake_dup', 'RSAC081': 'Collinsville_dup', 'CHDMC006': 'CVP_Intake_dup',
        'RSAC092': 'Emmaton_dup', 'RSAN018': 'Jersey_Point_dup', 'RSAC075': 'Mallard_Island_dup'}
    for rki in duplication_dict:
        new_name = duplication_dict[rki]
        df_model_ec[new_name] = df_model_ec[rki]

    # print('before error: columns='+str(df_model_ec.columns))
    # now add model output ec near CCC intakes
    b_part_dss_filename_dict = {'ROLD034': model_ec_file, 'SLRCK005': model_ec_file}
    b_part_c_part_dict = {'ROLD034': 'EC', 'SLRCK005': 'EC'}
    df_model_ec_2 = get_dss_data(b_part_dss_filename_dict, 'b_part', b_part_c_part_dict)
    df_model_ec = pd.merge(df_model_ec, df_model_ec_2, how='outer', left_index=True, right_index=True)
    # print('before error: columns='+str(df_model_ec.columns))

    # now add a copy of victoria intake ec
    df_model_ec['CHVCT000_dup'] = df_model_ec['CHVCT000']
    # now add another copy of Mtz ec
    df_model_ec['Martinez_input'] = df_model_ec['RSAC054']
        
    # now rename some of the columns
    col_rename_dict = {'CHDMC006': 'CHDMC006-CVP INTAKE', 'CHSWP003': 'CHSWP003-CCFB_INTAKE', 'CHVCT000': 'CHVCT000-VICTORIA INTAKE',
        'OLD_MID': 'OLD_MID-OLD RIVER NEAR MIDDLE RIVER', 'ROLD024': 'ROLD024-OLD RIVER AT BACON ISLAND', 
        'ROLD059': 'ROLD059-OLD RIVER AT TRACY BLVD', 'RSAC064': 'RSAC064-SACRAMENTO R AT PORT CHICAGO', 'RSAC075': 'RSAC075-MALLARDISLAND',
        'RSAC081': 'RSAC081-COLLINSVILLE', 'RSAC092': 'RSAC092-EMMATON', 'RSAC101': 'RSAC101-SACRAMENTO R AT RIO VISTA', 
        'RSAN007': 'RSAN007-ANTIOCH', 'RSAN018': 'RSAN018-JERSEYPOINT', 'RSAN032': 'RSAN032-SACRAMENTO R AT SAN ANDREAS LANDING',
        'RSAN037': 'RSAN037-SAN JOAQUIN R AT PRISONERS POINT', 'RSAN058': 'RSAN058-ROUGH AND READY ISLAND', 
        'RSAN072': 'RSAN072-SAN JOAQUIN R AT BRANDT BRIDGE', 'RSMKL008': 'RSMKL008-S FORK MOKELUMNE AT TERMINOUS',
        'SLCBN002': 'SLCBN002-CHADBOURNE SLOUGH NR SUNRISE DUCK CLUB', 'SLDUT007': 'SLDUT007-DUTCH SLOUGH', 
        'SLMZU011': 'SLMZU011-MONTEZUMA SL AT BELDONS LANDING', 'SLMZU025': 'SLMZU025-MONTEZUMA SL AT NATIONAL STEEL',
        'SLSUS012': 'SLSUS012-SUISUN SL NEAR VOLANTI SL', 'SLTRM004': 'SLTRM004-THREE MILE SLOUGH NR SAN JOAQUIN R', 'SSS': 'SSS-STEAMBOAT SL',
        'ROLD034': 'Old_River_Hwy_4', 'SLRCK005': 'CCWD_Rock', 'CHVCT000_dup': 'CCWD_Victoria_dup',
        'RSAC054': 'Martinez_input_dup'}

    df_model_ec.rename(columns=col_rename_dict, inplace=True)
    df_model_ec.to_csv(output_folder + '/df_model_ec.csv')

def csv_to_ann_xlsx(csv_dir, xlsx_filepath):
    csv_to_xlsx_dict = { 'base_ec_output':['df_model_ec.csv',None,None],
                        'sac_ec':['df_sac_ec.csv','RSAC139','sac_greens_ec'],
                        'sjr_vernalis_ec':['df_sjr_ec.csv','RSAN112','sjr_vernalis_ec'],
                        'mtz_tidal_nrg':['df_mtz_tidal_energy.csv','tidal_energy','daily_nrg'],
                        'net_delta_cu': ['df_cu_total.csv','cu_total','div+seep-drain_dcd+smcd'],
                        'dxc_gate_fraction':['dcc_gate_op.csv','gate_pos',
                                             'gate_pos'],
                        'suisun_gate_fraction':['suisun_gate_op.csv','gate_pos',
                                             's_gate_pos'],
                        'exports':['df_exports_flow.csv','exports','exports'],
                        'sjr_flow':['df_sjr_flow.csv','RSAN112','sjr_flow'],
                        'northern_flow':['df_northern_flow.csv', 'northern_flow', 'northern_flow']
    }

    ordered_sheets = ['northern_flow','sjr_flow','exports','dxc_gate_fraction','suisun_gate_fraction','net_delta_cu',
                      'mtz_tidal_nrg','sjr_vernalis_ec','sac_ec','base_ec_output']
    
    with pd.ExcelWriter(xlsx_filepath) as writer:
        for sheet in ordered_sheets:
            df = pd.read_csv(os.path.join(csv_dir,csv_to_xlsx_dict[sheet][0]))
            if csv_to_xlsx_dict[sheet][1] is not None:
                df = df.loc[:,['Unnamed: 0', csv_to_xlsx_dict[sheet][1]]]
                df = df.set_axis(['Time',csv_to_xlsx_dict[sheet][2]], axis=1)
                df.to_excel(writer, sheet_name=sheet, index=False)
            else:
                df.to_excel(writer, sheet_name=sheet, index=False)

def run_ann_input(in_fname, e_part_dict=None):

    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
    
    dsp_home = inputs.get('dsp_home')
    experiment = inputs.get('experiment')

    base_study_folder = inputs.get('base_study_folder').format(**locals())
    model_input_folder = inputs.get('model_input_folder').format(**locals())
    model_folder = inputs.get('model_folder').format(**locals())
    model_output_folder = inputs.get('model_output_folder').format(**locals())
    smcd_dss_file = inputs.get('smcd_dss_file').format(**locals())
    
    if 'cases' in inputs.keys():
        cases = inputs.get('cases')
        for case in cases:
            case_num = case.get('case_num')
            hist_dss_file = inputs.get('hist_dss_file').format(**locals())
            gate_dss_file = inputs.get('gate_dss_file').format(**locals())
            model_ec_file = inputs.get('model_ec_file').format(**locals())
            output_folder = inputs.get('output_folder').format(**locals())
            xlsx_filepath = inputs.get('xlsx_filepath').format(**locals())
            dcd_dss_file = inputs.get('dcd_dss_file').format(**locals())

            # generate aggregated ANN inputs from DSM2 outputs
            generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder,
                                b_part_e_part_dict=e_part_dict)

            # combine ANN input csv files into xlsx file (TODO: make this unecessary?)
            csv_to_ann_xlsx(output_folder, xlsx_filepath)

    else:
        hist_dss_file = inputs.get('hist_dss_file').format(**locals())
        gate_dss_file = inputs.get('gate_dss_file').format(**locals())
        model_ec_file = inputs.get('model_ec_file').format(**locals())
        output_folder = inputs.get('output_folder').format(**locals())
        xlsx_filepath = inputs.get('xlsx_filepath').format(**locals())
        dcd_dss_file = inputs.get('dcd_dss_file').format(**locals())
        
    # generate aggregated ANN inputs from DSM2 outputs
    generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder,
                        b_part_e_part_dict=e_part_dict)

    # combine ANN input csv files into xlsx file (TODO: make this unecessary?)
    csv_to_ann_xlsx(output_folder, xlsx_filepath)

if __name__ == '__main__':

    from schimpy import schism_yaml
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "./input/ann_config_lathypcub_v3_dsm2.yaml"
    # e_part_dict={'CHDMC006': '1DAY', 'CHSWP003': '1DAY', 'CHVCT000': '1DAY', 'OLD_MID': '1DAY', 
    #              'ROLD024': '1DAY', 'ROLD059': '1DAY', 'RSAC064': '1DAY', 'RSAC075': '1DAY', 
    #              'RSAC081': '1DAY', 'RSAC092': '1DAY', 'RSAC101': '1DAY', 'RSAN007': '1DAY', 
    #              'RSAN018': '1DAY', 'RSAN032': '1DAY', 'RSAN037': '1DAY', 'RSAN058':'1DAY',
    #              'RSAN072':'1DAY','RSMKL008':'1DAY','SLCBN002':'1DAY','SLDUT007':'1DAY',
    #              'SLMZU011':'1DAY', 'SLMZU025': '1DAY','SLSUS012': '1DAY', 'SLTRM004': '1DAY', 
    #              'SSS': '1DAY', 'RSAC054': '1HOUR'} # LAT: I don't remember why I did this...
    e_part_dict = None


    run_ann_input(in_fname, e_part_dict=e_part_dict)