# script by Lily Tomkovic to create DSM2 boundary data from the the meta-latinhypercube cases' datasets
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces boundary inputs for DSM2 within existing folders

import pandas as pd
import numpy as np
from schimpy import schism_yaml
from schimpy.prepare_schism import item_exist

from pydelmod.create_ann_inputs import get_dss_data
import pyhecdss

import config_python as config
from mrzecest.ec_boundary import ec_est
import string

import datetime as dt
import os


class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'


def build_dict(node):
    od = {}
    if isinstance(node, dict):
        for key in node.keys():
            if isinstance(node[key], dict):
                od[key] = build_dict(node[key])
            else:
                od[key] = node[key]
    elif isinstance(node, list):
        for d in node:
            od.update(build_dict(d))
    return od


def calc_dcu(dcd_filename, write_outfile=None):
    """calculate Delta Consumptive Use from DICULFILE inputs from DSM2 setup
    
    Parameters
    ----------
    dcd_filename: str | Path 
    
        name of boundary flow DSM2 input file in DSS format

    write_outfile: str | Path | None
    
        name of output file to write Net Delta Ouflow (df_ndo) out to

    Returns
    -------
    cu_total_dcd : pd.DataFrame
        Data with calculation of Delta Consumptive Use from boundary inputs
    
    """
    
    div_seep_dcd_c_part_dss_filename_dict = {'DIV-FLOW': dcd_filename, 'SEEP-FLOW': dcd_filename}
    drain_dcd_c_part_dss_filename_dict = {'DRAIN-FLOW': dcd_filename}

    df_div_seep_dcd = get_dss_data(div_seep_dcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True)
    df_drain_dcd = get_dss_data(drain_dcd_c_part_dss_filename_dict, 'c_part', filter_b_part_numeric=True)

    df_div_seep_dcd['dcd_divseep_total']=df_div_seep_dcd[df_div_seep_dcd.columns].sum(axis=1)
    df_drain_dcd['dcd_drain_total']=df_drain_dcd[df_drain_dcd.columns].sum(axis=1)

    cu_total_dcd = pd.merge(df_div_seep_dcd, df_drain_dcd, how='left', left_index=True, right_index=True)

    cu_total_dcd['cu_total'] = cu_total_dcd['dcd_divseep_total'] - cu_total_dcd['dcd_drain_total']

    cu_total_dcd = cu_total_dcd[['cu_total']]

    if write_outfile:
        cu_total_dcd.to_csv(write_outfile, index=True, float_format="%.2f")

    return cu_total_dcd


def calc_ndo(flow_filename, dcd_filename, write_outfile=None):
    """calculate Net Delta Outflow using the BNDRYINPUT and DICULFILE inputs from DSM2 setup
    
    Parameters
    ----------
    flow_filename: str | Path 
    
        name of boundary flow DSM2 input file in DSS format

    dcd_filename: str | Path 
    
        name of boundary flow DSM2 input file in DSS format

    write_outfile: str | Path | None
    
        name of output file to write Net Delta Ouflow (df_ndo) out to

    Returns
    -------
    df_ndo : pd.DataFrame
        Data with calculation of Net Delta Outflow from boundary inputs
    
    """

    cu_total_dcd = calc_dcu(dcd_filename)

    b_part_dss_filename_dict = {'RSAC155': flow_filename, 
                                'RSAN112': flow_filename,
                                'BYOLO040': flow_filename, 
                                'RMKL070': flow_filename, 
                                'RCSM075': flow_filename, 
                                'RCAL009': flow_filename, 
                                'SLBAR002': flow_filename, 
                                'CHSWP003': flow_filename, 
                                'CHDMC004': flow_filename, 
                                'CHVCT001': flow_filename, 
                                'ROLD034': flow_filename, 
                                'CHCCC006': flow_filename}
    df_ndo = get_dss_data(b_part_dss_filename_dict, 'b_part')
    df_ndo = pd.merge(df_ndo, cu_total_dcd, how='left', left_index=True, right_index=True)

    positive_flows = ['RSAC155', 'BYOLO040', 'RMKL070', 'RCSM075', 'RCAL009']
    negative_flows = ['SLBAR002', 'CHSWP003', 'CHDMC004', 'CHVCT001', 'ROLD034', 'CHCCC006', 'cu_total']

    df_ndo['ndo'] = df_ndo[positive_flows].sum(axis=1) - df_ndo[negative_flows].sum(axis=1)
    df_ndo = df_ndo[['ndo']]
    
    if write_outfile:
        df_ndo.to_csv(write_outfile, index=True, float_format="%.2f")

    return df_ndo

def create_perturbations(row, pert_vars):
    perturbations = [f"{var} {row[var].lower()}" for var in pert_vars if 'regular' not in row[var].lower()]
    return perturbations


def create_dsm2_ec_est(in_fname, dsm2_config_fname, write_dss=True, skip=0):
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
        output_dir = inputs['output_dir']

    if item_exist(inputs, 'cases'):
        print(f'Handling cases:')
        case_items = inputs.get('cases')        
        if 'filename' in case_items[0].keys():
            
            # setup cases with csv file and not yaml
            case_file = pd.read_csv(case_items[0]['filename'])
            pert_vars = case_items[0]['pert_vars']
            case_file['start'] = pd.to_datetime(case_file['start'])
            case_file['end'] = pd.to_datetime(case_file['end'])

            case_file['name'] = case_file['case']
            case_file.rename(columns={'start': 'case_start',
                                      'end':'case_end'}, inplace=True)
            case_file['perturbations'] = case_file.apply(create_perturbations, axis=1, args=(pert_vars,))
        
            case_items = case_file.iloc[skip:].copy() # this behaves the same as if the cases are defined in a yaml

    with open(dsm2_config_fname, 'r') as f:
        dsm2_inputs = schism_yaml.load(f)
    dsm2_config = dsm2_inputs.get('model_config')
    config_infile = dsm2_inputs.get('config_file')

    # create output dir
    out_dss_dir = dsm2_inputs.get('out_dss_dir')
    training_set = dsm2_inputs.get('training_set')
    inp_dir = os.path.abspath(os.path.join(dsm2_inputs.get('out_dss_dir'), training_set))
    dsm2_dir = os.path.abspath(os.path.join(out_dss_dir, f'{training_set}/timeseries/'))

    os.chdir(os.path.join(out_dss_dir, training_set))

    # retrieve and write out cases
    print('Handling cases:')
    for index, case in case_items.iterrows():

        cname = case.get('name')
        print(f'\t- {cname}')
        f_part_out = f'DSP_{cname}'
        
        config_filename = os.path.basename(config_infile).replace('CASE',str(''.join(filter(str.isdigit, cname)))) # replace with the proper filename
        config_case_filename = os.path.join(inp_dir, config_filename)

        config.setConfigVars(config_case_filename) # load the configuration parameters

        start = pd.to_datetime(f'{config.getAttr("START_DATE")} 0000', format="%d%b%Y %H%M")
        end = pd.to_datetime(f'{config.getAttr("END_DATE")} 0000', format="%d%b%Y %H%M")

        print(f"\t\tCalculating boundary salinity for the period {start.strftime('%m-%d-%Y')} to {end.strftime('%m-%d-%Y')}")
        
        ec_est_filename = os.path.join(dsm2_dir,f"{cname}_ec_est.dss")
        stage_bdy_filename = os.path.normpath(os.path.join(os.getcwd(),config.getAttr('STAGE_SOURCE_FILE'))).replace("\\", "/")
        flow_filename = os.path.normpath(os.path.join(os.getcwd(),config.getAttr('BNDRYINPUT'))).replace("\\", "/")
        dcd_filename = os.path.normpath(os.path.join(os.getcwd(),config.getAttr('DICUFILE'))).replace("\\", "/")

        df_ndo = calc_ndo(flow_filename, dcd_filename)
        # df_ndo.index = pd.to_datetime(df_ndo.index)

        if not os.path.exists(ec_est_filename) and write_dss:
            dumdss = pyhecdss.DSSFile(ec_est_filename, create_new=True) # create the file if writing out to DSS

        b_part_dss_filename_dict={'RSAC054': stage_bdy_filename}
        b_part_c_part_dict={'RSAC054': 'STAGE'}
        df_mtz_stage = get_dss_data(b_part_dss_filename_dict, 'b_part', primary_part_c_part_dict=b_part_c_part_dict, daily_avg=False)

        # parameters from estimation
        log10beta = 10.217 # x[0] from ec_boundary_fit_gee.py printout
        npow = 0.461 # x[1] from ec_boundary_fit_gee.py printout
        area_coef = -6127433509.04 # x[2] from ec_boundary_fit_gee.py printout
        energy_coef = 1495.91 # x[3] from ec_boundary_fit_gee.py printout
        beta0 = 1.6828 # from const coef result 
        beta1 = -23.0735 * 1e-3 # from gnpow coef result 
        filter_k0 = 6 # from fitting_config.yaml
        filt_coefs = np.array([0.111, 0.896, -0.606, 0.678, -0.745, -0.416, -0.046, 
                            1.161, 0.321, -1.069, 0.515, -0.965, 0.576]) * 1e-3 # z{n} from output coefs
        filter_dt = pd.Timedelta('3h') # from fitting_config.yaml
        so = 20000. # hardwired in ec_boundary_fit_gee.py
        sb = 200. # hardwired in ec_boundary_fit_gee.py

        mrzecest = ec_est(df_ndo, df_mtz_stage, start, end,
                      area_coef,energy_coef,
                      log10beta,
                      beta0, beta1, npow, filter_k0,
                      filt_coefs, filter_dt,
                      so, sb)
        
        unit_part = "uS/CM"
        pathname_out = f"/FILL+CHAN/RSAC054/EC/{start.strftime('%d%b%Y')}/1HOUR/{f_part_out}/"
        ptype = 'INST-VAL'
        mrzecest.index = mrzecest.index.to_timestamp()
        
        # write to DSS file
        with pyhecdss.DSSFile(ec_est_filename) as d_out:
            # print(f'Writing out {pathname_out}')
            # need to shift by one day because of DSS writing timestamp issues
            # mrzecest.index = mrzecest.index + pd.Timedelta(days=1)
            if mrzecest.index.freq is None:
                mrzecest.index.freq = pd.infer_freq(mrzecest.index)
            # write regular output to DSS file
            d_out.write_rts(pathname_out, mrzecest, unit_part, ptype)

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # model_dir = r"D:\projects\delta_salinity\model\schism\dsp_202311_baseline"
    in_fname = "./input/lathypcub_v4_setup.yaml"

    # cases = create_cases()
    dsm2_config_fname = "./input/lathypcub_v4_dsm2_config.yaml"
    # in_fname = "../../../../model/schism/dsp_202311_baseline/dsp_baseline_bay_delta.yaml"

    # args = Namespace(main_inputfile=in_fname)

    create_dsm2_ec_est(in_fname, dsm2_config_fname, skip=105)