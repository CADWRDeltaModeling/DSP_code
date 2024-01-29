# script by Lily Tomkovic to create boundary datasets for the meta-latinhypercube cases
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces case folders with any modified boundary inputs

import pandas as pd
from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist
from schimpy.schism_setup import check_and_suggest
from argparse import Namespace
import datetime as dt
import os

def build_dict(node):
    od = {}
    if isinstance(node, dict):
        for key in node.keys():
            if isinstance(node[key], dict):
                od[key] = build_dict(node[key])
            else:
                od[key] = node[key]
    return od



def apply_method(pert_dict, ssvar=None):
    method = pert_dict['method']
    file_in = pert_dict['args']['file']
    if method == 'read_dcd':
        dat_in = pd.read_csv(file_in.format(**locals()), parse_dates=[0], index_col=[0])
        dat_out = dat_in.copy()
    else:
        dat_in = pd.read_csv(file_in, parse_dates=[0], index_col=[0], header=None)

    if method == 'read':
        dat_out = dat_in.copy()
    elif method == 'shift':
        dat_out = dat_in.copy()
        dat_out.index = dat_out.index - pd.Timedelta(days=pert_dict['args']['shift_forward'])
    elif method == 'set':
        dat_out = dat_in.copy()
        dat_out[1] = float(pert_dict['args']['set_value'])
    elif method == 'scale':
        dat_out = dat_in.copy()
        dat_out[1] = float(pert_dict['args']['scale_factor'])*dat_out[1]
    elif method == 'read_dcd':
        pass
    else:
        raise ValueError(f'Perturbation method "{method}" is not defined in the code at the moment.')
    return dat_out

# read in yaml
def create_cases(in_fname):
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
        output_dir = process_output_dir(inputs)

    # create output folder
    if not os.path.exists(output_dir):
        print(f'Creating output directory {output_dir}')
        os.mkdir(output_dir)

    keys_top_level = ["training_set", "output_dir", "perturbations", "cases"] \
        + schism_yaml.include_keywords
    # check_and_suggest(list(inputs.keys()), keys_top_level)

    # check_nested_match(inputs)

    if item_exist(inputs, 'training_set'):
        training_set = inputs['training_set']

    if item_exist(inputs, 'perturbations'):
        perturb_items = inputs.get('perturbations')
        keys_perturb_section = ['name']

        # update perturbations
        perturbs = {}
        print(f'Storing perturbations:')
        for perturb in perturb_items:
            pname = perturb.get('name')
            print(f'\t- {pname}')
            perturbs[pname] = build_dict({k: perturb[k] for k in set(list(perturb.keys())) - set(['name'])})

    if item_exist(inputs, 'cases'):
        print(f'Handling cases:')
        case_items = inputs.get('cases')
        keys_case_section = ['name']

        # retrieve and write out cases
        for case in case_items:

            cname = case.get('name')
            print(f'\t- {cname}')
            # create output folder
            case_dir = os.path.join(output_dir, cname)
            if not os.path.exists(case_dir):
                print(f'\t\t Creating output directory {case_dir}')
                os.mkdir(case_dir)

            # crange = [dt.datetime.strftime(case.get('case_start'), format='%Y-%m-%d 00:00'),
            #           dt.datetime.strftime(case.get('case_end')+dt.timedelta(days=1), format='%Y-%m-%d 00:00')]
            start_date = case.get('case_start')
            end_date = case.get('case_end')
            crange = [dt.datetime.strftime(start_date, format='%Y%m%d0000'),
                      dt.datetime.strftime(end_date + dt.timedelta(days=1), format='%Y%m%d0000')]

            # get perturbations
            
            cperts = case.get('perturbations')
            if cperts is not None:
                for cp in cperts:
                    try:
                        pdict = perturbs[cp]
                    except:
                        raise ValueError(f"The perturbation {cp} needs to be defined in the perturbations section of the yaml file")

                    if 'components' in pdict.keys():
                        cdict = pdict.get('components')

                        # create output folder
                        subout_dir = os.path.join(case_dir, cp)
                        if not os.path.exists(subout_dir):
                            print(f'\t\t\t Creating sub-output directory {subout_dir}')
                            os.mkdir(subout_dir)
                        
                        for comp in cdict:

                            # loop through and apply methods
                            dat_out = apply_method(comp)
                            dat_out = dat_out[dat_out.index.to_series().between(pd.to_datetime(start_date),pd.to_datetime(end_date))]
                            dat_out.to_csv(os.path.join(subout_dir, f'{comp["model_input"]}_{comp["method"]}_{crange[0]}-{crange[1]}.csv'),
                                           header=None)
                    elif pdict.get('method') == 'read_dcd':
                        for ssvar in ['source','sink']:
                            dat_out = apply_method(pdict, ssvar=ssvar)
                            dat_out = dat_out[dat_out.index.to_series().between(pd.to_datetime(start_date),pd.to_datetime(end_date))]
                            dat_out.to_csv(os.path.join(case_dir, f'{pdict["model_input"]}_{pdict["method"]}_{ssvar}_{crange[0]}-{crange[1]}.csv'),
                                            header=None)
                    else:
                        dat_out = apply_method(pdict)
                        dat_out = dat_out[dat_out.index.to_series().between(pd.to_datetime(start_date),pd.to_datetime(end_date))]
                        dat_out.to_csv(os.path.join(case_dir, f'{pdict["model_input"]}_{pdict["method"]}_{crange[0]}-{crange[1]}.csv'),
                                           header=None)

    

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # model_dir = r"D:\projects\delta_salinity\model\schism\dsp_202311_baseline"
    in_fname = "../../data/lathypcub_v2_setup.yaml"
    # in_fname = "../../../../model/schism/dsp_202311_baseline/dsp_baseline_bay_delta.yaml"

    # args = Namespace(main_inputfile=in_fname)

    create_cases(in_fname)