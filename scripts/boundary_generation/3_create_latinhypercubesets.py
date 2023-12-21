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

def apply_method(pert_dict):
    method = pert_dict.get('method')
    if method == 'read':
        print('\tread it')
    elif method == 'shift':
        print('\tshift it')
    elif method == 'set':
        print('\tset it')
    elif method == 'scale':
        print('\tscale it')
    else:
        raise ValueError(f'Perturbation method "{method}" is not defined in the code at the moment.')

# read in yaml
def create_cases(in_fname):
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
        output_dir = process_output_dir(inputs)

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
        case_items = inputs.get('cases')
        keys_case_section = ['name']

        # retrieve and write out cases
        for case in case_items:
            cname = case.get('name')
            print(cname)
            crange = [dt.datetime.strftime(case.get('case_start'), format='%Y-%m-%d'),
                      dt.datetime.strftime(case.get('case_end'), format='%Y-%m-%d')]

            # get perturbations
            
            cperts = case.get('perturbations')
            if cperts is not None:
                for cp in cperts:
                    try:
                        pdict = perturbs[cp]
                    except:
                        raise ValueError(f"The perturbation {cp} needs to be defined in the yaml file")

                    if 'components' in pdict.keys():
                        cdict = pdict.get('components')
                        for comp in cdict:
                            print('\tTODO: start here and apply method per sub-component of perturbation')
                    else:
                        apply_method(pdict)
                        
                    print('\t\twrite out file(s) to folder')

    

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # model_dir = r"D:\projects\delta_salinity\model\schism\dsp_202311_baseline"
    in_fname = "../../data/lathypcub_v1_setup.yaml"
    # in_fname = "../../../../model/schism/dsp_202311_baseline/dsp_baseline_bay_delta.yaml"

    # args = Namespace(main_inputfile=in_fname)

    create_cases(in_fname)