# script by Lily Tomkovic to create SCHISM boundary data from the the meta-latinhypercube cases' datasets
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces boundary inputs for SCHISM within existing folders

import pandas as pd
from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist

from pydelmod.create_ann_inputs import get_dss_data
import pyhecdss

import datetime as dt
import os

# functions

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

# Get pathnames function
def get_pathname(dss_filename, b_part, c_part, e_part=None, f_part=None, filter_b_part_numeric=None):
    with pyhecdss.DSSFile(dss_filename) as d:
        catdf = d.read_catalog()
        dss_file_parts = dss_filename.split('/')
        dfilename = dss_file_parts[len(dss_file_parts)-1]
        filtered_df = None
        if b_part is not None:
            filtered_df = filtered_df[(catdf.B==b_part)] if filtered_df is not None else catdf[(catdf.B==b_part)]
        if c_part is not None:
            filtered_df = filtered_df[(catdf.C==c_part)] if filtered_df is not None else catdf[(catdf.C==c_part)]
        if e_part is not None:
            filtered_df = filtered_df[(catdf.E==e_part)] if filtered_df is not None else catdf[(catdf.E==e_part)]
        if f_part is not None:
            filtered_df = filtered_df[(catdf.F==f_part)] if filtered_df is not None else catdf[(catdf.F==f_part)]
        if filter_b_part_numeric:
            filtered_df = filtered_df[(catdf.B.str.isnumeric())]
        path_list = d.get_pathnames(filtered_df)

    return path_list

#### class definition
# TODO make the ModelBCGen object able to take in DSM2 input yamls and output DSM2 BCs
class ModelBCGen(object):

    def __init__(self, yml_fname, model_type):
        # read input yaml ---------------------------------------------------------------------------
        self.yml_fname = yml_fname
        with open(yml_fname, 'r') as f:
            self.inputs = schism_yaml.load(f)
        
        # assign env vars to format strings with ----------------------------------------------------
        self.env_vars = build_dict(self.inputs.get('env_vars'))
        for env_var in self.env_vars:
            self.env_vars[env_var] = self.env_vars[env_var].format(**self.env_vars)

        # get case bundling data --------------------------------------------------------------------
        self.case_setup = self.inputs.get('case_setup').format(**self.env_vars)
        with open(self.case_setup, 'r') as f:
            self.case_inputs = schism_yaml.load(f)
        self.case_data_dir = self.case_inputs.get('output_dir')
        perturb_items = self.case_inputs.get('perturbations')

        ## update perturbations
        perturbs = {}
        print(f'Storing perturbations:')
        for perturb in perturb_items:
            pname = perturb.get('name')
            print(f'\t- {pname}')
            perturbs[pname] = build_dict({k: perturb[k] for k in set(list(perturb.keys())) - set(['name'])})
        self.perturbs = perturbs

        ## update cases
        self.case_items = self.case_inputs.get('cases')

        # retrieve and write out cases
        print('Handling cases:')
        for case in self.case_items:

            cname = case.get('name')
            print(f'\t- {cname}')
            # create output folder
            case_dir = os.path.join(self.case_data_dir, cname)

            crange = [dt.datetime.strftime(case.get('case_start'), format='%Y%m%d0000'),
                    dt.datetime.strftime(case.get('case_end')+dt.timedelta(days=1), format='%Y%m%d0000')]
            
            case_perts = case_dir['perturbations']
            
            if model_type.lower() == 'schism':

                self.setup_schism_case(case_dir, crange, case_perts)

    def set_schism_vars(self):
        self.param_tropic_base = self.inputs.get('param_tropic_base').format(**self.env_vars)
        self.param_clinic_base = self.inputs.get('param_clinic_base').format(**self.env_vars)
        self.flux_file_in = self.inputs.get('flux_file_in').format(**self.env_vars)
        self.bc_tropic_in = self.inputs.get('bc_tropic_in').format(**self.env_vars)
        self.bc_clinic_in = self.inputs.get('bc_clinic_in').format(**self.env_vars)
        self.tropic_ocean_bc = self.inputs.get('tropic_ocean_bc').format(**self.env_vars)      
        self.th_repo = self.inputs.get('th_repo').format(**self.env_vars) 
        self.dcd_daily = [dcd.format(**self.env_vars) for dcd in self.inputs.get('dcd_daily')]
        
    def setup_schism_case(self, case_dir, crange, case_perts):

        self.set_schism_vars()

        # get param.nmls


    

## Read in param.tropic and param.clinic and modify

yml_fname = "./input/schism_lathypcub_v2.yaml"

mbc = ModelBCGen(yml_fname, 'schism')
mbc.set_schism_vars()


print('hi')