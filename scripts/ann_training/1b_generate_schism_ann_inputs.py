#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to gather input/output data from SCHISM for training ANN


"""

import os
from schimpy.station import *
import datetime
import string
import re
        
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

def get_start_date_from_param(param_in):

    with open(param_in, 'r') as param:
        for line in param.readlines():
            if 'start_year' in line:
                sy =  int(re.findall(r'\b\d+\b', line)[0])
            elif 'start_month' in line:
                sm =  int(re.findall(r'\b\d+\b', line)[0])
            elif 'start_day' in line:
                sd =  int(re.findall(r'\b\d+\b', line)[0])
                
    start_date = datetime.datetime(sy, sm, sd)
    
    return start_date

#### class definition ===============================================================================
# TODO make the ANNBCECGen object able to take in DSM2 input yamls and output DSM2 inputs/outputs for ANN
class ANNBCECGen(object):

    def __init__(self, yml_fname, model_type):
        # read input yaml ---------------------------------------------------------------------------
        self.yml_fname = yml_fname
        with open(yml_fname, 'r') as f:
            self.inputs = schism_yaml.load(f)
        
        # assign env vars to format strings with ----------------------------------------------------
        self.env_vars = build_dict(self.inputs.get('env_vars'))
        for env_var in self.env_vars:
            self.env_vars[env_var] = string.Formatter().vformat(self.env_vars[env_var],(),SafeDict((self.env_vars)))

        # define in/out vars
        self.in_vars = build_dict(self.inputs.get('in_vars'))
        self.out_vars = build_dict(self.inputs.get('out_vars'))

        # define header values
        with open(self.env_vars['th_header'], 'r') as thh:
            head = thh.readline()
            headers = head.split(' ')
        
        if 'meshes' in self.inputs.keys():
            meshes = self.inputs.get('meshes')
        else:
            raise ValueError("Need to define the meshes using key 'meshes' in the yaml file")
    
        if 'cases' in self.inputs.keys():
            cases = self.inputs.get('cases')
            for case in cases:
                case_num = case.get('case_num')
                for mesh in meshes:
                    case_dir = string.Formatter().vformat(self.env_vars['case_dir'],(),
                                                          SafeDict(({**self.env_vars, 
                                                                     **locals()})))
                    # EC outputs
                    station_fpath = string.Formatter().vformat(self.out_vars['station_in'],(),
                                                               SafeDict(({**self.env_vars, 
                                                                          **locals()})))
                    outputs_fpath = string.Formatter().vformat(self.out_vars['station_output'],(),
                                                               SafeDict(({**self.env_vars, 
                                                                          **locals()})))
                    param_fpath = string.Formatter().vformat(self.out_vars['param_clinic'],(),
                                                               SafeDict(({**self.env_vars, 
                                                                          **locals()})))
                    time_basis = get_start_date_from_param(param_fpath)
                    # all_ts = read_staout(outputs_fpath, station_fpath, time_basis)

                    # for ec_loc in self.out_vars['ec_locs']:
                    #     print('hi')
                        # look for upper or default (no lower)
                    
                    # Inputs

                    for invar in self.inpuin_vars:
                        in_name = invar['name']
                        if 'th_file' in invar.keys():
                            th_file = invar['th_file'].format({**self.env_vars, 
                                                               **locals()})
                            for inp in invar['inputs']:
                                if inp[0] in ['-']:
                                    # negative flow!
                                    mult = -1
                                else:
                                    mult = 1
                                print('hi')

                # hist_dss_file = inputs.get('hist_dss_file').format(**locals())
                # gate_dss_file = inputs.get('gate_dss_file').format(**locals())
                # model_ec_file = inputs.get('model_ec_file').format(**locals())
                # output_folder = inputs.get('output_folder').format(**locals())
                # xlsx_filepath = inputs.get('xlsx_filepath').format(**locals())
                # dcd_dss_file = inputs.get('dcd_dss_file').format(**locals())

                # # generate aggregated ANN inputs from DSM2 outputs
                # generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder,
                #                     b_part_e_part_dict=e_part_dict)

                # # combine ANN input csv files into xlsx file (TODO: make this unecessary?)
                # csv_to_ann_xlsx(output_folder, xlsx_filepath)
        else:

            raise ValueError("Need to define the casess using key 'cases' in the yaml file")

        # else:
            # hist_dss_file = inputs.get('hist_dss_file').format(**locals())
            # gate_dss_file = inputs.get('gate_dss_file').format(**locals())
            # model_ec_file = inputs.get('model_ec_file').format(**locals())
            # output_folder = inputs.get('output_folder').format(**locals())
            # xlsx_filepath = inputs.get('xlsx_filepath').format(**locals())
            # dcd_dss_file = inputs.get('dcd_dss_file').format(**locals())
            
        # # generate aggregated ANN inputs from DSM2 outputs
        # generate_ann_inputs(hist_dss_file, gate_dss_file, dcd_dss_file, smcd_dss_file, model_ec_file, output_folder,
        #                     b_part_e_part_dict=e_part_dict)

        # # combine ANN input csv files into xlsx file (TODO: make this unecessary?)
        # csv_to_ann_xlsx(output_folder, xlsx_filepath)

if __name__ == '__main__':

    from schimpy import schism_yaml
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "./input/ann_config_lathypcub_v3_schism.yaml"
    ANNBCECGen(in_fname, model_type="SCHISM")