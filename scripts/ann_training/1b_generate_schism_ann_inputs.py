#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to gather input/output data from SCHISM for training ANN


"""

import os
from schimpy.station import *
import datetime

# {'CHDMC006': 'EC', 'CHSWP003': 'EC',\
#         'CHVCT000': 'EC', 'OLD_MID': 'EC', 'ROLD024': 'EC',
#         'ROLD059': 'EC', 'RSAC064': 'EC', 'RSAC075': 'EC',
#         'RSAC081': 'EC', 'RSAC092': 'EC', 'RSAC101': 'EC',
#         'RSAN007': 'EC', 'RSAN018': 'EC', 'RSAN032': 'EC',
#         'RSAN037': 'EC', 'RSAN058': 'EC', 'RSAN072': 'EC',
#         'RSMKL008': 'EC', 'SLCBN002': 'EC', 'SLDUT007': 'EC',
#         'SLMZU011': 'EC', 'SLMZU025': 'EC', 'SLSUS012': 'EC',
#         'SLTRM004': 'EC', 'SSS': 'EC', 'RSAC054': 'EC'}

        
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
            self.env_vars[env_var] = self.env_vars[env_var].format(**self.env_vars)
    
        if 'cases' in self.inputs.keys():
            cases = self.inputs.get('cases')
            for case in cases:
                case_num = case.get('case_num')
                outputs_fpath = "../outputs/staout_1" # TODO: fix
                station_fpath = '../station.in' # TODO: fix
                time_basis = datetime.datetime(2009, 2, 10) # TODO: fix
                read_staout(outputs_fpath, station_fpath, time_basis)
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
    run_ann_input(in_fname)