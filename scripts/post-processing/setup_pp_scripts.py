# script by Lily Tomkovic to create SCHISM boundary data from the the meta-latinhypercube cases' datasets
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces boundary inputs for SCHISM within existing folders

import pandas as pd
from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist
from schimpy.model_time import file_to_elapsed

from pydelmod.create_ann_inputs import get_dss_data
from dms_datastore.read_ts import read_noaa, read_ts
import pyhecdss

import shutil
import datetime as dt
import numpy as np
import re
import os

import contextlib
import io
import sys
from collections import defaultdict
import string


# functions ================================================================

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

# Format file


def fmt_string_file(fn_in, fn_out, str_dict, method='format_map'):
    with open(fn_in, 'r') as f:
        fdata = f.read()

    if method == 'format_map':
        fdata = string.Formatter().vformat(fdata, (), SafeDict((str_dict)))
    elif method == 'replace':
        for key in str_dict.keys():
            fdata = fdata.replace(key, str_dict[key])

    with open(fn_out, 'w') as fout:
        fout.write(fdata)

# class definition ===============================================================================


class ModelPPSetup(object):

    def __init__(self, yml_fname, model_type, write_meta_bash=False):
        # read input yaml ---------------------------------------------------------------------------
        self.yml_fname = yml_fname
        with open(yml_fname, 'r') as f:
            self.inputs = schism_yaml.load(f)

        # assign env vars to format strings with ----------------------------------------------------
        self.env_vars = build_dict(self.inputs.get('env_vars'))
        for env_var in self.env_vars:
            self.env_vars[env_var] = self.env_vars[env_var].format(
                **self.env_vars)

        # get case bundling data --------------------------------------------------------------------
        self.case_setup = self.inputs.get('case_setup').format(**self.env_vars)
        with open(self.case_setup, 'r') as f:
            self.case_inputs = schism_yaml.load(f)
        self.case_data_dir = self.case_inputs.get('output_dir')
        self.simulation_dir = os.path.join(self.env_vars['simulation_dir'])
        self.model_dir = os.path.join(self.env_vars['exp_dir'])

        # update cases
        self.case_items = self.case_inputs.get('cases')
        self.meshes = self.inputs.get('meshes')
        self.cases_exclude = self.inputs.get('cases_exclude')
        self.meshes_exclude = self.inputs.get('meshes_exclude')
        self.meshcases_exclude = self.inputs.get('meshcases_exclude')

        # get setup files
        self.pp_yml_files = self.inputs.get('pp_yml_files')
        self.pp_yml_dir = string.Formatter().vformat(self.inputs.get('pp_yml_dir'), (), SafeDict(
            ({**self.env_vars, **locals()})))

        # get meschase_dir files to copy into each meshcase_dir
        self.meshcase_dir_files = self.inputs.get('meshcase_dir_files')

        # retrieve and write out cases
        print('Handling cases:')
        for case in self.case_items:
            cname = case.get('name')
            if cname not in self.cases_exclude:

                print(
                    f'CASE: {cname} ==================================================================')

                # specify model output directory
                modcase_dir = os.path.join(self.model_dir, cname)

                crange = [case.get('case_start'),
                          case.get('case_end')]
                case_start = crange[0]
                case_end = crange[1]
                runtimedays = (case_end - case_start).days - 1
                rndays_plusone = runtimedays + 1

                for mesh_info in self.meshes:
                    meshname = mesh_info['name']
                    meshcase = f'{meshname}_{cname}'
                    if meshcase not in self.meshcases_exclude and meshname not in self.meshes_exclude:

                        # Create the mesh_case directory where the simulations will go
                        print(
                            f"\tMesh: {meshname.upper()} -------------------------------------")
                        # This is where simulations will be run
                        meshcase_dir = os.path.join(
                            self.simulation_dir, f'{meshname}_{cname}')

                        rel_study_dir = f"{os.path.basename(self.env_vars['exp_dir'])}/{os.path.relpath(meshcase_dir,self.env_vars['exp_dir'])}"

                        # Copy case directory files
                        for meschase_file in self.meshcase_dir_files:

                            shutil.copy(string.Formatter().vformat(meschase_file, (), SafeDict(
                                ({**self.env_vars, **locals()}))), meshcase_dir)

                        for yml_file in self.pp_yml_files:
                            yml_file_in = string.Formatter().vformat(yml_file, (), SafeDict(
                                ({**self.env_vars, **locals()})))
                            pp_yml_file = os.path.join(self.pp_yml_dir, os.path.basename(
                                yml_file_in).replace("MESHNAME", meshname).replace("CASENAME", cname))

                            fmt_string_file(yml_file_in,
                                            pp_yml_file,
                                            {**locals()},
                                            method='format_map')


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    yml_fname = "./input/schism_x2_evap_setup.yaml"

    mpps = ModelPPSetup(yml_fname, 'schism')
