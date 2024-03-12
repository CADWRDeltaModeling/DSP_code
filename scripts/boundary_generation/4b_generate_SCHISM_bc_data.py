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

# Format file
def fmt_string_file(fn_in, fn_out, str_dict):
    with open(fn_in, 'r') as f:
            fdata = f.read()
            
    fdata.format(str_dict)

    with open(fn_out, 'w') as fout:
        fout.write(fdata)

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

        ## get model configurations
        print(f'Storing model configurations:')
        self.model_config = self.inputs.get('model_config')
        self.configs = {}
        for config in self.model_config:
            mname = config.get('model_input')
            print(f'\t- {mname}')

            self.configs[mname] = build_dict({k: config[k] for k in set(list(config.keys())) - set(['name'])})

        # retrieve and write out cases
        print('Handling cases:')
        for case in self.case_items:

            cname = case.get('name')
            print(f'\t- {cname}')
            # create output folder
            case_dir = os.path.join(self.case_data_dir, cname)
            if not os.path.exists(case_dir):
                raise ValueError(f'Need to have bundled data saved to {case_dir} from previous effort')

            crange = [case.get('case_start'),
                      case.get('case_end')]
            
            case_perts = case['perturbations']
            
            if model_type.lower() == 'schism':

                self.setup_schism_case(cname, case_dir, crange, case_perts)

    def set_schism_vars(self):
        self.meshes = [m.format(**self.env_vars) for m in self.inputs.get('meshes')]
        self.noaa_download = build_dict(self.inputs.get('noaa_download'))
        self.noaa_pt_reyes = self.noaa_download.get('pt_reyes').format(**self.env_vars)
        self.noaa_monterey = self.noaa_download.get('monterey').format(**self.env_vars)
        self.param_tropic_base = self.inputs.get('param_tropic_base').format(**self.env_vars)
        self.param_clinic_base = self.inputs.get('param_clinic_base').format(**self.env_vars)
        self.bash_tropic = self.inputs.get('bash_tropic').format(**self.env_vars)
        self.bash_clinic = self.inputs.get('bash_clinic').format(**self.env_vars)
        self.flux_file_in = self.inputs.get('flux_file_in').format(**self.env_vars)
        self.bc_tropic_in = self.inputs.get('bc_tropic_in').format(**self.env_vars)
        self.bc_clinic_in = self.inputs.get('bc_clinic_in').format(**self.env_vars)
        self.tropic_ocean_bc = self.inputs.get('tropic_ocean_bc').format(**self.env_vars)
        

    def handle_schism_model_input(self, mod_dir, cp, case_dir, crange, pdict, method, th_file):
        # Handle the different cases
        if 'tide' in cp:
            print('\t\t\t Modifying tidal boundary')
            if pdict['method'] == 'shift':
                if 'shift_forward' in pdict['args'].keys():
                    timedelt = pdict['args']['shift_forward']
                if 'shift_backward' in pdict['args'].keys():
                    timedelt = -pdict['args']['shift_backward']

                pt_reyes = read_noaa(self.noaa_pt_reyes,force_regular=True)
                pt_reyes.index = pt_reyes.index - pd.Timedelta(days=timedelt)
                pt_reyes.to_csv(os.path.join(mod_dir,'pt_reyes.csv'))

                monterey = read_noaa(self.noaa_monterey,force_regular=True)
                monterey.index = monterey.index - pd.Timedelta(days=timedelt)
                monterey.to_csv(os.path.join(mod_dir,'monterey.csv'))
            else:
                raise ValueError(f"Unknown method to modify tidal boundary: {method}")
            
        elif any(gate in cp for gate in ['dcc','suisun']):
            print('\t\t\t Modifying Gate Operations')
            # get modified operation. th_file is the model/framework to go off of. mod_file is the data to draw from
            mod_file = os.path.join(case_dir, 
                                    f"{pdict['model_input']}_{pdict['method']}_{dt.datetime.strftime(crange[0], format='%Y%m%d0000')}-{dt.datetime.strftime(crange[1]+dt.timedelta(days=1), format='%Y%m%d0000')}.csv")
            dat_in = pd.read_csv(mod_file, parse_dates=[0], index_col=[0], header=None)
            dat_in.index = dat_in.index.strftime('%Y-%m-%dT%H:%M')
            dat_in.index.name = 'datetime'

            th_head = pd.read_csv(th_file, delim_whitespace=True, nrows=1, header=0)
            schema = {'datetime':str,'install':int,'ndup':int,'op_down':np.float64, 
                        'op_up':np.float64,'elev':np.float64,'width':np.float64,'height':int}
            df = pd.DataFrame(columns=schema.keys()).astype(schema)
            df['datetime'] = dat_in.index
            for header in th_head.columns[1:-1]:
                df[header] = th_head[header][0]

            if pdict['method'] == 'set':
                dat_out = dat_in.copy()
                if float(pdict['args']['set_value']) == 2.0:
                    dat_out[1] = 10 # schism has height of 0/10
                elif float(pdict['args']['set_value']) == 0:
                    dat_out[1] = 0 # schism has height of 0/10
            elif pdict['method'] == 'read':
                dat_out = dat_in.copy()

            df['height'] = dat_out[1].to_list()

            # write out dataframe
            dated_fn = os.path.join(mod_dir,os.path.basename(th_file).replace('.th','_dated.th'))
            df.to_csv(dated_fn, sep=' ', index=False)

            # clip to model period
            clip_fn = os.path.join(mod_dir,os.path.basename(th_file))
            file_to_elapsed(dated_fn, 
                            outpath=clip_fn,
                            start=dt.datetime(year=crange[0].year,
                                                        month=crange[0].month,
                                                        day=crange[0].day))


        elif 'dcd' in cp:
            print('hi')


    def setup_schism_case(self, cname, case_dir, crange, cperts):

        self.set_schism_vars()

        # specify model output directory
        mod_dir = os.path.join(self.env_vars['case_dir'], cname)
        if not os.path.exists(mod_dir):
            os.mkdir(mod_dir)

        # DEAL WITH TH BOUNDARIES ===========================================================================

        case_start = crange[0]
        case_end = crange[1]
        runtimedays = (case_end - case_start).days
        
        # If the tide isn't perturbed then copy the noaa files to this directory
        if not any('tide' in cp for cp in cperts):
            shutil.copyfile(self.noaa_pt_reyes, os.path.join(mod_dir, 'pt_reyes.csv'))
            shutil.copyfile(self.noaa_monterey, os.path.join(mod_dir, 'monterey.csv'))

        # Go through perturbations in this case
        for cp in cperts:
            try:
                pdict = self.perturbs[cp]
            except:
                raise ValueError(f"The perturbation {cp} needs to be defined in the perturbations section of the yaml file")

            print(f'\t\t - {cp}')
            if 'components' in pdict.keys(): #perturbation is comprised of multiple inputs
                cdict = pdict.get('components')
                subout_dir = os.path.join(case_dir, cp)
                for comp in cdict:
                    model_setup = self.configs[comp['model_input']]
                    has_version = re.search(r'v\d',cp)

                    if has_version:
                        th_file = model_setup.get('th_file').format(**self.env_vars, **{'version':has_version.group(0)})
                    else:
                        th_file = model_setup.get('th_file').format(**self.env_vars)
                    method = model_setup.get('method')
                    
                    self.handle_schism_model_input(mod_dir, cp, case_dir, crange, pdict, method, th_file)
            else: # Single-element perturbation
                model_setup = self.configs[pdict['model_input']]
                th_file = model_setup.get('th_file').format(**self.env_vars)
                method = model_setup.get('method')

                self.handle_schism_model_input(mod_dir, cp, case_dir, crange, pdict, method, th_file)

        # DEAL WITH SPATIAL SETUP ===========================================================================

        for meshname in self.meshes:

            print(f"{meshname.upper()}-------------------")
            meshcase_dir = os.path.join(case_dir,meshname)
            if not os.path.exists(meshcase_dir):
                os.mkdir(meshcase_dir)

            # TROPIC ----------------------------------------------
            print(f"Handling the tropic inputs")
            
            run_time_dict = {'runtimedays':runtimedays,
                             'year_start':case_start.year,
                             'month_start':case_start.month,
                             'day_start':case_start.day}
            print(f"\t param.nml: {self.param_tropic_base}")
            fmt_string_file(self.param_tropic_base, os.path.join(meshcase_dir,'param.nml.tropic'), run_time_dict)
                


            print(f"\t tropic.sh: {self.bash_tropic}")

            # print(f"Handling the clinic inputs")
            # with open(self.param_clinic_base, 'r') as f:
            #     self.clinic_param = schism_yaml.load(f)
        


    
if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    ## Read in param.tropic and param.clinic and modify

    yml_fname = "./input/schism_lathypcub_v2.yaml"

    mbc = ModelBCGen(yml_fname, 'schism')


    print('hi')