# script by Lily Tomkovic to create DSM2 boundary data from the the meta-latinhypercube cases' datasets
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces boundary inputs for DSM2 within existing folders

import pandas as pd
import numpy as np

from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist

from pydelmod.create_ann_inputs import get_dss_data
import pyhecdss

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

def create_perturbations(row, pert_vars):
    perturbations = [f"{var} {row[var].lower()}" for var in pert_vars if 'regular' not in row[var].lower()]
    return perturbations


def create_dsm2_inps(in_fname, dsm2_config_fname):
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
        
            case_items = case_file # this behaves the same as if the cases are defined in a yaml
        else:
            # make into pd.DataFrame
            case_df = pd.DataFrame(columns=list(case_items[0].keys()))
            for case in case_items:
                case_df = pd.concat([case_df, pd.DataFrame.from_dict([case])], ignore_index=True)
            case_items = case_df


    with open(dsm2_config_fname, 'r') as f:
        dsm2_inputs = schism_yaml.load(f)
    dsm2_config = dsm2_inputs.get('model_config')
    # create output dir
    inp_dir = os.path.join(dsm2_inputs.get('out_dss_dir'), dsm2_inputs.get('training_set'))

    config_infile = dsm2_inputs.get('config_file')
    hydro_infile = dsm2_inputs.get('hydro_file')
    qual_infile = dsm2_inputs.get('qual_file')
    qual_x2_infile = dsm2_inputs.get('qual_x2_file')

    treatments = dsm2_inputs.get('treatments')

    mod_spinup_days = dsm2_inputs.get('mod_spinup_days')

    # retrieve and write out cases
    print('Handling cases:')
    for index, case in case_items.iterrows():

        cname = case.get('name')
        print(f'\t- {cname}')

        config_dict = {'sim_start': dt.datetime.strftime(case.get('case_start')-pd.Timedelta(days=mod_spinup_days), format='%d%b%Y').upper(),
                       'qual_start': dt.datetime.strftime(case.get('case_start')-pd.Timedelta(days=mod_spinup_days-1), format='%d%b%Y').upper(),
                       'sim_end': dt.datetime.strftime(case.get('case_end'), format='%d%b%Y').upper(),
                       'case_name': cname}

        perts = case.get('perturbations')
        if perts is np.nan:
            perts = [] 
            
        # config vars:
        for treat in treatments:
            if isinstance(treat['search'], str):
                # If it's a single string, check if it's in any of the elements of `perts`
                result = any(treat['search'] in pert for pert in perts)
            elif isinstance(treat['search'], list):
                # If it's a list of strings, check if any string in the list is in any of the elements of `perts`
                result = any(any(search in pert for search in treat['search']) for pert in perts)
            else:
                result = False
            if result:
                config_dict[f"{treat['name']}_treatment"] = f'DSP_{cname.upper()}'
            else:
                default_treat = treat['default'] if "{" not in treat['default'] else f"${treat['default']}"
                config_dict[f"{treat['name']}_treatment"] = default_treat
        
        config_filename = os.path.basename(config_infile).replace('CASE',str(''.join(filter(str.isdigit, cname)))) # replace with the proper filename
        hydro_filename = os.path.basename(hydro_infile).replace('CASE',str(''.join(filter(str.isdigit, cname)))) # replace with the proper filename
        qual_filename = os.path.basename(qual_infile).replace('CASE',str(''.join(filter(str.isdigit, cname)))) # replace with the proper filename
        qual_x2_filename = os.path.basename(qual_x2_infile).replace('CASE',str(''.join(filter(str.isdigit, cname)))) # replace with the proper filename

        fmt_string_file(config_infile, os.path.join(inp_dir, config_filename), config_dict, method='format_map')
        fmt_string_file(hydro_infile, os.path.join(inp_dir, hydro_filename), SafeDict(({**config_dict,
                                                                                        **locals()})), method='format_map')
        fmt_string_file(qual_infile, os.path.join(inp_dir, qual_filename), SafeDict(({**config_dict,
                                                                                      **locals()})), method='format_map')
        fmt_string_file(qual_x2_infile, os.path.join(inp_dir, qual_x2_filename), SafeDict(({**config_dict,
                                                                                           **locals()})), method='format_map')


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # in_fname = "./input/lathypcub_v4_setup.yaml"
    # dsm2_config_fname = "./input/lathypcub_v4_dsm2_config.yaml"
    
    in_fname = "./input/lathypcub_v3_setup.yaml"
    dsm2_config_fname = "./input/lathypcub_v3_dsm2_config.yaml"
    
    create_dsm2_inps(in_fname, dsm2_config_fname)