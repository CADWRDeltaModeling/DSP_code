# script by Lily Tomkovic to create boundary datasets for the meta-latinhypercube cases
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces case folders with any modified boundary inputs

import pandas as pd
from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist
from schimpy.schism_setup import check_and_suggest
from argparse import Namespace
import shutil
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

def write_dat_out(cdict, start_date, end_date, crange, out_dir):
    
    dat_out, header = apply_method(cdict)
    try:
        dat_out.index = dat_out.index.to_pydatetime()
    except:
        None
    dat_out = dat_out[dat_out.index.to_series().between(pd.to_datetime(start_date),pd.to_datetime(end_date))]
    dat_out.to_csv(os.path.join(out_dir, f'{cdict["model_input"]}_{cdict["method"]}_{crange[0]}-{crange[1]}.csv'),
                    header=header)

def apply_method(pert_dict):
    method = pert_dict['method']
    file_in = pert_dict['args']['file']
    header = None
    if method in ['read_dcd', 'read_suisun']:
        dat_in = pd.read_csv(file_in, parse_dates=[0], index_col=[0])
        dat_out = dat_in.copy()
        header = True
    else:
        dat_in = pd.read_csv(file_in, parse_dates=[0], index_col=[0], header=None)

    if method == 'read':
        dat_out = dat_in.copy()
    elif method == 'shift':
        dat_out = dat_in.copy()
        dat_out.index = dat_out.index + pd.Timedelta(days=pert_dict['args']['shift'])
    elif method == 'set':
        dat_out = dat_in.copy()
        dat_out[1] = float(pert_dict['args']['set_value'])
    elif method == 'scale':
        dat_out = dat_in.copy()
        dat_out[1] = float(pert_dict['args']['scale_factor'])*dat_out[1]
    elif method in ['read_dcd', 'read_suisun']:
        pass
    else:
        raise ValueError(f'Perturbation method "{method}" is not defined in the code at the moment.')
    return dat_out, header

# read in yaml
def create_cases(in_fname, skip=0):
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
        output_dir = inputs['output_dir']

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

        if 'filename' in case_items[0].keys():
            
            # setup cases with csv file and not yaml
            case_file = pd.read_csv(case_items[0]['filename'])
            pert_vars = case_items[0]['pert_vars']
            case_file['start'] = pd.to_datetime(case_file['start'])
            case_file['end'] = pd.to_datetime(case_file['end'])
            case_file['name'] = case_file['case']

            case_file = case_file.iloc[skip:].copy()


            # loop through cases
            for index, row in case_file.iterrows():
                
                cname = row['name']
                print(f'\t- {cname}')

                # create output folder
                case_dir = os.path.join(output_dir, cname)
                if not os.path.exists(case_dir):
                    print(f'\t\t Creating output directory {case_dir}')
                    os.mkdir(case_dir)

                start_date = row['start']
                end_date = row['end']
                crange = [dt.datetime.strftime(row['start'], format='%Y%m%d0000'),
                        dt.datetime.strftime(row['end'] + dt.timedelta(days=1), format='%Y%m%d0000')]

                # get perturbations
                for pv in pert_vars:

                    cp = f'{pv} {row[pv]}'.lower()
                    if "regular" in cp:
                        continue
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
                            write_dat_out(comp, start_date, end_date, crange, subout_dir)
                    else:
                        write_dat_out(pdict, start_date, end_date, crange, case_dir)
                        
                        if 'copy_files' in pdict['args'].keys():
                            for cfile in pdict['args'].get('copy_files'):
                                # so far this is for SCHISM which used the th files to generate the perturbation 
                                # so it was simpler to write the SCHISM th files out in 2b_perturb_binary
                                cbase = os.path.basename(cfile).split('.') # just in case we want to add something to this filename later
                                shutil.copy2(cfile,os.path.join(case_dir, f'{cbase[0]}.{cbase[1]}'))

        else:
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
                                write_dat_out(comp, start_date, end_date, crange, subout_dir)
                        else:
                            write_dat_out(pdict, start_date, end_date, crange, case_dir)
                            
                            if 'copy_files' in pdict['args'].keys():
                                for cfile in pdict['args'].get('copy_files'):
                                    # so far this is for SCHISM which used the th files to generate the perturbation 
                                    # so it was simpler to write the SCHISM th files out in 2b_perturb_binary
                                    cbase = os.path.basename(cfile).split('.') # just in case we want to add something to this filename later
                                    shutil.copy2(cfile,os.path.join(case_dir, f'{cbase[0]}.{cbase[1]}'))


    

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # model_dir = r"D:\projects\delta_salinity\model\schism\dsp_202311_baseline"
    in_fname = "./input/lathypcub_v4_setup.yaml"
    # in_fname = "../../../../model/schism/dsp_202311_baseline/dsp_baseline_bay_delta.yaml"

    # args = Namespace(main_inputfile=in_fname)

    create_cases(in_fname, skip=91)