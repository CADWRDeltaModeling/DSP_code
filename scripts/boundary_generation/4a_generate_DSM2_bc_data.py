# script by Lily Tomkovic to create DSM2 boundary data from the the meta-latinhypercube cases' datasets
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces boundary inputs for DSM2 within existing folders

import pandas as pd
from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist

from pydelmod.create_ann_inputs import get_dss_data

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

def create_dsm2_bcs(in_fname, dsm2_config_fname):
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
        output_dir = process_output_dir(inputs)

    # create output folder
    if not os.path.exists(output_dir):
        ValueError(f'Need to create the input data and store in {output_dir}')

    if item_exist(inputs, 'cases'):
        print(f'Handling cases:')
        case_items = inputs.get('cases')
        
    if item_exist(inputs, 'perturbations'):
        perturb_items = inputs.get('perturbations')

        # update perturbations
        perturbs = {}
        print(f'Storing perturbations:')
        for perturb in perturb_items:
            pname = perturb.get('name')
            print(f'\t- {pname}')
            perturbs[pname] = build_dict({k: perturb[k] for k in set(list(perturb.keys())) - set(['name'])})

    with open(dsm2_config_fname, 'r') as f:
        dsm2_inputs = schism_yaml.load(f)
    dsm2_config = dsm2_inputs.get('model_config')

    dsp_home = dsm2_inputs.get('dsp_home')
    in_dss_dir = os.path.join(dsp_home, dsm2_inputs.get('in_dss_dir'))
    hist_dss = dsm2_inputs.get('hist_dss')
    gates_dss = dsm2_inputs.get('gates_dss')
    hist_dss_file = os.path.join(in_dss_dir, hist_dss)
    gates_dss_file = os.path.join(in_dss_dir, gates_dss)
    
    # update perturbations
    configs = {}
    primary_pathname_part_dss_filename_dict = {}
    primary_part_c_part_dict = {}
    unit_part_dict = {}
    primary_pathname_part = 'b_part'

    print(f'Storing configurations:')
    for config in dsm2_config:
        mname = config.get('model_input')
        print(f'\t- {mname}')
        primary_part = config['primary_part']

        if config['dss_file'] == 'hist_dss':
            primary_pathname_part_dss_filename_dict[primary_part] = hist_dss_file
        elif config['dss_file'] == 'gates_dss':
            primary_pathname_part_dss_filename_dict[primary_part] = gates_dss_file
        primary_part_c_part_dict[primary_part] = config['part_c']
        unit_part_dict[primary_part] = config['unit_part']

        configs[mname] = build_dict({k: config[k] for k in set(list(config.keys())) - set(['name'])})

    print('Getting input DSS files')
    df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part, 
                            primary_part_c_part_dict=primary_part_c_part_dict,
                            primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)

    # retrieve and write out cases
    print('Handling cases:')
    for case in case_items:

        cname = case.get('name')
        print(f'\t- {cname}')
        # create output folder
        case_dir = os.path.join(output_dir, cname)
        os.mkdir(os.path.join(case_dir,'DSM2'))

        # crange = [dt.datetime.strftime(case.get('case_start'), format='%Y%m%d0000'),
        #           dt.datetime.strftime(case.get('case_end')+dt.timedelta(days=1), format='%Y%m%d0000')]
        
        ######################################################################################################################
        ######################################################################################################################
        ## Generate BC inputs
        ######################################################################################################################
        cperts = case.get('perturbations')
        if cperts is not None:
                for cp in cperts:
                    try:
                        pdict = perturbs[cp]
                    except:
                        raise ValueError(f"The perturbation {cp} needs to be defined in the perturbations section of the yaml file")

                    if 'components' in pdict.keys():
                        cdict = pdict.get('components')
                    else:
                        model_input = pdict['model_input']
                        print('hi')



# def get_change(configs, model_input):
#     # Load DSS data
#     # def load(data)
#     df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part, 
#                             primary_part_c_part_dict=primary_part_c_part_dict,
#                             primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)



# Load modifications
csv_files = list(filter(lambda f: f.endswith('.csv'), os.listdir(mod_dir)))

for experiment, mod_list in mod_dict_s.items():
    
    print(f'Creating outputs for {experiment}')
    yr_select = mod_list[0]
    mod_dict = mod_list[1]
    d_part_replace = f'01JAN{yr_select} - 31DEC{yr_select}'
    paths_out = []

    # Experiment setup
    f_part_out = f'DSP_{experiment}_202309'
    out_dss = {hist_dss: f'{model_dir}/timeseries/{experiment}_hist.dss',
               gates_dss: f'{model_dir}/timeseries/{experiment}_gates.dss'}
    # create DSS files
    for key, out_file in out_dss.items():
        if not os.path.exists(out_file):
            dumdss = pyhecdss.DSSFile(out_file.format(experiment), create_new=True)

    for mod_file in mod_dict.keys():
         
        mod_key = mod_file.format(yr_select)
        print(f'\t{mod_key}')

        # read in modified/augmented data
        mod_in_file = os.path.join(mod_dir,[c for c in csv_files if mod_key in c][0])
        mod_in_ts = pd.read_csv(mod_in_file, header=None, index_col=0, names=['Modified'])
        mod_in_ts.set_index(pd.DatetimeIndex(mod_in_ts.index).to_period(freq='D'), inplace=True)
        mod_in_ts.Modified = mod_in_ts.Modified.astype(float)
            
        if mod_dict[mod_file] == 'close' or mod_dict[mod_file] == 'open':
            out_ts = mod_in_ts.copy()
            b_part = 'RSAC128' 
            if mod_dict[mod_file] == 'close':
                out_ts.Modified = 0
            elif mod_dict[mod_file] == 'open':
                out_ts.Modified = 2         
            out_ts.Modified = out_ts.Modified.astype(float)
            if 'Timestamp' not in str(type(out_ts.index[0])):
                out_ts.index = out_ts.index.to_timestamp()

            in_dss = primary_pathname_part_dss_filename_dict[b_part]
            pathname = get_pathname(in_dss, b_part, primary_part_c_part_dict[b_part])[0]
            paths_out.append(pathname)
            d_part = pathname.split('/')[1:7][3]
            e_part = pathname.split('/')[1:7][4]
            f_part = pathname.split('/')[1:7][5]
            pathname_out = pathname.replace(f_part,f_part_out)
            pathname_out = pathname_out.replace(d_part,d_part_replace)
            with pyhecdss.DSSFile(in_dss) as d_in:
                df = None
                units = None
                ptype = None
                if d_in.parse_pathname_epart(pathname).startswith('IR-'):
                    df, units, ptype = d_in.read_its(pathname)
                else:
                    df, units, ptype = d_in.read_rts(pathname)
            # add historical timeseries to beginning and end using margin_ts (copy of df)
            margin_ts = df.copy()
            margin_ts = margin_ts[(margin_ts.index < out_ts.index[0]) | (margin_ts.index > out_ts.index[-1])]
            margin_ts.columns = ['Modified']
            out_ts = pd.concat([out_ts, margin_ts])
            out_ts.sort_index(inplace=True)

        # DCC Gate gets NaN-ed from the IR-DAY averaging in create_ann_inputs
        elif (len(mod_dict[mod_file]) == 1) or (mod_dict[mod_file] == 'plus_ten'):
            if mod_dict[mod_file] == 'plus_ten':
                b_part = 'RSAC155'
            else:
                b_part = mod_dict[mod_file][0]

            in_dss = primary_pathname_part_dss_filename_dict[b_part]
            pathname = get_pathname(in_dss, b_part, primary_part_c_part_dict[b_part])[0]
            paths_out.append(pathname)
            d_part = pathname.split('/')[1:7][3]
            e_part = pathname.split('/')[1:7][4]
            f_part = pathname.split('/')[1:7][5]
            pathname_out = pathname.replace(f_part,f_part_out)
            pathname_out = pathname_out.replace(d_part,d_part_replace)
            with pyhecdss.DSSFile(in_dss) as d_in:
                df = None
                units = None
                ptype = None
                if d_in.parse_pathname_epart(pathname).startswith('IR-'):
                    df, units, ptype = d_in.read_its(pathname)
                else:
                    df, units, ptype = d_in.read_rts(pathname)

            if mod_dict[mod_file] == 'plus_ten':
                # add 10% of flow to the input ts for sacramento
                out_ts = df.copy()
                out_ts.columns = ['Modified']
                out_ts.Modified = out_ts.Modified * 1.1
            elif 'Tidal_Amp' in mod_file:
                # manully shift timeseries by tide_shift days
                out_ts = df.copy()
                out_ts.index = out_ts.index + pd.Timedelta(days=tide_shift)
            else:
                out_ts = mod_in_ts.copy()
                # add historical timeseries to beginning and end using margin_ts (copy of df)
                if 'Timestamp' not in str(type(out_ts.index[0])):
                    out_ts.index = out_ts.index.to_timestamp()
                margin_ts = df.copy()
                if 'Timestamp' not in str(type(margin_ts.index[0])):
                    margin_ts.index = margin_ts.index.to_timestamp()
                margin_ts = margin_ts[(margin_ts.index < out_ts.index[0]) | (margin_ts.index > out_ts.index[-1])]
                margin_ts.columns = ['Modified']
                out_ts = pd.concat([out_ts, margin_ts])
                out_ts.sort_index(inplace=True)

        if type(mod_dict[mod_file]) != list or (type(mod_dict[mod_file]) == list and len(mod_dict[mod_file]) == 1):
            if 'Timestamp' not in str(type(out_ts.index[0])):
                out_ts.index = out_ts.index.to_timestamp()
            
            dss_out_file = out_dss[os.path.basename(in_dss)]
            with pyhecdss.DSSFile(dss_out_file) as d_out:
                print(f'Writing out {pathname_out}')
                if e_part.startswith('IR-'):
                    d_out.write_its(pathname_out, out_ts, unit_part_dict[b_part], ptype) # write to output DSS file
                else:
                    out_ts.index = out_ts.index + pd.Timedelta(days=1) # need to shift by one day because of DSS writing timestamp issues
                    if out_ts.index.freq is None:
                        out_ts.index.freq = pd.infer_freq(out_ts.index)
                    d_out.write_rts(pathname_out, out_ts, unit_part_dict[b_part], ptype) # write regular output to DSS file
        else:
            # scale flows per day
            in_ts = df_input[mod_dict[mod_file]]
            scale_in_ts = in_ts.copy()
            scale_in_ts = scale_in_ts.apply(lambda x: x.div(x.sum()), axis=1)

            # Join the scaled dataframe with the modified input timeseries 
            merge = pd.merge(mod_in_ts,scale_in_ts, how='inner', left_index=True, right_index=True)
            merge.iloc[:, 1:] =  merge.iloc[:, 1:].multiply(merge.iloc[:,0], axis='index') # distribute the modified timeseries across columns
            merge.index = merge.index.to_timestamp()

            # Loop through components and write out to DSS
            for b_part in mod_dict[mod_file]:
                in_dss = primary_pathname_part_dss_filename_dict[b_part]
                pathname = get_pathname(in_dss, b_part, primary_part_c_part_dict[b_part])[0]
                paths_out.append(pathname)
                d_part = pathname.split('/')[1:7][3]
                e_part = pathname.split('/')[1:7][4]
                f_part = pathname.split('/')[1:7][5]
                pathname_out = pathname.replace(f_part,f_part_out)
                pathname_out = pathname_out.replace(d_part,d_part_replace)
                with pyhecdss.DSSFile(in_dss) as d_in:
                    df = None
                    units = None
                    ptype = None
                    df, units, ptype = d_in.read_rts(pathname)
                    dss_out_file = out_dss[os.path.basename(in_dss)]
                    out_ts = merge[[b_part]].copy()
                    out_ts.columns = ['Modified']
                    # add historical timeseries to beginning and end using margin_ts (copy of df)
                    if 'Timestamp' not in str(type(out_ts.index[0])):
                        out_ts.index = out_ts.index.to_timestamp()
                    margin_ts = df.copy()
                    if 'Timestamp' not in str(type(margin_ts.index[0])):
                        margin_ts.index = margin_ts.index.to_timestamp()
                    margin_ts = margin_ts[(margin_ts.index < out_ts.index[0]) | (margin_ts.index > out_ts.index[-1])]
                    margin_ts.columns = ['Modified']
                    out_ts = pd.concat([out_ts, margin_ts])
                    out_ts.sort_index(inplace=True)
                    out_ts.index = out_ts.index + pd.Timedelta(days=1) # need to shift by one day because of DSS writing timestamp issues
                    with pyhecdss.DSSFile(dss_out_file) as d:
                        print(f'Writing out {pathname_out}')
                        if out_ts.index.freq is None:
                            out_ts.index.freq = pd.infer_freq(out_ts.index)
                        d.write_rts(pathname_out, out_ts, unit_part_dict[b_part], ptype)

    ###### Copy Remaining paths into modified DSS files

    for dss_in in out_dss:
        with pyhecdss.DSSFile(os.path.join(in_dss_dir, dss_in)) as d_in:
            with pyhecdss.DSSFile(out_dss[dss_in]) as d_out:
                incat = d_in.read_catalog() # all the pathnames in the DSS file
                paths_in = d_in.get_pathnames(incat)

                missing_paths = list(set(paths_in) - set(paths_out)) # paths in the input DSS that haven't been copied out

                for p in missing_paths:
                    df = None
                    units = None
                    ptype = None
                    print(f'Writing out {p} units {units} pytype {ptype}')
                    if d_in.parse_pathname_epart(p).startswith('IR-'):
                        df, units, ptype = d_in.read_its(p)
                        if units == 'und': units = 'UNSPECIF'
                        d_out.write_its(p, df, units, ptype) # write to output DSS file
                    else:
                        df,units,ptype=d_in.read_rts(p)
                        if units == 'und': units = 'UNSPECIF'
                        d_out.write_rts(p, df, units, ptype) # write to output DSS file


        

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # model_dir = r"D:\projects\delta_salinity\model\schism\dsp_202311_baseline"
    in_fname = "../../data/lathypcub_v1_setup.yaml"
    
    cases = create_cases()
    dsm2_config_fname = "../../data/lathypcub_v1_dsm2_config.yaml"
    # in_fname = "../../../../model/schism/dsp_202311_baseline/dsp_baseline_bay_delta.yaml"

    # args = Namespace(main_inputfile=in_fname)

    create_dsm2_bcs(in_fname, dsm2_config_fname)