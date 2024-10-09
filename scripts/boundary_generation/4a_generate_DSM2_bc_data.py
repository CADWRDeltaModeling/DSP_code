# script by Lily Tomkovic to create DSM2 boundary data from the the meta-latinhypercube cases' datasets
# takes input yaml (format in progress to be somewhat forward compatible with eventual format of repo)
# produces boundary inputs for DSM2 within existing folders

import pandas as pd
from schimpy import schism_yaml
from schimpy.prepare_schism import process_output_dir, check_nested_match, item_exist

from pydelmod.create_ann_inputs import get_dss_data
import pyhecdss

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
            filtered_df = filtered_df[(
                catdf.B == b_part)] if filtered_df is not None else catdf[(catdf.B == b_part)]
        if c_part is not None:
            filtered_df = filtered_df[(
                catdf.C == c_part)] if filtered_df is not None else catdf[(catdf.C == c_part)]
        if e_part is not None:
            filtered_df = filtered_df[(
                catdf.E == e_part)] if filtered_df is not None else catdf[(catdf.E == e_part)]
        if f_part is not None:
            filtered_df = filtered_df[(
                catdf.F == f_part)] if filtered_df is not None else catdf[(catdf.F == f_part)]
        if filter_b_part_numeric:
            filtered_df = filtered_df[(catdf.B.str.isnumeric())]
        path_list = d.get_pathnames(filtered_df)

    return path_list


def update_DSS(pdict, case_dir, crange, in_dss, out_dss, pathname,
               paths_out, f_part_out, d_part_replace, unit_part, ts_df=None):
    # Take in information about the perturbation and update the ouptut DSS file

    model_input = pdict['model_input']
    print(f'\t Creating outputs for {model_input}')

    mod_file = os.path.join(
        case_dir, f'{pdict["model_input"]}_{pdict["method"]}_{crange[0]}-{crange[1]}.csv')

    if ts_df is None:
        # read in modified/augmented data
        ts_df = pd.read_csv(mod_file, header=None,
                            index_col=0, names=['Modified'])

    ts_df.set_index(pd.DatetimeIndex(ts_df.index), inplace=True)
    ts_df.Modified = ts_df.Modified.astype(float)
    if 'Timestamp' not in str(type(ts_df.index[0])):
        ts_df.index = ts_df.index.to_timestamp()

    paths_out.append(pathname)

    d_part = pathname.split('/')[1:7][3]
    e_part = pathname.split('/')[1:7][4]
    f_part = pathname.split('/')[1:7][5]
    pathname_out = pathname.replace(f_part, f_part_out)
    pathname_out = pathname_out.replace(d_part, d_part_replace)

    # Get rest of historical record to append to either side of the input dataframe
    with pyhecdss.DSSFile(in_dss) as d_in:
        df = None
        units = None
        ptype = None
        if d_in.parse_pathname_epart(pathname).startswith('IR-'):
            df, units, ptype = d_in.read_its(pathname)
        else:
            df, units, ptype = d_in.read_rts(pathname)
    if 'Timestamp' not in str(type(df.index[0])):
        df.index = df.index.to_timestamp()

    # add historical timeseries to beginning and end using margin_ts (copy of df)
    margin_ts = df.copy()
    margin_ts = margin_ts[(margin_ts.index < ts_df.index[0])
                          | (margin_ts.index > ts_df.index[-1])]
    margin_ts.columns = ['Modified']
    ts_df = pd.concat([ts_df, margin_ts])
    ts_df.sort_index(inplace=True)

    # Write out
    dss_out_file = out_dss[locals()['in_dss']]

    with pyhecdss.DSSFile(dss_out_file) as d_out:
        print(f'Writing out {pathname_out}')
        if e_part.startswith('IR-'):
            d_out.write_its(pathname_out, ts_df, unit_part,
                            ptype)  # write to output DSS file
        else:
            # need to shift by one day because of DSS writing timestamp issues
            ts_df.index = ts_df.index + pd.Timedelta(days=1)
            if ts_df.index.freq is None:
                ts_df.index.freq = pd.infer_freq(ts_df.index)
            # write regular output to DSS file
            d_out.write_rts(pathname_out, ts_df, unit_part, ptype)

    return paths_out


def create_dsm2_bcs(in_fname, dsm2_config_fname):
    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
        output_dir = inputs['output_dir']

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
            perturbs[pname] = build_dict(
                {k: perturb[k] for k in set(list(perturb.keys())) - set(['name'])})

    with open(dsm2_config_fname, 'r') as f:
        dsm2_inputs = schism_yaml.load(f)
    dsm2_config = dsm2_inputs.get('model_config')

    training_set = dsm2_inputs.get('training_set')
    dsp_home = dsm2_inputs.get('dsp_home')
    in_dss_dir = os.path.join(dsp_home, dsm2_inputs.get('in_dss_dir'))
    out_dss_dir = os.path.join(dsp_home, dsm2_inputs.get('out_dss_dir'))
    hist_dss = os.path.join(in_dss_dir, dsm2_inputs.get('hist_dss'))
    gates_dss = os.path.join(in_dss_dir,  dsm2_inputs.get('gates_dss'))
    dcd_dss = os.path.join(in_dss_dir,  dsm2_inputs.get('dcd_dss'))

    # check that the timeseries output folder exists
    dsm2_dir = os.path.join(out_dss_dir, f'{training_set}/timeseries/')
    if not os.path.exists(dsm2_dir):
        os.makedirs(dsm2_dir)

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

        if 'get_parts_from_csv' in config.keys():
            csv_in = pd.read_csv(config.get(
                'get_parts_from_csv'), parse_dates=[0], index_col=[0])
            for cpath in csv_in.columns:
                path_parts = cpath.split("/")
                primary_part = path_parts[2]
                primary_pathname_part_dss_filename_dict[primary_part] = locals()[
                    config.get('dss_file')]
                primary_part_c_part_dict[primary_part] = path_parts[3]
                unit_part_dict[primary_part] = config['unit_part']
        else:
            primary_part = config['primary_part']

            primary_pathname_part_dss_filename_dict[primary_part] = locals()[
                config.get('dss_file')]
            primary_part_c_part_dict[primary_part] = config['part_c']
            unit_part_dict[primary_part] = config['unit_part']

        configs[mname] = build_dict(
            {k: config[k] for k in set(list(config.keys())) - set(['name'])})

    print('Getting input DSS files')
    # df_input = get_dss_data(primary_pathname_part_dss_filename_dict, primary_pathname_part,
    #                         primary_part_c_part_dict=primary_part_c_part_dict,
    #                         primary_part_e_part_dict=None, primary_part_f_part_dict=None, daily_avg=True, filter_b_part_numeric=False)

    # retrieve and write out cases
    print('Handling cases:')
    for case in case_items:

        cname = case.get('name')
        print(f'\t- {cname}')
        # create output folder
        case_dir = os.path.join(output_dir, cname)

        crange = [dt.datetime.strftime(case.get('case_start'), format='%Y%m%d0000'),
                  dt.datetime.strftime(case.get('case_end')+dt.timedelta(days=1), format='%Y%m%d0000')]

        ######################################################################################################################
        ######################################################################################################################
        # Generate BC inputs
        ######################################################################################################################
        # Case setup
        f_part_out = f'DSP_{cname}'
        out_dss = {hist_dss: f'{dsm2_dir}/{cname}_hist.dss',
                   gates_dss: f'{dsm2_dir}/{cname}_gates.dss',
                   dcd_dss: f'{dsm2_dir}/{cname}_dcd.dss'}
        # create DSS files
        for key, out_file in out_dss.items():
            if not os.path.exists(out_file):
                dumdss = pyhecdss.DSSFile(out_file, create_new=True)

        cperts = case.get('perturbations')
        start_date = case.get('case_start')
        end_date = case.get('case_end')
        d_part_replace = f'{dt.datetime.strftime(start_date, format="%d%B%Y")} - {dt.datetime.strftime(end_date, format="%d%B%Y")}'

        # track written files
        paths_out = []

        if cperts is not None:
            for cp in cperts:
                try:
                    pdict = perturbs[cp]
                except:
                    raise ValueError(
                        f"The perturbation {cp} needs to be defined in the perturbations section of the yaml file")

                if 'components' in pdict.keys():
                    cdict = pdict.get('components')
                    subout_dir = os.path.join(case_dir, cp)

                    for comp in cdict:
                        model_input = configs[comp['model_input']]
                        if 'get_parts_from_csv' in model_input.keys():
                            # need to cycle through all inputs of the DCD files and write them out to the dcd output file
                            # csv_in = pd.read_csv(config.get('get_parts_from_csv'), parse_dates=[0], index_col=[0])
                            mod_file = os.path.join(
                                subout_dir, f'{comp["model_input"]}_{comp["method"]}_{crange[0]}-{crange[1]}.csv')
                            in_df = pd.read_csv(
                                mod_file, parse_dates=[0], index_col=[0])

                            for cpath in in_df.columns:
                                path_parts = cpath.split("/")
                                primary_part = path_parts[2]
                                in_dss = locals()[model_input.get('dss_file')]
                                pathname = get_pathname(
                                    in_dss, primary_part, path_parts[3], e_part=path_parts[5], f_part=path_parts[6])[0]
                                unit_part = model_input.get('unit_part')
                                ts_df = in_df[[cpath]].copy()
                                ts_df.columns = ['Modified']

                                paths_out = update_DSS(comp, subout_dir, crange, in_dss, out_dss, pathname,
                                                       paths_out, f_part_out, d_part_replace, unit_part, ts_df=ts_df)

                        else:
                            b_part = model_input['primary_part']
                            in_dss = primary_pathname_part_dss_filename_dict[b_part]
                            pathname = get_pathname(
                                in_dss, b_part, primary_part_c_part_dict[b_part])[0]
                            unit_part = unit_part_dict[b_part]

                            paths_out = update_DSS(comp, subout_dir, crange, in_dss, out_dss, pathname,
                                                   paths_out, f_part_out, d_part_replace, unit_part)
                else:
                    model_input = configs[pdict['model_input']]
                    if 'get_parts_from_csv' in model_input.keys():
                        in_df = pd.read_csv(
                            model_input['get_parts_from_csv'], parse_dates=[0], index_col=[0])

                        for cpath in in_df.columns:
                            path_parts = cpath.split("/")
                            primary_part = path_parts[2]
                            in_dss = locals()[model_input.get('dss_file')]
                            pathname = get_pathname(
                                in_dss, primary_part, path_parts[3], e_part=path_parts[5], f_part=path_parts[6])[0]
                            unit_part = model_input.get('unit_part')
                            ts_df = in_df[[cpath]].copy()
                            ts_df.columns = ['Modified']

                            paths_out = update_DSS(pdict, case_dir, crange, in_dss, out_dss, pathname,
                                                   paths_out, f_part_out, d_part_replace, unit_part, ts_df=ts_df)

                    else:
                        b_part = model_input['primary_part']
                        in_dss = primary_pathname_part_dss_filename_dict[b_part]
                        pathname = get_pathname(
                            in_dss, b_part, primary_part_c_part_dict[b_part])[0]
                        unit_part = unit_part_dict[b_part]

                        paths_out = update_DSS(pdict, case_dir, crange, in_dss, out_dss, pathname,
                                               paths_out, f_part_out, d_part_replace, unit_part)

        # Copy Remaining paths into modified DSS files

        for dss_in in out_dss:
            with pyhecdss.DSSFile(dss_in) as d_in:
                with pyhecdss.DSSFile(out_dss[dss_in]) as d_out:
                    incat = d_in.read_catalog()  # all the pathnames in the DSS file
                    paths_in = d_in.get_pathnames(incat)

                    # paths in the input DSS that haven't been copied out
                    missing_paths = list(set(paths_in) - set(paths_out))

                    for p in missing_paths:
                        df = None
                        units = None
                        ptype = None
                        print(f'Writing out {p} units {units} pytype {ptype}')
                        if d_in.parse_pathname_epart(p).startswith('IR-'):
                            df, units, ptype = d_in.read_its(p)
                            if units == 'und':
                                units = 'UNSPECIF'
                            # write to output DSS file
                            d_out.write_its(p, df, units, ptype)
                        else:
                            df, units, ptype = d_in.read_rts(p)
                            if units == 'und':
                                units = 'UNSPECIF'
                            # write to output DSS file
                            d_out.write_rts(p, df, units, ptype)


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # model_dir = r"D:\projects\delta_salinity\model\schism\dsp_202311_baseline"
    in_fname = "./input/lathypcub_v3_setup.yaml"

    # cases = create_cases()
    dsm2_config_fname = "./input/lathypcub_v3_dsm2_config.yaml"
    # in_fname = "../../../../model/schism/dsp_202311_baseline/dsp_baseline_bay_delta.yaml"

    # args = Namespace(main_inputfile=in_fname)

    create_dsm2_bcs(in_fname, dsm2_config_fname)
