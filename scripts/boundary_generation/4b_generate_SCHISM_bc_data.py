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

def clip_file_to_start(infile, outpath, start):
    with contextlib.redirect_stdout(io.StringIO()):
        file_to_elapsed(infile, outpath=outpath, start=start)
        sys.stdout.close()


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

# make links


def make_links(start_date, end_date, src_dir, src_dir_narr, link_dir):
    # Creates sflux links from start to end _dates using the src and narr dirs
    # the link_dir is where the links will be created
    delt = dt.timedelta(days=1)
    current = start_date
    if (current >= end_date):
        print(
            f'ERROR: Start date {start_date.strftime("%b %d, %Y")} is after end date {end_date.strftime("%b %d, %Y")}')
    nfile = 0
    while (current <= end_date):
        # Air data
        # Ours
        src_str_air = os.path.join(src_dir,
                                   "baydelta_schism_air_%s%02d%02d.nc" % (current.year, current.month, current.day))
        src_str_air = os.path.relpath(src_str_air, link_dir)
        # NARR
        src_str_rad = os.path.join(src_dir_narr,
                                   "%4d_%02d/narr_rad.%4d_%02d_%02d.nc" % (current.year, current.month, current.year, current.month, current.day))
        src_str_rad = os.path.relpath(src_str_rad, link_dir)
        src_str_prc = os.path.join(src_dir_narr,
                                   "%4d_%02d/narr_prc.%4d_%02d_%02d.nc" % (current.year, current.month, current.year, current.month, current.day))
        src_str_prc = os.path.relpath(src_str_prc, link_dir)
        nfile += 1
        link_str_air = os.path.join(link_dir, "sflux_air_1.%04d.nc" % (nfile))
        link_str_rad = os.path.join(link_dir, "sflux_rad_1.%04d.nc" % (nfile))
        link_str_prc = os.path.join(link_dir, "sflux_prc_1.%04d.nc" % (nfile))
        if not os.path.islink(link_str_air):
            os.symlink(src_str_air, link_str_air)
        else:
            os.remove(link_str_air)
            os.symlink(src_str_air, link_str_air)
        if not os.path.islink(link_str_rad):
            os.symlink(src_str_rad, link_str_rad)
        else:
            os.remove(link_str_rad)
            os.symlink(src_str_rad, link_str_rad)
        if not os.path.islink(link_str_prc):
            os.symlink(src_str_prc, link_str_prc)
        else:
            os.remove(link_str_prc)
            os.symlink(src_str_prc, link_str_prc)
        current += delt


def write_noaa(infile, in_df, outfile, repl_col=1):
    with open(infile, 'r') as in_f:
        with open(outfile, 'w') as out_f:
            inlines = in_f.readlines()
            for inl in inlines:
                if re.match("[0-9]{4}-[0-9]{2}-[0-9]{2}", inl[:10]):
                    break
                else:
                    out_f.write(inl)
            for dateindex, row in in_df.iterrows():
                out_f.write(
                    f'{dateindex},{row["Water Level"]},999,0,0,0,0,v\n')

# class definition ===============================================================================
# TODO make the ModelBCGen object able to take in DSM2 input yamls and output DSM2 BCs


class ModelBCGen(object):

    def __init__(self, yml_fname, model_type, machine, write_meta_bash=False):
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
        perturb_items = self.case_inputs.get('perturbations')

        # update perturbations
        perturbs = {}
        print(f'Storing perturbations:')
        for perturb in perturb_items:
            pname = perturb.get('name')
            print(f'\t- {pname}')
            perturbs[pname] = build_dict(
                {k: perturb[k] for k in set(list(perturb.keys())) - set(['name'])})
        self.perturbs = perturbs

        # update cases
        self.case_items = self.case_inputs.get('cases')
        self.cases_exclude = self.inputs.get('cases_exclude')
        self.meshes_exclude = self.inputs.get('meshes_exclude')
        self.meshcases_exclude = self.inputs.get('meshcases_exclude')

        # get model configurations
        print(f'Storing model configurations:')
        self.model_config = self.inputs.get('model_config')
        self.configs = {}
        for config in self.model_config:
            mname = config.get('model_input')
            print(f'\t- {mname}')

            self.configs[mname] = build_dict(
                {k: config[k] for k in set(list(config.keys())) - set(['name'])})

        self.machine = machine

        if model_type.lower() == 'schism':
            self.setup_schism_model()
            # {new_list:[] for new_list in ['tropic','clinic']}
            self.bash_file_dict = defaultdict(list)

        # retrieve and write out cases
        print('Handling cases:')
        for case in self.case_items:
            cname = case.get('name')
            if cname not in self.cases_exclude:

                print(
                    f'CASE: {cname} ==================================================================')
                # define bundle output dir
                case_data_dir = os.path.join(self.case_data_dir, cname)
                if not os.path.exists(case_data_dir):
                    raise ValueError(
                        f'Need to have bundled data saved to {case_data_dir} from previous effort')

                crange = [case.get('case_start'),
                          case.get('case_end')]

                if 'perturbations' in case.keys():
                    case_perts = case['perturbations']
                else:
                    case_perts = {}

                if model_type.lower() == 'schism':

                    self.setup_schism_case(
                        cname, case_data_dir, crange, case_perts, case.get('model_year'))

        if model_type.lower() == 'schism' and write_meta_bash:
            # write meta-bash files
            for bash_cat in self.bash_file_dict:
                meta_bash = ''
                # create_combined bash output files
                self.bash_file_dict[bash_cat]
                for bash_file in self.bash_file_dict[bash_cat]:
                    bash_line = os.path.relpath(bash_file, self.model_dir)
                    # link any th files from the base directory
                    meta_bash += f'bash {bash_line}\n'
                with open(os.path.join(self.model_dir, f'{bash_cat}_meta_bash.sh'), 'w') as mbf:
                    mbf.write(meta_bash)

    def setup_schism_model(self):
        self.set_schism_vars()
        # list of files needed to be linked in .sh file
        self.linked_th_files = self.inputs.get('linked_th_files')
        # list of source/sink files needed to be linked in .sh file
        self.linked_ss_th_files = self.inputs.get('linked_ss_th_files')
        self.linked_tracer_th_files = self.inputs.get(
            'linked_tracer_th_files')  # list of tracer files

        # specify model output directory
        self.model_dir = os.path.join(self.env_vars['exp_dir'])
        self.simulations_dir = self.env_vars['simulation_dir']
        if not os.path.exists(self.simulations_dir):
            os.mkdir(self.simulations_dir)

    def setup_schism_case(self, cname, case_data_dir, crange, cperts, case_year):

        perturbed_th = {}

        # specify model output directory
        modcase_dir = os.path.join(self.model_dir, cname)
        if not os.path.exists(modcase_dir):
            os.mkdir(modcase_dir)

        # DEAL WITH TH BOUNDARIES ===========================================================================

        case_start = crange[0]
        case_end = crange[1]
        runtimedays = (case_end - case_start).days - 1
        write_flux = False
        fluxes_out = {}

        # If the tide isn't perturbed then copy the noaa files to this directory
        if not any('tide' in cp for cp in cperts):
            shutil.copyfile(self.noaa_pt_reyes, os.path.join(
                modcase_dir, 'pt_reyes.csv'))
            shutil.copyfile(self.noaa_monterey, os.path.join(
                modcase_dir, 'monterey.csv'))
            tidal_pert = ''

        # Go through perturbations in this case
        for cp in cperts:
            try:
                pdict = self.perturbs[cp]
            except:
                raise ValueError(
                    f"The perturbation {cp} needs to be defined in the perturbations section of the yaml file")

            print(f'\t - {cp}')
            if 'components' in pdict.keys():  # perturbation is comprised of multiple inputs
                cdict = pdict.get('components')
                subout_dir = os.path.join(case_data_dir, cp)
                for comp in cdict:
                    model_setup = self.configs[comp['model_input']]
                    has_version = re.search(r'v\d', cp)

                    if has_version:
                        th_file = string.Formatter().vformat(model_setup.get('th_file'), (),
                                                             SafeDict(({**self.env_vars,
                                                                        **locals(),
                                                                        **{'version': has_version.group(0)}})))
                    else:
                        th_file = string.Formatter().vformat(model_setup.get('th_file'), (), SafeDict(({**self.env_vars,
                                                                                                        **locals()})))
                    method = model_setup.get('method')
                    if 'flux_col' in model_setup.keys():
                        write_flux = True
                        sign_change = model_setup.get('sign_change')
                        perturbed_th, fluxes_out = self.format_schism_th(subout_dir, comp['model_input'], modcase_dir,
                                                                         crange, comp, method, th_file, cname, perturbed_th,
                                                                         fluxes_out=fluxes_out, fcol=model_setup[
                                                                             'flux_col'],
                                                                         enc_cp=cp, sign_change=sign_change)
                    else:
                        if 'tide' in cp:
                            perturbed_th, tidal_pert = self.format_schism_th(subout_dir, comp['model_input'], modcase_dir,
                                                                             crange, cdict, method, th_file, cname, perturbed_th,
                                                                             enc_cp=cp)
                        else:
                            perturbed_th = self.format_schism_th(subout_dir, comp['model_input'], modcase_dir,
                                                                 crange, cdict, method, th_file, cname, perturbed_th,
                                                                 enc_cp=cp)

            else:  # Single-element perturbation
                model_setup = self.configs[pdict['model_input']]
                th_file = string.Formatter().vformat(
                    model_setup.get('th_file'), (), SafeDict((self.env_vars)))
                method = model_setup.get('method')

                if 'flux_col' in model_setup.keys():
                    # need to handle this differently since it's all flux inputs in one file
                    write_flux = True
                    sign_change = model_setup.get('sign_change')
                    perturbed_th, fluxes_out = self.format_schism_th(case_data_dir, cp, modcase_dir,
                                                                     crange, pdict, method, th_file, cname, perturbed_th,
                                                                     fluxes_out=fluxes_out, fcol=model_setup[
                                                                         'flux_col'],
                                                                     sign_change=sign_change)
                else:
                    if 'tide' in cp:
                        perturbed_th, tidal_pert = self.format_schism_th(case_data_dir, cp, modcase_dir,
                                                                         crange, pdict, method, th_file, cname, perturbed_th)
                    else:
                        perturbed_th = self.format_schism_th(case_data_dir, cp, modcase_dir,
                                                             crange, pdict, method, th_file, cname, perturbed_th, enc_cp=cp)

        # creating and writing out flux file
        if write_flux:
            print('\t Writing modified flux file out....')
            flux_in = pd.read_csv(
                # Parse the datetime column
                self.flux_file_in, sep=' ', parse_dates=['datetime'],
                index_col='datetime'  # Set the datetime column as the index
            )

            for fcol in fluxes_out.keys():
                fdf = fluxes_out[fcol]
                fdf.columns = [fcol]
                flux_in[fcol] = np.nan
                fdf.index = pd.to_datetime(fdf.index)
                # fdf.index = fdf.index.strftime('%Y-%m-%d 00:00:00')
                # flux_in.loc['2012-01-01 00:00:00',:]
                flux_in.loc[flux_in.index.isin(fdf.index), fcol] = fdf[fcol]
            flux_in = flux_in.interpolate(axis=0, limit_direction='both')

            mod_fn = os.path.join(modcase_dir, f'flux.{cname}.dated.th')
            flux_in.index = flux_in.index.strftime('%Y-%m-%dT%H:%M')
            flux_in.to_csv(mod_fn, sep=' ', float_format="%.2f", index=True)
            clip_fn = mod_fn.replace('.dated.th', '.th')
            clip_file_to_start(mod_fn,
                               outpath=clip_fn,
                               start=dt.datetime(year=crange[0].year,
                                                 month=crange[0].month,
                                                 day=crange[0].day))
            perturbed_th[os.path.basename(self.flux_file_in)] = clip_fn

        print("Updating input files -----------------------------------------------")

        # Create param.nml files ------------------------------------------------------------------------------

        run_time_dict = {'{runtimedays}': str(runtimedays),
                         '{year_start}': str(case_start.year),
                         '{month_start}': str(case_start.month),
                         '{day_start}': str(case_start.day)}

        print(f"\t param.nml.tropic: {self.param_tropic_base}")
        fmt_string_file(self.param_tropic_base, os.path.join(
            modcase_dir, 'param.nml.tropic'), run_time_dict, method='replace')

        print(f"\t param.nml.clinic: {self.param_clinic_base}")
        fmt_string_file(self.param_clinic_base, os.path.join(
            modcase_dir, 'param.nml.clinic'), run_time_dict, method='replace')

        # copy interpolate_variables.in
        print(f"\t interpolate_variables.in")
        fmt_string_file(self.int_vars, os.path.join(
            modcase_dir, 'interpolate_variables.in'), run_time_dict, method='replace')

        # Create TH files ------------------------------------------------------------------------------------

        print('\t copying th files')
        th_files = []
        mod_th_files = dict()
        linked_th_file_strings = ''
        for th in perturbed_th.keys():
            th_files.append(perturbed_th[th])  # account for any modified files
            # link any modified files
            linked_th_file_strings += f'ln -sf {os.path.basename(perturbed_th[th])} {th}\n'
            mod_th_files[''.join(['{', f'{th}', '}'])] = perturbed_th[th]

        for fn in list(set(self.linked_th_files) - set(perturbed_th.keys())):
            mod_th_files[''.join(['{', f'{fn}', '}'])
                         ] = os.path.join(modcase_dir, fn)
            clip_fn = os.path.join(modcase_dir, fn)
            clip_file_to_start(f'{self.th_repo}/{fn}', outpath=clip_fn, start=dt.datetime(year=crange[0].year,
                                                                                          month=crange[0].month,
                                                                                          day=crange[0].day))
            th_files.append(clip_fn)  # account for any th files

        for fn in list(set(self.linked_ss_th_files) - set(perturbed_th.keys())):
            mod_th_files[''.join(['{', f'{fn}', '}'])
                         ] = os.path.join(modcase_dir, fn)
            clip_fn = os.path.join(modcase_dir, fn)
            clip_file_to_start(f'{self.dcd_repo}/{fn.replace(".th","_dated.th")}',  outpath=clip_fn, start=dt.datetime(year=crange[0].year,
                                                                                                                       month=crange[0].month,
                                                                                                                       day=crange[0].day))
            th_files.append(clip_fn)  # account for source/sink files

        for fn in self.versioned_th_files:
            linked_source = string.Formatter().vformat(
                fn['linked_fn'], (), SafeDict(({**self.env_vars, **locals()})))
            linked_fn = os.path.basename(linked_source)
            th_file = os.path.join(modcase_dir, linked_fn)
            shutil.copyfile(linked_source, th_file)
            th = fn['key']
            th_files.append(th_file)  # account for any modified files
            # link any modified files
            linked_th_file_strings += f'ln -sf {linked_fn} {th}\n'
            mod_th_files[''.join(['{', f'{th}', '}'])] = linked_fn

        # TROPIC BASH
        bash_tropic_dict = {**run_time_dict, **{'{year_end}': str(case_end.year),
                                                '{month_end}': str(case_end.month),
                                                '{day_end}': str(case_end.day),
                                                '{cname}': cname,
                                                '{linked_th_file_strings}': linked_th_file_strings},
                            '{tidal_pert}': tidal_pert}

        bash_case_tropic = os.path.join(modcase_dir, os.path.basename(
            self.bash_tropic).replace("CASENAME", cname))
        slurm_case_tropic = os.path.join(modcase_dir, os.path.basename(
            self.slurm_tropic).replace("CASENAME", cname))

        # CLINIC BASH
        bash_case_clinic = os.path.join(modcase_dir, os.path.basename(
            self.bash_clinic).replace("CASENAME", cname))
        slurm_case_clinic = os.path.join(modcase_dir, os.path.basename(
            self.slurm_clinic).replace("CASENAME", cname))

        # FINAL PACKAGING ===========================================================================

        for mesh_info in self.meshes:
            meshname = mesh_info['name']
            meshcase = f'{meshname}_{cname}'
            if meshcase not in self.meshcases_exclude and meshname not in self.meshes_exclude:

                # Create the mesh_case directory where the simulations will go
                print(
                    f"\tMesh: {meshname.upper()} -------------------------------------")
                # This is where simulations will be run
                meshcase_dir = os.path.join(
                    self.simulations_dir, f'{meshname}_{cname}')
                # This is where simulations will be run
                basemeshcase_dir = os.path.join(
                    self.simulations_dir, f'baseline_{cname}')
                if not os.path.exists(meshcase_dir):
                    os.mkdir(meshcase_dir)
                if not os.path.exists(os.path.join(meshcase_dir, 'outputs')):
                    os.mkdir(os.path.join(meshcase_dir, 'outputs'))

                for env_var in self.env_vars:
                    self.env_vars[env_var] = self.env_vars[env_var].format(
                        **self.env_vars)

                mesh_input_dir = string.Formatter().vformat(
                    mesh_info["indir"], (), SafeDict((self.env_vars)))

                # Copy case directory files
                shutil.copytree(modcase_dir, meshcase_dir, dirs_exist_ok=True)

                # make sflux links
                print('\t making sflux links')
                sflux_dir = os.path.join(meshcase_dir, 'sflux')
                if not os.path.exists(sflux_dir):
                    os.mkdir(sflux_dir)
                if self.machine.lower() == 'hpc5':
                    make_links(crange[0], crange[1], self.env_vars['sflux_src_dir'],
                               self.env_vars['sflux_narr_dir'], sflux_dir)
                elif self.machine.lower() == 'azure':
                    mk_links_file = self.inputs.get(
                        "make_links").format_map(self.env_vars)
                    fmt_string_file(mk_links_file,
                                    os.path.join(
                                        meshcase_dir, f'sflux/{os.path.basename(mk_links_file)}'),
                                    {**{"{start_date}": crange[0].strftime('%Y, %m, %d'),
                                        "{end_date}": crange[1].strftime('%Y, %m, %d'),
                                        "{src_dir}": os.path.relpath(self.env_vars['sflux_src_dir'], os.path.join(meshcase_dir, 'sflux')),
                                        "{src_dir_narr}": os.path.relpath(self.env_vars['sflux_narr_dir'], os.path.join(meshcase_dir, 'sflux')),
                                        "{link_dir}": './'}},
                                    method='replace')
                shutil.copyfile(os.path.join(self.env_vars['exp_dir'], 'sflux_inputs.txt'), os.path.join(
                    meshcase_dir, 'sflux/sflux_inputs.txt'))

                # copy spatial files to this dir
                print('\t copy spatial files')
                for gf in self.geometry_files:
                    infile = os.path.join(mesh_input_dir.format(**self.env_vars),
                                          gf.format_map({'case_year': str(case_year)}))
                    if not os.path.exists(infile):
                        raise ValueError(
                            f"infile: {infile} does not exist!\n Check geometry_files list in yaml")
                    else:
                        inbase = os.path.basename(infile)
                        shutil.copyfile(infile,
                                        os.path.join(meshcase_dir,
                                                     inbase))

                # copy th files to this dir
                print('\t copy time history files')
                for th in th_files:
                    shutil.copyfile(th, os.path.join(
                        meshcase_dir, os.path.basename(th)))

                # copy param files
                print('\t copy param files')
                shutil.copyfile(os.path.join(modcase_dir, 'param.nml.tropic'), os.path.join(
                    meshcase_dir, 'param.nml.tropic'))
                shutil.copyfile(os.path.join(modcase_dir, 'param.nml.clinic'), os.path.join(
                    meshcase_dir, 'param.nml.clinic'))

                # copy common files
                print('\t copy common files')
                for cf in self.common_files:
                    cff = string.Formatter().vformat(
                        cf, (), SafeDict(({**self.env_vars, **locals()})))
                    shutil.copyfile(cff, os.path.join(
                        meshcase_dir, os.path.basename(cff)))

                # TROPIC ----------------------------------------------
                print(f"\t\t Handling the tropic inputs")

                bash_meshcase_tropic = os.path.join(meshcase_dir, os.path.basename(
                    bash_case_tropic).replace("MESHNAME", meshname))
                fmt_string_file(self.bash_tropic,
                                bash_meshcase_tropic,
                                {**bash_tropic_dict, **{"{meshname}": meshname}},
                                method='replace')
                print(f"\t\t\t tropic.sh: {bash_meshcase_tropic}")
                self.bash_file_dict['tropic'].append(bash_meshcase_tropic)

                slurm_tropic_dict = {**locals(), **{'job_name': self.job_name,
                                                    'baro': 'tropic',
                                                    'output_log_file_base': self.output_log_file_base.format_map(locals())}}

                slurm_meshcase_tropic = os.path.join(meshcase_dir, os.path.basename(
                    slurm_case_tropic).replace("MESHNAME", meshname))
                if self.machine.lower() == 'hpc5':
                    fmt_string_file(self.slurm_tropic,
                                    slurm_meshcase_tropic,
                                    slurm_tropic_dict,
                                    method='format_map')

                # write out a yaml of modified th files
                mod_th_case = os.path.join(meshcase_dir, os.path.basename(
                    self.mod_th_dict).replace("CASENAME", cname))
                fmt_string_file(self.mod_th_dict,
                                mod_th_case,
                                mod_th_files,
                                method='replace')

                # CLINIC ----------------------------------------------
                print(f"\t\t Handling the clinic inputs")

                bash_clinic_dict = {**bash_tropic_dict, **{'{year_end}': str(case_end.year),
                                                           '{month_end}': str(case_end.month),
                                                           '{day_end}': str(case_end.day),
                                                           '{cname}': cname,
                                                           '{linked_th_file_strings}': linked_th_file_strings,
                                                           '{mesh_input_dir}': os.path.relpath(mesh_input_dir, meshcase_dir),
                                                           '{case_year}': str(case_year),
                                                           '{meshname}': meshname}}
                bash_meshcase_clinic = os.path.join(meshcase_dir, os.path.basename(
                    bash_case_clinic).replace("MESHNAME", meshname))
                fmt_string_file(self.bash_clinic,
                                bash_meshcase_clinic,
                                bash_clinic_dict,
                                method='replace')
                print(f"\t\t\t clinic.sh: {bash_meshcase_clinic}")
                self.bash_file_dict['clinic'].append(bash_meshcase_clinic)

                if self.machine.lower() == 'hpc5':
                    slurm_clinic_dict = {**locals(), **{'job_name': self.job_name,
                                                        'baro': 'clinic',
                                                        'output_log_file_base': self.output_log_file_base.format_map(locals())}}

                    slurm_meshcase_clinic = os.path.join(meshcase_dir, os.path.basename(
                        slurm_case_clinic).replace("MESHNAME", meshname))
                    fmt_string_file(self.slurm_clinic,
                                    slurm_meshcase_clinic,
                                    slurm_clinic_dict,
                                    method='format_map')
                elif self.machine.lower() == 'azure':
                    az_dict = {**locals(),
                               **self.env_vars,
                               **{'rel_study_dir': f"{os.path.basename(self.env_vars['exp_dir'])}/{os.path.relpath(meshcase_dir,self.env_vars['exp_dir'])}",
                                  'rel_study_base_dir': f"{os.path.basename(self.env_vars['exp_dir'])}/{os.path.relpath(basemeshcase_dir,self.env_vars['exp_dir'])}",
                                  'job_name': self.job_name}}
                    az_yml_meshcase = os.path.join(self.az_yml_dir,
                                                   os.path.basename(self.az_yml_file).replace("MESHNAME", meshname))
                    az_yml_meshcase = az_yml_meshcase.replace(
                        "CASENAME", cname)
                    fmt_string_file(self.az_yml_file,
                                    az_yml_meshcase,
                                    az_dict,
                                    method='format_map')
                else:
                    raise ValueError(
                        f"Machine needs to be defined as 'hpc5' or 'azure' to properly format bash files.")

    def set_schism_vars(self):
        self.meshes = self.inputs.get('meshes')
        self.noaa_download = build_dict(self.inputs.get('noaa_download'))
        self.noaa_pt_reyes = self.noaa_download.get(
            'pt_reyes').format(**self.env_vars)
        self.noaa_monterey = self.noaa_download.get(
            'monterey').format(**self.env_vars)
        self.param_tropic_base = self.inputs.get(
            'param_tropic_base').format(**self.env_vars)
        self.param_clinic_base = self.inputs.get(
            'param_clinic_base').format(**self.env_vars)
        self.bash_tropic = self.inputs.get(
            'bash_tropic').format(**self.env_vars)
        self.bash_clinic = self.inputs.get(
            'bash_clinic').format(**self.env_vars)
        self.int_vars = self.inputs.get('int_vars').format(**self.env_vars)
        self.flux_file_in = self.inputs.get(
            'flux_file_in').format(**self.env_vars)
        self.th_repo = self.env_vars['th_repo']
        self.dcd_repo = self.env_vars['dcd_repo']
        self.slurm_tropic = self.inputs.get(
            'slurm_tropic').format(**self.env_vars)
        self.slurm_clinic = self.inputs.get(
            'slurm_clinic').format(**self.env_vars)
        if self.machine.lower() == 'azure':
            self.az_yml_file = self.inputs.get(
                'az_yml_file').format(**self.env_vars)
            self.az_yml_dir = self.inputs.get('az_yml_dir').format(**self.env_vars)
        self.job_name = self.inputs.get('job_name')
        self.output_log_file_base = self.inputs.get('output_log_file_base')
        self.geometry_files = self.inputs.get('geometry_files')
        self.common_files = [string.Formatter().vformat(cf, (), SafeDict(
            (self.env_vars))) for cf in self.inputs.get('common_files')]
        self.mod_th_dict = string.Formatter().vformat(
            self.inputs.get('mod_th_dict'), (), SafeDict((self.env_vars)))
        self.versioned_th_files = self.inputs.get('versioned_th_files')

    def format_schism_th(self, case_data_dir, cp, modcase_dir, crange, pdict, method, th_file, cname, perturbed_th,
                         fluxes_out=None, fcol=None, enc_cp=None, sign_change=None):
        tidal_pert = ''
        # Handle the different types of perturbation
        if 'tide' in cp:
            print('\t\t Modifying tidal boundary')
            if pdict['method'] == 'shift':
                if 'shift_forward' in pdict['args'].keys():
                    timedelt = -pdict['args']['shift_forward']
                    tidal_pert = '_shift_forward'
                elif 'shift_backward' in pdict['args'].keys():
                    timedelt = pdict['args']['shift_backward']
                    tidal_pert = '_shift_backward'

                pt_reyes = read_noaa(self.noaa_pt_reyes, force_regular=True)
                pt_reyes.index = pt_reyes.index - pd.Timedelta(days=timedelt)
                write_noaa(self.noaa_pt_reyes, pt_reyes, os.path.join(
                    modcase_dir, f'pt_reyes{tidal_pert}.csv'))

                monterey = read_noaa(self.noaa_monterey, force_regular=True)
                monterey.index = monterey.index - pd.Timedelta(days=timedelt)
                write_noaa(self.noaa_monterey, monterey, os.path.join(
                    modcase_dir, f'monterey{tidal_pert}.csv'))
            else:
                raise ValueError(
                    f"Unknown method to modify tidal boundary: {method}")

            return perturbed_th, tidal_pert

        elif any(gate in cp for gate in ['dcc']):
            print('\t\t Modifying Gate Operations')
            # get modified operation. th_file is the model/framework to go off of. mod_file is the data to draw from
            mod_file = os.path.join(case_data_dir,
                                    (f"{pdict['model_input']}_{pdict['method']}"
                                     f"_{dt.datetime.strftime(crange[0], format='%Y%m%d0000')}"
                                     f"-{dt.datetime.strftime(crange[1]+dt.timedelta(days=1), format='%Y%m%d0000')}.csv"))
            dat_in = pd.read_csv(mod_file, parse_dates=[
                                 0], index_col=[0], header=None)
            dat_in.index = dat_in.index.strftime('%Y-%m-%dT%H:%M')
            dat_in.index.name = 'datetime'

            th_head = pd.read_csv(th_file, sep='\s+', nrows=1, header=0)
            schema = {'datetime': str, 'install': int, 'ndup': int, 'op_down': np.float64,
                      'op_up': np.float64, 'elev': np.float64, 'width': np.float64, 'height': int}
            df = pd.DataFrame(columns=schema.keys()).astype(schema)
            df['datetime'] = dat_in.index
            for header in th_head.columns[1:-1]:
                df[header] = th_head[header][0]

            if pdict['method'] == 'set':
                dat_out = dat_in.copy()
                if float(pdict['args']['set_value']) == 2.0:
                    dat_out[1] = 10  # schism has height of 0/10
                elif float(pdict['args']['set_value']) == 0:
                    dat_out[1] = 0  # schism has height of 0/10
            elif pdict['method'] == 'read':
                dat_out = dat_in.copy()*10  # the input is just 0 -> 1 and it needs to be 0 ->10

            df['height'] = dat_out[1].to_list()

            # write out dataframe
            dated_fn = os.path.join(modcase_dir, os.path.basename(
                th_file).replace('.th', f'.{cname}.dated.th'))
            df.to_csv(dated_fn, sep=' ', index=False, float_format="%.2f")

            # clip to model period
            clip_fn = dated_fn.replace('.dated.th', '.th')
            clip_file_to_start(dated_fn,
                               outpath=clip_fn,
                               start=dt.datetime(year=crange[0].year,
                                                 month=crange[0].month,
                                                 day=crange[0].day))
            # add to list of files that are perturbed
            perturbed_th[os.path.basename(th_file)] = clip_fn

            return perturbed_th

        elif 'suis' in cp:
            print(f'\t\t Copying Suisun gate operations')
            # Determine version
            gate_ver = enc_cp[-1]
            # Suisun needs to copy the three gates (radial, boatlock, and flashboard)
            for gate in ['radial', 'boat_lock', 'flash']:
                th_file_in = th_file.format_map(
                    locals())  # uses gate and gate_ver
                shutil.copyfile(th_file_in, os.path.join(
                    modcase_dir, os.path.basename(th_file_in)))
                clip_fn = os.path.join(
                    modcase_dir, f'{os.path.basename(th_file_in).split("_")[0]}_{gate}.{cname}.th')
                # print(f'{th_file_in} to {clip_fn}')
                clip_file_to_start(th_file_in,
                                   outpath=clip_fn,
                                   start=dt.datetime(year=crange[0].year,
                                                     month=crange[0].month,
                                                     day=crange[0].day))
                # add to list of files that are perturbed
                perturbed_th[f'{os.path.basename(th_file_in).split("_")[0]}_{gate}.th'] = clip_fn

            return perturbed_th

        elif 'dcd' in cp:
            print(
                f'\t\t Copying consumptive use file: {os.path.basename(th_file)}')
            shutil.copyfile(th_file, os.path.join(
                modcase_dir, os.path.basename(th_file)))
            clip_fn = os.path.join(
                modcase_dir, f'{os.path.basename(th_file).split("_")[0]}.{cname}.th')
            clip_file_to_start(th_file,
                               outpath=clip_fn,
                               start=dt.datetime(year=crange[0].year,
                                                 month=crange[0].month,
                                                 day=crange[0].day))
            # add to list of files that are perturbed
            perturbed_th[f'{os.path.basename(th_file).split("_")[0]}.th'] = clip_fn

            return perturbed_th

        elif fcol:
            print('\t\t Adding to flux file')
            mod_file = os.path.join(case_data_dir,
                                    (f"{pdict['model_input']}_{pdict['method']}"
                                     f"_{dt.datetime.strftime(crange[0], format='%Y%m%d0000')}"
                                     f"-{dt.datetime.strftime(crange[1]+dt.timedelta(days=1), format='%Y%m%d0000')}.csv"))
            dat_in = pd.read_csv(mod_file, parse_dates=[
                                 0], index_col=[0], header=None)
            dat_in.index = dat_in.index.strftime('%Y-%m-%dT%H:%M')
            dat_in.index.name = 'datetime'
            # sign_change: for inflow (eg Sac or SJR) this is negative, for exports (eg ccr) this is positive. Defined in yml
            dat_in[1] = sign_change * \
                dat_in[1].div(35.3147)  # convert to metric
            fluxes_out[fcol] = dat_in

            return perturbed_th, fluxes_out

        else:
            raise ValueError(f"cp method is not defined: {cp}")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # Read in param.tropic and param.clinic and modify
    if False:  # run for azure
        yml_fname = "./input/schism_lathypcub_v3_azure.yaml"

        mbc = ModelBCGen(yml_fname, 'schism', machine='azure')

    elif True:  # run for azure SLR
        yml_fname = "./input/schism_slr_lathypcub_v3_azure.yaml"

        mbc = ModelBCGen(yml_fname, 'schism', machine='hpc5')

    elif False:  # run for azure
        yml_fname = "./input/schism_lathypcub_v3_azure_from_base.yaml"

        mbc = ModelBCGen(yml_fname, 'schism', machine='azure')

    else:
        yml_fname = "./input/schism_lathypcub_v3_hpc5.yaml"

        mbc = ModelBCGen(yml_fname, 'schism', machine='hpc5')
