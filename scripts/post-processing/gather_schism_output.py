#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to gather input/output data from SCHISM for training ANN


"""

import os
from schimpy.station import *
from schimpy.model_time import file_to_timestamp
from vtools.functions.filter import cosine_lanczos
import datetime
import warnings
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
def read_th(infile):

    in_df = pd.read_table(infile,
                          sep='\s+',
                          index_col='datetime',
                          comment="#")
    in_df.index = pd.to_datetime(in_df.index,format="%Y-%m-%dT%H:%M")

    # monotonic increase check
    mon_inc = all(x<y for x, y in zip(in_df.index, in_df.index[1:])) # True for monotonically-increasing data
    if not mon_inc: print(in_df.loc[in_df.index.to_series().diff() < pd.to_timedelta('0 seconds')]) # prints the row(s) where monotonicity is broken

    in_df = in_df.reindex(
        pd.date_range(start=in_df.index.min(),
                      end=in_df.index.max(),
                      freq='D'),
                      method='ffill')
    
    print(f"TH FILE ==== {infile} min: {in_df.index.min()} max:{in_df.index.max()}")
    return in_df

def get_ts_from_th(infile, start, elapsed_unit="s"):
    if elapsed_unit == "s":
        elapsed_fac = 1.
    elif elapsed_unit == "d":
        elapsed_fac = 24 * 3600
    else:
        raise ValueError("elapsed_unit must be 's' or 'd'")
    if type(start) == str:
        start = datetime.datetime(*list(map(int, re.split('[^\d]', start))))

    if type(infile) is str:
        in_df = pd.read_table(infile,
                            sep='\s+',
                            comment="#",
                            index_col=0,
                            header=None)
    else:
        in_df = infile
        
    out_df = in_df.copy()
    out_df.index = pd.to_datetime(start) + pd.to_timedelta(round(in_df.index.to_series() * elapsed_fac), unit='s')

    return out_df

# define header values
def get_headers(infile, no_index=True):
    with open(infile, 'r') as headin:
        head = headin.readline().strip()
        headers = head.split()

    if no_index:
        headers = headers[1:]

    return headers

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
        self.in_vars = self.inputs.get('in_vars')
        self.comb_in_vars = self.inputs.get('comb_in_vars')
        self.out_vars = self.inputs.get('out_vars')

        # set station inputs and flow transect inputs
        # self.station_fpath = string.Formatter().vformat(self.inputs['station_in'],(),
        #                                                 SafeDict(({**self.env_vars, 
        #                                                            **locals()})))
        # self.station_df = read_station_in(self.station_fpath)

        flux_fpath = string.Formatter().vformat(self.inputs['flux_in'],(),
                                                SafeDict(({**self.env_vars, 
                                                           **locals()})))
        self.flux_df = flux_stations_from_yaml(flux_fpath)


        self.out_dir = string.Formatter().vformat(self.inputs.get('out_dir'),(),SafeDict((self.env_vars)))
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)

        # get names of th files
        self.mod_th_dict = self.inputs.get('mod_th_dict')
        
        if 'meshes' in self.inputs.keys():
            self.meshes = self.inputs.get('meshes')
        else:
            raise ValueError("Need to define the meshes using key 'meshes' in the yaml file")
    
    def run_all_cases(self):
        if 'cases' in self.inputs.keys():
            cases = self.inputs.get('cases')
            
            for case in cases:
                case_num = case.get('case_num')
                # cname = case.get('name')

                for mesh in self.meshes:

                    self.get_meshcase_inouts(self, mesh, case_num)

        else:

            raise ValueError("Need to define the cases using key 'cases' in the yaml file")

    def get_meshcase_inouts(self, mesh, case_num):
        meshcase_dir = string.Formatter().vformat(self.env_vars['mc_dir'],(),
                                                  SafeDict(({**self.env_vars, 
                                                             **locals()})))
        with open(string.Formatter().vformat(self.mod_th_dict,(),
                                                        SafeDict(({**self.env_vars, 
                                                                    **locals()}))), 'r') as f:
                mesh_mod_th_dict = schism_yaml.load(f) 
                mesh_mod_th_dict = build_dict(mesh_mod_th_dict)
        
        outputs_fpath = string.Formatter().vformat(self.inputs['station_output'],(),
                                                   SafeDict(({**self.env_vars, 
                                                              **locals()})))
        param_fpath = string.Formatter().vformat(self.inputs['param_clinic'],(),
                                                 SafeDict(({**self.env_vars, 
                                                            **locals()})))
        station_inpath = os.path.join(meshcase_dir,'station.in')
        time_basis = get_start_date_from_param(param_fpath)

        # Inputs ---------------------------------------------------------------------------------------
        invar_df = pd.DataFrame()
        for invar in self.in_vars:
            
            in_name = invar['name']
            if invar['method'] == "read_multiple_column_th": # th_file, time_basis, th_header, inputs 
                invar_df[in_name] = self.read_multiple_column_th(invar['th_files'][0].format_map({**self.env_vars,
                                                                                    **locals(),
                                                                                    **{"case_dir": meshcase_dir},
                                                                                    **mesh_mod_th_dict}),
                                                        time_basis,
                                                        invar['th_header'].format_map({**self.env_vars,
                                                                                    **locals(),
                                                                                    **{"case_dir": meshcase_dir},
                                                                                    **mesh_mod_th_dict}), 
                                                        invar['inputs'])
                
            elif invar['method'] == "read_single_column_th": # th_file, time_basis, th_header, in_col
                
                invar_df[in_name] = self.read_single_column_th(invar['th_files'][0].format_map({**self.env_vars,
                                                                                    **locals(),
                                                                                    **{"meshcase_dir": meshcase_dir},
                                                                                    **mesh_mod_th_dict}),
                                                        time_basis,
                                                        invar['th_header'].format_map({**self.env_vars,
                                                                                    **locals(),
                                                                                    **{"meshcase_dir": meshcase_dir},
                                                                                    **mesh_mod_th_dict}), 
                                                        invar['inputs'])
                
            elif invar['method'] == "calc_dcu":  # th_files, time_basis
                
                invar_df[in_name] = self.calc_dcu([invar['th_files'][0].format_map({**self.env_vars,
                                                                                    **locals(),
                                                                                    **{"meshcase_dir": meshcase_dir},
                                                                                    **mesh_mod_th_dict}),
                                                    invar['th_files'][1].format_map({**self.env_vars,
                                                                                    **locals(),
                                                                                    **{"meshcase_dir": meshcase_dir},
                                                                                    **mesh_mod_th_dict})],
                                                    time_basis)

            elif invar['method'] == "calc_tidal_energy": # outputs_fpath, time_basis, loc
                
                invar_df[in_name] = self.calc_tidal_energy(invar['station_output'].format_map({**self.env_vars,
                                                                                                **locals(),
                                                                                                **{"meshcase_dir": meshcase_dir},
                                                                                                **mesh_mod_th_dict}),
                                                            station_inpath,
                                                            time_basis,
                                                            invar['loc'])

            else:
                raise ValueError(f'There is no current process/function for {invar["method"]}')

        print('Done collecting inputs')

        # Combined Inputs ---------------------------------
        combinvar_df = pd.DataFrame()
        for combinvar in self.comb_in_vars:
            
            in_name = combinvar['name']
            combtemp = invar_df.loc[:,combinvar['vars']].copy()
            combinvar_df[in_name] = combtemp.sum(axis=1)

        # Outputs -----------------------------------------
        outvar_df = pd.DataFrame()
        for outvar in self.out_vars:

            in_name = outvar['name']
            if outvar['method'] == "assign_multiple_ec":
                
                ec_out_df = self.get_multiple_staout(outputs_fpath, station_inpath, time_basis, outvar['locs'])

                for col in ec_out_df.columns:
                    outvar_df[f"{col} EC"] = ec_out_df.loc[:,col]
                
            elif outvar['method'] == "assign_multiple_wse":
                
                wse_out_df = self.get_multiple_staout(outputs_fpath, station_inpath, time_basis, outvar['locs'])

                for col in wse_out_df.columns:
                    outvar_df[f"{col} WSE"] = wse_out_df.loc[:,col]

            elif outvar['method'] == "read_single_flux": # outputs_fpath, time_basis, loc TODO: check inputs
                
                outvar_df[in_name] = self.read_single_flux(outvar['flux_output'].format_map({**self.env_vars,
                                                                                                **locals(),
                                                                                                **{"meshcase_dir": meshcase_dir},
                                                                                                **mesh_mod_th_dict}),
                                                            time_basis,
                                                            outvar['loc'])
                
            else:
                raise ValueError(f'There is no current process/function for {outvar["method"]}')

        # Total output -------------------------------------
        tot_df = pd.merge(invar_df, combinvar_df, left_index=True, right_index=True, how='inner')
        tot_df = pd.merge(tot_df, outvar_df, left_index=True, right_index=True, how='inner')
        tot_df.to_csv(os.path.join(self.out_dir,f'{mesh}_lhc_{case_num}.csv'), index=True)

    def read_multiple_column_th(self, th_file, time_basis, th_header, inputs):
        
        # read single th file and add specified columns together
        th_df = get_ts_from_th(th_file, time_basis)
        th_df.columns = get_headers(th_header, no_index=True)
        in_cols = [i.split('-')[-1] for i in inputs]
        mults = [-1 if i.split('-')[0]=='' else 1 for i in inputs]

        if not set(in_cols).issubset(th_df.columns):
            raise ValueError("All inputs need to be found in the th_header to index outputs")
        
        th_df['total'] = 0

        for in_col, mult in zip(in_cols, mults):
            th_df['total'] = th_df['total'] + mult * th_df[in_col]

        out_df = th_df['total']

        return out_df
    
    def read_single_column_th(self, th_file, time_basis, th_header, in_col):
        
        # read single th file and a single column
        th_df = get_ts_from_th(th_file, time_basis)
        th_df.columns = get_headers(th_header)
        if isinstance(in_col, list) and len(in_col) == 1:
            in_col = in_col[0]

        col = in_col.split('-')[-1]
        mult = -1 if in_col.split('-')[0]=='' else 1

        if not col in th_df.columns:
            raise ValueError("in_col needs to be found in the th_header to index output")

        out_df = mult * th_df[col]
        out_df = out_df.resample('15min').ffill()

        return out_df
    
    def calc_dcu(self, th_files, time_basis):
        warnings.warn("'calc_dcu' method assumes an order to th_files of 1) source 2) sink")
        src = get_ts_from_th(th_files[0], time_basis)
        src = src.sum(axis=1)
        sink = get_ts_from_th(th_files[1], time_basis)
        sink = sink.sum(axis=1)

        net = src + sink
        net = net.resample('15min').ffill()

        return net
    
    def calc_tidal_energy(self, outputs_fpath, station_inpath, time_basis, loc):

        wse_ts = self.get_single_staout(outputs_fpath, station_inpath, time_basis, loc)

        if isinstance(wse_ts.index, pd.core.indexes.datetimes.DatetimeIndex):
            print('timeseries is inst-val, converting to per-aver')
            wse_ts.index = wse_ts.index.to_period()

        df_nrg = cosine_lanczos((wse_ts-cosine_lanczos(wse_ts.copy(), 
                                                       cutoff_period ='40D', padtype='odd'))**2, 
                                cutoff_period ='40D', padtype='odd') # = < (z- <z>)^2 >
        df_tidal_energy = df_nrg.resample('D', closed='right').mean()
        df_tidal_energy.columns=['tidal_energy']
        df_tidal_energy.index = df_tidal_energy.index.to_timestamp()
        df_tidal_energy = df_tidal_energy.resample('15min').ffill()

        return df_tidal_energy
    
    def get_single_staout(self, outputs_fpath, station_inpath, time_basis, loc, depth='upper'):

        all_ts = read_staout(outputs_fpath, station_inpath, time_basis)
        station_df = read_station_in(station_inpath)

        if isinstance(loc, list):
            loc = loc[0]

        if depth == 'upper':
            loc_outs = [stc for stc in station_df.index if loc in stc]
            if len(loc_outs) == 0:
                raise ValueError(f"There is no station output for {loc}")
            elif len([stc for stc in loc_outs if 'upper' in stc]) > 0:
                col = f'{loc}_upper'
            else:
                col = f'{loc}_default'
        else:
            raise ValueError('Currently only programmed to get the upper layer by default')

        out_ts = all_ts[col]

        return out_ts

    def get_multiple_staout(self, outputs_fpath, station_inpath, time_basis, locs):
        
        station_df = read_station_in(station_inpath)
        all_ts = read_staout(outputs_fpath, station_df, time_basis)

        cols = []
        for loc in locs:
            loc_outs = [stc for stc in station_df.index if loc in stc]
            if len(loc_outs) == 0:
                raise ValueError(f"There is no station output for {loc}")
            elif len([stc for stc in loc_outs if 'upper' in stc]) > 0:
                cols.append(f'{loc}_upper')
            else:
                cols.append(f'{loc}_default')

        out_ts = all_ts[cols] # this is the output df for all requested locs

        return out_ts

    def read_single_flux(self, outputs_fpath, time_basis, loc):

        if not loc[0] in self.flux_df:
            raise ValueError(f"Flux out location '{loc[0]}' needs to be found in the flow_station_xsects yaml file and is not")
        
        flux_ts = read_flux_out(outputs_fpath, self.flux_df, time_basis)
        flux_ts = flux_ts.loc[:,loc]
        flux_ts = flux_ts.resample('15min').ffill()

        return flux_ts

if __name__ == '__main__':

    from schimpy import schism_yaml
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "./input/pull_output_lathypcub_v3_schism.yaml"
    annbc = ANNBCECGen(in_fname, model_type="SCHISM")
    annbc.get_meshcase_inouts('suisun_basebc', 4)