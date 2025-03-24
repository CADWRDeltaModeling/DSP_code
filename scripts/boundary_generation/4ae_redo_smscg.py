import pandas as pd
import os
import pyhecdss

def get_pathname(dss_filename, b_part, c_part, e_part=None, f_part=None, filter_b_part_numeric=None):
    with pyhecdss.DSSFile(dss_filename) as d:
        catdf = d.read_catalog()
        dss_file_parts = dss_filename.split('/')
        dfilename = dss_file_parts[len(dss_file_parts)-1]
        filtered_df = None
        if b_part is not None:
            filtered_df = filtered_df[(
                filtered_df.B == b_part)] if filtered_df is not None else catdf[(catdf.B == b_part)]
        if c_part is not None:
            filtered_df = filtered_df[(
                filtered_df.C == c_part)] if filtered_df is not None else catdf[(catdf.C == c_part)]
        if e_part is not None:
            filtered_df = filtered_df[(
                filtered_df.E == e_part)] if filtered_df is not None else catdf[(catdf.E == e_part)]
        if f_part is not None:
            filtered_df = filtered_df[(
                filtered_df.F == f_part)] if filtered_df is not None else catdf[(catdf.F == f_part)]
        if filter_b_part_numeric:
            filtered_df = filtered_df[(filtered_df.B.str.isnumeric())]
        path_list = d.get_pathnames(filtered_df)

    return path_list

def schism_to_dsm2(cases, th_in, dss_fn_in, dss_fn_out):
    for case in cases:
        cname = f"lhc_{case}"
        print(cname)

        for sch_gate, dsm2_gate in zip(['boat_lock','flash','radial'], ['BOATLOCK','FLASHBOARD','RADIAL']):
            print(sch_gate)
            th_read = th_in.format(**locals())
            in_df = pd.read_table(th_read,
                                sep='\s+',
                                comment="#",
                                index_col=0,
                                parse_dates=['datetime'])


            dss_fn_in_fmt = dss_fn_in.format(**locals())
            dss_fn_out_fmt = dss_fn_out.format(**locals())
            if not os.path.exists(dss_fn_out_fmt):
                dumdss = pyhecdss.DSSFile(dss_fn_out_fmt, create_new=True)
            pathname = get_pathname(dss_fn_in_fmt, 'MTZSL', c_part=f"{dsm2_gate}_OP")[0]
            
            paths_out = []

            if dsm2_gate == "RADIAL":
                # also change radial_fract_to and radial_fract_from
                path_from = pathname.replace('RADIAL_OP','RADIAL_FRACT_FROM') # op_up
                pathname_from = get_pathname(dss_fn_in_fmt, 'MTZSL', c_part='RADIAL_FRACT_FROM')
                pathname_from = [p for p in pathname_from if "DSP_LHC" not in p][0]
                paths_out.append(pathname_from)
                path_to = pathname.replace('RADIAL_OP','RADIAL_FRACT_TO') # op_down
                pathname_to = get_pathname(dss_fn_in_fmt, 'MTZSL', c_part='RADIAL_FRACT_TO')
                pathname_to = [p for p in pathname_to if "DSP_LHC" not in p][0]
                paths_out.append(pathname_to)

                with pyhecdss.DSSFile(dss_fn_in_fmt) as d:
                    out_df, units, ptype = d.read_its(pathname_from)
                    # add historical timeseries to beginning and end using margin_ts (copy of df)
                    margin_ts = out_df.copy()
                    margin_ts = margin_ts[(margin_ts.index < in_df.index[0])
                                        | (margin_ts.index > in_df.index[-1])]
                    margin_ts.columns = ['op_up']
                    out_df = pd.concat([in_df['op_up'], margin_ts['op_up']]) # out_df.loc[in_df.index[0]:in_df.index[-1]]
                    out_df.sort_index(inplace=True)

                    to_df, units, ptype = d.read_its(pathname_to)
                    # add historical timeseries to beginning and end using margin_ts (copy of df)
                    margin_ts = to_df.copy()
                    margin_ts = margin_ts[(margin_ts.index < in_df.index[0])
                                        | (margin_ts.index > in_df.index[-1])]
                    margin_ts.columns = ['op_down']
                    to_df = pd.concat([in_df['op_down'], margin_ts['op_down']])
                    to_df.sort_index(inplace=True)

                with pyhecdss.DSSFile(dss_fn_out_fmt) as d:
                    d.write_its(path_from, out_df, units, ptype)
                    paths_out.append(path_from)
                    d.write_its(path_to, to_df, units, ptype)
                    paths_out.append(path_to)

                in_df['op_up'] = in_df['op_up'].where(in_df['op_up'] == 1, -10)

            with pyhecdss.DSSFile(dss_fn_in_fmt) as d:
                op_df, units, ptype = d.read_its(pathname)
                # add historical timeseries to beginning and end using margin_ts (copy of df)
                margin_ts = op_df.copy()
                margin_ts = margin_ts[(margin_ts.index < in_df.index[0])
                                    | (margin_ts.index > in_df.index[-1])]
                margin_ts.columns = ['op_up']
                op_df = pd.concat([in_df['op_up'], margin_ts['op_up']]) # op_df.loc[in_df.index[0]:in_df.index[-1]]
                op_df.sort_index(inplace=True)
            
            with pyhecdss.DSSFile(dss_fn_out_fmt) as d:
                d.write_its(pathname, op_df, units, ptype)
                paths_out.append(pathname)

        # write remaining data to DSS files
        with pyhecdss.DSSFile(dss_fn_in_fmt) as d_in:
            with pyhecdss.DSSFile(dss_fn_out_fmt) as d_out:
                incat = d_in.read_catalog()  # all the pathnames in the DSS file
                paths_in = d_in.get_pathnames(incat)

                # paths in the input DSS that haven't been copied out
                missing_paths = [p for p in paths_in if "MTZSL" not in p]

                for p in missing_paths:
                    df = None
                    units = None
                    ptype = None
                    # print(f'Writing out {p} units {units} pytype {ptype}')
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

def dsm2_smscg_from_to(cases, dss_fn_in, dss_fn_out):
    for case in cases:
        cname = f"lhc_{case}"
        print(cname)

        dsm2_gate = "RADIAL"

        dss_fn_in_fmt = dss_fn_in.format(**locals())
        dss_fn_out_fmt = dss_fn_out.format(**locals())
        if not os.path.exists(dss_fn_out_fmt):
            dumdss = pyhecdss.DSSFile(dss_fn_out_fmt, create_new=True)
        pathname = get_pathname(dss_fn_in_fmt, 'MTZSL', c_part=f"{dsm2_gate}_OP")[0]
    
        # change radial_fract_to and radial_fract_from
        path_from = pathname.replace('RADIAL_OP','RADIAL_FRACT_FROM') # op_up
        pathname_from = get_pathname(dss_fn_in_fmt, 'MTZSL', c_part='RADIAL_FRACT_FROM')
        pathname_from = [p for p in pathname_from if "DSP_LHC" not in p][0]
        
        path_to = pathname.replace('RADIAL_OP','RADIAL_FRACT_TO') # op_down
        pathname_to = get_pathname(dss_fn_in_fmt, 'MTZSL', c_part='RADIAL_FRACT_TO')
        pathname_to = [p for p in pathname_to if "DSP_LHC" not in p][0]
        
        with pyhecdss.DSSFile(dss_fn_in_fmt) as d:
            in_df, units, ptype = d.read_its(pathname)
            ts_df, units, ptype = d.read_its(pathname_from)
            # make sure RADIAL_FRACT_FROM is 0 when RADIAL_OP is -10 and 1 when 1
            from_df = in_df.copy()
            from_df = in_df.where(in_df > 0, 0) # keep value if greater than 0 else, change -10 to 0

            to_df, units, ptype = d.read_its(pathname_to)
            # just make RADIAL_FRACT_TO always = 1
            to_df[:] = 1
        
        with pyhecdss.DSSFile(dss_fn_out_fmt) as d:
            d.write_its(path_to, to_df, units, ptype)
            d.write_its(path_from, from_df, units, ptype)

        # write remaining data to DSS files
        with pyhecdss.DSSFile(dss_fn_in_fmt) as d_in:
            with pyhecdss.DSSFile(dss_fn_out_fmt) as d_out:
                incat = d_in.read_catalog()  # all the pathnames in the DSS file
                paths_in = d_in.get_pathnames(incat)

                # paths in the input DSS that haven't been copied out
                keywords = ["MTZSL", "RADIAL_FRACT_TO"]
                filtered_paths = [p for p in paths_in if not all(k in p for k in keywords)]
                keywords = ["MTZSL", "RADIAL_FRACT_FROM"]
                missing_paths = [p for p in filtered_paths if not all(k in p for k in keywords)]

                for p in missing_paths:
                    df = None
                    units = None
                    ptype = None
                    # print(f'Writing out {p} units {units} pytype {ptype}')
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

os.chdir(os.path.dirname(os.path.abspath(__file__)))


## SCHISM to DSM2 modification
# cases = range(1,8)
# dss_fn_in = "../../model/dsm2/DSP_DSM2_202307/latinhypercube_v3/timeseries/old_gates/{cname}_gates.dss"
# dss_fn_out = "../../model/dsm2/DSP_DSM2_202307/latinhypercube_v3/timeseries/{cname}_gates.dss"
# th_in = "../../model/schism/azure_dsp_2024_lhc_v3/{cname}/montezuma_{sch_gate}_{cname}.th"

# schism_to_dsm2(cases, th_in, dss_fn_in, dss_fn_out)

# ## DSM2 - RADIAL FROM-TO fix
# cases = range(1,108)
# dss_fn_in = "../../model/dsm2/DSP_DSM2_202412/latinhypercube_v4/timeseries/old_gates/{cname}_gates.dss"
# dss_fn_out = "../../model/dsm2/DSP_DSM2_202412/latinhypercube_v4/timeseries/{cname}_gates.dss"
# dsm2_smscg_from_to(cases, dss_fn_in, dss_fn_out)

## DSM2 - RADIAL FROM-TO fix
cases = range(1,8)
dss_fn_in = "../../model/dsm2/DSP_DSM2_202307/latinhypercube_v3/timeseries/old_gates/{cname}_gates.dss"
dss_fn_out = "../../model/dsm2/DSP_DSM2_202307/latinhypercube_v3/timeseries/{cname}_gates.dss"
dsm2_smscg_from_to(cases, dss_fn_in, dss_fn_out)

