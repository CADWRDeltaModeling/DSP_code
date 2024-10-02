# Taking the SCHISM csv files and converting them to the DSM2-like xlsx files
import pandas as pd
import os
import string
        
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

def schism_to_ann_xlsx(in_fname, mesh):

    with open(in_fname, 'r') as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
    
    # assign env vars to format strings with ----------------------------------------------------
    env_vars = build_dict(inputs.get('env_vars'))
    for env_var in env_vars:
        env_vars[env_var] = string.Formatter().vformat(env_vars[env_var],(),SafeDict((env_vars)))
    csv_fmt = inputs.get('csv_fmt')
    vars_map = inputs.get('vars_map')
    out_dir = env_vars['out_dir']
    if not os.path.exists(out_dir): os.makedirs(out_dir)

    ordered_sheets = ['northern_flow','sjr_flow','exports','dxc_gate_fraction','suisun_gate_fraction','net_delta_cu',
                      'mtz_tidal_nrg','sjr_vernalis_ec','sac_ec','base_ec_output']

    if 'cases' in inputs.keys():
        cases = inputs.get('cases')
        for case in cases:
            case_num = case.get('case_num')
            csv_file = string.Formatter().vformat(csv_fmt,(),
                                                SafeDict(({**env_vars, 
                                                           **locals()})))
            input_df = pd.read_table(csv_file,
                                    sep=',',
                                    index_col=0)
            # create daily averages of each variable
            input_df.index = pd.to_datetime(input_df.index)
            input_df_daily = input_df.resample('D', closed='right').mean()
            input_df_daily = input_df_daily.rename_axis('Time')

            # go through each sheet and create csvs
            for varmap in vars_map:
                sheet_name = varmap['sheet_name']
                sheet_colnames = varmap['sheet_colnames']
                csv_headers = varmap['csv_headers']

                if len(sheet_colnames) == 1:
                    # just one column
                    var_df = input_df_daily.loc[:,csv_headers[0]]
                    var_df.columns = sheet_colnames
                    var_df.to_csv(os.path.join(out_dir, f'{sheet_name}.csv'), float_format="%.2f")
                    
                else:
                    # multiple columns (base_ec_output)
                    var_df = input_df_daily.loc[:,csv_headers[1:]]
                    var_df = var_df.rename(columns=dict(zip(csv_headers[1:], sheet_colnames[1:])))
                    var_df = var_df.rename_axis(sheet_colnames[0])
                    var_df.to_csv(os.path.join(out_dir, f'{sheet_name}.csv'), float_format="%.2f")

            # combine ANN input csv files into xlsx file (TODO: make this unecessary?)
            xlsx_filepath = os.path.join(out_dir, f'{os.path.splitext(os.path.basename(csv_file))[0]}.xlsx')
            with pd.ExcelWriter(xlsx_filepath) as writer:
                for sheet_name in ordered_sheets:
                    df = pd.read_csv(os.path.join(out_dir, f'{sheet_name}.csv'))
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

    else:
        print("Needs cases")

if __name__ == '__main__':

    from schimpy import schism_yaml
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    in_fname = "./input/ann_config_lathypcub_v3_schism.yaml"
    mesh = 'baseline'

    schism_to_ann_xlsx(in_fname, mesh)