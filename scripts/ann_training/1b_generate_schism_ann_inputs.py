# Taking the SCHISM csv files and converting them to the DSM2-like xlsx files
import pandas as pd
import os
import string
import re
import numpy as np

from vtools.functions.unit_conversions import psu_ec_25c


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


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

    with open(in_fname, "r") as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)

    # assign env vars to format strings with ----------------------------------------------------
    env_vars = build_dict(inputs.get("env_vars"))
    for env_var in env_vars:
        env_vars[env_var] = string.Formatter().vformat(
            env_vars[env_var], (), SafeDict((env_vars))
        )
    csv_fmt = inputs.get("csv_fmt")
    in_vars = inputs.get("in_vars")
    out_dir = env_vars["out_dir"]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ordered_sheets = [
        "northern_flow",
        "sjr_flow",
        "exports",
        "dxc_gate_fraction",
        "suisun_gate_fraction",
        "net_delta_cu",
        "mtz_tidal_nrg",
        "sjr_vernalis_ec",
        "sac_ec",
        "base_ec_output",
    ]

    if "cases" in inputs.keys():
        cases = inputs.get("cases")
        for case in cases:
            case_num = case.get("case_num")
            csv_file = string.Formatter().vformat(
                csv_fmt, (), SafeDict(({**env_vars, **locals()}))
            )
            input_df = pd.read_table(csv_file, sep=",", index_col=0)
            # create daily averages of each variable
            input_df.index = pd.to_datetime(input_df.index)
            input_df_daily = input_df.resample("D", closed="right").mean()
            input_df_daily = input_df_daily.rename_axis("Time")

            # go through each sheet and create csvs
            for varmap in in_vars:
                sheet_name = varmap["sheet_name"]
                sheet_colnames = varmap["sheet_colnames"]
                csv_headers = varmap["csv_headers"]

                if len(sheet_colnames) == 1:
                    # just one column
                    var_df = input_df_daily.loc[:, csv_headers[0]]
                    var_df.columns = sheet_colnames
                    var_df.to_csv(
                        os.path.join(out_dir, f"{sheet_name}.csv"), float_format="%.2f"
                    )

                else:
                    # multiple columns (base_ec_output)
                    var_df = input_df_daily.loc[:, csv_headers[1:]]
                    var_df = var_df.rename(
                        columns=dict(zip(csv_headers[1:], sheet_colnames[1:]))
                    )
                    var_df = var_df.rename_axis(sheet_colnames[0])
                    var_df.to_csv(
                        os.path.join(out_dir, f"{sheet_name}.csv"), float_format="%.2f"
                    )

            # combine ANN input csv files into xlsx file (TODO: make this unecessary?)
            xlsx_filepath = os.path.join(
                out_dir, f"{os.path.splitext(os.path.basename(csv_file))[0]}.xlsx"
            )
            with pd.ExcelWriter(xlsx_filepath) as writer:
                for sheet_name in ordered_sheets:
                    df = pd.read_csv(os.path.join(out_dir, f"{sheet_name}.csv"))
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

    else:
        print("Needs cases")


def get_csv_df(csv_file):
    var_df = pd.read_table(csv_file, sep=",", index_col=0)
    # create daily averages of each variable
    var_df.index = pd.to_datetime(var_df.index)
    var_df_daily = var_df.resample("D", closed="right").mean()
    var_df_daily = var_df_daily.rename_axis("Time")

    return var_df_daily


def schism_to_ann_csv(in_fname, mesh, mesh_outname):

    with open(in_fname, "r") as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)

    # assign env vars to format strings with ---------------------------------------------------
    env_vars = build_dict(inputs.get("env_vars"))
    for env_var in env_vars:
        env_vars[env_var] = string.Formatter().vformat(
            env_vars[env_var], (), SafeDict((env_vars))
        )

    in_vars = inputs.get("in_vars")
    comb_vars = inputs.get("comb_vars")
    ann_var_names = [var["ann_colname"] for var in in_vars]
    ann_var_names = ann_var_names + [var["ann_colname"] for var in comb_vars]
    out_ec_locs = build_dict(inputs.get("out_ec_locs"))
    delete_vars = inputs.get("delete_vars")
    out_dir = env_vars["out_dir"]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if "cases" in inputs.keys():
        cases = inputs.get("cases")
        for case_num in cases:
            case_name = case_num  # f'lhc_{case_num}'
            print(case_name)
            case_number = re.findall(r"\d+", case_num)[0]
            # uses this csv just to get the datetime indices for the other variables and to get EC outputs
            csv_indx_fmt = string.Formatter().vformat(
                inputs.get("csv_indx_fmt"), (), SafeDict(({**env_vars, **locals()}))
            )
            input_df_daily = get_csv_df(csv_indx_fmt)
            output_df_daily = pd.DataFrame(
                index=input_df_daily.index,
                columns=["model", "scene", "case"]
                + ann_var_names
                + out_ec_locs["ann_colnames"],
            )  # creates an empty dataframe with spaces for each ANN variable
            output_df_daily["model"] = "SCHISM"
            output_df_daily["scene"] = mesh_outname
            output_df_daily["case"] = case_number

            # go through each sheet and store inputs
            for varmap in in_vars:
                ann_colname = varmap["ann_colname"]
                unit_conv = varmap["unit_conv"]
                csv_header = varmap["csv_header"]
                csv_file = varmap["csv_file"]

                var_df_daily = get_csv_df(
                    string.Formatter().vformat(
                        csv_file, (), SafeDict(({**env_vars, **locals()}))
                    )
                )
                var_df = var_df_daily.loc[:, csv_header]
                if isinstance(unit_conv, list):
                    var_df.loc[:] = np.interp(var_df.values, unit_conv[0], unit_conv[1])
                else:
                    var_df = var_df * unit_conv
                var_df.columns = ann_colname
                output_df_daily[ann_colname] = var_df

            # write combination variables
            for varmap in comb_vars:
                ann_colname = varmap["ann_colname"]
                vars = varmap["vars"]
                mult = varmap["mult"]

                var_df = (output_df_daily[vars].fillna(0) * mult).sum(axis=1)
                var_df.columns = ann_colname
                output_df_daily[ann_colname] = var_df

            # delete the columns that are not needed
            for var in delete_vars:
                if var in output_df_daily.columns:
                    del output_df_daily[var]

            # add EC ouptuts
            for ann_ec, col_ec in zip(
                out_ec_locs["ann_colnames"], out_ec_locs["csv_headers"]
            ):
                output_df_daily[ann_ec] = psu_ec_25c(input_df_daily.loc[:, col_ec])

            # write out csv
            out_fn = os.path.join(out_dir, f"schism_{mesh_outname}_{case_number}.csv")
            # clean up the start and end rows
            output_df_daily.index.name = "datetime"
            output_df_daily.dropna(inplace=True)
            output_df_daily.to_csv(out_fn, float_format="%.2f")
    else:
        print("Needs cases")


if __name__ == "__main__":

    from schimpy import schism_yaml

    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    in_fname = "./input/ann_csv_config_lathypcub_v3_schism.yaml"
    # mesh = "baseline"
    # mesh_outname = "base"
    # schism_to_ann_csv(in_fname, mesh, mesh_outname)
    # mesh = "suisun"
    # mesh_outname = "suisun"
    # schism_to_ann_csv(in_fname, mesh, mesh_outname)
    # mesh = "franks"
    # mesh_outname = "franks"
    # schism_to_ann_csv(in_fname, mesh, mesh_outname)
    mesh = "cache"
    mesh_outname = "cache"
    schism_to_ann_csv(in_fname, mesh, mesh_outname)

    # in_fname = "./input/ann_csv_config_lathypcub_v3_mss_schism.yaml"
    # mesh = "baseline"
    # mesh_outname = "base"
    # schism_to_ann_csv(in_fname, mesh, mesh_outname)

    # in_fname = "./input/ann_csv_config_lathypcub_v3_schism_slr.yaml"
    # mesh = "baseline"
    # mesh_outname = "slr_base"
    # schism_to_ann_csv(in_fname, mesh, mesh_outname)
