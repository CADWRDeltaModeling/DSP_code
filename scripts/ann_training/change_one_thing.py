import pandas as pd
import datetime as dt

from vtools.functions.merge import ts_splice
from vtools.functions.unit_conversions import m_to_ft

import pyhecdss
from pydelmod.create_ann_inputs import get_dss_data

from schimpy import schism_yaml

import os
import re
import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))

col_order = [
    "model",
    "scene",
    "case",
    "northern_flow",
    "sac_flow",
    "sjr_flow",
    "exports",
    "cu_total",
    "cu_delta",
    "cu_suisun",
    "ndo",
    "dcc",
    "smscg",
    "vern_ec",
    "mrz_tidal_energy",
    "mrz_tidal_filter",
    "sf_tidal_energy",
    "sf_tidal_filter",
    "anc",
    "anh",
    "bac",
    "bdl",
    "bdt",
    "bet",
    "cll",
    "cse",
    "dsj",
    "emm2",
    "frk",
    "god",
    "gys",
    "gzl",
    "hll",
    "hol2",
    "ibs",
    "jer",
    "mal",
    "mrz",
    "mtz",
    "nsl2",
    "obi",
    "oh4",
    "old",
    "pct",
    "ppt",
    "rri2",
    "rsl",
    "sal",
    "snc",
    "srv",
    "sss",
    "tms",
    "trp",
    "tss",
    "uni",
    "vcu",
    "vol",
    "wci",
    "x2",
]


def process_gate_data(
    dss_filename,
    b_part,
    c_part,
    map_zero_one=["df==0", "df==1"],
    startDateStr=None,
    endDateStr=None,
):
    """
    Read delta cross-channel gate operation data
    Create daily time series indicating fraction of maximum gate opening (100% means both gates open all day).
    """

    with pyhecdss.DSSFile(dss_filename) as d:
        fdname, generated = d._check_condensed_catalog_file_and_recatalog(
            condensed=True
        )
        catdf = pyhecdss.DSSFile._read_catalog_dsd(fdname)

        filtered_df = catdf[(catdf.B == b_part) & (catdf.C == c_part)]
        p = d.get_pathnames(filtered_df)[0]

        if d.parse_pathname_epart(p).startswith("IR-"):
            df, units, ptype = d.read_its(
                p, startDateStr=startDateStr, endDateStr=endDateStr
            )
        else:
            df, units, ptype = d.read_rts(p)

    df[eval(map_zero_one[0])] = 0
    df[eval(map_zero_one[1])] = 1

    # resample to 1 minute, then fill forward (with last value)
    df_1min = df.resample("T", closed="right").ffill()
    # now find daily averages of one minute data
    df_daily_avg = df_1min.resample("D", closed="right").mean()
    df_daily_avg = df_daily_avg.rename(columns={df_daily_avg.columns[0]: "gate_pos"})

    return df_daily_avg


def flip_exports(lhc_fn, case_nums):
    # rewrite data into dsm2 dataframes
    lhc_df = pd.read_csv(lhc_fn)
    for index, row in lhc_df.iterrows():
        case_num = re.search(r"(\d+)$", row["case"]).group(1)
        if int(case_num) in case_nums:
            print(row["case"])
            casanntra_casefile = os.path.join(
                cas_dir, csv_fmt.format(case_num=case_num)
            )
            cdf = pd.read_csv(casanntra_casefile, parse_dates=[0], index_col=0)

            cdf["exports"] = -cdf["exports"]

            if isinstance(cdf.index, pd.DatetimeIndex):
                cdf.index = cdf.index.to_period("D")
            cdf.to_csv(
                casanntra_casefile,
                float_format="%.2f",
                index=True,
                index_label="datetime",
            )


def add_mrz(in_fname, lhc_fn, case_nums, col_order, ec_locs=["mrz"]):

    with open(in_fname, "r") as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
    lhc_df = pd.read_csv(lhc_fn)
    for index, row in lhc_df.iterrows():
        case_num = re.search(r"(\d+)$", row["case"]).group(1)
        if int(case_num) in case_nums:
            casanntra_casefile = os.path.join(
                cas_dir, csv_fmt.format(case_num=case_num)
            )
            if int(case_num) > 1000:
                mod_case_num = True
                case_num = int(case_num) - 1000
            print(row["case"])
            cdf = pd.read_csv(casanntra_casefile, parse_dates=[0], index_col=0)
            experiment = inputs.get("experiment")

            base_study_folder = inputs.get("base_study_folder").format(**locals())
            model_input_folder = inputs.get("model_input_folder").format(**locals())
            model_folder = inputs.get("model_folder").format(**locals())
            model_output_folder = inputs.get("model_output_folder").format(**locals())
            model_ec_file = inputs.get("model_ec_file").format(**locals())

            b_part_dss_filename_dict = {name.upper(): model_ec_file for name in ec_locs}
            b_part_c_part_dict = {name.upper(): "EC" for name in ec_locs}
            b_part_e_part_dict = {name.upper(): "1HOUR" for name in ec_locs}

            df_model_ec = get_dss_data(
                b_part_dss_filename_dict,
                "b_part",
                primary_part_c_part_dict=b_part_c_part_dict,
                primary_part_e_part_dict=b_part_e_part_dict,
            )
            df_model_ec.columns = df_model_ec.columns.str.lower()
            df_model_ec.index = df_model_ec.index.to_timestamp()

            merge_df = pd.merge(
                cdf, df_model_ec, how="left", left_index=True, right_index=True
            )
            merge_df = merge_df[[col for col in col_order if col in merge_df.columns]]

            if isinstance(merge_df.index, pd.DatetimeIndex):
                merge_df.index = merge_df.index.to_period("D")
            merge_df.to_csv(
                casanntra_casefile,
                float_format="%.2f",
                index=True,
                index_label="datetime",
            )


def fix_smscg(in_fname, lhc_fn, case_nums, col_order, ec_locs=["mrz"]):

    with open(in_fname, "r") as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
    lhc_df = pd.read_csv(lhc_fn)
    for index, row in lhc_df.iterrows():
        case_num = re.search(r"(\d+)$", row["case"]).group(1)
        if int(case_num) in case_nums:
            casanntra_casefile = os.path.join(
                cas_dir, csv_fmt.format(case_num=case_num)
            )
            if int(case_num) > 1000:
                mod_case_num = True
                case_num = int(case_num) - 1000
            print(row["case"])
            cdf = pd.read_csv(casanntra_casefile, parse_dates=[0], index_col=0)
            experiment = inputs.get("experiment")

            base_study_folder = inputs.get("base_study_folder").format(**locals())
            model_input_folder = inputs.get("model_input_folder").format(**locals())
            model_folder = inputs.get("model_folder").format(**locals())
            model_output_folder = inputs.get("model_output_folder").format(**locals())
            gate_dss_file = inputs.get("gate_dss_file").format(**locals())

            b_part = "MTZSL"
            c_part = "RADIAL_OP"
            suisun_op_df = process_gate_data(
                gate_dss_file,
                b_part,
                c_part,
                map_zero_one=["df!=-10", "df==-10"],
                startDateStr="01JAN1953",
                endDateStr="01JAN2030",
            )

            cdf["smscg"] = suisun_op_df["gate_pos"]
            if isinstance(cdf.index, pd.DatetimeIndex):
                cdf.index = cdf.index.to_period("D")
            cdf.to_csv(
                casanntra_casefile,
                float_format="%.2f",
                index=True,
                index_label="datetime",
            )


def split_dsm2_cu(inputs, cdf, case_num=None):
    experiment = inputs.get("experiment")

    base_study_folder = inputs.get("base_study_folder").format(**locals())
    model_input_folder = inputs.get("model_input_folder").format(**locals())
    model_folder = inputs.get("model_folder").format(**locals())
    dcd_dss_file = inputs.get("dcd_dss_file").format(**locals())
    smcd_dss_file = inputs.get("smcd_dss_file").format(**locals())

    div_seep_dcd_c_part_dss_filename_dict = {
        "DIV-FLOW": dcd_dss_file,
        "SEEP-FLOW": dcd_dss_file,
    }
    div_seep_smcd_c_part_dss_filename_dict = {
        "DIV-FLOW": smcd_dss_file,
        "SEEP-FLOW": smcd_dss_file,
    }
    drain_dcd_c_part_dss_filename_dict = {"DRAIN-FLOW": dcd_dss_file}
    drain_smcd_c_part_dss_filename_dict = {"DRAIN-FLOW": smcd_dss_file}

    df_div_seep_dcd = get_dss_data(
        div_seep_dcd_c_part_dss_filename_dict, "c_part", filter_b_part_numeric=True
    )
    df_div_seep_smcd = get_dss_data(
        div_seep_smcd_c_part_dss_filename_dict, "c_part", filter_b_part_numeric=True
    )
    df_drain_dcd = get_dss_data(
        drain_dcd_c_part_dss_filename_dict, "c_part", filter_b_part_numeric=True
    )
    df_drain_smcd = get_dss_data(
        drain_smcd_c_part_dss_filename_dict, "c_part", filter_b_part_numeric=True
    )

    df_div_seep_dcd["dcd_divseep_total"] = df_div_seep_dcd[df_div_seep_dcd.columns].sum(
        axis=1
    )
    df_div_seep_smcd["smcd_divseep_total"] = df_div_seep_smcd[
        df_div_seep_smcd.columns
    ].sum(axis=1)

    df_drain_dcd["dcd_drain_total"] = df_drain_dcd[df_drain_dcd.columns].sum(axis=1)
    df_drain_smcd["smcd_drain_total"] = df_drain_smcd[df_drain_smcd.columns].sum(axis=1)

    cu_total_dcd = pd.merge(
        df_div_seep_dcd, df_drain_dcd, how="left", left_index=True, right_index=True
    )
    cu_total_smcd = pd.merge(
        df_div_seep_smcd, df_drain_smcd, how="left", left_index=True, right_index=True
    )
    cu_total = pd.merge(
        cu_total_dcd, cu_total_smcd, how="left", left_index=True, right_index=True
    )

    cu_total["cu_total"] = (
        cu_total["dcd_divseep_total"]
        + cu_total["smcd_divseep_total"]
        - cu_total["dcd_drain_total"]
        - cu_total["smcd_drain_total"]
    )
    cu_total["cu_delta"] = cu_total["dcd_divseep_total"] - cu_total["dcd_drain_total"]
    cu_total["cu_suisun"] = (
        +cu_total["smcd_divseep_total"] - cu_total["smcd_drain_total"]
    )
    out_df = cu_total[["cu_total", "cu_delta", "cu_suisun"]]
    out_df.index = out_df.index.to_timestamp()

    merge_df = pd.merge(cdf, out_df, how="left", left_index=True, right_index=True)
    merge_df = merge_df[[col for col in col_order if col in merge_df.columns]]

    if isinstance(merge_df.index, pd.DatetimeIndex):
        merge_df.index = merge_df.index.to_period("D")

    return merge_df


def split_dsm2_cu_cases(in_fname, lhc_fn, case_nums, col_order, pseudo_case=None):
    with open(in_fname, "r") as f:
        # loader = RawLoader(stream)
        inputs = schism_yaml.load(f)
    if pseudo_case is not None:
        casanntra_casefile = os.path.join(cas_dir, csv_fmt.format(case_num=pseudo_case))
        cdf = pd.read_csv(casanntra_casefile, parse_dates=[0], index_col=0)
        merge_df = split_dsm2_cu(inputs, cdf)
        merge_df.to_csv(
            casanntra_casefile,
            float_format="%.2f",
            index=True,
            index_label="datetime",
        )
    else:
        lhc_df = pd.read_csv(lhc_fn)
        for index, row in lhc_df.iterrows():
            case_num = re.search(r"(\d+)$", row["case"]).group(1)
            if int(case_num) in case_nums:
                if int(case_num) > 1000:
                    mod_case_num = True
                    case_num = int(case_num) - 1000
                casanntra_casefile = os.path.join(
                    cas_dir, csv_fmt.format(case_num=case_num)
                )
                print(row["case"])
                cdf = pd.read_csv(casanntra_casefile, parse_dates=[0], index_col=0)

                merge_df = split_dsm2_cu(inputs, cdf, case_num=case_num)
                merge_df.to_csv(
                    casanntra_casefile,
                    float_format="%.2f",
                    index=True,
                    index_label="datetime",
                )


# casanntra dir
cas_dir = "../../../scripts/casanntra/data"
csv_fmt = "dsm2_base_{case_num}.csv"

### FLIP EXPORTS -----------------------------------------------------------------------------------------------
# # load lhc_v3
# lhc_fn = "../boundary_generation/data_out/lhc_v3.csv"
# case_nums = range(1001, 1008)

# flip_exports(lhc_fn, case_nums)

# # load lhc_v4
# lhc_fn = "../boundary_generation/data_out/lhc_v4.csv"
# case_nums = range(1, 108)

# flip_exports(lhc_fn, case_nums)


### Add MRZ ----------------------------------------------------------------------------------------------------
# in_fname = "./input/ann_config_lathypcub_v3_dsm2.yaml"
# lhc_fn = "../boundary_generation/data_out/lhc_v3.csv"
# case_nums = range(1001, 1008)
# add_mrz(in_fname, lhc_fn, case_nums, col_order)

# in_fname = "./input/ann_config_lathypcub_v4_dsm2.yaml"
# lhc_fn = "../boundary_generation/data_out/lhc_v4.csv"
# case_nums = range(1, 108)
# add_mrz(in_fname, lhc_fn, case_nums, col_order)


# ### Fix SMSCG gates --------------------------------------------------------------------------------------------
# in_fname = "./input/ann_config_lathypcub_v3_dsm2.yaml"
# lhc_fn = "../boundary_generation/data_out/lhc_v3.csv"
# case_nums = range(1001, 1008)
# fix_smscg(in_fname, lhc_fn, case_nums, col_order)

### Split CU ----------------------------------------------------------------------------------------------------
in_fname = "./input/ann_config_lathypcub_v3_dsm2.yaml"
lhc_fn = "../boundary_generation/data_out/lhc_v3.csv"
case_nums = range(1001, 1008)
split_dsm2_cu_cases(in_fname, lhc_fn, case_nums, col_order)

# in_fname = "./input/ann_config_lathypcub_v4_dsm2.yaml"
# lhc_fn = "../boundary_generation/data_out/lhc_v4.csv"
# case_nums = range(1, 108)
# split_dsm2_cu_cases(in_fname, lhc_fn, case_nums, col_order)

# in_fname = "./input/ann_config_dsm2_historical_2000-2020.yml"
# split_dsm2_cu_cases(in_fname, None, None, col_order, pseudo_case="historical")
