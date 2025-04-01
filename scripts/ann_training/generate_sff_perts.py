# DSM2-output-specific code to generate the SFFPX tidal energy and filter data needed for the ANN

import pandas as pd
import datetime as dt

from vtools.functions.filter import cosine_lanczos
from vtools.functions.merge import ts_splice
from vtools.functions.unit_conversions import m_to_ft

from pydelmod.create_ann_inputs import get_dss_data
from sklearn.linear_model import LinearRegression

from schimpy import schism_yaml

import os
import re

import matplotlib.pyplot as plt

import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))


# shifted subtide back to instantaneous data (not used here)
def shift_subtide(df, shift_val, out_col):
    filt_df = cosine_lanczos(df.copy(), cutoff_period="40H", padtype="odd")
    shift_df = filt_df.copy()
    shift_df.index = shift_df.index + pd.Timedelta(days=shift_val)
    # shift_df = shift_df.to_frame()
    shift_df.columns = ["shift_filter"]

    shift_df["inst"] = df
    shift_df["eg_filter"] = filt_df

    shift_df[out_col] = (
        shift_df["inst"] - shift_df["eg_filter"] + shift_df["shift_filter"]
    )

    return shift_df[[out_col]]


# shifted instantaneous (not used here)
def shift_ts(dfs, shift_val, out_col):
    shift_dfs = []
    for df in dfs:
        shift_df = df.copy()
        shift_df.index = shift_df.index + pd.Timedelta(days=shift_val)
        # shift_df = shift_df.to_frame()
        shift_df.columns = [out_col]
        shift_dfs.append(shift_df[[out_col]])

    return shift_dfs


# Function: calculate filter on instantaneous data return single df with filter with daily resample
def calc_filt(df):
    filt_df = cosine_lanczos(df.copy(), cutoff_period="40H", padtype="odd")  # <z>
    filt_df = filt_df.resample("D", closed="right").mean()

    return filt_df


# Function: calculate energy on instantaneous data return single df with energy with daily resample
def calc_nrg(df):
    filt_df = cosine_lanczos(df.copy(), cutoff_period="40H", padtype="odd")  # <z>
    nrg_df = cosine_lanczos(
        (df - filt_df) ** 2, cutoff_period="40H", padtype="odd"
    )  # = < (z- <z>)^2 >
    nrg_df = nrg_df.resample("D", closed="right").mean()

    return nrg_df


# clean up the dataframe so you only get a continuous block of data
def clean_to_date_range(df, start=None, stop=None, max_nan_hrs=3):

    # first clip dataframe
    if start is not None:
        if stop is not None:
            df = df.loc[start:stop]
        else:
            df = df.loc[start:]
    else:
        df = df.loc[:stop]

    # Next fill datetime
    freq = pd.infer_freq(df.index[0:10])
    if not freq[0].isdigit():
        freq = f"1{freq}"
    df = df.resample(freq).asfreq()

    # Drop leading NaNs if they exist
    if df.first_valid_index() is not None:
        df = df.loc[df.first_valid_index() :]

    nan_idx = df.index[df.isna().any(axis=1)]

    first_nan_idx = df.index[-1]
    if len(nan_idx) > 0:
        for n, ni in enumerate(nan_idx):
            start_nan = ni
            int_nan = ni
            end_nan = None
            for next_ni in nan_idx[n + 1 :]:
                if (next_ni - int_nan) == pd.Timedelta(freq):
                    # it's a streak, see how long
                    end_nan = next_ni
                    int_nan = end_nan
                    # check if over limit
                    nan_len = end_nan - start_nan
                    if nan_len >= pd.Timedelta(hours=max_nan_hrs):
                        break
                elif end_nan is not None:
                    nan_len = end_nan - start_nan
                    break
                else:
                    nan_len = pd.Timedelta(hours=0)
                    break

            if nan_len >= pd.Timedelta(hours=max_nan_hrs):
                first_nan_idx = start_nan
                break
            else:
                next

    # Keep only the part before the first NaN
    df = df.loc[: first_nan_idx - pd.Timedelta("1s")]

    return df


sf_1hr_fn = "../boundary_generation/input/9414290_gageheight.txt"
sf_1hr_raw_ts = pd.read_csv(
    sf_1hr_fn,
    index_col=0,
    parse_dates=[0],
    usecols=["Date Time", "Water Level"],
    skipinitialspace=True,
    dtype={"Water Level": str},
    sep=",",
    comment="#",
)
sf_1hr_clean_inst = clean_to_date_range(sf_1hr_raw_ts)
sf_1hr_clean_inst["Water Level"] = pd.to_numeric(
    sf_1hr_clean_inst["Water Level"], errors="coerce"
)
sf_1hr_clean_inst = sf_1hr_clean_inst.ffill()
sf_1hr_clean_inst.index.freq = pd.infer_freq(sf_1hr_clean_inst.index)

sf_6min_fn = "../boundary_generation/input/noaa_download/noaa_sffpx_9414290_water_level_1980_2026_ft.csv"
sf_6min_raw_ts = pd.read_csv(
    sf_6min_fn,
    index_col=0,
    parse_dates=[0],
    usecols=["Date Time", "Water Level"],
    skipinitialspace=True,
    dtype={"Water Level": str},
    sep=",",
    comment="#",
)
sf_6min_clean_inst = clean_to_date_range(
    sf_6min_raw_ts, start=pd.to_datetime("2012-9-1")
)
sf_6min_clean_inst["Water Level"] = pd.to_numeric(
    sf_6min_clean_inst["Water Level"], errors="coerce"
)
sf_6min_clean_inst = sf_6min_clean_inst.ffill()
sf_6min_clean_inst.index.freq = pd.infer_freq(sf_6min_clean_inst.index)

# load harmonic data for tidal energy calc
sf_harmonic_fn = (
    "../boundary_generation/input/noaa_sffpx_9414290_prediction_1990_2025.csv"
)
sf_harm_raw_ts = pd.read_csv(sf_harmonic_fn, index_col=0, parse_dates=[0])
sf_harm_clean_inst = clean_to_date_range(
    sf_harm_raw_ts
)  # unecessary because it's all predicted anyway
sf_harm_clean_inst = m_to_ft(sf_harm_clean_inst)

# load mrz data
mrz_15min_fn = "../boundary_generation/input/dsm2_mrz_stage_ft_15min_clean.csv"
mrz_15min_raw_ts = pd.read_csv(mrz_15min_fn, index_col=0, parse_dates=[0])
mrz_15min_clean_inst = clean_to_date_range(mrz_15min_raw_ts)

# Regress Filtered data  -------------------------------------------------------------
# filter the data
sf_1hr_filt = calc_filt(sf_1hr_clean_inst)
sf_6min_filt = calc_filt(sf_6min_clean_inst)
mrz_15min_filt = calc_filt(mrz_15min_clean_inst)
mrz_15min_filt.columns = ["mrz"]

# create regression for tidal filter of <sf> ~ <mrz>
print(f"..creating <sf> ~ <mrz> regression")
reg_filt_df = pd.merge(
    mrz_15min_filt,
    pd.concat([sf_1hr_filt, sf_6min_filt]),
    how="left",
    left_index=True,
    right_index=True,
)
reg_filt_df.columns = ["mrz", "sf"]
reg_filt_df = reg_filt_df.dropna()
model = LinearRegression()  # create linreg model
model.fit(reg_filt_df[["mrz"]], reg_filt_df[["sf"]])  # Fit Linear Regression
print("======Tidal Filter Linear Regression Statistics=========")
r_squared = model.score(reg_filt_df[["mrz"]], reg_filt_df[["sf"]])
print(f"R² value: {r_squared:.4f}")
print(f"intercept value: {model.intercept_[0]:.4f}")
print(f"Coefficient/Scaling value: {model.coef_[0][0]:.4f}")

pred_filt_df = pd.DataFrame(index=mrz_15min_filt.index)
pred_filt_df["sf_pred"] = model.predict(mrz_15min_filt[["mrz"]])

# Plot regression of filtered data
if False:
    plt.figure(figsize=(12, 6))
    plt.plot(
        mrz_15min_filt.index,
        mrz_15min_filt["mrz"],
        label="MRZ",
        linestyle="-",
        alpha=0.7,
    )
    plt.plot(
        reg_filt_df.index, reg_filt_df["sf"], label="SF", linestyle="--", alpha=0.7
    )
    plt.plot(
        pred_filt_df.index,
        pred_filt_df["sf_pred"],
        label="SF Predicted",
        linestyle=":",
        alpha=0.7,
    )
    plt.xlabel("")
    plt.ylabel("Filtered Tide <ft>")
    plt.legend()
    plt.grid(True)
    plt.show()

# splice the data together
pred_filt_df.columns = sf_1hr_filt.columns
pred_filt_df = pred_filt_df.loc[
    (sf_1hr_filt.index[-1] - pd.Timedelta(weeks=2)) : (
        sf_6min_filt.index[0] + pd.Timedelta(weeks=2)
    )
]
sf_filt_df = ts_splice((sf_1hr_filt, pred_filt_df), transition="prefer_first")
sf_filt_df = ts_splice((sf_filt_df, sf_6min_filt), transition="prefer_last")

# Plot splice of filtered data
if False:
    plt.figure(figsize=(12, 6))
    plt.plot(
        sf_1hr_filt.index,
        sf_1hr_filt["Water Level"],
        label="1 Hour Data",
        linestyle="-",
        alpha=0.7,
    )
    plt.plot(
        sf_6min_filt.index,
        sf_6min_filt["Water Level"],
        label="6 Minute Data",
        linestyle="-.",
        alpha=0.7,
    )
    plt.plot(
        pred_filt_df.index,
        pred_filt_df["Water Level"],
        label="SF Predicted",
        linestyle=":",
        alpha=0.7,
    )
    plt.plot(
        sf_filt_df.index,
        sf_filt_df["Water Level"],
        label="SF Spliced",
        linestyle="--",
        alpha=0.7,
    )
    plt.xlabel("")
    plt.ylabel("Filtered SF Tide <ft>")
    plt.legend()
    plt.grid(True)
    plt.show()

# Regress the Tidal Energy -------------------------------------------------------------
# filter data
sf_1hr_nrg = calc_nrg(sf_1hr_clean_inst)
sf_6min_nrg = calc_nrg(sf_6min_clean_inst)
sf_harm_nrg = calc_nrg(sf_harm_clean_inst)
sf_harm_nrg.columns = ["sf_harm"]

# create regression for tidal energy of  <(sf_obs- <sf_obs>)^2> ~ <(sf_harm- <sf_harm>)^2>
print(f"..creating <(sf_obs- <sf_obs>)^2> ~ <(sf_harm- <sf_harm>)^2> regression")
reg_nrg_df = pd.merge(
    sf_harm_nrg,
    pd.concat([sf_1hr_nrg, sf_6min_nrg]),
    how="left",
    left_index=True,
    right_index=True,
)
reg_nrg_df.columns = ["sf_harm", "sf_obs"]
reg_nrg_df = reg_nrg_df.dropna()
model = LinearRegression()  # create linreg model
model.fit(reg_nrg_df[["sf_harm"]], reg_nrg_df[["sf_obs"]])  # Fit Linear Regression
print("======Tidal Energy Linear Regression Statistics=========")
r_squared = model.score(reg_nrg_df[["sf_harm"]], reg_nrg_df[["sf_obs"]])
print(f"R² value: {r_squared:.4f}")
print(f"intercept value: {model.intercept_[0]:.4f}")
print(f"Coefficient/Scaling value: {model.coef_[0][0]:.4f}")

pred_nrg_df = pd.DataFrame(index=sf_harm_nrg.index)
pred_nrg_df["sf_pred"] = model.predict(sf_harm_nrg[["sf_harm"]])

# Plot regression of tidal energy
if False:
    plt.figure(figsize=(12, 6))
    plt.plot(
        sf_harm_nrg.index,
        sf_harm_nrg["sf_harm"],
        label="SF Harmonic",
        linestyle=":",
        alpha=0.7,
    )
    plt.plot(
        reg_nrg_df.index,
        reg_nrg_df["sf_obs"],
        label="SF Observed",
        linestyle=":",
        alpha=0.7,
    )
    plt.plot(
        pred_nrg_df.index,
        pred_nrg_df["sf_pred"],
        label="SF Predicted",
        linestyle="--",
        alpha=0.7,
    )
    plt.xlabel("")
    plt.ylabel("Tidal Energy  <(ft- <ft>)^2>")
    plt.legend()
    plt.grid(True)
    plt.show()

# splice the data together
pred_nrg_df.columns = sf_6min_nrg.columns
pred_nrg_df = pred_nrg_df.loc[
    (sf_1hr_nrg.index[-1] - pd.Timedelta(weeks=2)) : (
        sf_6min_nrg.index[0] + pd.Timedelta(weeks=2)
    )
]
sf_nrg_df = ts_splice((sf_1hr_nrg.iloc[:-2], pred_nrg_df), transition="prefer_first")
sf_nrg_df = ts_splice((sf_nrg_df, sf_6min_nrg.iloc[1:]), transition="prefer_last")

# Plot splice of tidal energy
if False:
    plt.figure(figsize=(12, 6))
    plt.plot(
        sf_1hr_nrg.index,
        sf_1hr_nrg["Water Level"],
        label="1 Hour Data",
        linestyle="-",
        alpha=0.7,
    )
    plt.plot(
        sf_6min_nrg.index,
        sf_6min_nrg["Water Level"],
        label="6 Minute Data",
        linestyle="-.",
        alpha=0.7,
    )
    plt.plot(
        pred_nrg_df.index,
        pred_nrg_df["Water Level"],
        label="SF Predicted",
        linestyle=":",
        alpha=0.7,
    )
    plt.plot(
        sf_nrg_df.index,
        sf_nrg_df["Water Level"],
        label="SF Spliced",
        linestyle="--",
        alpha=0.7,
    )
    plt.xlabel("")
    plt.ylabel("Tidal Energy  <(ft- <ft>)^2>")
    plt.legend()
    plt.grid(True)
    plt.show()

# Generate various tidal perturbations ======================================================================================
# Regular
reg_out_df = pd.merge(
    sf_filt_df, sf_nrg_df, left_index=True, right_index=True, how="outer"
)
reg_out_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Shifted - 50d
filt, nrg = shift_ts([sf_filt_df, sf_nrg_df], -50, "shift-50")
shift_minus50_df = pd.merge(filt, nrg, left_index=True, right_index=True, how="outer")
shift_minus50_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Shifted + 100d
filt, nrg = shift_ts([sf_filt_df, sf_nrg_df], 100, "shift+100")
shift_plus100_df = pd.merge(filt, nrg, left_index=True, right_index=True, how="outer")
shift_plus100_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Subtidal Pert
pert_fn = "../boundary_generation/data_out/perturb_historical_subtide_v1.csv"
mrz_filt_fn = (
    "../boundary_generation/data_out/dsm2_mrz_stage_ft_15min_clean_filtered.csv"
)
mrz_pert_df = (
    pd.read_csv(pert_fn, index_col=[0], parse_dates=[0], header=None)
    .resample("D", closed="right")
    .mean()
)
mrz_filt_df = (
    pd.read_csv(mrz_filt_fn, index_col=[0], parse_dates=[0])
    .resample("D", closed="right")
    .mean()
)
mrz_pert_df.columns = mrz_filt_df.columns
mrz_filt_diff_df = mrz_pert_df - mrz_filt_df

sf_filt_df.columns = mrz_filt_diff_df.columns
sf_subtide_pert = sf_filt_df + mrz_filt_diff_df
sf_subtide_pert = sf_subtide_pert.dropna()

subtide_pert_df = pd.merge(
    sf_subtide_pert, sf_nrg_df, left_index=True, right_index=True, how="outer"
)
subtide_pert_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Shifted Subtide + 7 days
[filt] = shift_ts([sf_filt_df], 7, "shift+7")
shift_subtide_plus7_df = pd.merge(
    filt, sf_nrg_df, left_index=True, right_index=True, how="outer"
)
shift_subtide_plus7_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Shifted Subtide - 7 days
[filt] = shift_ts([sf_filt_df], -7, "shift-7")
shift_subtide_minus7_df = pd.merge(
    filt, sf_nrg_df, left_index=True, right_index=True, how="outer"
)
shift_subtide_minus7_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Subtidal Pert + Shifted Subtide + 7 days
[filt] = shift_ts([sf_subtide_pert], 7, "shift+7")
subtide_pert_shift_plus7_df = pd.merge(
    filt, sf_nrg_df, left_index=True, right_index=True, how="outer"
)
subtide_pert_shift_plus7_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Subtidal Pert + Shifted Subtide - 7 days
[filt] = shift_ts([sf_subtide_pert], -7, "shift-7")
subtide_pert_shift_minus7_df = pd.merge(
    filt, sf_nrg_df, left_index=True, right_index=True, how="outer"
)
subtide_pert_shift_minus7_df.columns = ["sf_tidal_filter", "sf_tidal_energy"]

# Run on cases =================================================================================================

col_order = [
    "model",
    "scene",
    "case",
    "northern_flow",
    "sac_flow",
    "sjr_flow",
    "exports",
    "cu_flow",
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

# load lhc_v3
lhc_fn = "../boundary_generation/data_out/lhc_v3.csv"
case_nums = range(1001, 1008)

# # load lhc_v4
# lhc_fn = "../boundary_generation/data_out/lhc_v4.csv"
# case_nums = range(106, 108)

run_wait = False

# casanntra dir
cas_dir = "../../../casanntra/data"
csv_fmt = "dsm2_base_{case_num}.csv"

# write data into dsm2 dataframes
lhc_df = pd.read_csv(lhc_fn)
for index, row in lhc_df.iterrows():
    case_num = re.search(r"(\d+)$", row["case"]).group(1)
    if int(case_num) in case_nums:
        print(row["case"])
        casanntra_casefile = os.path.join(cas_dir, csv_fmt.format(case_num=case_num))

        # to run this in parallel with ongoing/overnight check if the model is finished running
        if run_wait:
            # print(f"index: {index}, case_nums[-1]: {case_nums[-2]}")
            if index < case_nums[-2]:
                next_ann_fn = casanntra_casefile.replace(
                    f"_{case_num}", f"_{int(case_num)+1}"
                )
                while not os.path.exists(next_ann_fn):
                    print(
                        f"Waiting for file {next_ann_fn} to appear so that case {case_num} is done..."
                    )
                    time.sleep(120)  # Wait for 30 seconds before checking again

                print(
                    f"File {next_ann_fn} is now available! Post-processing model results for case {case_num}"
                )

            else:
                print(f"Last case, waiting 5 minutes to finish")
                time.sleep(60 * 5)

        cdf = pd.read_csv(casanntra_casefile, parse_dates=[0], index_col=0)

        if "sf_tidal_filter" in cdf.columns:
            cdf.drop(columns=["sf_tidal_filter", "sf_tidal_energy"], inplace=True)

        if row["tide"] == "Shifted + 100d":
            merge_df = pd.merge(
                cdf,
                shift_plus100_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Subtidal Pert":
            merge_df = pd.merge(
                cdf,
                subtide_pert_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Regular":
            merge_df = pd.merge(
                cdf,
                reg_out_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Shifted - 50d":
            merge_df = pd.merge(
                cdf,
                shift_minus50_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Shifted Subtide + 7 days":
            merge_df = pd.merge(
                cdf,
                shift_subtide_minus7_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Shifted Subtide - 7 days":
            merge_df = pd.merge(
                cdf,
                shift_subtide_plus7_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Subtidal Pert + Shifted Subtide + 7 days":
            merge_df = pd.merge(
                cdf,
                subtide_pert_shift_plus7_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        elif row["tide"] == "Subtidal Pert + Shifted Subtide - 7 days":
            merge_df = pd.merge(
                cdf,
                subtide_pert_shift_minus7_df.loc[cdf.index[0] : cdf.index[-1]],
                left_index=True,
                right_index=True,
                how="outer",
            )

        else:
            print("NO METHOD")

        merge_df = merge_df[[col for col in col_order if col in merge_df.columns]]
        if isinstance(merge_df.index, pd.DatetimeIndex):
            merge_df.index = merge_df.index.to_period("D")
        merge_df.to_csv(
            casanntra_casefile, float_format="%.2f", index=True, index_label="datetime"
        )
