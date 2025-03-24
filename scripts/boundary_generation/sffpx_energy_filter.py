# DSM2-output-specific code to generate the SFFPX tidal energy and filter data needed for the ANN

import pandas as pd
import numpy as np
import datetime as dt

from vtools.functions.filter import cosine_lanczos
from vtools.functions.merge import ts_splice
from vtools.functions.unit_conversions import m_to_ft, ft_to_m

from pydelmod.create_ann_inputs import get_dss_data
from sklearn.linear_model import LinearRegression

from schimpy import schism_yaml

import os
import re

import matplotlib.pyplot as plt

import time

os.chdir(os.path.dirname(os.path.abspath(__file__)))


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


# sf_1hr_fn = "./input/9414290_gageheight.txt"
# sf_1hr_raw_ts = pd.read_csv(
#     sf_1hr_fn,
#     index_col=0,
#     parse_dates=[0],
#     usecols=["Date Time", "Water Level"],
#     skipinitialspace=True,
#     dtype={"Water Level": str},
#     sep=",",
#     comment="#",
# )

sf_1hr_fn = r"\\nasbdo\delta_mod\Share\PlanningTide\PlanningTide_2023\Download\sf_stage_noaa_screened.csv"
column_names = ["Date", "Time", "Water Level"]
# Read the file into a DataFrame
sf_1hr_raw_ts = pd.read_csv(
    sf_1hr_fn, delim_whitespace=True, header=None, names=column_names
)
sf_1hr_raw_ts.index = pd.to_datetime(
    sf_1hr_raw_ts[["Date", "Time"]].apply(lambda x: " ".join(x), axis=1)
)
sf_1hr_raw_ts = sf_1hr_raw_ts.drop(columns=["Date", "Time"])
sf_1hr_clean_inst = clean_to_date_range(sf_1hr_raw_ts)
sf_1hr_clean_inst["Water Level"] = pd.to_numeric(
    sf_1hr_clean_inst["Water Level"], errors="coerce"
)
sf_1hr_clean_inst = sf_1hr_clean_inst.ffill()
sf_1hr_clean_inst.index.freq = pd.infer_freq(sf_1hr_clean_inst.index)

sf_6min_fn = "./input/noaa_download/noaa_sffpx_9414290_water_level_1980_2026_ft.csv"
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
sf_harmonic_fn = "./input/noaa_sffpx_9414290_prediction_1990_2025.csv"
sf_harm_raw_ts = pd.read_csv(sf_harmonic_fn, index_col=0, parse_dates=[0])
sf_harm_clean_inst = clean_to_date_range(
    sf_harm_raw_ts
)  # unecessary because it's all predicted anyway
sf_harm_clean_inst = m_to_ft(sf_harm_clean_inst)

# load mrz data
mrz_15min_fn = "./input/dsm2_mrz_stage_ft_15min_clean.csv"
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
if True:
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
if True:
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

# Detrend the Tidal Filter -------------------------------------------------------------
# Pivot Point
pivot_date = pd.to_datetime("1998-5-1")
rate = m_to_ft(1.98e-3)  # m/yr -> ft/yr
print(f"Sea Level Rise trend rate at: {1.98} mm/yr")
sf_filt_df_detrend = sf_filt_df.copy()
sf_filt_df_detrend["Water Level"] = sf_filt_df.apply(
    lambda row: (
        row["Water Level"] + rate * (pivot_date - row.name).days / 365
        if row.name < pivot_date
        else row["Water Level"] - rate * (row.name - pivot_date).days / 365
    ),
    axis=1,
)


# Calculate trend lines
sf_filt_df["Time"] = (sf_filt_df.index - sf_filt_df.index[0]).days
sf_filt_df_detrend["Time"] = (
    sf_filt_df_detrend.index - sf_filt_df_detrend.index[0]
).days

trend_sf_filt = np.polyfit(sf_filt_df["Time"], sf_filt_df["Water Level"], 1)
trend_sf_filt_detrend = np.polyfit(
    sf_filt_df_detrend["Time"], sf_filt_df_detrend["Water Level"], 1
)

sf_filt_trend = sf_filt_df.copy()
sf_detrend_filt_trend = sf_filt_df_detrend.copy()
sf_filt_trend["Trend"] = np.polyval(trend_sf_filt, sf_filt_df["Time"])
sf_detrend_filt_trend["Trend"] = np.polyval(
    trend_sf_filt_detrend, sf_filt_df_detrend["Time"]
)
# Print regression values
print("======Trend Line Regression Values=========")
print(
    f"SF Spliced Trend: Slope = {ft_to_m(trend_sf_filt[0])*10e3:.4f} mm/yr, Intercept = {trend_sf_filt[1]:.4f}"
)
print(
    f"SF Detrended Trend: Slope = {ft_to_m(trend_sf_filt_detrend[0])*10e3:.4f} mm/yr, Intercept = {trend_sf_filt_detrend[1]:.4f}"
)


# Plot detrended filtered data
if True:
    plt.figure(figsize=(12, 6))
    plt.plot(
        sf_filt_df.index,
        sf_filt_df["Water Level"],
        label="SF Spliced",
        linestyle="-",
        alpha=0.7,
    )
    plt.plot(
        sf_filt_df_detrend.index,
        sf_filt_df_detrend["Water Level"],
        label="SF Detrended",
        linestyle="--",
        alpha=0.7,
    )
    plt.plot(
        sf_filt_trend.index,
        sf_filt_trend["Trend"],
        label="Trend SF Spliced",
        linestyle=":",
        alpha=0.7,
    )
    plt.plot(
        sf_detrend_filt_trend.index,
        sf_detrend_filt_trend["Trend"],
        label="Trend SF Detrended",
        linestyle=":",
        alpha=0.7,
    )
    plt.axvline(x=pivot_date, color="r", linestyle="--", label="Pivot Date")
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
if True:
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
if True:
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

if isinstance(reg_out_df.index, pd.DatetimeIndex):
    reg_out_df.index = reg_out_df.index.to_period("D")
reg_out_df.to_csv(
    "./data_out/sffpx_filt_nrg.csv",
    float_format="%.2f",
    index=True,
    index_label="datetime",
)
