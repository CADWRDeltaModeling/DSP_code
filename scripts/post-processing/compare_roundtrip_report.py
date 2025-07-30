import os
import pandas as pd
import numpy as np
import datetime as dt

from vtools.functions.interpolate import rhistinterp
from vtools.functions.read_dss import read_dss
from vtools.functions.unit_conversions import ec_psu_25c
from vtools.data.vtime import months

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import pickle

plt.rcParams["font.family"] = "serif"  # Example: set a generic serif font
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams.update({"font.size": 11})


def get_schism_results(schism_fn, col_names, dt="ME"):
    df = pd.read_csv(schism_fn, index_col=0, parse_dates=[0])
    df = df.loc[:, col_names]

    df = df.resample(dt).mean()

    return df


def get_calsim_results(
    calsim_dss, pathname, col_names, start_date=None, end_date=None, dt=months(1)
):
    data = read_dss(
        calsim_dss,
        pathname,
        start_date=start_date,
        end_date=end_date,
        dt=dt,
        p=999e99,
    )
    data.columns = [col_names]

    # CalSim's EC and X2 results are off by 1 month
    data.index = data.index - pd.DateOffset(months=1)

    return data


def add_sch_cal_matplotlib(
    ax,
    calsim_df,
    schism_df,
    color_map,
    ylabel,
    ann_suffix,
    altname,
    std_df=None,
    ymax=-9999,
    ymin=9999,
):
    # Handle both DataFrame and Series input
    calsim_vals = (
        calsim_df.iloc[:, 0] if isinstance(calsim_df, pd.DataFrame) else calsim_df
    )
    schism_vals = (
        schism_df.iloc[:, 0] if isinstance(schism_df, pd.DataFrame) else schism_df
    )

    # Set ylims
    ymin = min(ymin, min(calsim_vals.min(), schism_vals.min()) * 0.8)
    ymax = max(ymax, max(calsim_vals.max(), schism_vals.max()) * 1.03)

    # Plot CalSim Standard if provided
    if std_df is not None:
        ax.fill_between(
            std_df.index,
            std_df.iloc[:, 0],
            y2=ymax,
            color="gray",
            alpha=0.3,
            step="post",  # Make fill stepped
        )
        ax.plot(
            std_df.index,
            std_df.iloc[:, 0],
            color="gray",
            label="Compliance Limit",
            linestyle="--",
            linewidth=0.8,  # Slightly smaller than default
            drawstyle="steps-post",  # Stepped line
        )

    # Plot CalSim
    color = color_map.get(f"CalSim-{ann_suffix}", "gray")
    ax.plot(
        calsim_df.index,
        calsim_df.iloc[:, 0],
        label=f"CalSim {altname}",
        color=color,
        linewidth=0.8,  # Slightly smaller than default
        drawstyle="steps-post",  # Stepped line
    )

    # Plot SCHISM
    color = color_map.get(f"SCHISM-{ann_suffix}", "gray")
    ax.plot(
        schism_df.index,
        schism_df.iloc[:],
        label=f"SCHISM {altname}",
        color=color,
        linestyle="dashdot",
        linewidth=0.8,  # Slightly smaller than default
        drawstyle="steps-post",  # Stepped line
    )

    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", linestyle=":", linewidth=0.5)
    ax.set_ylim([ymin, ymax])

    return ymin, ymax


def get_legend_handles(ax):
    handles, labels = ax.get_legend_handles_labels()
    calsim_handles = [h for h, l in zip(handles, labels) if "CalSim" in l]
    calsim_labels = [l for l in labels if "CalSim" in l]
    schism_handles = [h for h, l in zip(handles, labels) if "SCHISM" in l]
    schism_labels = [l for l in labels if "SCHISM" in l]
    compliance_handles = [h for h, l in zip(handles, labels) if "Compliance" in l]
    compliance_labels = [l for l in labels if "Compliance" in l]

    final_handles = calsim_handles + schism_handles + compliance_handles
    final_labels = calsim_labels + schism_labels + compliance_labels

    return final_handles, final_labels


def export_legend(legend, filename="plots/legend.png", expand=[-5, -5, 5, 5]):
    fig = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array(expand)))
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))

    start_date = "2015-02-18"
    end_date = "2016-05-01"

    dsp_repo = "D:/projects/delta_salinity/DSP_code"  # Local: D:/projects/delta_salinity/DSP_code, Remote: /scratch/tomkovic/DSP_code
    casanntra_repo = "D:/projects/delta_salinity/scripts/casanntra"  # Local: D:/projects/delta_salinity/scripts/casanntra, Remote: /scratch/tomkovic/casanntra

    figsize = [6, 4]
    plt_fn = "./plots/roundtrip_{var}_{ann_base}.png"

    alternative_dict = {
        "Suisun": ["suisun-base", "suisun-suisun"],
        "SLR": ["slr-base", "slr-slr"],
    }

    altname_dict = {
        "suisun-base": "Base",
        "suisun-suisun": "Suisun Geometry",
        "slr-base": "Base",
        "slr-slr": "w/ Sea Level Rise",
    }
    meshname_dict = {
        "suisun-base": "base",
        "suisun-suisun": "suisun",
        "slr-base": "base",
        "slr-slr": "base",
    }
    ann_suffix_dict = {
        "suisun-base": "Base",
        "suisun-suisun": "Alt",
        "slr-base": "Base",
        "slr-slr": "Alt",
    }
    x2_path = "/CALSIM/X2_PRV_KM/X2-POSITION-PREV//1MON//"
    ec_path = "/CALSIM/{calsim_ec_loc}_EC_MONTH/SALINITY//1MON//"
    ec_std_path = "/CALSIM/{calsim_ec_loc}_EC_STD/SALINITY//1MON//"
    calsim_ec_loc_dict = {
        "bdl": "BD",
        # "wci": "CI",
        "cse": "CO",
        "emm2": "EM",
        "jer": "JP",
        "bac": "RS",
        # "vcu": "VI",
    }
    plot_ec_dict = {
        "bdl": "Beldon's Landing",
        # "wci": "Clifton Court Intake",
        "cse": "Collinsville",
        "emm2": "Emmaton",
        "jer": "Jersey Point",
        "bac": "Rock Slough (@ Bacon)",
        # "vcu": "Victoria Island",
    }
    # ec_order = ["bdl", "cse", "emm2", "jer", "bac", "vcu", "wci"]
    ec_order = ["bdl", "cse", "emm2", "jer", "bac"]
    ec_stds = ["cse", "emm2", "jer", "bac"]
    ec_vars_with_std = []
    for loc in ec_order:
        ec_vars_with_std.append(loc)
        if loc in ec_stds:
            ec_vars_with_std.append(f"{loc}_std")
    ec_subplot_titles = [
        f"{scenario} - {plot_ec_dict[loc]}"
        for scenario in alternative_dict.keys()
        for loc in ec_order
        if loc in plot_ec_dict
    ]
    # Plotting options ------------------------------
    color_map = {
        "CalSim-Base": "#1f77b4",
        "SCHISM-Base": "#2488a3",
        "CalSim-Alt": "#d62728",
        "SCHISM-Alt": "#d63627",
    }
    color_map["Observed"] = "gray"

    # Matplotlib plotting
    import matplotlib as mpl

    mpl.rcParams["figure.dpi"] = 300

    # =============================================================================================
    # Gather Data
    # Try to load alt_data from file
    alt_data_path = "alt_data.pkl"
    if os.path.exists(alt_data_path):
        print("Loading alt_data from pickle...")
        with open(alt_data_path, "rb") as f:
            alt_data = pickle.load(f)
    else:
        print("Gathering data")
        alt_data = {}
        for ann, alternatives in alternative_dict.items():
            for alternative in alternatives:
                print(f"\t{alternative}")
                # determine mesh used
                ann_suffix = ann_suffix_dict[alternative]
                mesh = meshname_dict[alternative]

                # set schism v calsim data paths
                calsim_dss = f"{dsp_repo}/model/calsim/rma_roundtrip/schism-{alternative}/DSS/output/DCR2023_DV_9.3.1_Danube_Adj_v1.8.dss"
                schism_fn = f"{casanntra_repo}/data/schism_{mesh}_{alternative}.csv"

                # Get all CalSim vars:
                cs_paths = [x2_path]
                for ec_var in ec_order:
                    cs_paths.append(
                        ec_path.format(calsim_ec_loc=calsim_ec_loc_dict[ec_var])
                    )
                    if ec_var in ec_stds:
                        cs_paths.append(
                            ec_std_path.format(calsim_ec_loc=calsim_ec_loc_dict[ec_var])
                        )
                schism_df = get_schism_results(schism_fn, ["x2"] + ec_order)
                calsim_df = get_calsim_results(
                    calsim_dss,
                    cs_paths,
                    ["x2"] + ec_vars_with_std,
                    start_date=schism_df.index[0],
                    end_date=schism_df.index[-1],
                )

                alt_data[alternative] = (schism_df, calsim_df)
            # Save alt_data for next time
        with open(alt_data_path, "wb") as f:
            pickle.dump(alt_data, f)
        print("alt_data saved to pickle.")

    # =============================================================================================
    # X2 plots (One for each ANN basis)
    for ann, alternatives in alternative_dict.items():
        print(f"Making X2 plots for ANN basis: {ann}...")
        fig, ax = plt.subplots(figsize=(figsize))
        out_fn = plt_fn.format(var="x2", ann_base=ann)
        ymin, ymax = 999, -999
        for alternative in alternatives:
            altname = altname_dict[alternative]
            ann_suffix = ann_suffix_dict[alternative]
            schism_df, calsim_df = alt_data[alternative]
            ymin, ymax = add_sch_cal_matplotlib(
                ax,
                calsim_df["x2"],
                schism_df["x2"],
                color_map,
                ylabel="X2 Location (km)",
                ann_suffix=ann_suffix,
                altname=altname,
                std_df=None,
                ymax=ymax,
            )
        ax.set_title("X2")
        ax.tick_params(axis="x", rotation=45)
        ax.set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
        ax.set_ylim([45, 95])
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        plt.tight_layout()
        print(f"\tSaving {out_fn}...")
        fig.savefig(out_fn)
        plt.close(fig)

    # =============================================================================================
    # EC plots (One for each ANN basis per loc)
    for var in ec_stds:
        print(f"Making plots for EC loc: {var}...")
        for ann, alternatives in alternative_dict.items():
            print(f"\tANN basis: {ann}...")
            fig, ax = plt.subplots(figsize=(figsize))
            out_fn = plt_fn.format(var=var, ann_base=ann)
            ymin, ymax = 999, -999
            plt_std = False
            for alternative in alternatives:
                altname = altname_dict[alternative]
                ann_suffix = ann_suffix_dict[alternative]
                schism_df, calsim_df = alt_data[alternative]
                std_df = ec_psu_25c(calsim_df[f"{var}_std"]) if plt_std else None
                cs_df = ec_psu_25c(calsim_df[var].copy())
                sch_df = ec_psu_25c(schism_df[var].copy())

                ymin, ymax = add_sch_cal_matplotlib(
                    ax,
                    cs_df,
                    sch_df,
                    color_map,
                    ylabel="Salinity (PSU)",
                    ann_suffix=ann_suffix,
                    altname=altname,
                    std_df=std_df,
                    ymax=ymax,
                )
                plt_std = True

            ax.set_title(plot_ec_dict[var])
            ax.tick_params(axis="x", rotation=45)
            ax.set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
            plt.tight_layout()
            print(f"\t\tSaving {out_fn}...")
            fig.savefig(out_fn)
            if not (var == ec_stds[-1]):
                plt.close(fig)
            else:
                # --- Legend handling ---
                handles, labels = get_legend_handles(ax)
                legend = plt.legend(
                    handles,
                    labels,
                    ncol=3,
                    loc="lower center",
                    bbox_to_anchor=(0.5, -0.5),
                    fontsize=9,
                    frameon=False,
                )
                export_legend(
                    legend,
                    filename=f"plots/legend_{ann}.png",
                )
                plt.close(fig)
                # plt.tight_layout(pad=1.0)
