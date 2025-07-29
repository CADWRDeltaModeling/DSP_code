import os
import pandas as pd
import numpy as np
import datetime as dt

from vtools.functions.interpolate import rhistinterp
from vtools.functions.read_dss import read_dss
from vtools.data.vtime import months

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from itertools import cycle
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

import plotly.io as pio

pio.renderers.default = "browser"


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


def add_plot_traces(
    fig, calsim_df, color, hover, row, y_axis_title, std_df=None, dtick="M12"
):
    ymax = max(calsim_df.iloc[:, 0])

    show_legend = row == 1
    # Add the CalSim Standard EC value
    if std_df is not None:
        y_fill = std_df.iloc[:, 0]

        # First trace: invisible line at y_max
        fill_trace_top = go.Scatter(
            x=std_df.index,
            y=[ymax * 1.05] * len(std_df.index),
            mode="lines",
            line=dict(width=0),
            showlegend=False,
            hoverinfo="skip",
            name="",
        )
        fig.add_trace(fill_trace_top, row=row, col=1)

        # Second trace: std_df line, fill to previous trace
        fill_trace = go.Scatter(
            x=std_df.index,
            y=y_fill,
            mode="lines",
            line=dict(width=0),
            line_shape="hv",
            fill="tonexty",
            fillcolor="rgba(200,200,200,0.5)",  # light grey with opacity
            showlegend=False,
            hoverinfo="skip",
            name="",
        )
        fig.add_trace(fill_trace, row=row, col=1)

        # Now add your std_df line trace as usual
        trace = go.Scatter(
            x=std_df.index,
            y=std_df.iloc[:, 0],
            name=f"CalSim Standard",
            line=dict(color="gray"),
            line_shape="hv",
            mode="lines",
            hovertemplate=f"CalSim Standard<br>%{{y:.2f}} {hover[1]}<extra></extra>",
            showlegend=show_legend,
            legendgroup=f"CalSim Standard",
        )
        fig.add_trace(trace, row=row, col=1)

    # Add the CalSim predicted value
    trace = go.Scatter(
        x=calsim_df.index,
        y=calsim_df.iloc[:, 0],
        name=f"CalSim",
        line=dict(color=color),
        line_shape="hv",
        mode="lines",
        hovertemplate=f"CalSim<br>%{{y:.2f}} {hover[1]}<extra></extra>",
        showlegend=show_legend,  # Hide the legend
        legendgroup=f"CalSim",
    )
    fig.add_trace(trace, row=row, col=1)

    # Change the maximum y axis value:
    if y_axis_title == "X2 (km)":
        yrange = [45, 95]
    else:
        yrange = [0, ymax * 1.05]
    fig.update_yaxes(
        range=yrange, row=row, col=1, title_text=y_axis_title, fixedrange=True
    )  # Set the ymax
    fig.update_xaxes(showticklabels=True, dtick=dtick, row=row, col=1)

    return fig


def set_plot_options(fig, height=4000, title_text="Roundtrip CalSim"):
    # Finalize plotting options
    fig.update_xaxes(
        showspikes=True,
        spikecolor="black",
        spikethickness=1,
        showgrid=True,
        gridcolor="lightgrey",  # Major gridlines
        zeroline=False,
        minor=dict(showgrid=True, gridcolor="gainsboro"),  # Minor gridlines
        title_text="",
    )
    # Update layout
    fig.update_layout(
        height=height,
        title_text=title_text,
        showlegend=True,
        hoversubplots="axis",
        hovermode="x unified",
        plot_bgcolor="#f5f5f5",  # <-- Pale/white background
    )

    return fig


if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    dsp_repo = "D:/projects/delta_salinity/DSP_code"  # Local: D:/projects/delta_salinity/DSP_code, Remote: /scratch/tomkovic/DSP_code
    html_name = f"./plots/roundtrip_yr_select.html"
    compliance_figname = f"./plots/roundtrip_compliance.png"

    # set schism v calsim data paths
    calsim_dss = f"{dsp_repo}/model/calsim/9.3.1_danube_adj/DSS/output/DCR2023_DV_9.3.1_v2a_Danube_Adj_v1.8.dss"

    x2_path = "/CALSIM/X2_PRV_KM/X2-POSITION-PREV//1MON//"
    ndo_path = "/CALSIM/NDO/FLOW-NDO//1MON//"
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
    subplot_titles = ["X2", "Net Delta Outflow"] + [
        f"{plot_ec_dict[loc]}" for loc in ec_order if loc in plot_ec_dict
    ]
    doc_subplot_titles = [
        f"{plot_ec_dict[loc]}" for loc in ec_stds if loc in plot_ec_dict
    ]
    # Plotting options ------------------------------
    color = "#1f77b4"
    # Create subplots
    n_rows = len(subplot_titles) + 2
    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        subplot_titles=subplot_titles,
    )
    doc_fig_data = []

    # Get all CalSim vars:
    cs_paths = [x2_path, ndo_path]
    for ec_var in ec_order:
        cs_paths.append(ec_path.format(calsim_ec_loc=calsim_ec_loc_dict[ec_var]))
        if ec_var in ec_stds:
            cs_paths.append(
                ec_std_path.format(calsim_ec_loc=calsim_ec_loc_dict[ec_var])
            )
    calsim_df = get_calsim_results(
        calsim_dss,
        cs_paths,
        ["x2", "ndo"] + ec_vars_with_std,
        start_date=pd.to_datetime("2008-1-1"),
    )

    row = 1  # Initialize row index for subplot
    doc_row = 1
    for var in ["x2", "ndo"] + ec_order:
        print(f"\t\tProcessing var: {var}")
        std_df = None

        if var in ec_order:
            hover_text = ["EC", "psu"]
            y_axis_title = "Salinity (psu)"
            if var in ec_stds:
                main_df = calsim_df[var].copy()
                std_df = calsim_df[f"{var}_std"].copy()
                title = plot_ec_dict[var]
                doc_fig_data.append((main_df, std_df, title))

        elif var == "x2":
            hover_text = ["X2", "km"]
            y_axis_title = "X2 (km)"
        elif var == "ndo":
            hover_text = ["NDO", "cfs"]
            y_axis_title = "NDO (cfs)"
        fig = add_plot_traces(
            fig,
            calsim_df[var],
            color,
            hover_text,
            row,
            y_axis_title,
            std_df=std_df,
        )
        row += 1
    fig = set_plot_options(fig)

    # Create compliance points documentation plot
    doc_fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(12, 8))
    for ax, (main_df, std_df, title) in zip(axes, doc_fig_data):
        ymax = max(main_df.iloc[:, 0]) * 1.05
        # Plot std fill (as a band)
        ax.fill_between(
            std_df.index,
            std_df.iloc[:, 0],
            ymax,
            color="gray",
            alpha=0.3,
            label="CalSim Standard",
        )
        # Plot main line
        ax.plot(main_df.index, main_df.iloc[:, 0], color="#1f77b4", label="CalSim")
        ax.set_ylim(0, ymax)
        ax.set_title(title)
        ax.grid(True)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    # Shared y-label
    doc_fig.text(
        0.04, 0.5, "Salinity (psu)", va="center", rotation="vertical", fontsize=14
    )
    plt.rcParams["font.family"] = "serif"  # Example: set a generic serif font
    plt.rcParams["font.serif"] = ["Times New Roman"]
    plt.tight_layout(rect=[0.05, 0, 1, 1])

    # Export to HTML
    if True:
        # print(f"Exporting plot to HTML:\n{os.path.abspath(html_name)}")
        # fig.write_html(html_name)
        print(f"Exporting plot to PNG:\n{os.path.abspath(compliance_figname)}")
        plt.savefig(compliance_figname, dpi=150)
    else:
        fig.show()
        plt.show()
