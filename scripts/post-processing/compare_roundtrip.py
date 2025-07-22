import os
import pandas as pd
import numpy as np
import datetime as dt

from vtools.functions.interpolate import rhistinterp
from vtools.functions.read_dss import read_dss
from vtools.data.vtime import months

from itertools import cycle
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px


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


def add_sch_cal_plot_traces(
    fig, calsim_df, schism_df, color_map, hover, row, ann_suffix, std_df=None
):

    axis_key = "yaxis" if row == 1 else f"yaxis{row}"
    axis_range = fig.layout[axis_key].range
    figmax = axis_range[1] if axis_range is not None else 0
    ymax = max(schism_df.max().max(), calsim_df.max().max(), figmax)

    # Add the CalSim Standard EC value
    if std_df is not None:
        color = "gray"
        y_fill = std_df.iloc[:, 0]

        # First trace: invisible line at y_max
        fill_trace_top = go.Scatter(
            x=std_df.index,
            y=[ymax * 5] * len(std_df.index),
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
            line=dict(color=color),
            line_shape="hv",
            mode="lines",
            hovertemplate=f"CalSim Standard<br>%{{y:.2f}} {hover[1]}<extra></extra>",
            showlegend=show_legend,
            legendgroup=f"CalSim Standard",
        )
        fig.add_trace(trace, row=row, col=1)

    # Add the CalSim predicted value
    color = color_map.get(f"CalSim-{ann_suffix}", "gray")
    trace = go.Scatter(
        x=calsim_df.index,
        y=calsim_df.iloc[:, 0],
        name=f"CalSim-{ann_suffix}",
        line=dict(color=color),
        line_shape="hv",
        mode="lines",
        hovertemplate=f"CalSim-{ann_suffix}<br>%{{y:.2f}} {hover[1]}<extra></extra>",
        showlegend=show_legend,  # Hide the legend
        legendgroup=f"CalSim-{ann_suffix}",
    )
    fig.add_trace(trace, row=row, col=1)

    # Add the SCHISM predicted value
    color = color_map.get(f"SCHISM-{ann_suffix}", "gray")
    trace = go.Scatter(
        x=schism_df.index,
        y=schism_df.iloc[:],
        name=f"SCHISM-{ann_suffix}",
        line=dict(color=color, dash="dashdot"),
        line_shape="hv",
        mode="lines",
        hovertemplate=f"SCHISM-{ann_suffix}<br>%{{y:.2f}} km<extra></extra>",
        showlegend=show_legend,  # Hide the legend
        legendgroup=f"SCHISM-{ann_suffix}",
    )
    fig.add_trace(trace, row=row, col=1)

    # Change the maximum y axis value:
    fig.update_yaxes(range=[0, ymax * 1.05], row=row, col=1)  # Set the ymax
    fig.update_xaxes(showticklabels=True, dtick="M1", row=row, col=1)

    return fig


def set_plot_options(fig, y_axis_title, height=4000):
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
    fig.update_yaxes(
        showspikes=True,
        spikecolor="black",
        spikethickness=1,
        showgrid=True,
        gridcolor="lightgrey",  # Major gridlines
        zeroline=True,
        minor=dict(showgrid=True, gridcolor="gainsboro"),  # Minor gridlines
        title_text=y_axis_title,
    )
    # Update layout
    fig.update_layout(
        height=height,
        title_text="Roundtrip CalSim-SCHISM",
        showlegend=True,
        hoversubplots="axis",
        hovermode="x unified",
        plot_bgcolor="#f5f5f5",  # <-- Pale/white background
    )

    return fig


if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))
    dsp_repo = "D:/projects/delta_salinity/DSP_code"  # Local: D:/projects/delta_salinity/DSP_code, Remote: /scratch/tomkovic/DSP_code
    casanntra_repo = "D:/projects/delta_salinity/scripts/casanntra"  # Local: D:/projects/delta_salinity/scripts/casanntra, Remote: /scratch/tomkovic/casanntra

    alternative_dict = {
        "Suisun": ["suisun-base", "suisun-suisun"],
        "SLR": ["slr-base", "slr-slr"],
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
        "SCHISM-Base": "#1fb4a0",
        "CalSim-Alt": "#d62728",
        "SCHISM-Alt": "#c94c1e",
    }
    color_map["Observed"] = "gray"
    # Create subplots
    x2fig = make_subplots(
        rows=len(alternative_dict.keys()),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        subplot_titles=list(alternative_dict.keys()),
    )
    x2_html_name = f"./plots/roundtrip_x2_plot.html"
    ecfig = make_subplots(
        rows=len(ec_subplot_titles),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.01,
        subplot_titles=ec_subplot_titles,
    )
    ec_html_name = f"./plots/roundtrip_ec_plot.html"

    row = 0  # Initialize row index for subplot

    # After creating ecfig and before exporting/showing
    n_rows = len(ec_subplot_titles)
    divider_row = n_rows // 2

    # yref "paper" goes from 0 (bottom) to 1 (top)
    divider_y = 1 - divider_row / n_rows

    ecfig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=0,
        x1=1,
        y0=divider_y,
        y1=divider_y,
        line=dict(color="black", width=6),
    )

    for ann, alternatives in alternative_dict.items():
        print(f"ANN base: {ann}")
        ec_plt_rows = [
            i + 1 for i, title in enumerate(ec_subplot_titles) if title.startswith(ann)
        ]
        row += 1

        for alternative in alternatives:
            print(f"\tANN subset: {alternative}")
            plot_std = False

            # determine mesh used
            if alternative.split("-")[-1] == "suisun":
                mesh = "suisun"
                ann_suffix = "Alt"
            elif alternative.split("-")[-1] == "base":
                ann_suffix = "Base"
                mesh = "base"
                plot_std = True
            else:
                ann_suffix = "Alt"
                mesh = "base"

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
            for var in ["x2"] + ec_order:
                print(f"\t\tProcessing var: {var}")
                std_df = None

                if var == "x2":
                    show_legend = row == 1
                    x2fig = add_sch_cal_plot_traces(
                        x2fig,
                        calsim_df[var],
                        schism_df[var],
                        color_map,
                        ["X2", "km"],
                        row,
                        ann_suffix,
                        std_df=std_df,
                    )
                else:
                    if var in ec_stds and plot_std:
                        std_df = calsim_df[f"{var}_std"]
                    ec_row = ec_plt_rows[ec_order.index(var)]
                    show_legend = ec_row == 1
                    ecfig = add_sch_cal_plot_traces(
                        ecfig,
                        calsim_df[var],
                        schism_df[var],
                        color_map,
                        ["EC", "psu"],
                        ec_row,
                        ann_suffix,
                        std_df=std_df,
                    )

    x2fig = set_plot_options(x2fig, "X2 Location (km)", height=1000)
    x2fig.update_yaxes(range=[45, 95])  # Set the ymax
    ecfig = set_plot_options(ecfig, "Salinity (psu)")

    # Export to HTML
    if True:
        print(f"Exporting plot to HTML:\n{os.path.abspath(x2_html_name)}")
        x2fig.write_html(x2_html_name)
        print(f"Exporting plot to HTML:\n{os.path.abspath(ec_html_name)}")
        ecfig.write_html(ec_html_name)
    else:
        import plotly.io as pio

        pio.renderers.default = "browser"
        x2fig.show()
        ecfig.show()
