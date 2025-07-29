import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from clustergram import Clustergram
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as mtransforms

OMP_NUM_THREADS = 1


# Functions --------------------------------------------------------------
# Custom reshape function
def reshape_df(df, monthstart=10, interval="daily"):
    # df is a dataframe with a "Time" column
    df = df.copy()
    # Ensure 'Time' column is datetime
    df["Time"] = pd.to_datetime(df["Time"])
    colnms = df.columns

    # Add year, month, and initialize water year
    df["year"] = df["Time"].dt.year
    df["month"] = df["Time"].dt.month
    df["wy"] = np.nan

    # Assign water year
    for year in df["year"].unique():
        # months >= monthstart belong to next water year
        mask1 = (df["year"] == year) & (df["month"] >= monthstart)
        mask2 = (df["year"] == year + 1) & (df["month"] < monthstart)
        df.loc[mask1, "wy"] = int(year + 1)
        df.loc[mask2, "wy"] = int(year + 1)

    # Remove rows without water year
    df = df[~df["wy"].isna()]
    df["wy"] = df["wy"].astype(int)

    if interval == "daily":
        # Julian day within water year
        df["wyjday"] = np.nan
        for year in df["year"].unique():
            wy_start = pd.Timestamp(year=year, month=monthstart, day=1)
            mask1 = (df["year"] == year) & (df["month"] >= monthstart)
            mask2 = (df["year"] == year + 1) & (df["month"] < monthstart)
            df.loc[mask1, "wyjday"] = 1 + (df.loc[mask1, "Time"] - wy_start).dt.days
            df.loc[mask2, "wyjday"] = 1 + (df.loc[mask2, "Time"] - wy_start).dt.days

        # Reshape: water year as rows, julian day as columns
        value_col = colnms[1]  # assumes value is in second column
        df_melt = df[["wyjday", "wy", value_col]].dropna()
        df_wide = df_melt.pivot(index="wy", columns="wyjday", values=value_col)
        # Remove the columns index name and ensure columns are integers
        df_wide.columns.name = None
        df_wide.columns = df_wide.columns.astype(int)
        # Remove leap day (366)
        if 366 in df_wide.columns:
            df_wide = df_wide.drop(columns=366)
        # Remove years with missing data
        df_wide = df_wide.dropna()
        return df_wide

    elif interval == "monthly":
        # Month-day index within water year
        df["wymoday"] = np.nan
        for year in df["year"].unique():
            wy_start = pd.Timestamp(year=year, month=monthstart, day=1)
            mask1 = (df["year"] == year) & (df["month"] >= monthstart)
            mask2 = (df["year"] == year + 1) & (df["month"] < monthstart)
            # Leap year adjustment for March onward
            adj_day = (
                1 if pd.Timestamp(year=year + 1, month=1, day=1).is_leap_year else 0
            )
            # Jan/Feb: no leap year adjustment
            mask2_janfeb = mask2 & (df["month"] < 3)
            mask2_marplus = mask2 & (df["month"] >= 3)
            df.loc[mask1, "wymoday"] = 1 + (df.loc[mask1, "Time"] - wy_start).dt.days
            df.loc[mask2_janfeb, "wymoday"] = (
                1 + (df.loc[mask2_janfeb, "Time"] - wy_start).dt.days
            )
            df.loc[mask2_marplus, "wymoday"] = (
                1 - adj_day + (df.loc[mask2_marplus, "Time"] - wy_start).dt.days
            )

        value_col = colnms[1]
        df_melt = df[["wymoday", "wy", value_col]].dropna()
        df_wide = df_melt.pivot(index="wy", columns="wymoday", values=value_col)
        df_wide = df_wide.dropna()
        return df_wide

    else:
        raise ValueError("interval must be 'daily' or 'monthly'")


# Normalize flow and export data so variation is equally weighted
def norm_biv_df(df_rshp, nbv=2):
    """
    Normalize bivariate reshaped dataframe as in R's norm.biv.df.
    Assumes df_rshp columns: 'wy', then time series for each variable.
    """
    df_norm = pd.DataFrame(index=df_rshp.index, columns=df_rshp.columns)

    ncols = df_rshp.shape[1]
    # Exclude 'wy' column for calculation
    ts_cols = [col for col in df_rshp.columns if col != "wy"]
    n_ts = len(ts_cols) // nbv

    for b in range(nbv):
        # Get column indices for this variable
        start = b * n_ts
        end = (b + 1) * n_ts
        colinds = df_rshp.columns[start:end]
        ts = df_rshp.loc[:, colinds]

        divave_ts = ts / ts.values.mean()
        log_ts = np.log(divave_ts + 1e-6)
        min_bv = log_ts.min().min()
        max_bv = log_ts.max().max()
        norm_ts = (log_ts - min_bv) / (max_bv - min_bv)

        df_norm.loc[:, colinds] = norm_ts

    return df_norm


def prepare_cluster_data(
    dsm2_filename,
    start_date=pd.to_datetime("2001-05-01"),
    end_date=pd.to_datetime("2001-11-30"),
    refdate=pd.to_datetime("2000-11-30"),
):
    # Read in data ------------------------------------------------------------
    northern_flow = pd.read_excel(dsm2_filename, sheet_name="northern_flow")
    exports = pd.read_excel(dsm2_filename, sheet_name="exports")

    # Pre-prepare Data --------------------------------------------------------
    # reshape exports
    exp_df_rshp = reshape_df(exports, monthstart=12)
    exp_per_rshp = exp_df_rshp.loc[
        :,
        exp_df_rshp.columns.isin(
            ["wy"]
            + list(range((start_date - refdate).days, (end_date - refdate).days + 1))
        ),
    ]

    # reshape northern flow from December
    nfl_df_rshp = reshape_df(northern_flow, monthstart=12)
    nfl_per_rshp = nfl_df_rshp.loc[
        :,
        nfl_df_rshp.columns.isin(
            ["wy"]
            + list(range((start_date - refdate).days, (end_date - refdate).days + 1))
        ),
    ]

    # Combine exports and flows (North Flow EXports: nfex)
    nfex_full_rshp = pd.merge(nfl_df_rshp, exp_df_rshp, on="wy")
    df_rshp = pd.merge(nfl_per_rshp, exp_per_rshp, on="wy")

    df_norm = norm_biv_df(df_rshp, 2)
    df_norm_stat = df_norm.loc[:, ~df_norm.columns.isin(["wy"])]

    return df_norm_stat, df_norm, df_rshp


def plot_clustergram(df_norm_stat, plt_out_fn):
    cgram = Clustergram(range(1, 12), random_state=42)
    cgram.fit(df_norm_stat)
    sns.set_theme(style="whitegrid")

    cgram.plot(
        size=0.3,
        linewidth=0.2,
        cluster_style={"color": "grey", "edgecolor": "black"},
        line_style={"color": "#303030"},
    )
    plt.xticks(ticks=range(1, 12), labels=range(1, 12))

    # Add anotations:
    # Selection of 4 clusters
    ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    rect = patches.Rectangle(
        (3.5, ylim[0] * 1.2),  # x, y (left edge of cluster 4)
        1,  # width (one cluster column)
        ylim[1] * 0.98 - ylim[0] * 1.2,  # height
        linewidth=2,
        edgecolor="red",
        facecolor="none",
    )
    ax.add_patch(rect)

    # annotation on side of plot
    x_arrow = xlim[1] + 0.5  # adjust as needed
    y_start = ylim[0]
    y_end = ylim[1]
    ax.annotate(
        "",
        xy=(x_arrow, y_end),
        xytext=(x_arrow, y_start),
        arrowprops=dict(arrowstyle="<->", color="blue", linewidth=3),
        annotation_clip=False,
    )
    ax.text(
        x_arrow + 0.1,
        (y_start + y_end) / 2,
        "Difference in characteristics\nbetween clustered groups",
        va="center",
        ha="left",
        rotation=90,
        color="black",
        fontsize=10,
    )
    # To add arrow within plot:
    # ax.set_xlim(xlim[0], xlim[1] + 2)

    # Change x and y tick font size
    plt.xticks(fontsize=8)  # or your desired size
    plt.yticks(fontsize=8)

    # Change axis label font size
    ax.set_xlabel("Number of clusters (k)", fontsize=10)
    ax.set_ylabel("PCA weighted mean of the clusters", fontsize=10)

    # Change size of plot
    plt.gcf().set_size_inches(6.5, 3.5)
    plt.tight_layout()

    # plt.show()
    plt.savefig(plt_out_fn)
    # print("hi")

    from bokeh.plotting import output_file, save

    fig = cgram.bokeh(
        cluster_style={"color": "grey", "line_color": "black"},
        line_style={"color": "#303030"},
    )
    output_file(plt_out_fn.replace(".png", ".html"))
    save(fig)


def select_cluster_ts(df_norm_stat, df_norm, df_rshp, nclust, random_state=42):
    # Run k-means clustering
    kmeans = KMeans(n_clusters=nclust, random_state=random_state)
    clusts = kmeans.fit_predict(df_norm_stat)

    # Assign clusters to each water year
    wy_metrics = pd.DataFrame(
        {
            "wy": df_rshp.index,
            "clust": clusts,
            "skill": np.nan,
            "rmse": np.nan,
            "mmskill": None,
            "mmrmse": None,
        }
    )

    df_norm["clust"] = clusts
    jds = [col for col in df_norm_stat.columns]

    # Calculate cluster means and metrics
    for c in np.sort(np.unique(clusts)):
        # Get average values for cluster
        mean_clust = df_norm_stat[df_norm["clust"] == c].mean(axis=0)
        for wy in wy_metrics.loc[wy_metrics["clust"] == c, "wy"]:
            wy_pred = df_norm_stat.loc[df_rshp.index == wy, jds].values.flatten()
            rmse = np.sqrt(np.mean((wy_pred - mean_clust.values) ** 2))
            skill = 1 - (
                np.sum((wy_pred - mean_clust.values) ** 2)
                / np.sum(
                    (
                        np.abs(wy_pred - mean_clust.mean())
                        + np.abs(mean_clust.values - mean_clust.mean())
                    )
                    ** 2
                )
            )
            wy_metrics.loc[wy_metrics["wy"] == wy, "rmse"] = rmse
            wy_metrics.loc[wy_metrics["wy"] == wy, "skill"] = skill

        # Compute max and min skill/rmse
        cluster_mask = wy_metrics["clust"] == c
        # Max skill (for years > 2007)
        skill_mask = cluster_mask & (wy_metrics["wy"] > 2007)
        if skill_mask.any():
            max_skill = wy_metrics.loc[skill_mask, "skill"].max()
            wy_metrics.loc[
                (wy_metrics["skill"] == max_skill) & skill_mask, "mmskill"
            ] = f"Max Skill cluster {c}"
        min_skill = wy_metrics.loc[cluster_mask, "skill"].min()
        wy_metrics.loc[(wy_metrics["skill"] == min_skill) & cluster_mask, "mmskill"] = (
            f"Min Skill cluster {c}"
        )

        max_rmse = wy_metrics.loc[cluster_mask, "rmse"].max()
        wy_metrics.loc[(wy_metrics["rmse"] == max_rmse) & cluster_mask, "mmrmse"] = (
            f"Max RMSE cluster {c}"
        )
        # Min RMSE (for years > 2007)
        rmse_mask = cluster_mask & (wy_metrics["wy"] > 2007)
        if rmse_mask.any():
            min_rmse = wy_metrics.loc[rmse_mask, "rmse"].min()
            wy_metrics.loc[(wy_metrics["rmse"] == min_rmse) & rmse_mask, "mmrmse"] = (
                f"Min RMSE cluster {c}"
            )

    # Sort and select water years with minimum RMSE
    wy_metrics = wy_metrics.sort_values(["clust", "rmse"])
    wy_slcts = wy_metrics.loc[
        wy_metrics["mmrmse"].str.contains("Min RMSE", na=False), "wy"
    ]
    wy_slcts = wy_slcts[~wy_slcts.isin([2019])]  # Exclude 2019 if needed
    wy_slcts = np.sort(wy_slcts)

    return wy_metrics, wy_slcts


# def plot_cluster_ts(
#     df_rshp,
#     wy_metrics,
#     wy_slcts,
#     plot_clust_fn,
#     cluster_colors=[
#         "navy",
#         "royalblue",
#         "mediumorchid",
#         "orangered",
#         "palegreen",
#         "goldenrod",
#         "green",
#         "paleturquoise",
#         "yellow",
#         "pink",
#         "slateblue",
#     ],

#     refdate=pd.to_datetime("2000-11-30"),
# ):
#     # --- Plot final selection ------------------------------------------------
#     col_clust = pd.DataFrame(
#         {"clust": np.arange(1, len(cluster_colors) + 1), "color": cluster_colors}
#     )
#     wy_to_clust = wy_metrics[["wy", "clust"]]
#     wycol = pd.merge(wy_to_clust, col_clust, on="clust", how="left")

#     # Prepare data for plotting (reshape as needed)
#     # Here, you would reconstruct a long-form DataFrame for plotting, similar to df.plt.rshp in R

#     # Example plot: highlight selected years in color, others in grey
#     plt.figure(figsize=(12.5, 7))
#     for wy in df_rshp.index:
#         color = "grey"
#         linewidth = 0.5
#         alpha = 0.5
#         if wy in wy_slcts:
#             color = wycol.loc[wycol["wy"] == wy, "color"].values[0]
#             linewidth = 1
#             alpha = 1
#         plt.plot(
#             df_rshp.columns,
#             df_rshp.loc[df_rshp.index == wy].values[0][:],
#             color=color,
#             linewidth=linewidth,
#             alpha=alpha,
#         )

#     plt.xlabel("Julian Day")
#     plt.ylabel("Discharge (cfs)")
#     # plt.title("Final Water Year Selection")
#     plt.tight_layout()
#     plt.savefig(plot_clust_fn)


if __name__ == "__main__":

    # Set working directory to this file's directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    dsm2_filename = "../../data/dsm2_ann_inputs_base.xlsx"
    plt_out_fn = "./plots/python_clustergram.png"

    df_norm_stat, df_norm, df_rshp = prepare_cluster_data(dsm2_filename)

    if False:
        plot_clustergram(df_norm_stat, plt_out_fn)

    nclust = 4  # Using 4 clusters as most appropriate

    if False:
        random_state = 0
        slct_yrs = True
        while slct_yrs:
            random_state += 1
            wy_metrics, wy_slcts = select_cluster_ts(
                df_norm_stat, df_norm, df_rshp, nclust, random_state=random_state
            )
            print(f"With random_state: {random_state}, produced: {wy_slcts}")

            if (len(wy_slcts) == 4) and (
                wy_metrics.loc[
                    wy_metrics["wy"].isin([2008, 2010, 2012, 2014]), "clust"
                ].nunique()
                == 4
            ):
                break
    else:
        wy_metrics, wy_slcts = select_cluster_ts(
            df_norm_stat, df_norm, df_rshp, nclust, random_state=10
        )
        wy_slcts = [2008, 2010, 2012, 2014]

    plot_clust_fn = f"./plots/python_cluster_ts_{nclust}clusters.png"

    plot_cluster_ts(df_rshp, wy_metrics, wy_slcts, plot_clust_fn)
