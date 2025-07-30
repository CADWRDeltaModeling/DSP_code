import os
import pandas as pd

from bdschism.parse_cu import calc_net_calsim, calc_net_schism, calc_net_dsm2

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.cm as cm

os.chdir(os.path.dirname(__file__))


alternative = "suisun-base"
sch_dir = f"../../model/schism/roundtrip/{alternative}"
bds_sch_dir = "/home/tomkovic/BayDeltaSCHISM/data/channel_depletion"
bds_vs_fmt = "{vs}_dated.th"
cs_dss_fn = f"../../model/calsim/rma_roundtrip/schism-{alternative}/DSS/output/DCR2023_DV_9.3.1_Danube_Adj_v1.8.dss"
dcd_dss_file = "../../model/dsm2/2021DSM2FP_202301/timeseries/DCD_hist_Lch5.dss"

alt_dict = {
    "suisun": ["suisun-base", "suisun-suisun"],
    "slr": ["slr-base", "slr-slr"],
}
start_date = "2015-02-18"
end_date = "2016-05-01"

vs_fmt = "{vs}.{alternative}.dated.th"

figname = "./plots/consumptive_use_calsim_schism.png"
sch_kwargs = {
    "vsource": bds_vs_fmt.format(vs="vsource"),
    "vsink": bds_vs_fmt.format(vs="vsink"),
}
# Historical net cu
hist_net = calc_net_schism(bds_sch_dir, **sch_kwargs)

# Plot alternatives!
color_map = plt.get_cmap("Dark2")
# fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12.7, 5))
fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(12.7, 5))
ax.plot(hist_net.index, hist_net.values, color="grey", label="Historical SCHISM")
ax.set_title("Net Consumptive Use")
# axes[0].plot(hist_net.index, hist_net.values, color="grey", label="Historical SCHISM")
# axes[0].set_title("Net Consumptive Use")
# axes[1].set_title("Difference From Historical")

print(f"Gathering {alternative}...")

# SCHISM
sch_kwargs = {
    "vsource": vs_fmt.format(alternative=alternative, vs="vsource"),
    "vsink": vs_fmt.format(alternative=alternative, vs="vsink"),
}
sch_net = calc_net_schism(sch_dir, **sch_kwargs)
sch_diff = sch_net.copy() - hist_net.copy()

# CalSim
cs_net = calc_net_calsim(
    cs_dss_fn, start_date=start_date, end_date=end_date, **{"p": 20}
)
cs_diff = cs_net.copy() - hist_net.copy()

# DSM2
dsm2_net = calc_net_dsm2(dcd_dss_file)
cs_dsm2_diff = cs_net.copy() - dsm2_net.copy()


print(f"Plotting {alternative}....")
# Add raw plots
ax.plot(
    sch_net.index,
    sch_net.values,
    color=color_map(0),
    label="Roundtrip SCHISM Input",
)
ax.plot(
    dsm2_net.index,
    dsm2_net.values,
    color=color_map(1),
    label="DSM2",
)
ax.plot(
    cs_net.index,
    cs_net.values,
    color=color_map(2),
    label="CalSim",
)
# # Add difference plots
# axes[1].plot(
#     sch_diff.index,
#     sch_diff.values,
#     color=color_map(0),
#     label="Roundtrip SCHISM - Historical SCHISM",
# )
# axes[1].plot(
#     cs_dsm2_diff.index,
#     cs_dsm2_diff.values,
#     color=color_map(1),
#     linestyle="--",
#     label="Roundtrip CalSim - Historical DSM2",
# )
# axes[1].plot(
#     cs_diff.index,
#     cs_diff.values,
#     color=color_map(2),
#     label="Roundtrip CalSim - Historical SCHISM",
# )
# for ax in axes:
ax.axhline(0, color="black", linewidth=1.5, zorder=0)

if False:
    sch_diff.to_csv(f"./data_out/cu_difference_roundtrip.csv")

leg_loc = ["upper left", "lower left"]
i = 0
# for i, ax in enumerate(axes):
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc=leg_loc[i])
ax.grid(True)
ax.xaxis.set_major_locator(mdates.MonthLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
ax.tick_params(axis="x", rotation=45)
ax.set_xlim(pd.to_datetime(start_date), pd.to_datetime(end_date))
ax.set_ylim([-4000, 7500])
# Shared y-label
fig.text(0.04, 0.5, "Flow (cfs)", va="center", rotation="vertical", fontsize=14)
plt.rcParams["font.family"] = "serif"  # Example: set a generic serif font
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.tight_layout(rect=[0.05, 0, 1, 1])
plt.savefig(figname, dpi=250)
# plt.show()
