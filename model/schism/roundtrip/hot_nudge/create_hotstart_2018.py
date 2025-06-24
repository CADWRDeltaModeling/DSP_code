# -*- coding: utf-8 -*-
"""
Hotstart example for a basic schism run with TEMP and SALT as tracers.
"""

import schimpy.schism_hotstart as sh
import matplotlib.pyplot as plt
import os

yaml_fn = "./hotstart_2018.yaml"
modules = ["TEM", "SAL"]
hotstart_fn = "hotstart.20180315.nc"  # output hotstart file

# create a hotstart file for SCHISM
h = sh.hotstart(yaml_fn)
h.create_hotstart()
hnc = h.nc_dataset
hnc.to_netcdf(os.path.join(h.out_dir, hotstart_fn))

# %% making a 2D surface plot
coll = h.mesh.plot_elems(hnc["tr_el"].values[:, 0, 0], clim=(14, 18))  # clim=[0,35])
cb = plt.colorbar(coll)
plt.axis("off")
plt.axis("equal")
plt.title("Regional Temperature")
plt.tight_layout(pad=1)

# %% converting hotstart file to schism output format so that it can be viewd by VisIt
# sh.hotstart_to_outputnc(
#     hotstart_fn,
#     str(h.date),
#     hgrid_fn="../rt_v1/hgrid.gr3",
#     vgrid_fn="../rt_v1/vgrid.in.3d",
#     vgrid_version=h.vgrid_version,
#     outname="schout_hotstart.nc",
# )
