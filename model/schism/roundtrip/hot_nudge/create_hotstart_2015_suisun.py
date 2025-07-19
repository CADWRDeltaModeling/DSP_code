# -*- coding: utf-8 -*-
"""
Hotstart example for a basic schism run with TEMP and SALT as tracers.
"""

import schimpy.schism_hotstart as sh
import matplotlib.pyplot as plt
import os

yaml_fn = "./hotstart_2015_suisun.yaml"
modules = ["TEM", "SAL"]
hotstart_fn = "hotstart.20150218.nc"  # output hotstart file

# create a hotstart file for SCHISM
h = sh.hotstart(yaml_fn)
h.create_hotstart()
hnc = h.nc_dataset
hnc.to_netcdf(os.path.join(h.out_dir, hotstart_fn))
