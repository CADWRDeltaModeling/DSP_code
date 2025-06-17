#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import datetime as dtm
from schimpy.model_time import multi_file_to_elapsed

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# User will need to change the start time and location of BayDeltaSCHISM
start = dtm.datetime(2018, 3, 15)
outdir = "/scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1"
mod_home = "/scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1"
repo_home = "/home/tomkovic/BayDeltaSCHISM/data/"


def custom_name_transform(input_name):
    # Extract the base name of the file
    base_name = os.path.basename(input_name)

    # Remove "_dated" if it exists at the end of the base name
    if "_dated" in base_name:
        base_name = base_name.replace("_dated", "")

    # Replace "repo" with "output" and append "_elapsed_20160427" before the extension
    transformed_name = base_name.replace("repo", "output").replace(
        ".th", "_elapsed_20180315.th"
    )
    return transformed_name


th_repo = os.path.join(mod_home, "*rt_v1.dated.th")
multi_file_to_elapsed(th_repo, outdir, start, name_transform=custom_name_transform)

th_files = [
    os.path.join(os.path.join(repo_home, "time_history"), fname)
    for fname in [
        "grantline_barrier.th",
        "grantline_culvert.th",
        "grantline_weir.th",
        "midr_culvert_l.th",
        "midr_culvert_r.th",
        "midr_weir.th",
        "oldr_head_barrier.th",
        "oldr_tracy_barrier.th",
        "oldr_tracy_culvert.th",
        "oldr_tracy_weir.th",
        "salt.th",
        "temp.th",
        "tom_paine_sl_culvert.th",
        "west_false_river_barrier_leakage.th",
    ]
]
th_files.append(os.path.join(repo_home, "channel_depletion/msource_dated.th"))
multi_file_to_elapsed(th_files, outdir, start, name_transform=custom_name_transform)
