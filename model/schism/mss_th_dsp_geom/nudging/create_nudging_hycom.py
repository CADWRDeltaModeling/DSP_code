#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from schimpy import nudging
from schimpy.schism_mesh import write_mesh, read_mesh
import xarray as xr
import datetime

# Beloe are needed for plotting
from schimpy.geo_tools import ll2utm
from shapely.ops import cascaded_union
from shapely.geometry import Polygon

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

yaml_fn = 'nudge_hycom.yaml'

# These three lines are all you need to create the nudging file
nudging = nudging.Nudging(yaml_fn,crs ='EPSG:26910')
nudging.read_yaml()
nudging.create_nudging()

# The rest is for visualization, highly recommended for QA/QC

