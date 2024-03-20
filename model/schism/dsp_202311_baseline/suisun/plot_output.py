import schimpy
import pandas as pd
import os

from schimpy.plot_default_formats import *
from schimpy.station import *
import matplotlib.pyplot as plt
import datetime

os.chdir(os.path.dirname(__file__))
station_fpath = 'station.in'
time_basis = datetime.datetime(2006, 11, 14)
set_color_cycle_dark2()

outputs_fpath = 'outputs_tropic/staout_1'

station_name = "bdl" # "Montezuma Slough near Beldons Landing"
sub_loc = "default"

all_ts = read_staout(outputs_fpath, station_fpath, time_basis)

all_ts[station_name + "_" + sub_loc].plot()
plt.show()