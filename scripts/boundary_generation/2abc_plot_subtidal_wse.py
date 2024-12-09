
import os
import pandas as pd
from vtools.functions.filter import cosine_lanczos
import matplotlib.pyplot as plt

from bokeh.io import show
from bokeh.layouts import column
from bokeh.models import RangeTool
from bokeh.plotting import figure, show, output_file, save
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis, Colorblind, Bokeh

def read_df(filename):
    wse_ts = pd.read_csv(filename, index_col='datetime')
    wse_ts.index = pd.to_datetime(wse_ts.index, format="%Y-%m-%d %H:%M:%S")

    return wse_ts

def calc_filtered(filename, unit='ft'):
    
    wse_ts = read_df(filename)

    if isinstance(wse_ts.index, pd.core.indexes.datetimes.DatetimeIndex):
        print('timeseries is inst-val, converting to per-aver')
        wse_ts.index = wse_ts.index.to_period()

    # convert to feet
    if unit=="m":
        wse_ts = wse_ts * 3.28084  # ft/m 

    df_filt = cosine_lanczos(wse_ts.copy(), cutoff_period='40H', padtype='odd')
    df_filt.columns = ['filtered']

    return df_filt

def calc_tidal_energy(filename, unit='ft'):

    wse_ts = read_df(filename)

    if isinstance(wse_ts.index, pd.core.indexes.datetimes.DatetimeIndex):
        print('timeseries is inst-val, converting to per-aver')
        wse_ts.index = wse_ts.index.to_period()

    # convert to feet
    if unit=="m":
        wse_ts = wse_ts * 3.28084  # ft/m

    df_nrg = cosine_lanczos((wse_ts-cosine_lanczos(wse_ts.copy(),
                                                    cutoff_period='40H', padtype='odd'))**2,
                            cutoff_period='40H', padtype='odd')  # = < (z- <z>)^2 >
    if not isinstance(df_nrg.index, pd.DatetimeIndex):
        df_nrg.index = df_nrg.index.to_timestamp()
    df_tidal_energy = df_nrg.resample('D', closed='right').mean()
    df_tidal_energy.columns = ['tidal_energy']
    # df_tidal_energy.index = df_tidal_energy.index.to_timestamp()
    df_tidal_energy = df_tidal_energy.resample('15min').ffill()

    return df_tidal_energy

def plot_energies(nrg_df_dict, filt_df_dict, inst_df_dict, save_dir=None, date_range=None):
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    # default colors same as sns.color_palette("tab10", n_colors=10).as_hex() # '#d62728' is red

    keys = list(nrg_df_dict.keys())
    
    if date_range is not None:
        for key, df in nrg_df_dict.items():
            nrg_df_dict[key] = df.loc[date_range[0]:date_range[1]]
        for key, df in filt_df_dict.items():
            filt_df_dict[key] = df.loc[date_range[0]:date_range[1]]
        for key, df in inst_df_dict.items():
            inst_df_dict[key] = df.loc[date_range[0]:date_range[1]]

    # get the first DateTime in the dataframe index
    first = nrg_df_dict[keys[0]].index[0]
    # get the last DateTime in the dataframe index
    last = nrg_df_dict[keys[0]].index[-1]

    nrg = figure(height=400, width=1200, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
            x_axis_type="datetime", x_axis_location="above",
            background_fill_color="#efefef", x_range=(first, last), title='MTZ Tidal Energy')

    for key, df in nrg_df_dict.items():
        nrg.line(df.index, df.iloc[:,0], color=colors[keys.index(key)], legend_label=f'{key} Data')
        
    select = figure(title="Drag the middle and edges of the selection box to change the range",
                    height=130, width=1200, y_range=nrg.y_range,
                    x_axis_type="datetime", y_axis_type=None,
                    tools="", toolbar_location=None, background_fill_color="#efefef")

    range_tool = RangeTool(x_range=nrg.x_range)
    range_tool.overlay.fill_color = "navy"
    range_tool.overlay.fill_alpha = 0.2

    for key, df in filt_df_dict.items():
        select.line(df.index, df.iloc[:,0], color=colors[keys.index(key)], legend_label=f'{key} Data')

    select.ygrid.grid_line_color = None
    select.add_tools(range_tool)

    filt = figure(height=400, width=1200, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
            x_axis_type="datetime", x_axis_location="above",
            background_fill_color="#efefef", x_range=nrg.x_range, title='MTZ Filtered Stage')

    for key, df in filt_df_dict.items():
        filt.line(df.index, df.iloc[:,0], color=colors[keys.index(key)], legend_label=f'{key} Data')

    inst = figure(height=400, width=1200, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
            x_axis_type="datetime", x_axis_location="above",
            background_fill_color="#efefef", x_range=nrg.x_range, title='MTZ Inst Stage')

    for key, df in inst_df_dict.items():
        inst.line(df.index, df.iloc[:,0], color=colors[keys.index(key)], legend_label=f'{key} Data')

    plot_final = column(nrg,select,filt,inst)
    if save_dir is not None:
        output_file(f'{save_dir}/mtz_subtidal_filter_v1.html', mode='inline')
        save(plot_final)
    else:
        show(plot_final)

os.chdir(os.path.dirname(os.path.abspath(__file__)))

original = "./input/dsm2_mtz_stage_ft_15min_clean.csv"
perturbed = "./data_out/perturb_historical_subtide_v1.csv"

# use this for re-creating subtide in R script
# calc_filtered(original).to_csv("./data_out/dsm2_mtz_stage_ft_15min_clean_filtered.csv", index=True, float_format="%.2f")


# calc and plot tidal energy and filter
nrg_df = {'Original':calc_tidal_energy(original), 'Perturbed':calc_tidal_energy(perturbed)}
filt_df = {'Original':calc_filtered(original), 'Perturbed':calc_filtered(perturbed)}
inst_df = {'Original':read_df(original), 'Perturbed':read_df(perturbed)}

plot_energies(nrg_df, filt_df, inst_df, save_dir='./plots',
               date_range=pd.to_datetime(['2000-01-01 00:00','2010-01-01 00:00']))