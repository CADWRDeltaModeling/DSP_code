# adapted from Sophie Munger by Lily Tomkovic

import shutil, os, sys, io
import pyhecdss
import pandas as pd
from vtools.functions.unit_conversions import FT2M

import contextlib
from schimpy.model_time import file_to_elapsed
import datetime as dt

import matplotlib.pyplot as plt

from bokeh.io import show
from bokeh.layouts import column
from bokeh.models import RangeTool
from bokeh.plotting import figure, show, output_file, save
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis, Colorblind, Bokeh

#%% Functions
def create_scenario(fname, template_dir, study_dir):
    
    case_dir = os.path.join(study_dir, fname.split('.')[0])
    
    if not os.path.exists(case_dir):
        shutil.copytree(template_dir, case_dir)
    else:
        pass
    return case_dir

def postprocess(case_dir, case_name, ccfb_ref, out_file, start_time, ndup=5):
    
    if not os.path.exists('../output_th/'):
        os.mkdir('../output_th/')
    out_dss = os.path.join(case_dir, f'output/{case_name}_CCF.dss')

    gate = list(pyhecdss.get_ts(out_dss,'//CCF_GATE/OP-FROM-NODE//15MIN//'))[0].data
    gate.rename(columns={gate.columns[0]: "op" }, inplace = True)
    gate['height_ft'] = list(pyhecdss.get_ts(out_dss,'//CCF_GATE/HEIGHT//15MIN//'))[0].data

    # Read gate information in standard ccfb_gate format
    df_ccfb_ref = pd.read_csv(ccfb_ref, index_col=0, header=0, sep="\s+", dtype=str)

    # specify first time slice as one before start time to be sure that we get the first 0ish hours
    dt_ind = (df_ccfb_ref.index>start_time).argmax()-1
    df_ccfb_case = df_ccfb_ref[dt_ind:].copy()
    df_ccfb_case.index = pd.to_datetime(df_ccfb_case.index)
    df_ccfb_case = df_ccfb_case.reindex(gate.index, method="ffill")
    df_ccfb_case = df_ccfb_case[gate.index[0]:]
    df_ccfb_case['height'] = gate['op'] * gate['height_ft'] * FT2M # sets height to 0 if operation is closed (=0)
    df_ccfb_case = df_ccfb_case.fillna(method="ffill")
    df_ccfb_case['ndup'] = ndup
    df_ccfb_case.index.name = "datetime"
    df_ccfb_case.index = df_ccfb_case.index.strftime('%Y-%m-%dT%H:%M')
    df_ccfb_case.to_csv(out_file, sep=" ", float_format="%.3f")

    print("postprocessing done.")
    return None

def clip_file_to_start(infile, outpath, start):
    with contextlib.redirect_stdout(io.StringIO()):
        file_to_elapsed(infile, outpath=outpath, start=start)
        sys.stdout.close()

def clip_th(dated_fn, year, month, day):
    # clip to model period
    
    clip_fn = dated_fn.replace('.th' ,'_clip.th')
    clip_file_to_start(dated_fn, 
                       outpath=clip_fn,
                       start=dt.datetime(year=int(year),
                                        month=int(month),
                                        day=int(day)))


def plot_compare(exg_th, mod_th, out_dss, case_name):
    # Read gate information in standard ccfb_gate format
    df_ccfb_exg = pd.read_csv(exg_th, index_col=0, header=0, sep="\s+", dtype=str)
    df_ccfb_mod = pd.read_csv(mod_th, index_col=0, header=0, sep="\s+", dtype=str)
    ccf_stage = list(pyhecdss.get_ts(out_dss,'//CCF_RES/STAGE//15MIN//'))[0].data

    gate = pd.merge(df_ccfb_exg, df_ccfb_mod, left_index=True, right_index=True, how='inner')
    gate.index = pd.to_datetime(gate.index)
    gate = gate.astype(float)

    # get the first DateTime in the dataframe index
    first = gate.index[0]
    # get the last DateTime in the dataframe index
    last = gate.index[-1]

    p = figure(height=250, width=900, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
               x_axis_type="datetime", x_axis_location="above",
               background_fill_color="#efefef", x_range=(first, last), title='CCFB compared historical to modified')

    p.line(gate.index, gate.height_x, color='black', legend_label='Historical')
    p.line(gate.index, gate.height_y, color='blue', legend_label='Modified')

    r = figure(height=250, width=900, tools = "xpan,wheel_zoom,box_zoom,reset,save,hover",
               x_axis_type="datetime", x_axis_location="above", x_range=p.x_range,
               background_fill_color="#efefef")

    r.line(ccf_stage.index, ccf_stage.iloc[:,0], color='blue', legend_label='Modified')

    select = figure(title="Drag the middle and edges of the selection box to change the range",
                    height=130, width=900, y_range=p.y_range,
                    x_axis_type="datetime", y_axis_type=None,
                    tools="", toolbar_location=None, background_fill_color="#efefef")

    range_tool = RangeTool(x_range=p.x_range)
    range_tool.overlay.fill_color = "navy"
    range_tool.overlay.fill_alpha = 0.2

    select.line(ccf_stage.index, ccf_stage.iloc[:,0], color='blue', legend_label='Modified')

    select.ygrid.grid_line_color = None
    select.add_tools(range_tool)

    plot_final = column(p,r, select)
    output_file(f'ccfb_{case_name}.html', mode='inline')
    save(plot_final)

#%% Run postprocessor

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # case_name = "lhc_1"
    # start_time = "2013-10-24"
    # case_name = "lhc_2"
    # start_time = "2011-10-18"
    # case_name = "lhc_3"
    # start_time = "2013-10-24"
    # case_name = "lhc_5"
    # start_time = "2006-11-14"
    # case_name = "lhc_6"
    # start_time = "2011-10-18"
    case_name = "lhc_7"
    start_time = "2009-10-29"
    case_dir = "../"

    # Reference ccfb th, whose height column will be replaced.
    ccfb_ref = './data/ccfb_gate.th'
    out_file = os.path.join('../output_th/', f'ccfb_{case_name}.th')

    #POSTPROCESS 
    postprocess(case_dir, case_name, ccfb_ref, out_file, start_time, ndup=5)

    # CHECK FILE BEFORE CONVERTING
    clip_th(out_file, start_time.split('-')[0], start_time.split('-')[1], start_time.split('-')[2])

    out_dss = os.path.join(case_dir, "output", case_name+'_CCF.dss')

    plot_compare(ccfb_ref, out_file, out_dss, case_name)