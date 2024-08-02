
import shutil, os, sys, io
import pyhecdss
import pandas as pd
from vtools.functions.unit_conversions import FT2M

import contextlib
from schimpy.model_time import file_to_elapsed
import datetime as dt

#%% Functions
def create_scenario(fname, template_dir, study_dir):
    
    case_dir = os.path.join(study_dir, fname.split('.')[0])
    
    if not os.path.exists(case_dir):
        shutil.copytree(template_dir, case_dir)
    else:
        pass
    return case_dir

def postprocess(case_dir, case_name, ccfb_ref, out_file, ndup=5):
    
    if not os.path.exists('../output_th/'):
        os.mkdir('../output_th/')
    out_dss = os.path.join(case_dir, f'output/{case_name}_CCF.dss')

    gate = list(pyhecdss.get_ts(out_dss,'//CCF_GATE/OP-FROM-NODE//15MIN//'))[0].data
    gate.rename(columns={gate.columns[0]: "op" }, inplace = True)
    gate['height_ft'] = list(pyhecdss.get_ts(out_dss,'//CCF_GATE/HEIGHT//15MIN//'))[0].data

    # Write gate information in standard ccfb_gate format
    df_ccfb_ref = pd.read_csv(ccfb_ref, index_col=0, header=0, sep="\s+", dtype=str)

    df_ccfb_case = df_ccfb_ref[start_time:].copy()
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

#%% Run postprocessor

if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # case_name = "lhc_5"
    # start_time = "2006-11-14"
    case_name = "lhc_3"
    start_time = "2013-10-24"
    case_dir = "D:/projects/delta_salinity/scripts/DSP_code/model/dsm2/DSP_DSM2_202307/CCF_new_gate_op_June2024-lily"

    # Reference ccfb th, whose height column will be replaced.
    ccfb_ref = './data/ccfb_gate.th'
    out_file = os.path.join('../output_th/', f'ccfb_{case_name}.th')

    #POSTPROCESS 
    postprocess(case_dir, case_name, ccfb_ref, out_file, ndup=5)

    # CHECK FILE BEFORE CONVERTING
    clip_th(out_file, start_time.split('-')[0], start_time.split('-')[1], start_time.split('-')[2])