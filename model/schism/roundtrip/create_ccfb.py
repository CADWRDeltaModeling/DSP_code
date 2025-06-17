from bdschism.ccf_gate_height import ccf_gate, sffpx_level
import tempfile
import string
import os
import csv

os.chdir(os.path.dirname(os.path.abspath(__file__)))

years = [2018]

start_dates = {2018: "2018-03-15"}
end_dates = {2018: "2019-04-01"}
out_dir = {2018: "rt_v1"}

for year in years:
    sdate = start_dates[year]
    edate = end_dates[year]
    dest = f"./{out_dir[year]}/ccfb_gate_syn.{out_dir[year]}.dated.th"
    astro_file = "/scratch/tomkovic/DSP_code/data/oh4/oh4_15min_predicted_10y_14_25.out"
    export_file = f"./{out_dir[year]}/flux.{out_dir[year]}.dated.th"
    sf_data_repo = "/scratch/tomkovic/DSP_code/data/noaa_sf/"

    sffpx_elev = sffpx_level(sdate, edate, sf_data_repo)

    ccf_gate(sdate, edate, dest, astro_file, export_file, sffpx_elev, False)
