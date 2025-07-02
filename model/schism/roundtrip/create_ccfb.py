from bdschism.ccf_gate_height import ccf_gate, sffpx_level
import tempfile
import string
import os
import csv

os.chdir(os.path.dirname(os.path.abspath(__file__)))
scenarios = [
    "slr-slr",
]  # ["rt_v1", "suisun-base", "suisun-suisun", "slr-base", "slr-slr"]
years = {
    "rt_v1": 2018,
    "suisun-base": 2015,
    "suisun-suisun": 2015,
    "slr-base": 2015,
    "slr-slr": 2015,
}

start_dates = {2018: "2018-03-15", 2015: "2015-02-18"}
end_dates = {2018: "2019-04-01", 2015: "2016-05-15"}

for scene in scenarios:
    year = years[scene]
    sdate = start_dates[year]
    edate = end_dates[year]
    dest = f"./{scene}/ccfb_gate_syn.{scene}.dated.th"
    astro_file = "/scratch/tomkovic/DSP_code/data/oh4/oh4_15min_predicted_10y_14_25.out"
    export_file = f"./{scene}/flux.{scene}.dated.th"
    sf_data_repo = "/scratch/tomkovic/DSP_code/data/noaa_sf/"

    sffpx_elev = sffpx_level(sdate, edate, sf_data_repo)

    ccf_gate(sdate, edate, dest, astro_file, export_file, sffpx_elev, False)
