# script to push the tidal boundary forward consistent with the model run

import pandas as pd
from dms_datastore.read_ts import read_noaa
import re, os

def write_noaa(infile, in_df, outfile, repl_col=1):
    with open(infile,'r') as in_f:
        with open(outfile,'w') as out_f:
            inlines = in_f.readlines()
            for inl in inlines:
                if re.match("[0-9]{4}-[0-9]{2}-[0-9]{2}",inl[:10]):
                    break
                else:
                    out_f.write(inl)
            for dateindex, row in in_df.iterrows():
                out_f.write(f'{dateindex},{row["Water Level"]},999,0,0,0,0,v\n')
                



if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    tidal_infile = 'noaa_download/noaa_sffpx_9414290_water_level_2006_2017.csv'
    timedelt = 100
    tidal_outfile = 'noaa_sffpx_shift_forward_100d.csv'

    tidal_df = read_noaa(tidal_infile,force_regular=True)
    tidal_df.index = tidal_df.index - pd.Timedelta(days=timedelt)
    write_noaa(tidal_infile, tidal_df, tidal_outfile)