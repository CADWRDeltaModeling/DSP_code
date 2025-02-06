from dms_datastore.read_ts import read_noaa
from vtools.functions.unit_conversions import m_to_ft
import os
import re

os.chdir(os.path.dirname(os.path.abspath(__file__)))

infile = "./noaa_sffpx_9414290_water_level_2017_2020.csv"
outfile = "./noaa_sffpx_9414290_water_level_2017_2020_ft.csv"

with open(infile, "r") as inf, open(outfile, "w") as outf:
    for line in inf:
        # print(line)
        if re.match(r"^\d{4}-\d{2}-\d{2}", line):  # Check if line starts with a date
            parts = line.strip().split(",")  # Split by comma
            parts[1] = f"{m_to_ft(float(parts[1])):.2f}"  # Convert 2nd item to feet
            parts = parts[:5]  # Keep only the first 5 columns
            lineout = ",".join(parts) + "\n"  # Reconstruct the line
        else:
            lineout = line
        outf.write(lineout)  # Write to output file

# in_df = read_noaa(infile)

# out_df = m_to_ft(in_df)
# out_df.to_csv(outfile, float_format="%.2f", index=True, index_label='datetime')