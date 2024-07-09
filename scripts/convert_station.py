from schimpy.station import convert_stations

input = "/scratch/tomkovic/gis/dsp_added_sta.shp"
output =  "/scratch/tomkovic/gis/dsp_added_sta.in"

convert_stations(input, output)

# input =  "/scratch/tomkovic/gis/dsp_added_sta.in"
# output =  "/scratch/tomkovic/gis/test.shp"

# convert_stations(input, output)