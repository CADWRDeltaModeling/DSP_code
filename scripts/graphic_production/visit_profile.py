import geopandas as gpd

pt_path = "D:/projects/delta_salinity/data/GIS/post-processing/visit_longitudinal_pts.shp"

gdf = gpd.read_file(pt_path)
gdf = gdf.sort_values(by='Id')

xy_list = [(point.x, point.y) for point in gdf.geometry]

print(xy_list)
