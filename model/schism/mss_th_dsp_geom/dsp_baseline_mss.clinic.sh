# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
#gen_elev2d --outfile baseline.mss.elev2D.th.nc --hgrid=hgrid.gr3 --stime 2021-1-01 --etime 2024-1-6 --slr 0.0 noaa_pryc1_9415020_elev_2020_2025.csv noaa_mtyc1_9413450_elev_2020_2025.csv
ln -sf baseline.mss.elev2D.th.nc elev2D.th.nc

# run script to create uv3d.th.nc
ulimit -s unlimited
							   
						   

# START !no-interp commands
[ "$1" != "no-interp" ] && 
# CHANGE OUTPUTS DIR -------------------------------------
mkdir -p outputs_tropic && 
mv outputs/* outputs_tropic && 
mkdir -p outputs &&

# CREATE OCEAN BOUNDARY ----------------------------------

cd outputs_tropic &&
rsync -avz ../interpolate_variables.in . &&

# link necessary files
ln -sf ../hgrid.gr3 bg.gr3 &&
ln -sf ../hgrid.gr3 fg.gr3 &&
ln -sf ../vgrid.in.2d vgrid.bg &&
ln -sf ../vgrid.in.3d vgrid.fg &&

interpolate_variables8 && # this takes quite a while
cp uv3D.th.nc ../uv3D.th.nc &&
cd ../
# END !no-interp commands

# Re link barotropic files -------------------------------

# CREATE CLIINIC SYMBOLIC LINKS ----------------------------

# add new links
ln -sf SAL_nu_obshycom.nc SAL_nu.nc
ln -sf TEM_nu_obshycom.nc TEM_nu.nc
ln -sf salinity_nudge_obshycom.gr3 SAL_nudge.gr3
ln -sf temperature_nudge_obshycom.gr3 TEM_nudge.gr3
ln -sf hotstart_2021.nc hotstart.nc

# change existing links
ln -sf bctides.in.3d bctides.in
ln -sf vgrid.in.3d vgrid.in
ln -sf param.nml.clinic param.nml
# CHECK OUTPUT DIRECTORY IS PRESENT ----------------------

mkdir -p outputs # only creates outputs if it's not preset

# model runs with yml file for Azure setup