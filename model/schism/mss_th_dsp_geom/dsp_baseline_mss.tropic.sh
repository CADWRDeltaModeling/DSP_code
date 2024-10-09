# CREATE OCEAN BOUNDARY ----------------------------------

## this has already been run in earliler iterations but is kept here to show where the data comes from
# download_noaa --syear 2006 --eyear 2016 --param water_level noaa_stations.txt 

# generate .th.nc file
gen_elev2d --outfile baseline.mss.elev2D.th.nc --hgrid=hgrid.gr3 --stime 2021-1-01 --etime 2024-1-6 --slr 0.0 noaa_pryc1_9415020_elev_2020_2025.csv noaa_mtyc1_9413450_elev_2020_2025.csv
ln -sf baseline.mss.elev2D.th.nc elev2D.th.nc

# CREATE OTHER SYMBOLIC LINKS ----------------------------
# shared spatial inputs
ln -sf bctides.in.2d bctides.in
ln -sf vgrid.in.2d vgrid.in

# add gr3 links
ln -sf salinity_nudge_roms.gr3 SAL_nudge.gr3
ln -sf temperature_nudge_roms.gr3 TEM_nudge.gr3

# modified TH inputs
# ln -sf ccfb_lhc_4.th ccfb_gate.th

# inputs specific to this setup
ln -sf param.nml.tropic param.nml
# ln -sf ../launch.tropic.pbs launch.pbs # only necessary with HPC4

# CHECK OUTPUT DIRECTORY IS PRESENT ----------------------

mkdir -p outputs # only creates outputs if it's not preset

# model runs with yml file for Azure setup