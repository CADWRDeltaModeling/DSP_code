# Go to folder and create simulation folder
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip;
mkdir -p rt_v1;
cd rt_v1;
</pre>

# Copy geometry files in
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
GEOM_DIR=/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/baseline_inputs
cp $GEOM_DIR/hgrid.gr3 \
  $GEOM_DIR/vgrid.in.2d \
  $GEOM_DIR/vgrid.in.3d \
  $GEOM_DIR/albedo.gr3 \
  $GEOM_DIR/diffmax.gr3 \
  $GEOM_DIR/diffmin.gr3 \
  $GEOM_DIR/estuary.gr3 \
  $GEOM_DIR/hgrid.gr3 \
  $GEOM_DIR/hgrid.ll \
  $GEOM_DIR/krvel.gr3 \
  $GEOM_DIR/manning.gr3 \
  $GEOM_DIR/nlayer.gr3 \
  $GEOM_DIR/nlayer_default.gr3 \
  $GEOM_DIR/rough.gr3 \
  $GEOM_DIR/sav_cd.gr3 \
  $GEOM_DIR/sav_D.gr3 \
  $GEOM_DIR/sav_h.gr3 \
  $GEOM_DIR/sav_N.gr3 \
  $GEOM_DIR/watertype.gr3 \
  $GEOM_DIR/windfactor.gr3 \
  $GEOM_DIR/windrot_geo2proj.gr3 \
  $GEOM_DIR/xlsc.gr3 \
  $GEOM_DIR/hydraulics.in \
  $GEOM_DIR/source_sink.in \
  $GEOM_DIR/fluxflag.prop \
  $GEOM_DIR/split_quad.prop \
  .
</pre>

Copy others:
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
OTHER_DIR=/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3
cp $OTHER_DIR/bctides.in.2d \
  $OTHER_DIR/bctides.in.3d \
  $OTHER_DIR/fluxflag.prop \
  $OTHER_DIR/station.in \
  $OTHER_DIR/vgrid.in.2d \
  .
</pre>

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
cp /scratch/projects/summer_x2_2025/preprocess/common/2016/check_progress.sh ./
</pre>

# Time-related

## Determine model run period
Start date: 2018-03-15
End date: 2019-04-01
rndays: 382

## Create ocean boundary files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip;
conda activate schism;
mkdir -p noaa_download;
cd noaa_download;
download_noaa --syear 2018 --eyear 2019 --param water_level --dest ./ noaa_stations.txt;
cd ../rt_v1;
gen_elev2d --outfile elev2D.th.nc --hgrid=hgrid.gr3 --stime=2018-3-15 --etime=2019-04-15 --slr 0.0 ../noaa_download/noaa_pryc1_9415020_water_level_2018_2020.csv ../noaa_download/noaa_mtyc1_9413450_water_level_2018_2020.csv;
</pre>

## Create param.nml files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
cp /home/tomkovic/BayDeltaSCHISM/templates/bay_delta/param.nml.tropic ./;
cp /home/tomkovic/BayDeltaSCHISM/templates/bay_delta/param.nml.clinic ./;
</pre>

change start and rndays on both

## Create slurm file
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/slurm_MESHNAME_CASENAME.tropic.sh ./slurm.sh;
cp /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/schism.sh ./;
</pre>

Change parameters like job name etc.

## Create Nudging files

### Gather data
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
hot_nudge_data --start_date 2018-3-15 --nudge_len 383 --repo_dir /scratch/nasbdo/modeling_data/repo/continuous/screened/ --dest_dir 2018
</pre>

### Copy nudging yaml
Nudging will only be at the ocean boundary. Hotstart will be at the start of the model period(?)

Copy nudging files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
cp /home/tomkovic/BayDeltaSCHISM/examples/nudging/nudge_cencoos.yaml ./nudge_cencoos_2018.yaml;
</pre>

Edit start date, rndays, vgrid adn hgrid sources and data source to /scratch/nasbdo/modeling_data/coastal/roms/cencoos_ca_roms/processed/cencoos_hourly_pst

### Create nudging nc files and copy to simulation folder

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
create_nudging --input ./nudge_cencoos_2018.yaml 
cp ./TEM_*2018* ../rt_v1/
cp ./SAL_*2018* ../rt_v1/
</pre>

## Create Hotstart

From select_polaris_date.R 3-15-2018 has 35 stations - golden.
copy from \\nasbdo\Modeling_Data\usgs_cruise -> ./2018/usgs_cruise_2018.csv

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
cp /home/tomkovic/BayDeltaSCHISM/examples/hotstart/examples/basic/create_hotstart.py ./create_hotstart_2018.py
cp /home/tomkovic/BayDeltaSCHISM/examples/hotstart/examples/basic/hotstart.yaml ./hotstart_2018.yaml
</pre>

Edit necessary fields

Run create_hotstart_2018.py
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge
cp ./hotstart.20180315.nc ../rt_v1/
cp ./elev.ic ../rt_v1/
</pre>

## Sflux

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
mkdir -p sflux;
rsync -avz /scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/simulations/baseline_lhc_1/sflux/make_links_az.py ./sflux/make_links.py;
rsync -avz /scratch/projects/itp202411/simulations/itp2024_hist/sflux/*.txt ./sflux/
</pre>

change source and run make_links.py in rt_v1/sflux

## Create boundary files from CalSim
<pre>
cd /scratch/tomkovic/DSP_code/scripts/boundary_generation/input;
bds port_bc port_calsim_schism.yaml -- --sd 2018/3/1 --ed 2019/4/1
</pre>

produces *.rt_v1.dated.th files

## Run CCFB gate op estimator
run create_ccfb.py to create ccfb_gate_syn.rt_v1.dated.th

## Create elapsed .th files

Modify multi_clip_2018.py

## Link boundary files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
ln -sf flux.rt_v1.dated_elapsed_20180315.th flux.th;
ln -sf delta_cross_channel.rt_v1.dated_elapsed_20180315.th delta_cross_channel.th;
ln -sf vsource.rt_v1.dated_elapsed_20180315.th vsource.th;
ln -sf vsink.rt_v1.dated_elapsed_20180315.th vsink.th;
ln -sf montezuma_radial.rt_v1.dated_elapsed_20180315.th montezuma_radial.th;
ln -sf montezuma_boat_lock.rt_v1.dated_elapsed_20180315.th montezuma_boat_lock.th;
ln -sf montezuma_flash.rt_v1.dated_elapsed_20180315.th montezuma_flash.th;
ln -sf grantline_culvert_elapsed_20180315.th grantline_culvert.th;
ln -sf grantline_barrier_elapsed_20180315.th grantline_barrier.th;
ln -sf grantline_weir_elapsed_20180315.th grantline_weir.th;
ln -sf midr_culvert_l_elapsed_20180315.th midr_culvert_l.th;
ln -sf midr_culvert_r_elapsed_20180315.th midr_culvert_r.th;
ln -sf midr_weir_elapsed_20180315.th midr_weir.th;
ln -sf oldr_head_barrier_elapsed_20180315.th oldr_head_barrier.th;
ln -sf oldr_tracy_culvert_elapsed_20180315.th oldr_tracy_culvert.th;
ln -sf oldr_tracy_barrier_elapsed_20180315.th oldr_tracy_barrier.th;
ln -sf oldr_tracy_weir_elapsed_20180315.th oldr_tracy_weir.th;
ln -sf salt_elapsed_20180315.th salt.th;
ln -sf temp_elapsed_20180315.th temp.th;
ln -sf tom_paine_sl_culvert_elapsed_20180315.th tom_paine_sl_culvert.th;
ln -sf west_false_river_barrier_leakage_elapsed_20180315.th west_false_river_barrier_leakage.th;
ln -sf msource_elapsed_20180315.th msource.th
ln -sf ccfb_gate_syn.rt_v1.dated_elapsed_20180315.th ccfb_gate.th
</pre>

# Setup run
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1;
ln -sf param.nml.tropic param.nml;
ln -sf bctides.in.2d bctides.in;
ln -sf vgrid.in.2d vgrid.in;
bds set_nudge cencoos_2018;
ln -sf hotstart.20180315.nc hotstart.nc;
</pre>


### copying to hpc4
rsync -avz tomkovic@10.3.80.51:/scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1/* /scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/rt_v1/  --dry-run
