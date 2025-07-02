# Go to folder and create simulation folder
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip;
mkdir -p suisun-base;
cd suisun-base;
</pre>

# Copy geometry files in
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
GEOM_DIR=/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/baseline_inputs
cp $GEOM_DIR/hgrid.gr3 \
  $GEOM_DIR/vgrid.in.3d \
  $GEOM_DIR/albedo.gr3 \
  $GEOM_DIR/diffmax.gr3 \
  $GEOM_DIR/diffmin.gr3 \
  $GEOM_DIR/estuary.gr3 \
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
  $GEOM_DIR/tvd.prop \
  .
</pre>

Copy others:
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
OTHER_DIR=/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3
cp $OTHER_DIR/bctides.in.2d \
  $OTHER_DIR/bctides.in.3d \
  $OTHER_DIR/fluxflag.prop \
  $OTHER_DIR/station.in \
  $OTHER_DIR/vgrid.in.2d \
  $OTHER_DIR/interpolate_variables.in \
  .
</pre>

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
cp /scratch/projects/summer_x2_2025/preprocess/common/2016/check_progress.sh ./
</pre>

# Time-related

## Determine model run period
Start date: 2015-02-18
End date: 2016-05-01
rndays: 438

### Set number of days in interpolate_variables.in
Fix number of days in interpolate_variables.in

## Create ocean boundary files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip;
conda activate schism;
mkdir -p noaa_download;
cd noaa_download;
download_noaa --syear 2015 --eyear 2016 --param water_level --dest ./ noaa_stations.txt;
cd ../suisun-base;
gen_elev2d --outfile elev2D.th.nc --hgrid=hgrid.gr3 --stime=2015-2-18 --etime=2016-05-15 --slr 0.0 ../noaa_download/noaa_pryc1_9415020_water_level_2015_2017.csv ../noaa_download/noaa_mtyc1_9413450_water_level_2015_2017.csv;
</pre>

## Create param.nml files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
cp ../rt_v1/param.nml.tropic ./;
cp ../rt_v1/param.nml.clinic ./;
</pre>

change start and rndays on both

## Create slurm file
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
cp /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1/slurm.sh ./slurm.sh;
cp /scratch/tomkovic/DSP_code/model/schism/roundtrip/rt_v1/schism.sh ./;
</pre>

Change parameters like job name etc.

## Create Nudging files

### Gather data
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
hot_nudge_data --start_date 2015-2-18 --nudge_len 383 --repo_dir /scratch/nasbdo/modeling_data/repo/continuous/screened/ --dest_dir 2015
</pre>

### Copy nudging yaml
Nudging will only be at the ocean boundary. Hotstart will be at the start of the model period(?)

Copy nudging files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
cp ./nudge_cencoos_2018.yaml ./nudge_cencoos_2015.yaml;
</pre>

Edit start date, rndays, vgrid adn hgrid sources and data source to /scratch/nasbdo/modeling_data/coastal/roms/cencoos_ca_roms/processed/cencoos_hourly_pst

### Create nudging nc files and copy to simulation folder

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
create_nudging --input ./nudge_cencoos_2015.yaml 
cp ./TEM_*2015* ../suisun-base/
cp ./SAL_*2015* ../suisun-base/
</pre>

## Create Hotstart

From select_polaris_date.R 2-18-2015 has 33 stations - golden.
copy from \\nasbdo\Modeling_Data\usgs_cruise -> ./2015/usgs_cruise_2015.csv

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge;
cp ./create_hotstart_2018.py ./create_hotstart_2015.py
cp ./hotstart_2018.yaml ./hotstart_2015.yaml
</pre>

Edit necessary fields

Run create_hotstart_2015.py
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge
cp ./2015/hotstart.20150218.nc ../suisun-base/
cp ./2015/elev.ic ../suisun-base/
</pre>

## Sflux

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
mkdir -p sflux;
rsync -avz ../rt_v1/sflux/make_links.py ./sflux/make_links.py;
rsync -avz ../rt_v1/sflux/*.txt ./sflux/
</pre>

change source and run make_links.py in suisun-base/sflux

## Create boundary files from CalSim
<pre>
cd /scratch/tomkovic/DSP_code/scripts/boundary_generation/input;
bds port_bc port_calsim_schism.suisun-base.yaml -- --sd 2015/2/18 --ed 2016/5/15
</pre>

produces *.suisun-base.dated.th files

## Run CCFB gate op estimator
run /scratch/tomkovic/DSP_code/model/schism/roundtrip/create_ccfb.py to create ccfb_gate_syn.suisun-base.dated.th

## Create elapsed .th files

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/;
cp ./multi_clip_2018.py ./multi_clip.suisun-base.py
</pre>

Modify multi_clip.suisun-base.py and run

## Link boundary files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
ln -sf flux.suisun-base.dated_elapsed_20150218.th flux.th;
ln -sf delta_cross_channel.suisun-base.dated_elapsed_20150218.th delta_cross_channel.th;
ln -sf vsource.suisun-base.dated_elapsed_20150218.th vsource.th;
ln -sf vsink.suisun-base.dated_elapsed_20150218.th vsink.th;
ln -sf montezuma_radial.suisun-base.dated_elapsed_20150218.th montezuma_radial.th;
ln -sf montezuma_boat_lock.suisun-base.dated_elapsed_20150218.th montezuma_boat_lock.th;
ln -sf montezuma_flash.suisun-base.dated_elapsed_20150218.th montezuma_flash.th;
ln -sf grantline_culvert_elapsed_20150218.th grantline_culvert.th;
ln -sf grantline_barrier_elapsed_20150218.th grantline_barrier.th;
ln -sf grantline_weir_elapsed_20150218.th grantline_weir.th;
ln -sf midr_culvert_l_elapsed_20150218.th midr_culvert_l.th;
ln -sf midr_culvert_r_elapsed_20150218.th midr_culvert_r.th;
ln -sf midr_weir_elapsed_20150218.th midr_weir.th;
ln -sf oldr_head_barrier_elapsed_20150218.th oldr_head_barrier.th;
ln -sf oldr_tracy_culvert_elapsed_20150218.th oldr_tracy_culvert.th;
ln -sf oldr_tracy_barrier_elapsed_20150218.th oldr_tracy_barrier.th;
ln -sf oldr_tracy_weir_elapsed_20150218.th oldr_tracy_weir.th;
ln -sf salt_elapsed_20150218.th SAL_1.th;
ln -sf temp_elapsed_20150218.th TEM_1.th;
ln -sf tom_paine_sl_culvert_elapsed_20150218.th tom_paine_sl_culvert.th;
ln -sf west_false_river_barrier_leakage_elapsed_20150218.th west_false_river_barrier_leakage.th;
ln -sf msource_elapsed_20150218.th msource.th
ln -sf ccfb_gate_syn.suisun-base.dated_elapsed_20150218.th ccfb_gate.th
</pre>

# Setup Barotropic Run
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
mkdir -p outputs;
ln -sf param.nml.tropic param.nml;
ln -sf bctides.in.2d bctides.in;
ln -sf vgrid.in.2d vgrid.in;
bds set_nudge cencoos_2015 --param param.nml.clinic;
ln -sf hotstart.20150218.nc hotstart.nc;
</pre>


<!-- ### copying to hpc4
rsync -avz tomkovic@10.3.80.51:/scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base/* /scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/suisun-base/  --dry-run -->

# Setup Baroclinic Run
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-base;
<!-- cd /scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/suisun-base; -->
mkdir -p outputs_tropic;
mv outputs/* outputs_tropic; 
mkdir -p outputs;
cd outputs_tropic;
rsync -avz ../interpolate_variables.in ./;
# link necessary files
ln -sf ../hgrid.gr3 bg.gr3;
ln -sf ../hgrid.gr3 fg.gr3;
ln -sf ../vgrid.in.2d vgrid.bg;
ln -sf ../vgrid.in.3d vgrid.fg;
ulimit -s unlimited;
interpolate_variables8; # this takes quite a while
cp uv3D.th.nc ../uv3D.th.nc;
cd ../;
ln -sf param.nml.clinic param.nml;
ln -sf bctides.in.3d bctides.in;
ln -sf vgrid.in.3d vgrid.in;
bds set_nudge cencoos_2015;
ln -sf hotstart.20150218.nc hotstart.nc;
</pre>