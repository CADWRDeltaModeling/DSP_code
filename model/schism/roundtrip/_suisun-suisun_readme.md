# Go to folder and create simulation folder
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip;
mkdir -p suisun-suisun;
cd suisun-suisun;
</pre>

# Copy geometry files in
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun;
GEOM_DIR=/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/suisun_inputs
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
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun;
BASE_DIR=../suisun-base
cp $BASE_DIR/bctides.in.2d \
  $BASE_DIR/bctides.in.3d \
  $BASE_DIR/fluxflag.prop \
  $BASE_DIR/station.in \
  $BASE_DIR/vgrid.in.2d \
  $BASE_DIR/interpolate_variables.in \
  $BASE_DIR/elev2D.th.nc \
  $BASE_DIR/check_progress.sh \
  $BASE_DIR/param.nml.tropic \
  $BASE_DIR/param.nml.clinic \
  $BASE_DIR/slurm.sh \
  $BASE_DIR/schism.sh \
  .
HOTNUDGE_DIR=/scratch/tomkovic/DSP_code/model/schism/roundtrip/hot_nudge
cp $HOTNUDGE_DIR/TEM_*cencoos_2015_suisun* \
  $HOTNUDGE_DIR/SAL_*cencoos_2015_suisun* \
  $HOTNUDGE_DIR/2015_suisun/hotstart.20150218.nc \
  $HOTNUDGE_DIR/2015_suisun/elev.ic \
  .
</pre>

## Sflux

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun;
BASE_DIR=../suisun-base
mkdir -p sflux;
rsync -avz $BASE_DIR/sflux/make_links.py ./sflux/make_links.py;
rsync -avz $BASE_DIR/sflux/*.txt ./sflux/
cd sflux;
python make_links.py;
cd ../;
</pre>

## Create boundary files from CalSim
<pre>
cd /scratch/tomkovic/DSP_code/scripts/boundary_generation/input;
cp ./port_calsim_schism.suisun-base.yaml ./port_calsim_schism.suisun-suisun.yaml;
</pre>

Change the version to suisun-suisun

<pre>
cd /scratch/tomkovic/DSP_code/scripts/boundary_generation/input;
bds port_bc port_calsim_schism.suisun-suisun.yaml -- --sd 2015/2/18 --ed 2016/5/15
</pre>

produces *.suisun-suisun.dated.th files

## Run CCFB gate op estimator
run /scratch/tomkovic/DSP_code/model/schism/roundtrip/create_ccfb.py to create ccfb_gate_syn.suisun-suisun.dated.th

## Create elapsed .th files

<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/;
cp ./multi_clip.suisun-base.py ./multi_clip.suisun-suisun.py
</pre>

Modify multi_clip.suisun-suisun.py by changing the version and run

## Link boundary files
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun;
ln -sf flux.suisun-suisun.dated_elapsed_20150218.th flux.th;
ln -sf delta_cross_channel.suisun-suisun.dated_elapsed_20150218.th delta_cross_channel.th;
ln -sf vsource.suisun-suisun.dated_elapsed_20150218.th vsource.th;
ln -sf vsink.suisun-suisun.dated_elapsed_20150218.th vsink.th;
ln -sf montezuma_radial.suisun-suisun.dated_elapsed_20150218.th montezuma_radial.th;
ln -sf montezuma_boat_lock.suisun-suisun.dated_elapsed_20150218.th montezuma_boat_lock.th;
ln -sf montezuma_flash.suisun-suisun.dated_elapsed_20150218.th montezuma_flash.th;
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
ln -sf ccfb_gate_syn.suisun-suisun.dated_elapsed_20150218.th ccfb_gate.th
</pre>

# Setup Barotropic Run
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun;
mkdir -p outputs;
ln -sf param.nml.tropic param.nml;
ln -sf bctides.in.2d bctides.in;
ln -sf vgrid.in.2d vgrid.in;
bds set_nudge cencoos_2015 --param param.nml.clinic;
ln -sf hotstart.20150218.nc hotstart.nc;
</pre>


<!-- ### copying to hpc4
rsync -avz tomkovic@10.3.80.51:/scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun/* /scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun/  --dry-run -->

# Setup Baroclinic Run
<pre>
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun;
<!-- cd /scratch/dms/tomkovic/DSP_code/model/schism/roundtrip/suisun-suisun; -->

<!-- # some rigamarole so that I can create uv3d from a separate run -->
rsync -avz ../interpolate_variables.in ./;
# link necessary files
ln -sf ../../suisun-base/hgrid.gr3 bg.gr3;
ln -sf ../hgrid.gr3 fg.gr3;
ln -sf ../../suisun-base/vgrid.in.2d vgrid.bg;
ln -sf ../vgrid.in.3d vgrid.fg;
ulimit -s unlimited;
interpolate_variables8; # this takes quite a while
cp uv3D.th.nc ../uv3D.th.nc;
cd ../;
ln -sf param.nml.clinic param.nml;
ln -sf bctides.in.3d bctides.in;
ln -sf vgrid.in.3d vgrid.in;
bds set_nudge cencoos_2015_suisun;
ln -sf hotstart.20150218.nc hotstart.nc;
</pre>