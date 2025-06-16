# Go to folder and create simulation folder
cd /scratch/tomkovic/DSP_code/model/schism/roundtrip;
mkdir -p rt_v1;
cd rt_v1;

# Copy geometry files in
GEOM_DIR=/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3/baseline_inputs
cp $GEOM_DIR/hgrid.gr3 \
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
  $GEOM_DIR/tvd.prop \
  .

