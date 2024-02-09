#!/bin/sh
# Purge any modules that are already loaded.
module purge
# Load intel compilers
module use /opt/intel/oneapi/modulefiles
module load compiler/2022.1.0
module load mpi/2021.5.1

# Load dependencies: HDF5, NetCDF
module load hdf5/1.12.2_intel2022.1
module load netcdf-c/4.9.0_intel2022.1
module load netcdf-fortran/4.6.0_intel2022.1
ulimit -s unlimited
/home/tomkovic/schism_build/build_20220927/bin/pschism_PREC_EVAP_GOTM_TVD-VL 10