#!/bin/bash
module purge
module load intel/2024.0 openmpi/5.0.2 hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1 schism/5.11.0
 
ulimit -s unlimited
NSCRIBES=10
SCHISM_BIN=pschism_PREC_EVAP_GOTM_TVD-VL
$SCHISM_BIN $NSCRIBES