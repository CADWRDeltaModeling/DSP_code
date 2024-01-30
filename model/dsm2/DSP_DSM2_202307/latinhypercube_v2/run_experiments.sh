#!/bin/sh
# purge modules and load intel
module purge
source /opt/intel/parallel_studio_xe_2019_update1/parallel_studio_xe_2019.1.053/bin/psxevars.sh ia32

ulimit -s unlimited
RUNHYDRO="/home/tomkovic/dsm2_repo/dsm2/BUILD/hydro"
RUNQUAL="/home/tomkovic/dsm2_repo/dsm2/BUILD/qual"
EXPERIMENT="latinhypercube"
NUMEXPS=7

for NUMBER in $(seq 1 $NUMEXPS); # all experiments
#for NUMBER in {5,7}; # just the selected experiments that didn't run the first time
do
    echo "${RUNHYDRO} hydro_${EXPERIMENT}_${NUMBER}.inp"
	${RUNHYDRO} hydro_${EXPERIMENT}_${NUMBER}.inp
    echo "${RUNQUAL} qual_ec_${EXPERIMENT}_${NUMBER}.inp"
	${RUNQUAL} qual_ec_${EXPERIMENT}_${NUMBER}.inp
done
