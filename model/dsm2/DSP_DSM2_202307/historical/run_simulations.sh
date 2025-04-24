#!/bin/sh
# purge modules and load intel
#module purge
#source /opt/intel/parallel_studio_xe_2019_update1/parallel_studio_xe_2019.1.053/bin/psxevars.sh ia32

ulimit -s unlimited
RUNHYDRO=D:/dsm2/DSM2v822/bin/hydro
RUNQUAL=D:/dsm2/DSM2v822/bin/qual

echo "${RUNHYDRO} hydro.inp"
${RUNHYDRO} hydro.inp
echo "${RUNQUAL} qual_ec.inp"
${RUNQUAL} qual_ec.inp
echo "${RUNQUAL} qual_ec_x2.inp"
${RUNQUAL} qual_ec_x2.inp
done
