
ulimit -s unlimited
RUNHYDRO=hydro #"D:\dsm2\DSM2v821\bin\hydro.exe"
RUNQUAL=qual #"D:\dsm2\DSM2v821\bin\qual.exe"
EXPERIMENT="latinhypercube"
NUMBER=$1

echo "${RUNHYDRO} hydro_${EXPERIMENT}_${NUMBER}.inp"
${RUNHYDRO} hydro_${EXPERIMENT}_${NUMBER}.inp
echo "DONE"