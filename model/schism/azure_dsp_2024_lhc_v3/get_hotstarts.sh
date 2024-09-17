if [ "$#" -ne 2 ]; then
    echo "Error: This script requires exactly two inputs."
    echo "Usage: get_hotstarts.sh dir baro"
    exit 1
fi

dir=$1 # example 'simulations/baseline_lhc_2'
baro=$2 # either clinic or tropic

parent_dir="${dir%/*}"
sim_dir="${dr##*/}"

az_dir="/scratch/tomkovic/DSP_code/model/schism/azure_dsp_2024_lhc_v3"

mod_dir="${az_dir}/${dir}"

if [[ ! -d "${mod_dir}" ]]; then
    echo "Error: The model directory doesn't exist: ${mod_dir}"
    exit 1
fi

echo "Model directory: ${mod_dir}"
sleep 1.5

# download mirror.out
export study_dir="azure_dsp_2024_lhc_v3/${dir}"
azcopy copy "${AZLINK}${study_dir}?${sas}" "${az_dir}/${parent_dir}" --recursive --preserve-symlinks --include-path="outputs/mirror.out" --exclude-path="sflux"
echo "Downloaded mirror.out"
sleep 1.5

# copy restart_from_hotstart.py to directory
cp -f "${az_dir}/restart_from_hotstart.py" "${mod_dir}/restart_from_hotstart.py"
cd "${mod_dir}"

# get last hotstart file timestep from mirror.out
echo "python restart_from_hotstart.py --mod_dir "${mod_dir}" --baro "${baro}""
last_restart=$(python restart_from_hotstart.py --mod_dir "${mod_dir}" --baro "${baro}")
echo "Creating hotstart from timestep: ${last_restart}"
sleep 1

# download the hotstart files
azcopy copy "${AZLINK}${study_dir}?${sas}" "${az_dir}/${parent_dir}" --recursive --preserve-symlinks --include-pattern="hotstart_*_${last_restart}.nc;local_to_global*" --exclude-path="sflux"

# check that size of hotstart_000000_${last_restart}.nc is over 0.5 MB
size=$(stat -c%s "outputs/hotstart_000000_${last_restart}.nc")
if [ "${size}" -gt 500000 ]; then
    echo "hotstart.nc file is sufficient size"
else
    echo "ERROR: The first hotstart file size is too small to be able to create a combined hotstart"
    exit 1
fi

# combine hotstarts
echo "python restart_from_hotstart.py --mod_dir "${mod_dir}" --baro "${baro}" --combine"
python restart_from_hotstart.py --mod_dir "${mod_dir}" --baro "${baro}" --combine

# copy hotstart_it and param.nml back to Azure
export study_dir="azure_dsp_2024_lhc_v3/${parent_dir}"
azcopy copy "${mod_dir}" "${AZLINK}${study_dir}?${sas}"  --recursive --preserve-symlinks --include-regex="hotstart_it=*;param.nml.${baro}.hot*"

# print out ln -sf lines
echo -e "Copy the following lines into the run_*.yml file \n"
echo "  ln -sf hotstart_it=${last_restart}.nc hotstart.nc;
  ln -sf param.nml.${baro}.hot${last_restart} param.nml;"