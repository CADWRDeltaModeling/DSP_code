# Get number of days from param.nml
# To run: bash check_progress.sh simulations/baseline_lhc_3
param="$(readlink -f "$1/param.nml")"
runtime=$(cat $param | grep 'rnday = ' | grep -Eo '[0-9]*')

# Get highest out2d file in outputs
highest=-1
for file in $1/outputs/out2d_*.nc
do
  if [[ $file =~ out2d_([0-9]+).* ]]
  then
    [[ "${BASH_REMATCH[1]}" -gt "$highest" ]] && highest=${BASH_REMATCH[1]}
  fi
done

# Print result in console
printf "\n\t\tThe simulation $1 has run $highest out of $runtime days.\n"