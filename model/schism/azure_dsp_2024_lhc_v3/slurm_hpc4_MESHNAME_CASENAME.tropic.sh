#!/bin/bash
#SBATCH --job-name={job_name}_{cname}_{meshname}_{baro}     # Job name
#SBATCH --partition=work        # Partition name to submit the job
#SBATCH --mail-type=NONE            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lily.tomkovic@water.ca.gov       # Where to send mail
#SBATCH --ntasks=144                      # Number of MPI ranks (or cores)
#SBATCH --nodes=6                       # Number of nodes
#SBATCH --ntasks-per-node=24             # How many tasks on each node
#SBATCH --output={output_log_file_base}_%j.log     # Standard output and error log
 
module load intel/2024.1
module load hmpt/2.30
 
mpiexec bash schism_hpc4.sh