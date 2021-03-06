#!/bin/bash

# s02_run_rmd.sh
# Alexey Larionov, 12Sep2020

#SBATCH -J s02_run_rmd
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake-himem
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s02_run_rmd.log
#SBATCH --qos=INTR

## Modules section (required, do not remove)
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4

## Set initial working folder
cd "${SLURM_SUBMIT_DIR}"

## Report settings and run the job
echo "Job id: ${SLURM_JOB_ID}"
echo "Allocated node: $(hostname)"
echo "$(date)"
echo ""
echo "Job name: ${SLURM_JOB_NAME}"
echo ""
echo "Initial working folder:"
echo "${SLURM_SUBMIT_DIR}"
echo ""
echo " ------------------ Job progress ------------------ "
echo ""

# Stop at runtime errors
set -e

# Start message
echo "Started s02_run_rmd: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
scripts_folder="${project_folder}/scripts/s03_hard_filters"
cd "${scripts_folder}"

# Tools
module load pandoc # required by rmarkdown for html rendering
module load R/3.6 # make sure its the same R version that was used in rmd script

# Execute rmd script and render html log
R -e "library('rmarkdown'); render('s02_check_data.Rmd')"

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
