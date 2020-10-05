#!/bin/bash

# s01_extract_with_QUAL.sh
# Alexey Larionov, 20Sep2020

#SBATCH -J s01_extract_with_QUAL
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s01_extract_with_QUAL.log
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
echo "Started s01_extract_with_QUAL: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
data_folder="${project_folder}/data/s00_source_vcfs"
scripts_folder="${project_folder}/scripts/s00_explore_source_vcfs"
cd "${scripts_folder}"

# Files
vcf="${data_folder}/CCLG_GL_hg38.bwa.vcf.gz"

# Tools
tools_folder="${base_folder}/tools"
bcftools_bin="${tools_folder}/bcftools/bcftools-1.10.2/bin"
export PATH="${bcftools_bin}:$PATH"

# Progress report
echo "vcf: ${vcf}"
echo ""

# Specific annotations
annotations="${vcf%.vcf.gz}_annotations_QUAL.txt"
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\n' "${vcf}" > "${annotations}"

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
