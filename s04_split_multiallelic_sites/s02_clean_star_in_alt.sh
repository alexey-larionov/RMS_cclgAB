#!/bin/bash

# s02_clean_star_in_alt.sh
# Alexey Larionov, 22Sep2020

#SBATCH -J s02_clean_star_in_alt
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s02_clean_star_in_alt.log
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
echo "Started s02_clean_star_in_alt: $(date +%d%b%Y_%H:%M:%S)"
echo ""

base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
data_folder="${project_folder}/data/s04_split_multiallelic_sites"
scripts_folder="${project_folder}/scripts/s04_split_multiallelic_sites"
cd "${scripts_folder}"

# Files
source_vcf="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma.vcf"
clean_vcf="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln.vcf.gz"
clean_log="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln.log"

# Tools
tools_folder="${base_folder}/tools"
bcftools="${tools_folder}/bcftools/bcftools-1.10.2/bin/bcftools"

# Count variants
echo "Variant counts before cleaning:"
echo ""
"${bcftools}" +counts "${source_vcf}"
echo ""

# Filter variants
echo "Removing variants with * in ALT ..."
echo "(appeared after splitting multiallelic sitrs ?)"
echo ""

"${bcftools}" view \
  "${source_vcf}" \
  --exclude 'ALT="*"' \
  --output-type z \
  --threads 4 \
  --output-file "${clean_vcf}" \
  &> "${clean_log}"

"${bcftools}" index "${clean_vcf}"

echo "Variant counts after removing * in ALT: "
echo ""
"${bcftools}" +counts "${clean_vcf}"
echo ""

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
