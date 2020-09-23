#!/bin/bash

# s01_add_tags_check_ac.sh
# Add tags, check for variants with AC==0 (remove, if any)
# Alexey Larionov, 22Sep2020

#SBATCH -J s01_add_tags_check_ac
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s01_add_tags_check_ac.log
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
echo "Started s01_add_tags_check_ac: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"

scripts_folder="${project_folder}/scripts/s05_annotate"
cd "${scripts_folder}"

source_data_folder="${project_folder}/data/s04_split_multiallelic_sites"

results_data_folder="${project_folder}/data/s05_annotate"
rm -fr "${results_data_folder}"
mkdir -p "${results_data_folder}"

# Files
source_vcf="${source_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln.vcf.gz"

tagged_vcf="${results_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag.vcf.gz"
tagging_log="${results_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag.log"

ac_vcf="${results_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag_ac.vcf.gz"
ac_log="${results_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag_ac.log"

# Tools
tools_folder="${base_folder}/tools"
bcftools="${tools_folder}/bcftools/bcftools-1.10.2/bin/bcftools"

# Count variants
echo "Variant counts:"
echo ""
"${bcftools}" +counts "${source_vcf}"
echo ""

# Add tags
echo "Adding tags ..."
echo ""

"${bcftools}" +fill-tags "${source_vcf}" \
  --output "${tagged_vcf}" \
  --output-type z \
  --threads 4 \
  &> "${tagging_log}"

# Index output vcf
"${bcftools}" index "${tagged_vcf}"

echo "Removing variants with AC==0 ..."
echo ""
"${bcftools}" view "${tagged_vcf}" \
  --exclude 'AC=0' \
  --output-type z \
  --threads 4 \
  --output-file "${ac_vcf}" \
  &> "${ac_log}"

"${bcftools}" index "${ac_vcf}"

echo "Counts after removing variants with AC=0: "
echo ""
"${bcftools}" +counts "${ac_vcf}"
echo ""

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
