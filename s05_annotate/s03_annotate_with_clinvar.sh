#!/bin/bash

# s03_annotate_with_clinvar.sh
# Alexey Larionov, 23Sep2020

# Check consistency of chromosome naming in the data and ClinVar

#SBATCH -J s03_annotate_with_clinvar
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s03_annotate_with_clinvar.log
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
echo "Started s03_annotate_with_clinvar: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
clinvar_folder="${project_folder}/resources/clinvar"
data_folder="${project_folder}/data/s05_annotate"
scripts_folder="${project_folder}/scripts/s05_annotate"
cd "${scripts_folder}"

# Tools
tools_folder="${base_folder}/tools"
bcftools="${tools_folder}/bcftools/bcftools-1.10.2/bin/bcftools"

# Files
source_vcf="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag_ac_id.vcf.gz"
output_vcf="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag_ac_id_clinvar.vcf.gz"
clinvar_log="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma_cln_tag_ac_id_clinvar.log"
clinvar_vcf="${clinvar_folder}/clinvar_20200905_chr.vcf.gz"

# Progress report
echo "--- Input and output files ---"
echo ""
echo "source_vcf: ${source_vcf}"
echo "output_vcf: ${output_vcf}"
echo "clinvar_vcf: ${clinvar_vcf}"
echo ""
echo "Working..."

# Annotate using bcftools
"${bcftools}" annotate "${source_vcf}" \
  --annotations "${clinvar_vcf}" \
  --columns INFO \
  --output "${output_vcf}" \
  --output-type z \
  --threads 4 \
  &> "${clinvar_log}"

# Index output vcf
"${bcftools}" index "${output_vcf}"

# Summary of INFO fields

echo ""
echo "Number of INFO fields in source vcf:"
"${bcftools}" view -h "${source_vcf}" | grep ^##INFO | wc -l
echo ""

echo ""
echo "Number of INFO fields in ClinVar-annotated vcf:"
"${bcftools}" view -h "${output_vcf}" | grep ^##INFO | wc -l
echo ""

echo "List of INFO fields in ClinVar-annotated vcf:"
"${bcftools}" view -h "${output_vcf}" | grep ^##INFO
echo ""

# Progress report
echo ""
echo "selected lines from output_vcf:"
echo "" # print some lines by --query to extract fields from bcf/vcf into user-frienlty format #-i include #-f format
"${bcftools}" query \
  -i 'ALLELEID != "."' \
  -f '%ID\t%ALLELEID\t%CLNSIG\t%CLNDN\n' "${output_vcf}" | head
echo ""

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
