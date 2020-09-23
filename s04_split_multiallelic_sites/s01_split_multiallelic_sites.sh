#!/bin/bash

# s01_split_multialelic_sites.sh
# Alexey Larionov, 22Sep2020

#SBATCH -J s01_split_multialelic_sites
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s01_split_multialelic_sites.log
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
echo "Started s01_split: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"

scripts_folder="${project_folder}/scripts/s04_split_multiallelic_sites"
cd "${scripts_folder}"

source_data_folder="${project_folder}/data/s03_hard_filters"

results_data_folder="${project_folder}/data/s04_split_multiallelic_sites"
rm -fr "${results_data_folder}"
mkdir -p "${results_data_folder}"

# Tools
tools_folder="${base_folder}/tools"
java="${tools_folder}/java/jre1.8.0_40/bin/java"
gatk="${tools_folder}/gatk/gatk-3.7-0/GenomeAnalysisTK.jar"
bcftools="${tools_folder}/bcftools/bcftools-1.10.2/bin/bcftools"

# Resources
ref_genome="${project_folder}/resources/ref_genome/hg38.bwa.fa"
targets_intervals="${project_folder}/resources/targets/nexterarapidcapture_exome_targetedregions_v1.2_hg38.bed"

# Files
source_vcf="${source_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp.vcf"
split_vcf="${results_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma.vcf"
split_log="${results_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp_sma.log"

# Count variants
echo "Variant counts before splitting:"
echo ""
"${bcftools}" +counts "${source_vcf}"
echo ""

# Split variants
echo "Splitting multiallelic variants ..."
# Split ma sites
"${java}" -Xmx60g -jar "${gatk}" \
  -T LeftAlignAndTrimVariants \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -o "${split_vcf}" \
  --splitMultiallelics &> "${split_log}"
echo ""

# Count variants
echo "Variant counts after splitting: "
echo ""
"${bcftools}" +counts "${split_vcf}"
echo ""

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
