#!/bin/bash

# s01_vqsr_add_filter.sh
# Alexey Larionov, 21Sep2020

#SBATCH -J s01_vqsr_add_filter
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake-himem
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:00:00
#SBATCH --output=s01_vqsr_add_filter.log
#SBATCH --qos=INTR

## Modules section (required, do not remove)
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load texlive/2015 # for pdf plots by gatk?

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
echo "Started s01_vqsr_add_filter: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
source_vcf_folder="${project_folder}/data/s01_reduce_padding"

output_folder="${project_folder}/data/s02_vqsr"
rm -fr "${output_folder}"
mkdir "${output_folder}"

vqsr_data_folder="${output_folder}/vqsr_data"
mkdir "${vqsr_data_folder}"

logs_folder="${output_folder}/logs"
mkdir "${logs_folder}"

scripts_folder="${project_folder}/scripts/s02_vqsr"
cd "${scripts_folder}"

# Tools
tools_folder="${base_folder}/tools"
java="${tools_folder}/java/jre1.8.0_40/bin/java"
gatk="${tools_folder}/gatk/gatk-3.7-0/GenomeAnalysisTK.jar"
bcftools="${tools_folder}/bcftools/bcftools-1.10.2/bin/bcftools"

r_bin_folder="${tools_folder}/r/R-3.2.0/bin"
PATH="${r_bin_folder}:${PATH}" # for GATK plots

# Resources
ref_genome="${project_folder}/resources/ref_genome/hg38.bwa.fa"
targets_intervals="${project_folder}/resources/targets/nexterarapidcapture_exome_targetedregions_v1.2_hg38.bed"

# b38 VQSR Resources
b38_bundle_folder="${base_folder}/resources/gatk_bundle/b38/hg38bundle"
hapmap="${b38_bundle_folder}/hapmap_3.3.hg38.vcf.gz"
omni="${b38_bundle_folder}/1000G_omni2.5.hg38.vcf.gz"
phase1_1k_hc="${b38_bundle_folder}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
dbsnp_144="${b38_bundle_folder}/dbsnp_144.hg38.vcf.gz" # The current dbsnp build in Sep2020 is 153
mills="${b38_bundle_folder}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# VQSR TS settings
SNP_TS="97.0"
INDEL_TS="95.0"

# --- Train vqsr snp model --- #

# Progress report
echo "Started training vqsr snp model"

# File names
dataset="CCLG_GL_hg38.bwa_10bp"

source_vcf="${source_vcf_folder}/${dataset}.vcf"
recal_snp="${vqsr_data_folder}/${dataset}_snp.recal"
plots_snp="${vqsr_data_folder}/${dataset}_snp_plots.R"
tranches_snp="${vqsr_data_folder}/${dataset}_snp.tranches"
log_train_snp="${logs_folder}/${dataset}_snp_train.log"

# Train vqsr snp model
"${java}" -Xmx80g -jar "${gatk}" \
  -T VariantRecalibrator \
  -input "${source_vcf}" \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${hapmap}" \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 "${omni}" \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 "${phase1_1k_hc}" \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${dbsnp_144}" \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
  -recalFile "${recal_snp}" \
  -tranchesFile "${tranches_snp}" \
  -rscriptFile "${plots_snp}" \
  --target_titv 3.2 \
  -mode SNP \
  -tranche 100.0 -tranche 99.0 -tranche 97.0 -tranche 95.0  -tranche 90.0 \
  -nt 8 &>  "${log_train_snp}"

# Progress report
echo "Completed training vqsr snp model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Apply vqsr snp model --- #

# Progress report
echo "Started applying vqsr snp model"

# File names
vqsr_snp_vcf="${vqsr_data_folder}/${dataset}_snp_vqsr.vcf"
log_apply_snp="${logs_folder}/${dataset}_snp_apply.log"

# Apply vqsr snp model
"${java}" -Xmx80g -jar "${gatk}" \
  -T ApplyRecalibration \
  -input "${source_vcf}" \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -recalFile "${recal_snp}" \
  -tranchesFile "${tranches_snp}" \
  --ts_filter_level "${SNP_TS}" \
  -mode SNP \
  -o "${vqsr_snp_vcf}" \
  -nt 8 &>  "${log_apply_snp}"

# Progress report
echo "Completed applying vqsr snp model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Train vqsr indel model --- #

# Progress report
echo "Started training vqsr indel model"

# File names
recal_indel="${vqsr_data_folder}/${dataset}_indel.recal"
plots_indel="${vqsr_data_folder}/${dataset}_indel_plots.R"
tranches_indel="${vqsr_data_folder}/${dataset}_indel.tranches"
log_train_indel="${logs_folder}/${dataset}_indel_train.log"

# Train vqsr indel model
"${java}" -Xmx80g -jar "${gatk}" \
  -T VariantRecalibrator \
  -input "${vqsr_snp_vcf}" \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -resource:mills,known=false,training=true,truth=true,prior=12.0 "${mills}" \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "${dbsnp_144}" \
  -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
  -recalFile "${recal_indel}" \
  -tranchesFile "${tranches_indel}" \
  -rscriptFile "${plots_indel}" \
  -tranche 100.0 -tranche 99.0 -tranche 97.0 -tranche 95.0 -tranche 90.0 \
  --maxGaussians 4 \
  -mode INDEL \
  -nt 8 &>  "${log_train_indel}"

# Progress report
echo "Completed training vqsr indel model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# --- Apply vqsr indel model --- #

# Progress report
echo "Started applying vqsr indel model"

# File names
out_vcf="${output_folder}/${dataset}_vqsr.vcf"
log_apply_indel="${logs_folder}/${dataset}_indel_apply.log"

# Apply vqsr indel model
"${java}" -Xmx80g -jar "${gatk}" \
  -T ApplyRecalibration \
  -input "${vqsr_snp_vcf}" \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -recalFile "${recal_indel}" \
  -tranchesFile "${tranches_indel}" \
  --ts_filter_level "${INDEL_TS}" \
  -mode INDEL \
  -o "${out_vcf}" \
  -nt 8 &>  "${log_apply_indel}"

# Progress report
echo "Completed applying vqsr indel model: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Count variants
echo "Variant counts:"
echo ""
"${bcftools}" +counts "${out_vcf}"

# Extract annotations
annotations="${out_vcf%.vcf}_annotations.txt"
"${bcftools}" query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%VQSLOD\n' "${out_vcf}" > "${annotations}"

echo ""
echo "--- FILTER ---"
echo ""
awk 'NR>1 {print $6}' "${annotations}" | sort |  uniq -c | sort -nr
echo ""

# Completion message
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
