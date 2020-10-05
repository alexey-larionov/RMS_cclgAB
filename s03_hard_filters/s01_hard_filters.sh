#!/bin/bash

# s01_hard_filters.sh
# Alexey Larionov, 21Sep2020

#SBATCH -J s01_hard_filters
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=s01_hard_filters.log
#SBATCH --qos=INTR

## Modules section (required, do not remove)
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load texlive/2015

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
echo "Started s01_hard_filters: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
source_data_folder="${project_folder}/data/s02_vqsr"

output_data_folder="${project_folder}/data/s03_hard_filters"
rm -fr "${output_data_folder}"
mkdir "${output_data_folder}"

stats_folder="${output_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp.stats"
mkdir "${stats_folder}"

scripts_folder="${project_folder}/scripts/s03_hard_filters"
cd "${scripts_folder}"

# Files
source_vcf="${source_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_filt.vcf"

qual_vcf="${output_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual.vcf"
qual_log="${output_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual.log"

dp_vcf="${output_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp.vcf"
dp_log="${output_data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp.log"

stats_file="${stats_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_qual_dp.vchk"

# Tools
tools_folder="${base_folder}/tools"

java="${tools_folder}/java/jre1.8.0_40/bin/java"
gatk="${tools_folder}/gatk/gatk-3.7-0/GenomeAnalysisTK.jar"

bcftools_bin="${tools_folder}/bcftools/bcftools-1.10.2/bin"
python_bin="${tools_folder}/python/python_3.8.3/bin" # Includes mathplotlib for bcftools vcfcheck
export PATH="${bcftools_bin}:${python_bin}:$PATH"

# Resources
ref_genome="${project_folder}/resources/ref_genome/hg38.bwa.fa"
targets_intervals="${project_folder}/resources/targets/nexterarapidcapture_exome_targetedregions_v1.2_hg38.bed"

# Count variants
echo "Variant counts in the source vcf:"
echo ""
bcftools +counts "${source_vcf}"

# Filter by QUAL
echo ""
echo "Filtering by QUAL >= 150 ..."
"${java}" -Xmx40g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -o "${qual_vcf}" \
  -select "QUAL >= 150" \
  &>  "${qual_log}"

# Count variants
echo ""
echo "Variant counts after applying QUAL filter:"
echo ""
bcftools +counts "${qual_vcf}"

# Filter by DP
echo ""
echo "Filtering by DP >= 840 (10 DP x 84 samples) ..."
"${java}" -Xmx40g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${qual_vcf}" \
  -o "${dp_vcf}" \
  -select "DP >= 840" \
  &>  "${dp_log}"

# Count variants
echo ""
echo "Variant counts after applying DP filter:"
echo ""
bcftools +counts "${dp_vcf}"

# Extract annotations
annotations="${dp_vcf%.vcf}_annotations.txt"
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%DP\t%VQSLOD\n' "${dp_vcf}" > "${annotations}"

# Calculate bcfstats
echo ""
echo "Calculating bcfstats..."
echo ""
bcftools stats -s- "${dp_vcf}" > "${stats_file}"

# Plot the stats (PDF requires texlive)
echo "Making plots..."
plot-vcfstats -p "${stats_folder}" "${stats_file}"
echo ""

# Completion message
echo ""
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
