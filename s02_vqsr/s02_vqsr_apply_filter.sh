#!/bin/bash

# s02_vqsr_apply_filter.sh
# Alexey Larionov, 21Sep2020

#SBATCH -J s02_vqsr_apply_filter
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:00:00
#SBATCH --output=s02_vqsr_apply_filter.log
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
echo "Started s02_vqsr_apply_filter: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept/RMS_cclgAB"
data_folder="${project_folder}/data/s02_vqsr"
stats_folder="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_filt.stats"
rm -fr "${stats_folder}"
mkdir "${stats_folder}"
scripts_folder="${project_folder}/scripts/s02_vqsr"
cd "${scripts_folder}"

# Files
source_vcf="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr.vcf"
output_vcf="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_filt.vcf"
log="${data_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_filt.log"
stats_file="${stats_folder}/CCLG_GL_hg38.bwa_10bp_vqsr_filt.vchk"

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

# Remove filtered variants
echo ""
echo "Selecting variants ..."
"${java}" -Xmx40g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  --excludeFiltered \
  -V "${source_vcf}" \
  -o "${output_vcf}" \
  -nt 8 &>  "${log}"

# Count variants
echo ""
echo "Variant counts after applying vqsr filters:"
echo ""
bcftools +counts "${output_vcf}"

# Extract annotations
annotations="${output_vcf%.vcf}_annotations.txt"
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t%QUAL\t%DP\t%VQSLOD\n' "${output_vcf}" > "${annotations}"

echo ""
echo "--- FILTER ---"
echo ""
awk 'NR>1 {print $6}' "${annotations}" | sort |  uniq -c | sort -nr
echo ""

# Calculate bcfstats
echo ""
echo "Calculating bcfstats..."
echo ""
bcftools stats -s- "${output_vcf}" > "${stats_file}"

# Plot the stats (PDF requires texlive)
echo "Making plots..."
plot-vcfstats -p "${stats_folder}" "${stats_file}"
echo ""

# Completion message
echo ""
echo "Done: $(date +%d%b%Y_%H:%M:%S)"
echo ""
