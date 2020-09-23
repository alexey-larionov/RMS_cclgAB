#!/bin/bash

# s01_reduce_padding.sh
# Alexey Larionov, 21Sep2020

#SBATCH -J s01_reduce_padding
#SBATCH -A TISCHKOWITZ-SL2-CPU
#SBATCH -p skylake
#SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:00:00
#SBATCH --output=s01_reduce_padding.log
#SBATCH --qos=INTR

## Modules section (required, do not remove)
. /etc/profile.d/modules.sh
module purge
module load rhel7/default-peta4
module load texlive/2015 # for ploting pdf summary for vcfstats

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
echo "Started s01_reduce_padding: $(date +%d%b%Y_%H:%M:%S)"
echo ""

# Folders
base_folder="/rds/project/erf33/rds-erf33-medgen"
project_folder="${base_folder}/users/mae/RMS_2020Sept"
source_folder="${project_folder}/data/s00_source_vcfs"

output_folder="${project_folder}/data/s01_reduce_padding"
rm -fr "${output_folder}"
mkdir "${output_folder}"

stats_folder="${output_folder}/CCLG_GL_hg38.bwa_10bp.stats"
rm -fr "${stats_folder}"
mkdir "${stats_folder}"

scripts_folder="${project_folder}/scripts/s01_reduce_padding"
cd "${scripts_folder}"

# Files
source_vcf="${source_folder}/CCLG_GL_hg38.bwa.vcf"
output_vcf="${output_folder}/CCLG_GL_hg38.bwa_10bp.vcf"
log="${output_folder}/CCLG_GL_hg38.bwa_10bp.log"
stats_file="${stats_folder}/CCLG_GL_hg38.bwa_10bp.vchk"

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
echo "Variant counts in the source vcf with 100bp padding:"
echo ""
bcftools +counts "${source_vcf}"

# Reduce padding
echo ""
echo "Selecting variants ..."
"${java}" -Xmx40g -jar "${gatk}" \
  -T SelectVariants \
  -R "${ref_genome}" \
  -L "${targets_intervals}" -ip 10 \
  -V "${source_vcf}" \
  -o "${output_vcf}" \
  -nt 8 &>  "${log}"

# Count variants
echo ""
echo "Variant counts after changing padding to 10 bp:"
echo ""
bcftools +counts "${output_vcf}"

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
