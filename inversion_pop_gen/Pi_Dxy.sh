#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=1000GB
#SBATCH --job-name=pixy
#SBATCH --time=5:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/pixy.txt 

set -x  # Print each command as it executes
set -e  # Exit on any error

source ~/miniconda3/etc/profile.d/conda.sh

# Activate environment with error checking
conda activate pixy || { echo "Failed to activate pixy environment"; exit 1; }

# Verify activation
which pixy
conda info --envs

module load htslib
module load bcftools
echo "bcftools loaded"

# directory
dir="/QRISdata/Q6656/chapter_II/new_inv_discovery/alt"
# vcf file including monomorphic and variant sites 
gzvcf="/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/pi_dxy/sf7_noD1_reheader.vcf.gz"
# list of inversions 
inversions="/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/inversions.csv"

# Convert potential Windows line endings and clean the file
tr -d '\r' < "$inversions" > "${inversions}.tmp"

# Debug: print first few lines
echo "First few lines of inversions file:"
head -n 3 "${inversions}.tmp"

while IFS=, read -r inv population scaffold || [[ -n "$inv" ]]; do
    # Debug: print current values
    echo "Processing: inv='$inv' population='$population' scaffold='$scaffold'"
    # Skip empty lines
    [[ -z "$inv" ]] && continue
    scaffold=$(echo ${scaffold} | tr -d '\n\r')
    bcftools view --threads 12 -O z \
        -S "${dir}/${population}_sample_list.txt" \
        -r "${scaffold}" \
        --output "${dir}/pi_dxy/${inv}_${population}_pixy.vcf.gz" \
        "${gzvcf}"

    echo "filtered vcf created"
    # index the filtered VCF with tabix (tabix index is required for pixy)
    tabix -p vcf ${dir}/pi_dxy/${inv}_${population}_pixy.vcf.gz
    echo "vcf indexed"

    # create files for genotype lists
    POPS="${dir}/pi_dxy/${inv}_${population}_POPS.txt"
    # read in the sample genotype data for all inversions 
    genotypes="${dir}/genotypes.csv"

    # for each individual add them to the corre
    while IFS="," read -r pop sample inversion genotype; do
        if [ "$inversion" == "$inv" ]; then
            echo -e "$sample\t$genotype" >> "$POPS"
            else
            echo "inversion: $inversion not equal to $inv"
        fi
    done < "$genotypes"

    # Run pixy - calculate pi and dxy values
    pixy --stats pi dxy --vcf ${dir}/pi_dxy/${inv}_${population}_pixy.vcf.gz --populations ${POPS} --window_size 10000 --output_prefix ${inv}_${population} --output_folder ${dir}/pi_dxy --n_cores 12

done < "${inversions}.tmp"
  