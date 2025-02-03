#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180GB
#SBATCH --job-name=bam2fastq
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/bam2fastq.log

OUTPUT_DIR="/QRISdata/Q6656/long_read_data/2_combined_raw_hifi"

# Combine all BAM files into a single FASTQ file
BASE_DIR="/QRISdata/Q6656/long_read_data"
for folder in "$BASE_DIR"/1_raw_hifi/*; do
    if [ "$folder" = "9_D03" ] && [ -d "$folder" ]; then
        sample_name=$(basename "$folder")
        # Get list of BAM files
        BAMS=($folder/*.bam)
        # Check if there are any BAM files before running bam2fastq
        if [ ${#BAMS[@]} -gt 0 ]; then
            for bam_file in "${BAMS[@]}"; do
                echo "Indexing: $bam_file"
                /home/uqkmcla4/pbindex "$bam_file"
            done
            /home/uqkmcla4/bam2fastq -o "$OUTPUT_DIR/${sample_name}_combined" "${BAMS[@]}"
            echo "Processed sample: $sample_name"
        else
            echo "No BAM files found for sample: $sample_name"
        fi
    fi
done

echo "All samples processed!"
