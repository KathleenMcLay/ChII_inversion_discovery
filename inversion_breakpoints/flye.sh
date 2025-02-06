#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=450GB
#SBATCH --job-name=flye_H01_contigs
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe
#SBATCH --output=/home/uqkmcla4/scripts/flye_H01_contigs.log

IN_DIR="/QRISdata/Q6656/long_read_data/2_combined_raw_hifi"
OUT_DIR="/QRISdata/Q6656/long_read_data/3_denovo_contigs/H01"

# Run Flye on the combined PacBio HiFi reads for contig level denovo assembly
python /home/uqkmcla4/Flye/bin/flye --pacbio-hifi "$IN_DIR"/1_H01_combined.fastq.gz --out-dir "$OUT_DIR" --threads 48 --asm-coverage 40 --genome-size 2.8g
