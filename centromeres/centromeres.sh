#!/bin/bash -l
#PBS -N repeatOBserver
#PBS -l walltime=48:00:00
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l jobfs=400GB
#PBS -o job_output.out
#PBS -l storage=scratch/ht96+gdata/ht96

### Docker container: kathleenmclay/repeat_observer

module load singularity 

directory="/g/data/ht96/McLay_UQ/inversion_paper/centromeres"
species="seneciolautus" #no underscores in name
ref_genome="/g/data/ht96/McLay_UQ/inversion_paper/centromeres/2_D01_30.fasta"
haplotype="H0"
cpus="24"
mem="194560"

cd ${directory}

singularity exec --bind /g/data/ht96/McLay_UQ/inversion_paper/centromeres:/mnt /home/564/km6006/repeat_observer_latest.sif ${directory}/Setup_Run_Repeats.sh -i ${species} -f 2_D01_24.fasta -h ${haplotype} -c ${cpus} -m ${mem} -g FALSE


