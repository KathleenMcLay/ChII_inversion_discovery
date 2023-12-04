### Genotype each samples for structural variant sites found across all samples*
# GNU parallel by sample 

module load singularity

# set variables 
name="${1}_joint"
ref="/scratch/ht96/km6006/2_D01_30.fasta"
outdir="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints/joint_genotyped"
merged="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints/merged.sites.vcf.gz"
input="/g/data/ht96/McLay_UQ/inversion_paper/3_clean/${1}_PCRm_cln_srt.bam"

# genotype each sample for all SV sites 
singularity exec /home/564/km6006/smoove_latest.sif smoove genotype -d -x -p 1 --name ${name} --fasta ${ref} --outdir ${outdir} --vcf ${merged} ${input}