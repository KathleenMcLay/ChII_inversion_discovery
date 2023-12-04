### Create a vcf file containing sites across all samples with a structural variant called 

module load singularity

# set variables 
name="merged"
ref="/scratch/ht96/km6006/2_D01_30.fasta"
outdir="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints"
merged="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints/merged.sites.vcf.gz"
input="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints/*-smoove.genotyped.vcf.gz"

# create merged file ./merged.sites.vcf.gz
singularity exec /home/564/km6006/smoove_latest.sif smoove merge --name ${name} --fasta ${ref} --outdir ${outdir} ${input}