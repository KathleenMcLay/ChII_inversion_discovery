### Call structural variants for each sample
# GNU parallel by sample 

module load singularity

# set variables 
name="${1}"
ref="/scratch/ht96/km6006/2_D01_30.fasta"
outdir="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints"
merged="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints/merged.sites.vcf.gz"
input="/g/data/ht96/McLay_UQ/inversion_paper/3_clean/${1}_PCRm_cln_srt.bam"

# call structual variants 
singularity exec /home/564/km6006/smoove_latest.sif smoove call --name ${1} --outdir ${outdir} --fasta ${ref} -p 1 --genotype ${input}