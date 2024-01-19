### Combine the joint SV genotyped samples into a single vcf 

module load singularity

# set variables
name="cohort_final"
outdir="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints"
input="/g/data/ht96/McLay_UQ/inversion_paper/breakpoints/joint_genotyped/*.vcf.gz"

singularity exec /home/564/km6006/smoove_latest.sif smoove paste --name ${name} --outdir ${outdir} ${input}