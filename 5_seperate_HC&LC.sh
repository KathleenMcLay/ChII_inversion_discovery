### Seperate high coverage and low coverage samples, for seperate SNP filtering pipelines
### GNU parallel by populations containing HC and LC samples (H01, H05, D01, D04, A03,& A07) 

module load bcftools
module load gatk

# create directory variable
DIR="/g/data/ht96/McLay_UQ/inversion_paper/5_joint_calling"

# Extract HC sample
bcftools view --threads 12 -O z -s ${1} -o ${DIR}/HC/${2}_HC_jntcl.vcf.gz ${DIR}/${2}_jntcl.vcf.gz
bcftools index --threads 12 ${DIR}/HC/${2}_HC_jntcl.vcf.gz

# Exclude HC sample 
bcftools view --threads 12 -O z -s ^${1} -o ${DIR}/LC/${2}_LC_jntcl.vcf.gz ${DIR}/${2}_jntcl.vcf.gz
bcftools index --threads 12 ${DIR}/LC/${2}_LC_jntcl.vcf.gz
