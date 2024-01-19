### Merge the the population files from 'split_chroms.sh' for each chromosome, ready for SNP filtering

module load bcftools

DIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/pop_chrom"
ODIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering"

cov="HC"

### Merge the the population files from 'split_chroms.sh' for each chromosome, ready for SNP filtering

# index the population chromosome files 
for file in ${DIR}/*.vcf.gz; 
do
    bcftools index --threads 24 -t ${file} &
done

wait

# array of chrom numbers 
chromosomes=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")

# for each chromosome create a file containing a list of population files for that chromosome
for chrom in "${chromosomes[@]}"; do
    # Find files with the current prefix and store their names in a text file
    find "${DIR}" -maxdepth 1 -type f -name "${chrom}_*.vcf.gz" -exec readlink -f {} \; >> "${DIR}/${chrom}.fofn"
done

# merge the chromosomes files 
for chrm in ${DIR}/*.fofn;
do 
    chrm_name=$(basename "$chrm" .fofn)
    out_vcf="${ODIR}/${chrm_name}_${cov}allpops.vcf.gz"
    bcftools merge --threads 6 -O z --file-list ${chrm} -o ${out_vcf} &
done

wait 

# index the chromosome files 
for file in ${ODIR}/*.vcf.gz; 
do
    bcftools index --threads 24 -t ${file} &
done

wait