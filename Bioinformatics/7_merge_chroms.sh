### Merge the the population files from 'split_chroms.sh' for each chromosome, ready for SNP filtering

module load bcftools

ODIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering"

### Merge the the population files from 'split_chroms.sh' for each chromosome, ready for SNP filtering

# array of chrom numbers 
chromosomes=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20")

# for each chromosome create a file containing a list of population files for that chromosome
for chrom in ${ODIR}/"${chromosomes[@]}"; do
    # Find files with the current prefix and store their names in a text file
    find "$(pwd)" -maxdepth 1 -type f -name "${chrom}_*" -exec readlink -f {} \; >> "${chrom}_files.txt"
done

# merge the chromosomes files 
for chrm in ${ODIR}/*.fofn;
do 
    out_vcf="${ODIR}/LC/merge.${chrm}.vcf.gz"
    bcftools merge --threads 12 -O z --file-list ${chrm} -o ${out_vcf} & 
done

wait 

# index the chromosome files 
for file in ${ODIR}/*.vcf.gz; 
do
    bcftools index --threads 12 -t ${file} &
done

wait