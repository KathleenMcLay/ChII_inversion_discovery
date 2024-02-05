### Calculate Pi and Dxy for each putative inversion
### GNU parallel by chromosome/inverison

#activate the miniconda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

module load bcftools 

DIR="/g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/pixy"
FILE="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/all_data/AD_allpops_VO_0.8_final_PIXY.vcf.gz"

# filter to the chromosome
bcftools filter ${FILE} --regions ${1} --output ${DIR}/${2}.vcf 

wait

# Zip and index the inversion files 
bgzip -@ 12 ${DIR}/${2}.vcf
tabix ${DIR}/${2}.vcf.gz

wait

POPS="${DIR}/${2}_POPS.txt"
genotypes="/home/564/km6006/Scripts/inversion_paper/local_pca/inversion_genotypes.txt"

# for each inversion check the genotype, and assign the genotype value as the 'population' for that individual
while IFS=$'\t' read -r inversion name genotype population; do
    echo "current $inversion has "$genotype" genotype"
    if [ "$inversion" == "$2" ] && [ "$genotype" -eq 0 ]; then
        echo -e "$name\t$genotype" >> "$POPS"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 1 ]; then
        echo -e "$name\t$genotype" >> "$POPS"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 2 ]; then
        echo -e "$name\t$genotype" >> "$POPS"        
    else
        echo "no match"
    fi
done < "$genotypes"

# Run pixy - calculate pi and dxy values
pixy --stats pi dxy --vcf ${DIR}/${2}.vcf.gz --populations ${POPS} --window_size 10000 --output_prefix ${2} --output_folder ${DIR} --n_cores 12