### Calculate Fst for each putative inversion
### GNU parallel by chromosome/inversion

#activate the miniconda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

module load bcftools 

FILE="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only/HF_allpops_VO_0.8_final.vcf.gz"
DIR="/g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/inv_files"
OUTPUT="/g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/fst/${2}"

# Filter to the chromosome
bcftools filter ${FILE} --regions ${1} --output ${DIR}/${2}.vcf 

G0="${DIR}/${2}_G0.txt"
G2="${DIR}/${2}_G2.txt"

genotypes="/home/564/km6006/Scripts/inversion_paper/local_pca/inversion_genotypes.txt"

# For the current inversion add each homozygous sample to the correct list, based on their genotype
while IFS=$'\t' read -r inversion genotype name population; do
    echo "current $inversion has "$genotype" genotype"
    if [ "$inversion" == "$2" ] && [ "$genotype" -eq 0 ]; then
        echo -e "$name\t$population" >> "$G0"
    elif [ "$inversion" == "$2" ] && [ "$genotype" -eq 2 ]; then
        echo -e "$name\t$population" >> "$G2"
    else
        echo "het sample"
    fi
done < "$genotypes"

wait

# Calculate Weir and Cockerham’s Fst between the two homozygous genotypes 
/home/564/km6006/bin/vcftools \
    --vcf ${DIR}/${2}.vcf \
    --weir-fst-pop ${G0} \
    --weir-fst-pop ${G2} \
    --out ${OUTPUT}

