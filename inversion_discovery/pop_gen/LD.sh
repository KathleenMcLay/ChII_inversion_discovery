### Calculate LD for each putative inversion
### GNU parallel by chromosome/inversion

module load bcftools 
 
FILE="/g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/${2}.vcf"
DIR="/g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/LD"
OUTPUT="/g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/LD/${2}"

# calculate LD for all samples 
/home/564/km6006/bin/vcftools --vcf ${FILE} --geno-r2 --maf 0.05 --thin 1000 --out ${DIR}/${2}_all.geno.ld
cat ${DIR}/${2}_all.geno.ld | perl /home/564/km6006/Scripts/inversion_paper/local_pca/emerald2windowldcounts.pl > /g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/LD/${2}_all_ld.txt


# determine the most common homozygous genotype group
cmn_gen="/home/564/km6006/Scripts/inversion_paper/local_pca/${2}_cmn_gen.txt"
genotypes="/home/564/km6006/Scripts/inversion_paper/local_pca/inversion_genotypes.txt"

declare -A genotype_counts

while IFS=$'\t' read -r inversion name genotype population; do
    if [ "$inversion" == "$2" ] && [ "$genotype" == "0" -o "$genotype" == "2" ]; then
        ((genotype_counts[$genotype]++))
    else
        echo "wrong chromosome or wrong genotype"
    fi
done < "$genotypes"

max_count=0
max_genotype=""
for gen in "${!genotype_counts[@]}"; do
    count="${genotype_counts[$gen]}"
    if [ "$count" -gt "$max_count" ]; then
        max_count="$count"
        max_genotype="$gen"
    fi
done

while IFS=$'\t' read -r inversion name genotype population; do
    if [ "$inversion" == "$2" ] && [ "$genotype" == "$max_genotype" ]; then
        echo -e "$(echo -n "$population" | tr -d '[:space:]')\t$(echo -n "$name" | tr -d '[:space:]')" >> "$cmn_gen"
    else
        echo "wrong chromosome or not the most common genotype"
    fi
done < "$genotypes"

# subset the vcf file to just the samples from the most common homozygous genotype group
bcftools view --threads 12 --samples-file ${cmn_gen} -o ${DIR}/${2}_cmn_gen_only.vcf ${FILE}

# calculate LD for the most common homozygous group only, in each chromosome with an inversion
/home/564/km6006/bin/vcftools --vcf ${DIR}/${2}_cmn_gen.vcf --geno-r2 --maf 0.05 --thin 1000 --out ${DIR}/${2}_cmn_gen.geno.ld
cat ${DIR}/${2}_cmn_gen.geno.ld | perl /home/564/km6006/Scripts/inversion_paper/local_pca/emerald2windowldcounts.pl > /g/data/ht96/McLay_UQ/inversion_paper/inv_pop_gen/LD/${2}_cmn_gen_ld.txt
