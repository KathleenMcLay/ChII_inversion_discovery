### Calculate LD for each putative inversion
### GNU parallel by chromosome/inversion

module load bcftools
module load vcftools 
 
ld="/QRISdata/Q6656/new_lpca/LD_out"
directory="/QRISdata/Q6656/new_lpca/inversion_vcfs"

for file in $directory/*.vcf.gz; do 
    # calculate LD for all samples 
    filename=$(basename -- "$file")
    filename="${filename%.*.*}"
    current_inv="${filename#*_}"

    vcftools --gzvcf ${file} --geno-r2 --maf 0.05 --thin 1000 --out ${ld}/${filename}_all
    cat ${ld}/${filename}_all.geno.ld | perl /home/uqkmcla4/scripts/reformat/emerald2windowldcounts.pl > ${ld}/${filename}_all_ld.txt

    # determine the most common homozygous genotype group
    cmn_gen="${ld}/${filename}_cmn_gen.txt"
    genotypes="${ld}/${filename}_genotypelist.txt"
    echo $genotypes

    declare -A genotype_counts

    while IFS=$'\t' read -r name inversion genotype population; do 
        echo $current_inv
        if [ "$inversion" == "$current_inv" ] && [ "$genotype" == "0" -o "$genotype" == "2" ]; then
            echo $inversion
            echo $genotype
            ((genotype_counts[$genotype]++))
        else
            echo "not the current inversion"
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

    while IFS=$'\t' read -r name inversion genotype population; do
        if [ "$inversion" == "$current_inv" ] && [ "$genotype" == "$max_genotype" ]; then
            echo -e "$(echo -n "$population" | tr -d '[:space:]')\t$(echo -n "$name" | tr -d '[:space:]')" >> "$cmn_gen"
        else
            echo "not the current inversion"
        fi
    done < "$genotypes"

    # subset the vcf file to just the samples from the most common homozygous genotype group
    bcftools view --threads 12 --samples-file ${cmn_gen} -O z -o ${ld}/${filename}_cmn_gen.vcf.gz ${file}
    bcftools index -t ${ld}/${filename}_cmn_gen.vcf.gz

    # calculate LD for the most common homozygous group only, in each chromosome with an inversion
    vcftools --gzvcf ${ld}/${filename}_cmn_gen.vcf.gz --geno-r2 --maf 0.05 --thin 1000 --out ${ld}/${filename}_cmn_gen
    cat ${ld}/${filename}_cmn_gen.geno.ld | perl /home/uqkmcla4/scripts/reformat/emerald2windowldcounts.pl > ${ld}/${filename}_cmn_gen_ld.txt
done