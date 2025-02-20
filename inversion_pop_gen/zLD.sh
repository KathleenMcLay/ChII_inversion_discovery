### Calculate LD for each putative inversion

module load bcftools
module load vcftools 
 
ld="/QRISdata/Q6656/new_lpca/LD_out"
popdir="/QRISdata/Q6656/new_lpca/population_vcfs"

poplist="/home/uqkmcla4/scripts/current/pops.txt"

while IFS= read -r pop; do # for each population
    echo "$pop"
  
    gzvcf="${popdir}/${pop}.vcf.gz"

    for file in "$ld"/"${pop}"_*genotypelist.txt; do  # for each inversion
        # calculate LD for all samples 
        filename=$(basename -- "$file")
        filename="${filename%_*}"
        current_inv="${filename#*_}"

        # filter vcf to the chromosome
        scaffold="${current_inv%%:*}"
        # echo $scaffold
        # bcftools filter ${gzvcf} --regions ${scaffold} --output ${ld}/${filename}_${scaffold}.vcf

        # vcftools --vcf ${ld}/${filename}_${scaffold}.vcf --geno-r2 --maf 0.05 --thin 1000 --max-missing-count 0 --out ${ld}/${filename}_all
        # cat ${ld}/${filename}_all.geno.ld | perl /home/uqkmcla4/scripts/reformat/emerald2windowldcounts.pl > ${ld}/${filename}_all_ld.txt

        ### determine the most common homozygous genotype group
        
        # create a file to hold data 
        cmn_gen="${ld}/${filename}_cmn_gen.txt"
        
        # read in the sample genotype data for all inversions
        genotypes="$file" 

        declare -A genotype_counts

        # for the current inversion count each hom group 
        while IFS=$'\t' read -r name inversion genotype population; do 
            if [ "$inversion" == "$current_inv" ] && [ "$genotype" == "0" -o "$genotype" == "2" ]; then
                ((genotype_counts[$genotype]++))
            else
                echo "not the current inversion"
            fi
        done < "$genotypes"

        # determine which hom group has the most samples and set as the most common
        max_count=0
        max_genotype=""
        for gen in "${!genotype_counts[@]}"; do
            count="${genotype_counts[$gen]}"
            if [ "$count" -gt "$max_count" ]; then
                max_count="$count"
                max_genotype="$gen"
            fi
        done

        # write to file 
        while IFS=$'\t' read -r name inversion genotype population; do
            if [ "$inversion" == "$current_inv" ] && [ "$genotype" == "$max_genotype" ]; then
                echo -e "$(echo -n "$name" | tr -d '[:space:]')" >> "$cmn_gen"
            else
                echo "not the current inversion"
            fi
        done < "$genotypes"

        # subset the vcf file to just the samples from the most common homozygous genotype group
        bcftools view --threads 12 --samples-file ${cmn_gen} -O z -o ${ld}/${filename}_cmn_gen.vcf.gz ${ld}/${filename}_${scaffold}.vcf
        bcftools index -t ${ld}/${filename}_cmn_gen.vcf.gz

        # calculate LD for the most common homozygous group only, in each chromosome with an inversion
        vcftools --gzvcf ${ld}/${filename}_cmn_gen.vcf.gz --geno-r2 --maf 0.05 --thin 1000 --max-missing-count 0 --out ${ld}/${filename}_cmn_gen
        cat ${ld}/${filename}_cmn_gen.geno.ld | perl /home/uqkmcla4/scripts/reformat/emerald2windowldcounts.pl > ${ld}/${filename}_cmn_gen_ld.txt
    done
done < "$poplist"