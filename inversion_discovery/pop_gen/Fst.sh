### For each inversion, calculate Fst for the scaffold between the 0 and 2 genotype groups 

module load bcftools
module load vcftools 

fst="/QRISdata/Q6656/new_lpca/fst_out"
ld="/QRISdata/Q6656/new_lpca/LD_out"
popdir="/QRISdata/Q6656/new_lpca/population_vcfs"

poplist="/home/uqkmcla4/scripts/current/pops.txt"

while IFS= read -r pop; do # for each population
    echo "$pop"
    
    gzvcf="${popdir}/${pop}.vcf.gz"

    for file in "$ld"/"${pop}"_*genotypelist.txt; do  # for each inversion 
        filename=$(basename -- "$file")
        filename="${filename%_*}"
        current_inv="${filename#*_}"

        # filter vcf to the chromosome
        scaffold="${current_inv%%:*}"
        bcftools filter ${gzvcf} --regions ${scaffold} --output ${fst}/${filename}_${scaffold}.vcf 

        #create files for genotype lists
        G0="${fst}/${filename}_G0.txt"
        G2="${fst}/${filename}_G2.txt"

        # read in the sample genotype data for all inversions
        genotypes="$file" 

        # For the current inversion add each homozygous sample to the correct list, based on their genotype
        while IFS=$'\t' read -r name inversion genotype population; do
            echo "current $inversion has $genotype genotype"
            if [ "$inversion" == "$current_inv" ] && [ "$genotype" -eq 0 ]; then
                echo -e "$name\t$population" >> "$G0"
            elif [ "$inversion" == "$current_inv" ] && [ "$genotype" -eq 2 ]; then
                echo -e "$name\t$population" >> "$G2"
            else
                echo "het sample"
            fi
        done < "$genotypes"

        # Calculate Weir and Cockerham’s Fst between the two homozygous genotypes
        vcftools \
            --vcf ${fst}/${filename}_${scaffold}.vcf \
            --weir-fst-pop "$G0" \
            --weir-fst-pop "$G2" \
            --fst-window-size 10000 \
            --fst-window-step 10000 \
            --out "${fst}/${filename}_${scaffold}"
    done
done < "$poplist"