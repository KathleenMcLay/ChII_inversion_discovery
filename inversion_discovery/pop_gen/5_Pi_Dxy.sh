### Calculate Pi and Dxy for each putative inversion

#activate the miniconda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

module load bcftools 

pixy="/g/data/ht96/McLay_UQ/inversion_paper/pixy"
popdir="/g/data/ht96/McLay_UQ/inversion_paper/pixy"

poplist="/home/564/km6006/Scripts/inversion_paper/pops.txt"

while IFS= read -r pop; do # for each population
    echo "$pop"
  
    gzvcf="${popdir}/${pop}.vcf.gz"

    for file in "${pixy}"/"${pop}"_*genotypelist.txt; do  # for each inversion
        # calculate LD for all samples 
        filename=$(basename -- "$file")
        filename="${filename%_*}"
        current_inv="${filename#*_}"

        # filter vcf to the chromosome
        scaffold="${current_inv%%:*}"
        bcftools filter ${gzvcf} --regions ${scaffold} --output ${pixy}/${filename}.vcf

        # Zip and index the inversion files 
        bgzip -@ 12 ${pixy}/${filename}.vcf
        tabix ${pixy}/${filename}.vcf.gz

        POPS="${pixy}/${filename}_POPS.txt"
        genotypes="$file"

        # for each inversion check the genotype, and assign the genotype value as the 'population' for that individual
        while IFS=$'\t' read -r name inversion genotype population; do
            echo "current $inversion has "$genotype" genotype"
            if [ "$inversion" == "$current_inv" ] && [ "$genotype" -eq 0 ]; then
                echo -e "$name\t$genotype" >> "$POPS"
            elif [ "$inversion" == "$current_inv" ] && [ "$genotype" -eq 1 ]; then
                echo -e "$name\t$genotype" >> "$POPS"
            elif [ "$inversion" == "$current_inv" ] && [ "$genotype" -eq 2 ]; then
                echo -e "$name\t$genotype" >> "$POPS"        
            else
                echo "no match"
            fi
        done < "$genotypes"

        # Run pixy - calculate pi and dxy values
        pixy --stats pi dxy --vcf ${pixy}/${filename}.vcf.gz --populations ${POPS} --window_size 10000 --output_prefix ${filename} --output_folder ${pixy} --n_cores 12
    done
done < "$poplist"    