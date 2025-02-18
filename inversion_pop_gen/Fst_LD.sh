dir="/QRISdata/Q6656/chapter_II/new_inv_discovery/alt"
gzvcf="/QRISdata/Q6656/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz"

inversions="/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/4_regions_results_summary.csv"

while IFS="," read -r inv population scaffold; do
    echo "current inversion is ${inv}"

    # check if the files already exist
    if [ -f "${dir}/pop_vcfs/${inv}_${population}.vcf.gz" ]; then
        echo "vcf file already exists"
        continue
    fi

    # filter the vcf to the current population and scaffold 
    bcftools view --threads 12 -O z --regions ${scaffold} -S ${dir}/${population}_sample_list.txt --output ${dir}/pop_vcfs/${inv}_${population}.vcf.gz ${gzvcf}
    echo "filtered vcf created"
    # index the filtered VCF
    bcftools index --threads 24 ${dir}/pop_vcfs/${inv}_${population}.vcf.gz
    echo "vcf indexed"
    # create files for genotype lists
    G0="${dir}/${inv}_${population}_G0.txt"
    G2="${dir}/${inv}_${population}_G2.txt"

    # read in the sample genotype data for all inversions
    genotypes="${dir}/3_PCA_het_all_inversions_final.csv"
    awk -F',' '{print $1","$2","$3","$4}' "$genotypes" > "${dir}/${inv}_${population}_genotypes.csv"
    genotypes="${dir}/${inv}_${population}_genotypes.csv"
    head "$genotypes"

    # For the current inversion add each homozygous sample to the correct list, based on their genotype
    while IFS="," read -r pop sample inversion genotype; do
        echo "Raw genotype: '$genotype'"
        if [ "$inversion" == "$inv" ]; then
            if [ "$genotype" == "0" ]; then
                echo -e "$sample" >> "$G0"
            elif [ "$genotype" == "2" ]; then
                echo -e "$sample" >> "$G2"
            else
                echo "Het sample: $sample"
            fi
        fi
    done < "$genotypes"

    # Calculate Weir and Cockerham’s Fst between the two homozygous genotypes
    vcftools \
        --gzvcf ${dir}/pop_vcfs/${inv}_${population}.vcf.gz \
        --weir-fst-pop "$G0" \
        --weir-fst-pop "$G2" \
        --fst-window-size 10000 \
        --fst-window-step 10000 \
        --out "${dir}/fst/${inv}_${population}"

    # Calculate LD for all samples in the inversion
    vcftools --gzvcf ${dir}/pop_vcfs/${inv}_${population}.vcf.gz --geno-r2 --maf 0.05 --thin 1000 --max-missing-count 0 --out ${dir}/ld/${inv}_${population}

    # subset the vcf file to just the samples from each homozygous genotype group
    bcftools view --threads 12 --samples-file ${G0} -O z -o ${dir}/ld/${inv}_${population}_G0.vcf.gz ${dir}/pop_vcfs/${inv}_${population}.vcf.gz
    bcftools index -t ${dir}/ld/${inv}_${population}_G0.vcf.gz

    bcftools view --threads 12 --samples-file ${G2} -O z -o ${dir}/ld/${inv}_${population}_G2.vcf.gz ${dir}/pop_vcfs/${inv}_${population}.vcf.gz
    bcftools index -t ${dir}/ld/${inv}_${population}_G2.vcf.gz

    # Calculate LD for each homozygous genotype
    vcftools --gzvcf ${dir}/ld/${inv}_${population}_G0.vcf.gz --geno-r2 --maf 0.05 --thin 1000 --max-missing-count 0 --out ${dir}/ld/${inv}_${population}_G0
    cat ${dir}/ld/${inv}_${population}_G0.geno.ld | perl /home/uqkmcla4/scripts/reformat/emerald2windowldcounts.pl > ${dir}/ld/${inv}_${population}_G0_ld.txt

    vcftools --gzvcf ${dir}/ld/${inv}_${population}_G2.vcf.gz --geno-r2 --maf 0.05 --thin 1000 --max-missing-count 0 --out ${dir}/ld/${inv}_${population}_G2
    cat ${dir}/ld/${inv}_${population}_G2.geno.ld | perl /home/uqkmcla4/scripts/reformat/emerald2windowldcounts.pl > ${dir}/ld/${inv}_${population}_G2_ld.txt

done < "$inversions"
