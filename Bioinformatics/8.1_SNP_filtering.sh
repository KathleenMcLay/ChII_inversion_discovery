### SNP filtering part 2

module load bcftools 

# Create directory variables 
ODIR="6_SNP_filtering/variant_only"
LDIR="6_SNP_filtering/population_missing_lists"

# Create file naming variables
cov="LC"
ds="VO"

# Merge the chromosome files into a single file
file_list="${ODIR}/file_list.txt"
sort_list="${ODIR}/sorted_file_list.txt"

find "$ODIR" -type f -name "*LC*_2.vcf.gz" >> "$file_list"

if [ -s "$file_list" ]; then
    awk '{print $1}' "$file_list" | sort -t'_' -k1 -n > "$sort_list" 
    # Concatenate and index 
    bcftools concat --threads 6 --file-list "$sort_list" --output ${ODIR}/${cov}${ds}_allpops_snpfil_3.vcf.gz
    bcftools index --threads 6 -t ${ODIR}/${cov}${ds}_allpops_snpfil_3.vcf.gz

else
    echo "No match" 
fi 

# Generate missing data per sample
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}${ds}_all_pops_snpfil_5.vcf.gz \
    --missing-indv \
    --out ${ODIR}/${cov}${ds}_missing_inds

# Filter out individuals with too much missing data
mawk '$5 > 0.5' ${ODIR}/${cov}${ds}_missing_inds.imiss | cut -f1 > ${ODIR}/lowDP.indv
/home/564/km6006/bin/vcftools \
     --gvcf ${ODIR}/${cov}${ds}_all_pops_snpfil_3.vcf.gz \
     --remove ${ODIR}/lowDP.indv \
     --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}${ds}_all_pops_snpfil_4.vcf.gz

# Check missing site data - per population 
for value in ${LDIR}/populations.txt: do
    /home/564/km6006/bin/vcftools \
        --gvcf ${ODIR}/${cov}${ds}_all_pops_snpfil_4.vcf.gz \
        --keep ${LDIR}/${value}.txt \
        --missing-site \
        --out ${ODIR}/${cov}${ds}_${value}_snpfil_5_popmiss
done

# Concatenate the population lists
for file in ${ODIR}/*.lmiss; do
    if [ -f "$file" ]; then
        cat "$file" >> ${ODIR}/concat_LCVO_snpfil_6.txt
    fi
done

# Filter the list for those sites where the missing data for any population is > 40%  
awk '!/CHR/ && $6 > 0.4 {print $1, $2}' ${ODIR}/concat_LCVO_snpfil_6.txt > ${ODIR}/badloci.txt

# Filter out low data sites for populations 
/home/564/km6006/bin/vcftools \
    --gvcf ${ODIR}/${cov}${ds}_pops_snpfil_3.vcf.gz \
    --exclude-positions ${ODIR}/badloci.txt \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}${ds}_pops_snpfil_4.vcf.gz
 
# Visualise the distribution of mean depth per site 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}${ds}_all_pops_snpfil_4.vcf.gz \
    --site-mean-depth \
    --out ${ODIR}/${cov}${ds}_mean_depth

# Filter for max-meanDP (20 for LC, 60 for HC)
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}${ds}_all_pops_snpfil_4.vcf.gz \
    --max-meanDP 60 \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}${ds}_all_pops_snpfil_5.vcf.gz

# Filter for min-meanDP
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}${ds}_all_pops_snpfil_5.vcf.gz \
    --min-meanDP 3 \
    --recode --recode-INFO-all --out ${ODIR}/${cov}${ds}_all_pops_snpfil_6_final