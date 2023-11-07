module load bcftools 

# create directory variables 
ODIR="6_SNP_filtering/variant_only"
LDIR="6_SNP_filtering/population_missing_lists"

#create file naming variables
cov="LC"
ds="VO"

# Merge the chromosome files into a single file
file_list="${ODIR}/file_list.txt"
sort_list="${ODIR}/sorted_file_list.txt"

find "$ODIR" -type f -name "*LC*_2.vcf.gz" >> "$file_list"

if [ -s "$file_list" ]; then
    awk '{print $1}' "$file_list" | sort -t'_' -k1 -n > "$sort_list" 
    # Concatenate and index 
    bcftools concat --threads 6 --file-list "$sort_list" --output "${ODIR}/LCVO_allpops_snpfil_3.vcf.gz" 
    gatk IndexFeatureFile --input "${ODIR}/LCVO_allpops_snpfil_3.vcf.gz"  

else
    echo "No match" 
fi 

# Generate missing data per sample
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_5.vcf.gz \
    --missing-indv \
    --out ${ODIR}/${cov}_${ds}_missing_inds

# Filter out individuals with too much missing data
mawk '$5 > 0.5' ${ODIR}/${cov}_${ds}_missing_inds.imiss | cut -f1 > ${ODIR}/lowDP.indv
/home/564/km6006/bin/vcftools \
     --gvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_3.vcf.gz \
     --remove ${ODIR}/lowDP.indv \
     --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_all_pops_snpfil_4.vcf.gz

# Check missing data per site - per population  
/home/564/km6006/bin/vcftools \
    --gvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_4.vcf.gz \
    --keep ${LDIR}/${1}.pop \
    --missing-site \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_${1}_snpfil_5_popmiss

# Concatenate the population lists
for file in ${LDIR}: do 
    if [ -f "$file" ]; then
        concatenated_content+=$(cat ${file})
    fi
done

echo -e "$concatenated_content" > ${LDIR}/${1}_concat_pop_list.txt

# Filter the list for those sites where the missing data for any population is > 10%  
${LDIR}/${1}_concat_pop_list.txt | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci

# Filter out low data sites for populations 
/home/564/km6006/bin/vcftools \
    --gvcf ${ODIR}/${1}_${cov}_${ds}_pops_snpfil_3.vcf.gz \
    --exclude-positions badloci \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_pops_snpfil_4.vcf.gz
 
# Visualise the distribution of mean depth per site 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_4.vcf.gz \
    --site-mean-depth \
    --out ${ODIR}/${1}_${cov}_${ds}_mean_depth

# Filter for max-meanDP
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_4.vcf.gz \
    --max-meanDP 60 \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_5.vcf.gz

# Filter for min-meanDP
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_5.vcf.gz \
    --min-meanDP 3 \
    --recode --recode-INFO-all --out ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_6_final