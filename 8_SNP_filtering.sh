### SNP Filtering for data with only variant sites 

module load gatk

# create directory variables 
ODIR="6_SNP_filtering/variant_only"
LDIR="6_SNP_filtering/population_missing_lists"

#create file naming variables
cov="LC"
ds="VO"
input_file="6_SNP_filtering/${1}_LC_VO_all_pops.vcf.gz"

# Remove invariant sites
gatk SelectVariants \
    -R /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
    -V ${input_file} \
    --exclude-non-variants \
    -O ${ODIR}/${1}_${cov}_${ds}_all_pops.vcf.gz

# Gatk filtering - adds PASS to the filter field, otherwise, if failed adds the name/s of the failed filter
gatk VariantFiltration \
    -V ${ODIR}/${cov}_${ds}_all_pops.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    > /scratch/ht96/km6006/prog.out.${PBS_JOBID} 2> /scratch/ht96/km6006/prog.err.${PBS_JOBID}

# Read depth per individuals. Note: the --remove-filtered-all option removes all sites that have been filtered out by the previous step
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    --minDP 3 \
    --max-missing 0.8 \
    --remove-filtered-all \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_2.vcf.gz

# Check missing data per sample 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_2.vcf.gz \
    --missing-indv
    --out ${ODIR}/${1}_${cov}_${ds}_missing_inds

# Filter out individuals with too much missing data
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
/home/564/km6006/bin/vcftools \
     --gvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_2.vcf.gz \
     --remove lowDP.indv \
     --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_3.vcf.gz

# Check missing data per site - per population  
/home/564/km6006/bin/vcftools \
    --gvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_3.vcf.gz \
    --keep ${LDIR}/${1}.pop \
    --missing-site \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_${1}_snpfil_3_popmiss

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