### SNP Filtering for data with only variant sites 

module load gatk

# create directory variables 
ODIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only"
LDIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/population_missing_lists"

#create file naming variables
cov="LC"
ds="VO"
input_file="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/LC/merge.subset_vcf${1}.fofn.vcf.gz"

# Remove invariant sites
gatk SelectVariants \
    -R /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
    -V ${input_file} \
    --exclude-non-variants \
    -O ${ODIR}/${1}_${cov}_${ds}_all_pops.vcf.gz

# Extract annotation data for QD, MQ, FS, SOR, MQRankSum,ReadPosRankSum
zcat ${ODIR}/${1}_LC_VO_all_pops.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^QD=" | awk -F "=" '{print $2}' > ${ODIR}/ann_tables/${1}_LC_VO_all_pops_QD${1}.txt
zcat ${ODIR}/${1}_LC_VO_all_pops.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^MQ=" | awk -F "=" '{print $2}' > ${ODIR}/ann_tables/${1}_LC_VO_all_pops_MQ${1}.txt
zcat ${ODIR}/${1}_LC_VO_all_pops.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^FS=" | awk -F "=" '{print $2}' > ${ODIR}/ann_tables/${1}_LC_VO_all_pops_FS${1}.txt
zcat ${ODIR}/${1}_LC_VO_all_pops.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^SOR=" | awk -F "=" '{print $2}' > ${ODIR}/ann_tables/${1}_LC_VO_all_pops_SOR${1}.txt
zcat ${ODIR}/${1}_LC_VO_all_pops.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^MQRankSum=" | awk -F "=" '{print $2}' > ${ODIR}/ann_tables/${1}_LC_VO_all_pops_MQRankSum${1}.txt
zcat ${ODIR}/${1}_LC_VO_all_pops.vcf.gz | grep -v "^#" | awk '{print $8}' | sed 's/;/\n/g' | grep "^ReadPosRankSum=" | awk -F "=" '{print $2}' > ${ODIR}/ann_tables/${1}_LC_VO_all_pops_ReadPosRankSum${1}.txt

# combine the files for each annotation 
declare -a arr=("QD" "MQ" "FS" "SOR" "MQRankSum" "ReadPosRankSum")

all_files_exist=true

for i in "${arr[@]}"; do
    pattern="${ODIR}/*$i*.txt"
    files=($pattern)  
    if [ ${#files[@]} -ne 4 ]; then
        all_files_exist=false
        break
    fi
done

if [ "$all_files_exist" = true ]; then
    for i in "${arr[@]}"; do
        ann=$i
        for file in "${ODIR}"/*"$ann"*.txt; do
            if [ -f "$file" ]; then
                cat "$file" >> "${ODIR}/LC_VO_all_pops_${ann}.txt"
            fi
        done
    done
else
    echo "Skipping"
fi

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
    -O ${ODIR}/${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    > /scratch/ht96/km6006/prog.out.${PBS_JOBID} 2> /scratch/ht96/km6006/prog.err.${PBS_JOBID}

# optional - check VCF format 
gatk ValidateVariants \
    -R /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
    -V ${ODIR}/${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    --validation-type-to-exclude ALL

# Remove the variants that didn't pass filtering 
gatk SelectVariants \
    -R /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
    -V ${ODIR}/${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    --exclude-filtered \
    -O ${ODIR}/${cov}_${ds}_all_pops_snpfil_2.vcf.gz

# Read depth per individuals 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_2.vcf.gz \
    --minDP 1.5 \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_all_pops_snpfil_3.vcf.gz

# Missing data per site - filter for both 50 and 95 and see the difference in no. SNPs kept 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_3.vcf.gz \
    --max-missing 0.95 \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_all_pops_snpfil_4.vcf.gz

# Check missing data per sample 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_4.vcf.gz \
    --missing-indv
    --out ${ODIR}/${cov}_${ds}_missing_inds_95

# Filter out individuals with too much missing data
mawk '$5 > 0.5' out.imiss | cut -f1 > lowDP.indv
/home/564/km6006/bin/vcftools \
     --gvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_4.vcf.gz \
     --remove lowDP.indv \
     --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_all_pops_snpfil_5.vcf.gz

# # Check missing data per site - per population  
/home/564/km6006/bin/vcftools \
    --gvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_5.vcf.gz \
    --keep ${LDIR}/${1}.pop \
    --missing-site \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_${1}_snpfil_5_popmiss

# # Concatenate the population lists
for file in ${LDIR}: do 
    if [ -f "$file" ]; then
        concatenated_content+=$(cat ${file})
    fi
done

echo -e "$concatenated_content" > ${LDIR}/concat_pop_list.txt

# Filter the list for those sites where the missing data for any population is > 10%  
${DIR}/concat_pop_list.txt | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> badloci

# Filter low data sites for populations 
/home/564/km6006/bin/vcftools \
    --gvcf ${ODIR}/${cov}_${ds}_pops_snpfil_5.vcf.gz \
    --exclude-positions badloci \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_pops_snpfil_6.vcf.gz
 
# Visualise the distribution of mean depth per site 
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_6.vcf.gz \
    --site-mean-depth \
    --out ${ODIR}/${cov}_${ds}_mean_depth

# Max-meanDP
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_6.vcf.gz \
    --max-meanDP 60 \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${cov}_${ds}_all_pops_snpfil_7.vcf.gz

# Min-meanDP
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${cov}_${ds}_all_pops_snpfil_7.vcf.gz \
    --min-meanDP 3 \
    --recode --recode-INFO-all --out ${ODIR}/${cov}_${ds}_all_pops_snpfil_8_final