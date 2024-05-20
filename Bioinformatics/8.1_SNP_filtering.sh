### SNP filtering part 2

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio
module load bcftools 
module load gatk

# Create directory variables 
DIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering"
LDIR="/home/564/km6006/Scripts"

# Create file naming variables
cov="HC"
ds="VO"

# Merge the chromosome files into a single file
file_list="${DIR}/${cov}${ds}_file_list.txt"
sort_list="${DIR}/${cov}${ds}_sorted_file_list.txt"

find "$DIR" -type f -name "*${cov}${ds}*sf2.vcf.gz" >> "$file_list"

if [ -s "$file_list" ]; then
    awk '{print $1}' "$file_list" | sort -t'_' -k1 -n > "$sort_list" 
    # Concatenate and index 
    bcftools concat --threads 6 --file-list "$sort_list" --output ${DIR}/${cov}${ds}_sf3.vcf.gz
    bcftools index --threads 6 -t ${DIR}/${cov}${ds}_sf3.vcf.gz

else
    echo "No match" 
fi 

wait

# Visualise the distribution of mean depth per site 
/home/564/km6006/bin/vcftools \
    --gzvcf ${DIR}/${cov}${ds}_sf3.vcf.gz \
    --site-mean-depth \
    --out ${DIR}/${cov}${ds}_sf3

# Filter for max-meanDP and min-meanDP
/home/564/km6006/bin/vcftools \
    --gzvcf ${DIR}/${cov}${ds}_sf3.vcf.gz \
    --max-meanDP 60 \
    --min-meanDP 3 \
    --recode --recode-INFO-all --stdout | gzip -c > ${DIR}/${cov}${ds}_sf4.vcf.gz

# Generate list of missing data per sample
/home/564/km6006/bin/vcftools \
    --gzvcf ${DIR}/${cov}${ds}_sf4.vcf.gz \
    --missing-indv \
    --out ${DIR}/${cov}${ds}_sf4_missing_inds

# Filter out individuals with too much missing data
awk '$5 > 0.4' ${DIR}/${cov}${ds}_sf4_missing_inds.imiss | cut -f1 > ${DIR}/${cov}${ds}_sf4_lowDP.indv

/home/564/km6006/bin/vcftools \
     --gzvcf ${DIR}/${cov}${ds}_sf4.vcf.gz \
     --remove ${DIR}/${cov}${ds}_sf4_lowDP.indv \
     --recode --recode-INFO-all --out ${DIR}/${cov}${ds}_sf5

bgzip -@ 12 ${DIR}/HCVO_sf5.recode.vcf
bcftools index --threads 12 -t ${DIR}/HCVO_sf5.recode.vcf.gz

### merge LC & HC 
# Merge the chromosome files into a single file
merge_list="${DIR}/file_list_LCHC.txt"

find "$DIR" -type f -name "*sf5.vcf.gz" >> "$merge_list"

if [ -s "$merge_list" ]; then
    awk '{print $1}' "$merge_list" 
    # Concatenate and index 
    bcftools merge --threads 6 --file-list "$merge_list" --output ${DIR}/sf6.vcf.gz
    bcftools index --threads 6 -t ${DIR}/sf6.vcf.gz

else
    echo "No match" 
fi 

wait

# Check missing site data - per population 
populations="/home/564/km6006/Scripts/inversion_paper/populations.txt"

while IFS= read -r value; do
    /home/564/km6006/bin/vcftools \
        --gzvcf ${DIR}/sf6.vcf.gz \
        --keep ${LDIR}/${value}.txt \
        --missing-site \
        --out ${DIR}/sf6_${value}
done < "$populations"

wait

# Concatenate the population lists (omitting any populations that have < 5 samples)
touch "${DIR}/sf6_popmiss.txt"

for file in "${DIR}"/*.lmiss; do
    if [[ $file == *A03* || $file == *A07* || $file == *A10* || $file == *A11* || $file == *A14* || $file == *D32* || $file == *H15* ]]; then
        continue  
    else
        cat "$file" >> "${DIR}/sf6_popmiss.txt"
    fi
done

# Filter the list for those sites where the missing data for any population is > 50%  
awk '!/CHR/ && $6 > 0.5 {print $1, $2}' ${DIR}/sf6_popmiss.txt > ${DIR}/sf6_badloci.txt

# Filter out low data sites for populations and filter to biallelic sites only 
/home/564/km6006/bin/vcftools \
    --gzvcf ${DIR}/sf6.vcf.gz \
    --exclude-positions ${DIR}/sf6_badloci.txt \
    --max-missing 0.8 \
    --max-alleles 2 \
    --recode --recode-INFO-all --out ${DIR}/sf7

bgzip -@ 12 ${DIR}/sf7.recode.vcf
bcftools index --threads 12 -t ${DIR}/sf7.recode.vcf.gz
wait

#Filter invariant sites 
gatk SelectVariants \
    -R /scratch/ht96/km6006/2_D01_30.fasta \
    -V ${DIR}/sf7.recode.vcf.gz \
    --exclude-non-variants \
    -O ${DIR}/sf8_final.vcf.gz

# calculate missing data per sample for the final dataset 
/home/564/km6006/bin/vcftools \
    --gzvcf ${DIR}/sf8_final.vcf.gz \
    --missing-indv \
    --out ${DIR}/sf8_final

# create bcf file 
bcftools convert -O b ${DIR}/sf8_final.vcf.gz > ${DIR}/sf8_final.bcf
bcftools index -f ${DIR}/sf8_final.bcf

# remove singletons - singletons also removed from vcf
/home/564/km6006/bin/vcftools \
    --bcf ${DIR}/sf8_no_variants_noD1_reheader_final.bcf \
    --mac 2 \
    --out ${DIR}/sf8_no_variants_noD1_reheader_final_nsng