#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH --job-name=zoe_pairs
#SBATCH --output=zoe_pairs_%j.out
#SBATCH --time=8:00:00
#SBATCH --partition=general
#SBATCH --account=a_ortiz_barrientos_coe

module load parallel 
module load bcftools/1.15.1-gcc-11.3.0
module load gatk
module load vcftools
module load htslib

# set directory paths
#JC_DIR="/QRISdata/Q6656/chapter_II/1_bioinformatics/joint_calling/para_pairs"
DIR="/scratch/user/uqkmcla4"

### Subset each population vcf files by chromosome to allow parallel SNP filtering
# declare -a arr=("scaffold_1" "scaffold_2" "scaffold_3" "scaffold_4" "scaffold_5" "scaffold_6" "scaffold_7" "scaffold_8" "scaffold_9" "scaffold_10" "scaffold_11" "scaffold_12" "scaffold_13" "scaffold_14" "scaffold_15" "scaffold_16" "scaffold_17" "scaffold_18" "scaffold_19" "scaffold_20" "scaffold_21" "scaffold_22" "scaffold_23" "scaffold_24" "scaffold_25" "scaffold_26" "scaffold_27" "scaffold_28" "scaffold_29" "scaffold_30")

# for i in "${arr[@]}";
# do
#     if [ ${#i} -eq 10 ]; then
#         modified_i="${i:0:9}0${i:9:1}"
#     else
#         modified_i="${i}"
#     fi
#     chr="${modified_i: -2}"
#     bcftools view --threads 2 -O z --regions $i -o ${DIR}/${chr}_${1}.vcf.gz ${JC_DIR}/${1}_jntcl.vcf.gz &
# done 

# wait

### Merge the the population files by chromosome, ready for SNP filtering 

# index the population chromosome files 
# for file in ${DIR}/*.vcf.gz; 
# do
#     bcftools index --threads 24 -t ${file} &
# done

# wait

# array of chrom numbers 
# chromosomes=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30")

# # for each chromosome create a file containing a list of population files for that chromosome
# for chrom in "${chromosomes[@]}"; do
#     # Find files with the current prefix and store their names in a text file
#     find "${DIR}" -maxdepth 1 -type f -name "${chrom}_*.vcf.gz" -exec readlink -f {} \; >> "${DIR}/${chrom}.txt"
# done

#!/bin/bash

# # Check if DIR is set
# if [[ -z "$DIR" ]]; then
#     echo "Error: DIR variable not set"
#     exit 1
# fi

# # Check if directory exists
# if [[ ! -d "$DIR" ]]; then
#     echo "Error: Directory $DIR does not exist"
#     exit 1
# fi

# # Process each chromosome
# for chrm in ${DIR}/04.txt; do
#     if [[ ! -f "$chrm" ]]; then
#         echo "Warning: No .txt files found in $DIR"
#         continue
#     fi
    
#     chrm_name=$(basename "$chrm" .txt)

#     # Skip files that start with 01, 02, or 03
#     if [[ "$chrm_name" =~ ^0[123] ]]; then
#     echo "Skipping $chrm_name (starts with 01, 02, or 03)"
#     continue
#     fi
    
#     out_vcf="${DIR}/${chrm_name}_allpops.vcf.gz"
    
#     echo "Processing $chrm_name..."
#     echo "Input file: $chrm"
#     echo "Output file: $out_vcf"
    
#     # Check if input file is not empty
#     if [[ ! -s "$chrm" ]]; then
#         echo "Warning: $chrm is empty, skipping..."
#         continue
#     fi
    
#     # Run bcftools merge with error checking
#     if bcftools merge --threads 12 -O z --file-list "$chrm" -o "$out_vcf"; then
#         echo "Successfully created $out_vcf"
#     else
#         echo "Error: Failed to create $out_vcf"
#         echo "Check the contents of $chrm and verify all listed files exist"
#     fi
# done

# wait

# index the chromosome files 
# for file in ${DIR}/*allpops.vcf.gz; 
# do
#     bcftools index --threads 24 -t ${file} &
# done

# wait

## SNP filtering part 1
### GNU parallel by chromosome 

#input_file="${DIR}/${1}_allpops.vcf.gz"

# Gatk filtering - adds PASS to the filter field, otherwise, if failed adds the name/s of the failed filter
# gatk VariantFiltration \
#     -V ${input_file} \
#     -filter "QD < 2.0" --filter-name "QD2" \
#     -filter "QUAL < 30.0" --filter-name "QUAL30" \
#     -filter "SOR > 3.0" --filter-name "SOR3" \
#     -filter "FS > 60.0" --filter-name "FS60" \
#     -filter "MQ < 40.0" --filter-name "MQ40" \
#     -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#     -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#     -O ${DIR}/${1}_1_VF.vcf.gz \
#     > /scratch/user/uqkmcla4/prog.out.${PBS_JOBID} 2> /scratch/user/uqkmcla4/prog.err.${PBS_JOBID}

#max missing per site. Note: the --remove-filtered-all option removes all sites that have been filtered out by the previous step
# vcftools \
#     --gzvcf ${DIR}/${1}_1_VF.vcf.gz \
#     --max-missing 0.8 \
#     --remove-filtered-all \
#     --max-meanDP 60 \
#     --min-meanDP 2 \
#     --remove-indv D01208 \
#     --remove-indv D01209 \
#     --remove-indv D01288 \
#     --recode --recode-INFO-all --stdout | bgzip -c > ${DIR}/${1}_2_VF_0.8_meanDP.vcf.gz

# index the chromosome files 
# for file in ${DIR}/*_2_VF_0.8_meanDP.vcf.gz; 
# do
#     echo ${file}
#     bcftools index --threads 24 -t ${file} &
# done

# index the chromosome files 
# for file in ${DIR}/*_2_VF_0.8_meanDP.vcf.gz; do
#     echo "Processing ${file}"
#     bcftools view --threads 24 -S ${DIR}/sample_order.txt --force-samples \
#         -Oz -o ${file%.vcf.gz}_reordered.vcf.gz \
#         $file
# done

# index the chromosome files
# for file in ${DIR}/*_2_VF_0.8_meanDP_reordered.vcf.gz; do
#     echo "Processing ${file}"
#     bcftools index --threads 24 -t ${file} &
# done    

#Merge the chromosome files into a single file
# file_list="${DIR}/file_list.txt"
# sort_list="${DIR}/sorted_file_list.cleaned.txt"

# find "$DIR" -type f -name "*_2_VF_0.8_meanDP_reordered.vcf.gz" >> "$file_list"

# if [ -s "$sort_list" ]; then
#     #awk '{print $1}' "$file_list" | sort -t'_' -k1 -n > "$sort_list" 
#     # Concatenate and index 
#     bcftools concat --threads 24 --file-list "$sort_list" --output ${DIR}/3_VF_0.8_meanDP_reordered_allscaff.vcf.gz
#     bcftools index --threads 24 -t ${DIR}/3_VF_0.8_meanDP_reordered_allscaff.vcf.gz

# else
#     echo "No match" 
# fi 

# wait

#bcftools index --threads 12 -t ${DIR}/3_VF_0.8_meanDP_reordered_allscaff.vcf.gz

# # Filter to biallelic sites only 
vcftools \
--gzvcf ${DIR}/3_VF_0.8_meanDP_reordered_allscaff.vcf.gz \
--min-alleles 1 \
--max-alleles 2 \
--recode --recode-INFO-all --out ${DIR}/4_VF_0.8_meanDP_reordered_allscaff_biallelic
