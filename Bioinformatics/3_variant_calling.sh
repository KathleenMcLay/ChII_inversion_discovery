### Call variants with GATK HaplotypeCaller by chromosome and concatenate output files using bcftools concat 
### GNU parallel by sample 

module load gatk 
module load bcftools

# create directory variables 
CLDIR="/g/data/ht96/McLay_UQ/inversion_paper/3_clean"
SPDIR="/g/data/ht96/McLay_UQ/inversion_paper/4_variant_calling/split"
VCDIR="/g/data/ht96/McLay_UQ/inversion_paper/4_variant_calling"

# create an array of scaffold names to use as interval values with HaplotypeCaller
declare -a arr=("scaffold_1" "scaffold_2" "scaffold_3" "scaffold_4" "scaffold_5" "scaffold_6" "scaffold_7" "scaffold_8" "scaffold_9" "scaffold_10" "scaffold_11" "scaffold_12" "scaffold_13" "scaffold_14" "scaffold_15" "scaffold_16" "scaffold_17" "scaffold_18" "scaffold_19" "scaffold_20" "scaffold_21" "scaffold_22" "scaffold_23" "scaffold_24" "scaffold_25" "scaffold_26" "scaffold_27" "scaffold_28" "scaffold_29" "scaffold_30")

# run HaplotypeCaller 
for i in "${arr[@]}";
do
    # create a variable to isolate the scaffold no. for file naming  
    if [ ${#i} -eq 10 ]; then
        modified_i="${i:0:9}0${i:9:1}"
    else
        modified_i="${i}"
    fi
    chr="${modified_i: -2}"
    # run HaplotypeCaller per chromosome using -intervals
    gatk --java-options "-Xmx8g" HaplotypeCaller \
        --input ${CLDIR}/${1}_PCRm_cln_srt.bam \
        --output ${SPDIR}/${1}_${chr}_vrnt.g.vcf.gz \
        --reference /scratch/ht96/km6006/2_D01_30.fasta \
        --emit-ref-confidence GVCF \
        --intervals $i &
done

wait

# concatenate the scaffolds for each sample

#create list variables 
file_list="${VCDIR}/${1}_file_list.txt"
sort_list="${VCDIR}/${1}_sorted_file_list.txt"

# find all files for the current sample in the output directory 
find "$SPDIR" -type f -name "${1}_*.vcf.gz" >> "$file_list"

# print the file to the file_list variable, sort it in the sort list variable and then concatenate 
if [ -s "$file_list" ]; then
    awk '{print $1}' "$file_list" | sort -t'_' -k5 -n > "$sort_list" 
    # Concatenate and index 
    bcftools concat --threads 6 --file-list "$sort_list" --output "${VCDIR}/${1}_vrnt.g.vcf.gz" 
    bcftools index --threads 6 "${VCDIR}/${1}_vrnt.g.vcf.gz" 

else
    echo "No match for: $1" 
fi  
