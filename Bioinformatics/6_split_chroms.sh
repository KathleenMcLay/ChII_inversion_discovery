### Subset each population vcf files by chromosome to allow parallel SNP filtering
### GNU parallel by population

module load bcftools 

DIR="/g/data/ht96/McLay_UQ/inversion_paper/joint"
ODIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/pop_chrom"

#subset the population vcf files by chromosome
declare -a arr=("scaffold_1" "scaffold_2" "scaffold_3" "scaffold_4" "scaffold_5" "scaffold_6" "scaffold_7" "scaffold_8" "scaffold_9" "scaffold_10" "scaffold_11" "scaffold_12" "scaffold_13" "scaffold_14" "scaffold_15" "scaffold_16" "scaffold_17" "scaffold_18" "scaffold_19" "scaffold_20" "scaffold_21" "scaffold_22" "scaffold_23" "scaffold_24" "scaffold_25" "scaffold_26" "scaffold_27" "scaffold_28" "scaffold_29" "scaffold_30")

for i in "${arr[@]}";
do
    if [ ${#i} -eq 10 ]; then
        modified_i="${i:0:9}0${i:9:1}"
    else
        modified_i="${i}"
    fi
    chr="${modified_i: -2}"
    bcftools view --threads 6 -O z --regions $i -o ${ODIR}/${chr}_${1}.vcf.gz ${DIR}/${1}_jntcl.vcf.gz &
done 

wait