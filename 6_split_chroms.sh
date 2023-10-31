### Subset each population vcf files by chromosome to allow parallel SNP filtering
### GNU parallel by population

module load bcftools 

DIR="/g/data/ht96/McLay_UQ/inversion_paper/5_joint_calling"
ODIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering"

#subset the population vcf files by chromosome (35 pops total) 
declare -a arr=("SLv141Ch1" "SLv141Ch2" "SLv141Ch3" "SLv141Ch4" "SLv141Ch5" "SLv141Ch6" "SLv141Ch7" "SLv141Ch8" "SLv141Ch9" "SLv141Ch10" "SLv141Ch11" "SLv141Ch12" "SLv141Ch13" "SLv141Ch14" "SLv141Ch15" "SLv141Ch16" "SLv141Ch17" "SLv141Ch18" "SLv141Ch19" "SLv141Ch20")

for i in "${arr[@]}";
do
    if [ ${#i} -eq 9 ]; then
        modified_i="${i:0:8}0${i:8:1}"
    else
        modified_i="${i}"
    fi
    chr="${modified_i: -2}"
    bcftools view --threads 6 -O z --regions $i -o ${DIR}/${chr}_${1}_cohort.vcf.gz ${DIR}/${1}_cohort.vcf.gz &
done 

wait

# index the chromosome files 
for file in ${ODIR}/*.vcf.gz; 
do
    bcftools index --threads 12 -t ${file} &
done

wait