### Produce a table of HiC interactions per sample
### GNU parallel by sample

# set variables 
hic="/g/data/ht96/McLay_UQ/Honours/Raw_data/Hi_C"
dir="/g/data/ht96/McLay_UQ/inversion_paper/hic"
ref="/scratch/ht96/km6006/2_D01_30.fasta"

# trim, rename and zip files 
/home/564/km6006/homer/bin/homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 -q 15 ${hic}/${1}_HiC_${2}.fastq.gz

mv ${hic}/${1}_HiC_${2}.fastq.gz.trimmed ${dir}/1_${1}_HiC_${2}_trimmed.fastq
/home/564/km6006/bin/pigz/pigz ${dir}/1_${1}_HiC_${2}_trimmed.fastq

# align to reference genome 
ngm -r ${ref} -t 24 -q ${dir}/1_${1}_HiC_${2}_trimmed.fastq.gz -o ${dir}/2_${1}_HiC_${2}_aligned.sam; 

# make tag directory
/home/564/km6006/homer/bin/makeTagDirectory ${dir}/3_${1}_HiC-TAG ${dir}/2_${1}_HiC_R1_aligned.sam,${dir}/2_${1}_HiC_R2_aligned.sam -tbp 1 -mapq 10’ 

# array of chroms 
declare -a arr=("scaffold_1" "scaffold_2" "scaffold_3" "scaffold_4" "scaffold_5" "scaffold_6" "scaffold_7" "scaffold_8" "scaffold_9" "scaffold_10" "scaffold_11" "scaffold_12" "scaffold_13" "scaffold_14" "scaffold_15" "scaffold_16" "scaffold_17" "scaffold_18" "scaffold_19" "scaffold_20" "scaffold_21" "scaffold_22" "scaffold_23" "scaffold_24" "scaffold_25" "scaffold_26" "scaffold_27" "scaffold_28" "scaffold_29" "scaffold_30")

# create Hi-C interaction matrix and convert it to a table 
for i in "${arr[@]}";
do
    if [ ${#i} -eq 10 ]; then
        modified_i="${i:0:9}0${i:9:1}"
    else
        modified_i="${i}"
    fi
    chr="${modified_i: -2}"
    /home/564/km6006/homer/bin/analyzeHiC ${dir}/${1}_HiC-TAG -chr ${i} -res 10000 -coverageNorm > ${dir}/int_matrix/${1}_${i}_hic_matrix.txt
    cat ${dir}/${1}_${i}_hic_matrix.txt | perl /home/564/km6006/Scripts/inversion_paper/hic/hicmatrix2table.pl ${1} > ${dir}/${1}_${i}_hic_table.txt
done

wait