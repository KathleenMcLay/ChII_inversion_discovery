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

# make a matrix of interactions and convert it to a table
/home/564/km6006/homer/bin/analyzeHiC ${dir}/3_${1}_HiC-TAG -res 10000 -coverageNorm > ${dir}/4_${1}_hic_matrix.txt
cat ${dir}/4_${1}_hic_matrix.txt | perl /home/564/km6006/Scripts/inversion_paper/hic/hicmatrix2table.pl ${1} > ${dir}/4_${1}_hic_table.txt 