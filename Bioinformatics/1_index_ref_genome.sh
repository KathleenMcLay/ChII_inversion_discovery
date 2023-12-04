### Index the reference genome for use with downstream programs

module load samtools
module load bwa

#index the reference genome using 'bwa index'
bwa index -a bwtsw /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta

#index the reference genome with samtools for use with GATK
samtools faidx /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta

#creating a dictionary file of the reference genome for GATK 
java -jar /home/564/km6006/picard.jar CreateSequenceDictionary R=/scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta O=/scratch/ht96/km6006/SLv141Asm_Ch20RN.dict