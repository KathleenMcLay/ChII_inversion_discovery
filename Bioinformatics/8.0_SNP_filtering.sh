### SNP filtering part 1
### GNU parallel by chromosome 

module load gatk

#create file naming variables
cov="LC"
ds="VO"

# create directory variables 
DIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/"

input_file="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/${1}_chrom_${cov}allpops.vcf.gz"

#Remove indels
gatk SelectVariants \
    -R /scratch/ht96/km6006/2_D01_30.fasta \
    -V ${input_file} \
    --select-type-to-include SNP \
    -O ${DIR}/${1}_${cov}${ds}_sf0.vcf.gz

# Gatk filtering - adds PASS to the filter field, otherwise, if failed adds the name/s of the failed filter
gatk VariantFiltration \
    -V ${DIR}/${1}_${cov}${ds}_sf0.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${DIR}/${1}_${cov}${ds}_sf1.vcf.gz \
    > /scratch/ht96/km6006/prog.out.${PBS_JOBID} 2> /scratch/ht96/km6006/prog.err.${PBS_JOBID}

#max missing per site. Note: the --remove-filtered-all option removes all sites that have been filtered out by the previous step
/home/564/km6006/bin/vcftools \
    --gzvcf ${DIR}/${1}_${cov}${ds}_sf1.vcf.gz \
    --max-missing 0.8 \
    --remove-filtered-all \
    --recode --recode-INFO-all --stdout | gzip -c > ${DIR}/${1}_${cov}${ds}_sf2.vcf.gz
