### SNP filtering part 1
### GNU parallel by chromosome 
### for dataset including non-variant sites exclude non-variant site, QD and QUAL filters 

module load gatk

# create directory variables 
ODIR="6_SNP_filtering/variant_only"
LDIR="6_SNP_filtering/population_missing_lists"

#create file naming variables
cov="LC"
ds="VO"
input_file="6_SNP_filtering/${1}_LC_VO_all_pops.vcf.gz"

# Remove invariant sites
gatk SelectVariants \
    -R /scratch/ht96/km6006/SLv141Asm_Ch20RN.fasta \
    -V ${input_file} \
    --exclude-non-variants \
    -O ${ODIR}/${1}_${cov}_${ds}_all_pops.vcf.gz

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
    -O ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    > /scratch/ht96/km6006/prog.out.${PBS_JOBID} 2> /scratch/ht96/km6006/prog.err.${PBS_JOBID}

# Read depth and max missing per site. Note: the --remove-filtered-all option removes all sites that have been filtered out by the previous step
/home/564/km6006/bin/vcftools \
    --gzvcf ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_1.vcf.gz \
    --minDP 3 \
    --max-missing 0.8 \
    --remove-filtered-all \
    --recode --recode-INFO-all --stdout | gzip -c > ${ODIR}/${1}_${cov}_${ds}_all_pops_snpfil_2.vcf.gz
