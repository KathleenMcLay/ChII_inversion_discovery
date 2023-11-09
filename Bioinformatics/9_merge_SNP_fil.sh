
module load bcftools

DIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only"

# Zip the files and index
bgzip -@ 12 ${DIR}/HC/HC_VO_all_pops_snpfil_6_final.vcf
bcftools index --threads 12 -t ${DIR}/HC/HC_VO_all_pops_snpfil_6_final.vcf.gz

bgzip -@ 12 ${DIR}/LC/LCVO_allpops_snpfil_9_final.vcf
bcftools index --threads 12 -t ${DIR}/LC/LCVO_allpops_snpfil_9_final.vcf.gz


# Merge the seperately SNP filtered high coverage (HC) and low coverage (LC) data and filter sites that have > 20% missing data in the final dataset 
hcdata="${DIR}/HC/HCVO_allpops_snpfil_6_final.vcf.gz"
lcdata="${DIR}/LC/LCVO_allpops_snpfil_9_final.vcf.gz"

bcftools merge ${hcdata} ${lcdata} -o ${DIR}/allpops_VO.vcf

/home/564/km6006/bin/vcftools \
    --vcf ${DIR}/allpops_VO.vcf \
    --max-missing 0.8 \
    --recode --recode-INFO-all --out ${DIR}/allpops_VO_0.8_final.vcf

# Convert to bcf and index
bcftools convert -O b ${DIR}/allpops_VO_0.8_final.recode.vcf > ${DIR}/allpops_VO_0.8_final.bcf
bcftools index -f ${DIR}/sallpops_VO_0.8_final.bcf