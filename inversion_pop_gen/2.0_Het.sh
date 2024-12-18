### Calculate HET for each individual for all SNPs in the putative inversion

module load vcftools 

directory="/QRISdata/Q6656/new_lpca/inversion_vcfs"

for file in $directory/*.vcf.gz; do 
    filename=$(basename -- "$file")
    filename="${filename%.*}"
    vcftools --gzvcf "$file" --het --out /QRISdata/Q6656/new_lpca/het_out/"$filename"
    --vcf ${output} \
    --het \
    --out /QRISdata/Q6684/working_data/het/${1}