### Calculate HET for each individual for all SNPs in the putative inversion

module load vcftools 

directory="/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/inv_vcfs"

for file in $directory/*.vcf.gz; do 
    filename=$(basename -- "$file")
    filename="${filename%.*}"
    vcftools --gzvcf "$file" --het --out /QRISdata/Q6656/chapter_II/new_inv_discovery/alt/het/"$filename" 
done