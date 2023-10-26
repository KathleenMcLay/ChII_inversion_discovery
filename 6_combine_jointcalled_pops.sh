### Combine joint called populations for SNP filtering 

module load bcftools 

DIR="/g/data/ht96/McLay_UQ/inversion_paper/5_joint_calling"
ODIR="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering"

# create combined population subsets for low coverage WGS data for parallel SNP filtering 

# number of populations per subset 
n=7
# create file of filenames (fofn) for each subset of files 
ls -d -1 ${DIR}/LC/*.vcf.gz | split -l ${n} --additional-suffix '.fofn' - subset_vcfs

# merge the files in parallel loop using bcftools merge 
for sub in *.fofn;
do
    out_vcf="${ODIR}/LC/merge.${sub}.vcf.gz"
    bcftools merge --threads 12 -O z --file-list ${sub} -o ${out_vcf} & 
done

# index the files using bcftools index
wait

for file in ${ODIR}/*.vcf.gz; 
do
    bcftools index --threads 12 -t ${file} &
done

wait
