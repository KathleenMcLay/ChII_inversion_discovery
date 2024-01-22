### Create Phylogeny 

module load bcftools

INFILE="/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/sf7.vcf.gz"
PR_OUT="/g/data/ht96/McLay_UQ/inversion_paper/phy/sf7_pruned_10kb.vcf"
PHY_DIR="/g/data/ht96/McLay_UQ/inversion_paper/phy"
RAX_IN="/g/data/ht96/McLay_UQ/inversion_paper/phy/sf7_pruned_10kb.phy"

# Filter to 1 SNP per 100kb 
bcftools +prune -w 10kb --nsites-per-win 1 ${INFILE} -O v -o ${PR_OUT}

# convert to phy file 
/home/564/km6006/miniconda3/bin/python /home/564/km6006/vcf2phylip/vcf2phylip.py --input ${PR_OUT} --output-folder ${PHY_DIR}

# RAxML phylogeny 
/home/564/km6006/bin/raxmlHPC-PTHREADS-AVX2 -T 46 -f a -d -m GTRGAMMA -p 45687 -N 1000 -e 0.01 -x 56043 -B 0.03 --bootstop-perms 1000 -n ML -w ${PHY_DIR} -s ${RAX_IN}