ref_file="/QRISdata/Q6656/long_read_data/2_D01_30.fasta"
query_file="/QRISdata/Q6656/long_read_data/3_denovo_contigs/H01/assembly.fasta"
outdir="/scratch/user/uqkmcla4/breakpoints"
Q=$(basename ${query_file} .fasta)

# align denovo contigs to reference with minimap2
/home/uqkmcla4/minimap2 -ax map-hifi ${ref_file} ${query_file} > ${outdir}/${Q}_H01_flye_aln.sam