#align long read assembly to reference with minimap2, plot alignment

module load samtools

ref_file="/QRISdata/Q6656/reference_genome/2_D01_30.fasta"
query_file="/QRISdata/Q6656/long_read_data/AGRF_top_30_scaffold_assemblies/1_H01_30.fasta"
outdir="/QRISdata/Q6656/long_read_data"

# align long read assembly to reference with minimap2
/home/uqkmcla4/minimap2 -ax map-hifi ${ref_file} ${query_file} > ${outdir}/H01_hifi_wga.sam

# sort and index the alignment
samtools view -bS ${outdir}/H01_hifi_wga.sam > ${outdir}/H01_hifi_wga.bam
samtools sort ${outdir}/H01_hifi_wga.bam -o ${outdir}/H01_hifi_wga.sorted.bam
samtools index ${outdir}/H01_hifi_wga.sorted.bam

dir="/home/uqkmcla4/ChII_inversion_discovery/inversion_breakpoints"

# plot the alignment - ref seq and query seq name is the scaffold 
python ${dir}/plot_bam_alignment.py ${outdir}/H01_hifi_wga.sorted.bam ${outdir}/H01_hifi_wga scaffold_1 1-271379410 contig_1 ${REVERSE}

