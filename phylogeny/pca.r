###Create a 'neutral' PCA for the dataset
library(pals)
library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggrepel)
library(ggplot2)

samples <- read.table("/home/564/km6006/Scripts/inversion_paper/phylogeny/PCA_populationlist.tsv", header=TRUE)
vcf_fn <- ("/g/data/ht96/McLay_UQ/inversion_paper/phy/sf7_pruned_10kb_noD1.vcf")

# Convert VCF file 
snpgdsVCF2GDS(vcf_fn, "/g/data/ht96/McLay_UQ/inversion_paper/phy/senecio.gds") 
# Assign the new file to an object 
genofile <- snpgdsOpen("/g/data/ht96/McLay_UQ/inversion_paper/phy/senecio.gds")
# Run PCA 
senecio_pca <- snpgdsPCA(genofile, autosome.only=FALSE)

# get PVE 
pc1.percent <- head(round(senecio_pca$varprop[1]*100, 2))
pc2.percent <- head(round(senecio_pca$varprop[2]*100, 2))

# Extract PCs and make dataset with samples 
PC1 <- senecio_pca$eigenvect[,1]
PC2 <- senecio_pca$eigenvect[,2]
PCA <- cbind(PC1, PC2)
colnames(PCA) <- c("PC1", "PC2")
PCA <- cbind(samples, PCA)
write.table(PCA, "/g/data/ht96/McLay_UQ/inversion_paper/phy/sf7_pruned_PCA_10kb_noD1.txt", sep="\t")

#PCA <- read.table("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /2_phylogeny_pop_structure/data/sf7_pruned_PCA_10kb_noD1.txt")
#ecotypes <- unique(PCA$ecotype)

pca_theme <- theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
                 axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
                 axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
                 plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
                 axis.ticks = element_line(colour="black"), panel.background=element_rect(colour="black"), panel.border=element_rect(colour="black", size=0.75),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.title = element_text(size=8, face="bold"), legend.text=element_text(size=6, face="bold"), legend.spacing.y = unit(0.1, "cm")
                 ) 

tiff("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /2_phylogeny_pop_structure/sf7_pruned_PCA_10kb_noD1.tiff", units="in", width=6, height=4.5, res=500)
ggplot(PCA, aes(x=PC1, y=PC2, col=ecotype)) + geom_point(size = 2) + theme_bw() + 
  scale_colour_manual(values=c("#333399", "#66CCFF", "#009933", "#FF9933", "#6699CC", "#663300"), name="Ecotype") + 
  guides(color = guide_legend(override.aes = list(size = 5))) +
  xlab(paste("PC1 11.72%")) +
  ylab(paste("PC2 4.62%")) +
  pca_theme +
  guides(color = guide_legend(override.aes = list(size = 2)))
dev.off()

file.remove("/g/data/ht96/McLay_UQ/inversion_paper/phy/senecio.gds")