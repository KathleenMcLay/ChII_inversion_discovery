### Create PCAs for SNP filtered dataset 

library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggrepel)
library(ggplot2)

samples <- read_tsv("/home/564/km6006/Scripts/inversion_paper/populationlist.tsv")
vcf_fn <- ("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/sf8_final.vcf.gz")

# # Convert VCF file 
snpgdsVCF2GDS(vcf_fn, "/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/senecio.gds", method="biallelic.only") #, compress.annotation="LZMA_RA", compress.gen
# Assign the new file to an object 
genofile <- snpgdsOpen("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/senecio.gds")
# Run PCA 
senecio_pca <- snpgdsPCA(genofile, autosome.only=FALSE)

# # Extract PCs and make dataset with samples 
PC1 <- senecio_pca$eigenvect[,1]
PC2 <- senecio_pca$eigenvect[,2]
PCA <- cbind(PC1, PC2)
PCA <- cbind(samples, PCA)
write.table(PCA, "/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/sf8_final_PCA.txt", sep="\t")

# Plot the PCA 
tiff("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/sf8_final_PCA.tiff", units="in", width=9, height=7, res=300)
ggplot(PCA, aes(x=PC1, y=PC2, col=population)) + geom_point(size = 1) + theme_bw() +
  scale_colour_discrete(name="Populations") +
  #geom_label_repel(aes(label = sample), size=1, segment.size=0.25, fill = "NA", show.legend = FALSE, segment.color = "transparent") + 
  xlab("PC1") + ylab("PC2") +
  xlim(-0.05, 0.05) +
  ylim(-0.05, 0.05) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"), axis.title = element_text(size = 12),axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)), text = element_text(family = "Helvetica") ) 
dev.off()
file.remove("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/senecio.gds")