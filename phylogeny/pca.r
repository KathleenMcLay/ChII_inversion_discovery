###Create a 'neutral' PCA for the dataset
library(pals)
library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggrepel)
library(ggplot2)

samples <- read.table("/home/564/km6006/Scripts/inversion_paper/phylogeny/PCA_populationlist.tsv", header=TRUE)
vcf_fn <- ("/g/data/ht96/McLay_UQ/inversion_paper/phy/sf7_pruned_10kb.vcf")

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

# Plot the PCA 
tiff("/g/data/ht96/McLay_UQ/inversion_paper/phy/sf7_pruned_PCA_10kb.tiff", units="in", width=9, height=7, res=500)
ggplot(PCA, aes(x=PC1, y=PC2, col=population)) + geom_point(size = 2) + theme_bw() +
  scale_colour_manual(values=unname(polychrome()), name="Populations") + 
  #geom_label_repel(aes(label = sample), label.size = NA, fill = "NA", show.legend = FALSE, segment.color = "transparent") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  xlab(paste("PC1 (",pc1.percent, "%)")) +
  ylab(paste("PC1 (",pc2.percent, "%)")) +
  xlim(-0.08, -0.03) +
  ylim(-0.06, -0.03) +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"), 
        axis.title = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 20)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)), 
        text = element_text(family = "Helvetica") ) 
dev.off()
file.remove("/g/data/ht96/McLay_UQ/inversion_paper/phy/senecio.gds")