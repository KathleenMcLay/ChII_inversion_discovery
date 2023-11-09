### Create PCAs for each SNP filtered data set to compare

library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggrepel)
library(doParallel)
library(stringr)
library(ggplot2)

files <- list.files(path = "/home/564/km6006/Scripts/inversion_paper/", pattern="recode", all.files=FALSE, full.names = FALSE)

foreach(i=1:length(files)) %dopar% {
  samples <- read_tsv(paste("/home/564/km6006/Scripts/inversion_paper/populationlist_", str_split(files[i], "_", n = 2, simplify=TRUE), ".tsv", sep = ""))
  vcf_fn <- (paste("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only/", (files[i])))

  # Convert VCF file 
  snpgdsVCF2GDS(vcf_fn, "senecio.gds", method="copy.num.of.ref") 
  # Assign the new file to an object 
  genofile <- snpgdsOpen("senecio.gds")
  # Run PCA 
  senecio_pca <- snpgdsPCA(genofile, autosome.only=FALSE)

  # Extract PCs and make dataset with samples 
  PC1 <- senecio_pca$eigenvect[,1]
  PC2 <- senecio_pca$eigenvect[,2]
  PCA <- cbind(PC1, PC2)
  colnames(PCA) <- c("PC1", "PC2")
  PCA <- cbind(samples, PCA)

  # Plot the PCA 
  tiff(paste("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only/", str_split(files[i], "_", n = 2, simplify=TRUE), ".tiff"), units="in", width=9, height=7, res=300)
  ggplot(PCA, aes(x=PC1, y=PC2, col=population)) + geom_point(size = 1) + theme_bw() +
    scale_colour_discrete(name="Populations") +
    geom_label_repel(aes(label = sample), label.size = NA, fill = "NA", show.legend = FALSE, segment.color = "transparent") +
    xlab("PC1") + ylab("PC2") +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"), axis.title = element_text(size = 12),axis.title.x = element_text(margin = margin(t = 20)), axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)), text = element_text(family = "Helvetica") ) 
  dev.off()
  file.remove("/home/564/km6006/senecio.gds")
}