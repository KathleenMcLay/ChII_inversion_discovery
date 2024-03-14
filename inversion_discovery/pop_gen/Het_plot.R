### HET values are calculate in the LPC script. This script plots the heterozygosity values for each inversion.

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)

heterozygosity <- read.table("/.../results/lpca_tables/9_heterozygosity.txt")
heterozygosity$genotype <- as.character(heterozygosity$genotype)

new <- tibble(mds_coord=character(), name=character(), 
              PC1=numeric(), PC2=numeric(), genotype=character(), het=character())

pdf("/.../results/HET_plots.pdf", height=7, width=7)
for (i in 1:nrow(mds_info)) {
  for (j in 1:nrow(heterozygosity)) {
    if (mds_info$mds_coord[i] == heterozygosity$mds_coord[j]) {
      new <- rbind(new, heterozygosity[j,]) 
    }
  }
  het_plot <- new %>%
    ggplot(., aes(x=as.character(genotype), y=het, fill=as.character(genotype))) + geom_boxplot() + theme_bw() +
    ggtitle(paste(mds_info$mds_coord[i], " ", mds_info$chromosome[i])) +
    scale_fill_manual(name="Genotype", values=c("red","purple","blue")) +
    xlab("Genotype") + ylab("Heterozygosity") +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12))
  print(het_plot)
  new <- tibble(mds_coord=character(), name=character(), 
                  PC1=numeric(), PC2=numeric(), genotype=character(), het=character())
}
dev.off()