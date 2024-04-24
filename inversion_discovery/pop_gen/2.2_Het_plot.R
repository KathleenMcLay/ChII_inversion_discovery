### HET values are calculate in the LPC script. This script plots the heterozygosity values for each inversion.

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)

het_data <- read.csv("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/data/Het/het_data_all_inversions.csv")
het_data$genotype <- as.character(het_data$genotype)
het_data$genotype <- trimws(het_data$genotype, whitespace = '"')

unique_inversions <- unique(het_data$inversion) # get unique inversions

pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/figures/HET_plots.pdf", height=4.5, width=4.5)
for (inv in unique_inversions) {
  current_inv <- subset(het_data, inversion == inv) # subset the data to current inversion
  current_pop <- unique(current_inv$pop)
  het_plot <- current_inv %>%
    ggplot(., aes(x=as.character(genotype), y=het, fill=as.character(genotype))) + geom_boxplot() + theme_bw() +
    ggtitle(paste(current_pop, " ", inv)) +
    scale_fill_manual(name="Genotype", values=c("red","purple","blue")) +
    xlab("Genotype") + ylab("Heterozygosity") +
    theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12), legend.position = "none")
  print(het_plot)
}
dev.off()