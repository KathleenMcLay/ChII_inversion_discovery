###MDS PLOT###

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)

#read in files 
mds_info <- read.table("/.../results/lpca_tables/inversions.csv")
mds_info$mds <- tolower(mds_info$MDS)
mds_info$mds <- substr(mds_info$MDS, 1, 5)
outliers <- read.table("/.../results/lpca_tables/3_outlier_windows.txt")
win_regions <- read.table("/.../results/lpca_tables/1_win_regions.txt")
cm <- read.table("/.../centromeres.tsv", header = TRUE)

#Filter the outliers and mds data(win_regions) by mds_info which lists the outlier clusters along the genome for all MDS. Then plot.  
pdf("/.../results/MDS_plots.pdf", height=7, width=7)
for (i in 1:nrow(mds_info)) {
  outliers %>%
    filter(mds_coord == mds_info$MDS[i]) -> outs
  win_regions %>% select(1:4, mds_info$mds[i], 45) %>%
    filter(chrom == mds_info$scaffold[i]) %>%
    mutate(outlier=case_when(n %in% outs$n ~ "Outlier", 
                             TRUE ~ "Non-outlier")) -> new
  colnames(new) <- c("chrom", "start", "end", "mid", "mds", "n", "outlier")
  cm %>%
    mutate(centromere=as.double(centromere)) %>%
    filter(scaffold == mds_info$scaffold[i]) -> cent
  genome_plot <- new %>%
    ggplot(., aes(x=mid/1000000, y=mds, color=outlier)) + geom_point() + theme_bw() +
    geom_vline(xintercept = cent$centromere, colour="black", linetype = "longdash", linewidth=0.25) +
    scale_color_manual(values=c("gray90","green4")) +
    ggtitle(paste(mds_info$MDS[i], " ", mds_info$scaffold[i])) +
    xlab("Position (Mb)") + ylab(toupper(mds_info$mds[i])) +
    theme(legend.position="none", 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.ticks = element_line(colour="black"), panel.background=element_rect(colour="black"), panel.border=element_rect(colour="black", size=0.75))
  print(genome_plot)
}
dev.off()
