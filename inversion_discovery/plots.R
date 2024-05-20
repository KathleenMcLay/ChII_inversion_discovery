### Panel plots 

library(fANCOVA)
library(ggplot2)
library(tools)
library(stringr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggrepel)
library(tidyverse)

### READ IN DATA 

# MDS 
mds_data <- read.csv("/Users/kathleenmclay/.../MDSinfo_data_all.csv", header=TRUE)
out_data <- read.csv("/Users/kathleenmclay/.../MDS_outliers_data_all.csv", header=TRUE)
win_data <- read.csv("/Users/kathleenmclay/.../MDS_win_regions_data_all.csv", header=TRUE)
cm <- read.table("/Users/kathleenmclay/.../centromeres.tsv", header = TRUE)
cm$centromere <- as.double(cm$centromere)
# PCA 
pca_data <- read.csv("/Users/kathleenmclay/.../PCA_genotypes_all.csv", header=TRUE)
pca_stat_data <- read.csv("/Users/kathleenmclay/.../PCA_stats_all.csv", header=TRUE)
# Het 
het_data <- read.csv("/Users/kathleenmclay/.../het_data_all.csv")
het_data$genotype <- as.character(het_data$genotype)
het_data$genotype <- trimws(het_data$genotype, whitespace = '"')
# LD 
ld_data <- read.csv("/Users/kathleenmclay/.../new_LD.csv")
ld_data <- drop_na(ld_data) 
ld_data <- ld_data[2:11]
ld_cg_data <- read.csv("/Users/kathleenmclay/.../new_LD_CG.csv")
ld_cg_data <- drop_na(ld_cg_data)
ld_cg_data <- ld_cg_data[2:11]
# Fst
fst_data <- read.csv("/Users/kathleenmclay/.../fst_data_all.csv")
fst_data <- na.omit(fst_data)
fst_data$mid <- floor((fst_data$BIN_START + fst_data$BIN_END) / 2)
fst_data$WEIGHTED_FST[fst_data$WEIGHTED_FST < 0] <- 0
fst_data$MEAN_FST[fst_data$MEAN_FST < 0] <- 0
# Pi 
pi_data <- read.csv("/Users/kathleenmclay/.../pi_data_all.csv")
pi_data$pop <- as.character(pi_data$pop)
pi_data <- drop_na(pi_data) 
pi_data$mid <- floor((pi_data$window_pos_1 + pi_data$window_pos_2) / 2)
# Dxy 
dxy_data <- read.csv("/Users/kathleenmclay/.../dxy_data_all.csv")
dxy_data$pop1 <- as.character(dxy_data$pop1)
dxy_data$pop2 <- as.character(dxy_data$pop2)
dxy_data <- drop_na(dxy_data)
dxy_data$mid <- floor((dxy_data$window_pos_1 + dxy_data$window_pos_2) / 2)

### Get list of unique inversions 
unique_inversions <- unique(mds_data$inversion)

### Set plot theme variable
p_theme <- theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"),
                 axis.title.x = element_text(margin = margin(t = 20), size = 12), 
                 axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 12),
                 plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
                 axis.ticks = element_line(colour="black"), panel.background=element_rect(colour="black"), panel.border=element_rect(colour="black", size=0.75),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.position = "none") 

### Initialize PDF
#pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/figures/lpca_panel_plots_TEST.pdf", height=21, width=3)

### PLOT 
for (inv in unique_inversions) {
  
  current_mds <- subset(mds_data, inversion == inv)
  
  tiff(paste("/Users/kathleenmclay/.../figures/", inv, "_", current_mds$population, "_panel_plot.tiff"), width = 3, height = 21, units = 'in', res = 300)
  
  # get current inversion start and end
  split_string <- unlist(strsplit(inv, "[:-]"))
  inv_start <- as.numeric(split_string[2])
  inv_end <- as.numeric(split_string[3])
  
  # MDS
  current_outliers <- subset(out_data, mds_coord == current_mds$mds_coord & population == current_mds$population)

  selected_columns <- c(names(win_data)[1:4], tolower(current_mds$mds), names(win_data)[45:47])
  current_windows <- win_data[selected_columns]
  current_windows <- subset(current_windows, population == current_mds$population & chrom == current_mds$chromosome)
  current_windows <- transform(current_windows, outlier = ifelse(n %in% current_outliers$n, "Outlier", "Non-outlier"))
  colnames(current_windows)[5] <- "mds"
  current_pop <- unique(current_windows$population)

  cent <- subset(cm, scaffold == current_mds$chromosome)

  genome_plot <- ggplot(current_windows, aes(x=mid/1000000, y=mds, color=outlier)) +
    geom_point() + theme_bw() +
    geom_vline(xintercept = cent$centromere, colour="black", linetype = "longdash", linewidth=0.25) +
    scale_color_manual(values=c("gray90","purple")) +
    xlab("Position (Mb)") + ylab(toupper(current_mds$mds)) +
    p_theme
  
  # PCA
  current_pca <- subset(pca_data, inversion == inv)
  PC1_perc <- subset(pca_stat_data[,3], pca_stat_data$inversion == inv)
  PC2_perc <- subset(pca_stat_data[,4], pca_stat_data$inversion == inv)

  PCA_plot <- ggplot(current_pca, aes(x=PC1, y=PC2, col=genotype)) +
    geom_point(size = 2) + theme_bw() +
    scale_color_manual(name="Genotype", values=c("red","purple","blue")) +
    xlab(paste("PC1 (",PC1_perc, "%)")) + ylab(paste("PC2 (",PC2_perc, "%)")) +
    p_theme

  # Het
  current_het <- subset(het_data, inversion == inv)

  het_plot <- ggplot(current_het, aes(x=as.character(genotype), y=het, fill=as.character(genotype))) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(name="Genotype", values=c("red","purple","blue")) +
    xlab("Genotype") + ylab("Heterozygosity") +
    p_theme

  # LD
  ld_dat <- subset(ld_data, inversion == inv)
  ld_cg_dat <- subset(ld_cg_data, inversion == inv)

  ld_plot <- ggplot(ld_dat, aes(x = win1 / 1000000, y = win2 / 1000000)) +
    theme_classic() +
    geom_segment(mapping=aes(x=inv_start/1000000,xend=inv_end/1000000,y=-4,yend=-4), col="purple",size = 2.5) +
    geom_segment(mapping=aes(x=-4,xend=-4,y=inv_start/1000000,yend=inv_end/1000000), col="purple",size = 2.5) +
    geom_tile(aes(fill = std_mean_r2)) + #add all samples LD to upper
    scale_fill_gradientn(colours = c("grey95", "blue", "red"), values = c(0, 0.5, 1), name = "LD") +
    geom_tile(data=ld_cg_dat, aes(x = win2 / 1000000, y = win1 / 1000000, fill = std_mean_r2)) + #add common genotype samples LD to lower
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    coord_fixed(ratio = 1) +
    xlab("Position (Mb)") + ylab("Position (Mb)") +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 12),
          legend.position = "none") 

  # Fst
  current_fst <- subset(fst_data, inversion == inv)
  current_fst$index <- 1:nrow(current_fst)
  # determine optimal loess smoothing values
  fst.lo <- loess.as(current_fst$index, current_fst$WEIGHTED_FST, degree = 2, criterion ="aicc", user.span = NULL, plot = F)
  fst.lo.pred <- predict(fst.lo)
  # add loess smoothing values to the dataset
  fst.dat <- cbind(current_fst,fst.lo.pred)
  # plot
  fst_plot <- ggplot(current_fst, aes(x=mid/1000000, y=WEIGHTED_FST)) +
    geom_point(col="gray90") + theme_bw() +
    geom_line(data=fst.dat, aes(x=mid/1000000, y=fst.lo.pred),col="purple") +
    annotate("rect", fill = "#C1C4AD", alpha = 0.5,
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = 0, ymax = 1) +
    xlab("Positions (Mbp)") +ylab("Fst (weighted)") +
    p_theme

  # Pi
  pi_dat <- subset(pi_data, inversion == inv)
  #subset to population, calculate loess smoothing values
  pi_dat_0 <- subset(pi_dat, pop == "0")
  pi_dat_0$index <- 1:nrow(pi_dat_0)
  pi.lo <- loess.as(pi_dat_0$index, pi_dat_0$avg_pi, degree = 0, criterion ="aicc", user.span = NULL, plot = F)
  pi.lo.pred <- predict(pi.lo)
  pi_dat_0 <- cbind(pi_dat_0,pi.lo.pred)

  pi_dat_2 <- subset(pi_dat, pop == "2")
  pi_dat_2$index <- 1:nrow(pi_dat_2)
  pi.lo <- loess.as(pi_dat_2$index, pi_dat_2$avg_pi, degree = 0, criterion ="aicc", user.span = NULL, plot = F)
  pi.lo.pred <- predict(pi.lo)
  pi_dat_2 <- cbind(pi_dat_2,pi.lo.pred)

  if (any(pi_dat$pop == "1")) {
    pi_dat_1 <- subset(pi_dat, pop == "1")
    pi_dat_1$index <- 1:nrow(pi_dat_1)
    pi.lo <- loess.as(pi_dat_1$index, pi_dat_1$avg_pi, degree = 0, criterion = "aicc", user.span = NULL, plot = FALSE)
    pi.lo.pred <- predict(pi.lo)
    pi_dat_1 <- cbind(pi_dat_1, pi.lo.pred)

    rm(pi_dat)
    pi_dat <- rbind(pi_dat_0, pi_dat_1, pi_dat_2)
    pi_dat$pop <- as.character(pi_dat$pop)

  } else {
    rm(pi_dat)
    pi_dat <- rbind(pi_dat_0, pi_dat_2)
    pi_dat$pop <- as.character(pi_dat$pop)
  }

  ymin <- min(pi_dat$pi.lo.pred)
  ymax <- max(pi_dat$pi.lo.pred)

  # plot
  pi_plot <- ggplot(pi_dat, aes(x=mid/1000000, y=pi.lo.pred, col=pop, group=pop)) +
    geom_line() + theme_bw() +
    scale_colour_manual(name="Genotype", values=c("red","purple","blue")) +
    annotate("rect", fill = "#C1C4AD", alpha = 0.5,
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = ymin, ymax = ymax) +
    xlab("Positions (Mbp)") + ylab(expression(paste(pi))) +
    p_theme

  # Dxy
  dxy_dat <- subset(dxy_data, inversion == inv)
  # subset to population, calculate loess smoothing values
  dxy_dat <- dxy_dat[(dxy_dat$pop1 == "0" & dxy_dat$pop2 == "2") | (dxy_dat$pop1 == "2" & dxy_dat$pop2 == "0"), ]
  dxy_dat$index <- 1:nrow(dxy_dat)
  dxy.lo <- loess.as(dxy_dat$index, dxy_dat$avg_dxy, degree = 0, criterion ="aicc", user.span = NULL, plot = F)
  dxy.lo.pred <- predict(dxy.lo)
  dxy_dat <- cbind(dxy_dat,dxy.lo.pred)

  ymin <- min(dxy_dat$dxy.lo.pred)
  ymax <- max(dxy_dat$dxy.lo.pred)
  # plot
  dxy_plot <- ggplot(dxy_dat, aes(x=mid/1000000, y=dxy.lo.pred)) +
    geom_line(col="purple") + theme_bw() +
    annotate("rect", fill = "#C1C4AD", alpha = 0.5,
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = ymin, ymax = ymax) +
    xlab("Positions (Mbp)") + ylab("Dxy") +
    p_theme
  
  #format plot
  print(
    heights = 3,
    widths = 3,
    grid.arrange(genome_plot, PCA_plot, het_plot, ld_plot, fst_plot, pi_plot, dxy_plot,
      layout_matrix = cbind(c(1, 2, 3, 4, 5, 6, 7)),
      top = textGrob((paste(current_pop, " ", inv)) , gp=gpar(fontsize=10, font=1))
    )
  )
  dev.off()
}

#dev.off()