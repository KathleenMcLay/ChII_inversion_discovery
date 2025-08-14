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
install.packages("ggnewscale")
library(ggnewscale)

### READ IN DATA 

# MDS 
inv_data <- read.csv("4_inversion_discovery/final_inversion_list.csv", header=TRUE)
out_data <- read.csv("4_inversion_discovery/regions_windows_summary.csv", header=TRUE)
win_data <- read.csv("4_inversion_discovery/MDS_win_regions_data_all.csv", header=TRUE)
cm <- inv_data[, c(4, 13)]

# PCA & Heterozygosity data
pca_het_data <- read.csv("5_inversion_pop_gen/PCA_het.csv", header=TRUE)
pca_het_data$pca_genotype <- as.character(pca_het_data$pca_genotype)
pca_het_data$pca_genotype <- trimws(pca_het_data$pca_genotype, whitespace = '"')

# LD 

# LD for all samples togeather
ld_data <- read.csv("5_inversion_pop_gen/ld_concatenated_ALL.csv")
ld_data <- drop_na(ld_data) 
# Create LD common genotype dataset 
ld_g0 <- read.csv("5_inversion_pop_gen/ld_concatenated_G0.csv")
ld_g0 <- ld_g0[ld_g0$inversion %in% c("scaffold_1:121452175-140943131", "scaffold_1:144787767-146039534",
                                        "scaffold_10:112787027-121639654", "scaffold_10:98975439-103069101",
                                        "scaffold_11:3167182-7087909", "scaffold_13:36478957-45395630",
                                        "scaffold_2:176510351-185799375", "scaffold_2:189437492-196374760",
                                        "scaffold_4:167106062-174439192", "scaffold_5:161724490-166529893",
                                        "scaffold_5:17269039-25882675", "scaffold_5:195968395-201339920",
                                        "scaffold_6:167476923-170468194",
                                        "scaffold_8:127818797-139663425", "scaffold_8:140439606-145099475",
                                        "scaffold_8:141837386-150869394", "scaffold_8:26711179-27397213",
                                        "scaffold_8:9748739-19694111", "scaffold_9:71975361-89342649"), ]


ld_g2 <- read.csv("5_inversion_pop_gen/ld_concatenated_G2.csv")
ld_g2 <- ld_g2[ld_g2$inversion %in% c("scaffold_10:99316558-102000668", "scaffold_11:23336059-25636162",
                                        "scaffold_12:32242414-35864311", "scaffold_12:47567182-50302659",
                                        "scaffold_12:51459708-54083055", "scaffold_2:207418890-212464244",
                                        "scaffold_2:251165261-254898012", "scaffold_5:177368560-188528516",
                                        "scaffold_5:179145606-195409507", "scaffold_5:206101837-212778787",
                                        "scaffold_5:28034585-39379225", "scaffold_6:18238202-28839170",
                                        "scaffold_6:33674742-41444463", "scaffold_6:51049655-52640067",
                                        "scaffold_6:81485-2851736", "scaffold_7:28599735-33301123",
                                        "scaffold_7:28852469-32048012", "scaffold_8:119984128-124368086",
                                        "scaffold_8:126704823-135605926", "scaffold_9:71043689-94971448", 
                                        "scaffold_6:33664350-37600176"), ]

ld_cg_data <- rbind(ld_g0, ld_g2)

# Fst
fst_data <- read.csv("5_inversion_pop_gen/fst_data.csv")
fst_data <- na.omit(fst_data)
fst_data$mid <- floor((fst_data$BIN_START + fst_data$BIN_END) / 2)
fst_data$WEIGHTED_FST[fst_data$WEIGHTED_FST < 0] <- 0
fst_data$MEAN_FST[fst_data$MEAN_FST < 0] <- 0
# Pi 
pi_data <- read.csv("5_inversion_pop_gen/pi_data.csv")
pi_data$pop <- as.character(pi_data$pop)
pi_data <- drop_na(pi_data) 
pi_data$mid <- floor((pi_data$window_pos_1 + pi_data$window_pos_2) / 2)
# Dxy 
dxy_data <- read.csv("5_inversion_pop_gen/dxy_data.csv")
dxy_data$pop1 <- as.character(dxy_data$pop1)
dxy_data$pop2 <- as.character(dxy_data$pop2)
dxy_data <- drop_na(dxy_data)
dxy_data$mid <- floor((dxy_data$window_pos_1 + dxy_data$window_pos_2) / 2)

### Get list of unique inversions 
unique_inversions <- unique(inv_data$inversion)

### Set plot theme variable
p_theme <- theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"),
                 axis.title.x = element_text(margin = margin(t = 20), size = 10), 
                 axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 10),
                 plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
                 axis.ticks = element_line(colour="black"), panel.background=element_rect(colour="black"), panel.border=element_rect(colour="black", size=0.75),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.position = "none") 

### Initialize PDF
#pdf("lpca_panel_plots_TEST.pdf", height=21, width=3)

### PLOT 
for (inv in unique_inversions) {
  
  current_mds <- subset(inv_data, inversion == inv)
  inv_code <- current_mds$inversion_code
  
  tiff(paste(inv, "_", current_mds$discovery_population, "_panel_plot.tiff"), width = 3, height = 21, units = 'in', res = 400)
  
  # get current inversion start and end
  split_string <- unlist(strsplit(inv, "[:-]"))
  inv_start <- as.numeric(split_string[2])
  inv_end <- as.numeric(split_string[3])
  
  # MDS
  current_outliers <- subset(out_data, haploblock == current_mds$inversion & population == current_mds$discovery_population)
  
  selected_columns <- c(names(win_data)[1:4], tolower("mds01"), names(win_data)[45:46])
  current_windows <- win_data[selected_columns]
  current_windows <- subset(current_windows, population == current_mds$discovery_population & chrom == current_mds$scaffold)
  current_windows <- transform(current_windows, outlier = ifelse(n %in% current_outliers$n, "Outlier", "Non-outlier"))
  colnames(current_windows)[5] <- "mds"
  current_pop <- unique(current_windows$population)
  
  cent <- subset(cm, scaffold == current_mds$scaffold)
  cent <- cent[!duplicated(cent$hist_centromere_position_Mb), ]
  
  genome_plot <- ggplot(current_windows, aes(x=mid/1000000, y=mds, color=outlier)) +
    geom_point() + theme_bw() +
    geom_vline(xintercept = cent$hist_centromere_position_Mb, colour="black", linetype = "longdash", linewidth=0.25) +
    scale_color_manual(values=c("gray90","#440154FF")) +
    xlab("Position (Mb)") + ylab("MDS1") +
    p_theme

  # PCA
  current_pca_het <- subset(pca_het_data, inversion == inv)

  PC1_perc <- current_mds[,14]
  PC2_perc <- current_mds[,15]

  PCA_plot <- ggplot(current_pca_het, aes(x=PC1, y=PC2, col=pca_genotype)) +
    geom_point(size = 2) + theme_bw() +
    scale_color_manual(name="Genotype", values=c("#39568CFF","#1F968BFF","#95D840FF")) +
    xlab(paste("PC1 (",PC1_perc, "%)")) + ylab(paste("PC2 (",PC2_perc, "%)")) +
    p_theme

  # Het
  het_plot <- ggplot(current_pca_het, aes(x=as.character(pca_genotype), y=HET, fill=as.character(pca_genotype))) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(name="Genotype", values=c("#39568CFF","#1F968BFF","#95D840FF")) +
    xlab("Genotype") + ylab("Heterozygosity") +
    p_theme

  # LD
  ld_dat <- subset(ld_data, inversion == inv)
  ld_cg_dat <- subset(ld_cg_data, inversion == inv)

  ld_plot <- ggplot(ld_dat, aes(x = win1 / 1000000, y = win2 / 1000000)) +
    theme_classic() +
    geom_segment(mapping=aes(x=inv_start/1000000,xend=inv_end/1000000,y=-4,yend=-4), col="#440154FF",size = 2) +
    geom_segment(mapping=aes(x=-4,xend=-4,y=inv_start/1000000,yend=inv_end/1000000), col="#440154FF",size = 2) +

    # First geom_tile with "fill" aesthetic
    geom_tile(aes(fill = mean_r2)) +
    scale_fill_gradientn(colours = c("grey95", "#95D840FF", "#39568CFF"),
                         values = c(0, 0.5, 1),
                         name = "mean r2") +

    ggnewscale::new_scale_fill() +
    geom_tile(data=ld_cg_dat, aes(x = win2 / 1000000, y = win1 / 1000000, fill = mean_r2)) +
    scale_fill_gradientn(colours = c("grey95", "grey", "darkgrey"),
                         values = c(0, 0.5, 1),
                         name = "mean r2") +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    coord_fixed(ratio = 1) +
    xlab("Position (Mb)") + ylab("Position (Mb)") +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), size = 10),
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 10),
          legend.position = "right",  # Changed to show both legends
          legend.key.size = unit(0.2, "cm"),        # Size of legend keys (color boxes)
          legend.title = element_text(size = 5),    # Legend title text size
          legend.text = element_text(size = 4),     # Legend label text size
          legend.margin = margin(l = 5, r = 5),     # Reduce margin around legends
          legend.box.spacing = unit(0.1, "cm"))     # Reduce space between multiple legends

    
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
    geom_line(data=fst.dat, aes(x=mid/1000000, y=fst.lo.pred),col="#440154FF") +
    annotate("rect", fill = "#C1C4AD", alpha = 0.5,
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = 0, ymax = 1) +
    xlab("Position (Mb)") +ylab("Fst (weighted)") +
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
    scale_colour_manual(name="Genotype", values=c("#39568CFF","#1F968BFF","#95D840FF")) +
    annotate("rect", fill = "#C1C4AD", alpha = 0.5,
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = ymin, ymax = ymax) +
    xlab("Position (Mb)") + ylab(expression(paste(pi))) +
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
    geom_line(col="#440154FF") + theme_bw() +
    annotate("rect", fill = "#C1C4AD", alpha = 0.5,
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = ymin, ymax = ymax) +
    xlab("Position (Mb)") + ylab("Dxy") +
    p_theme

  #format plot
  print(
    heights = 3,
    widths = 3,
    grid.arrange(genome_plot, PCA_plot, het_plot, fst_plot, dxy_plot, pi_plot, ld_plot,
                 layout_matrix = cbind(c(1, 2, 3, 4, 5, 6, 7)),
                 top = textGrob((paste(current_pop, " ", inv_code)) , gp=gpar(fontsize=10, font=1))
    )
  )
  dev.off()
}

#dev.off()
