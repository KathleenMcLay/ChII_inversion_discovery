### Plot Fst with loess smoothing 

install.packages("fANCOVA")
library(fANCOVA)
library(ggplot2)
library(tools)
library(stringr)

fst_data <- read.csv("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/data/Fst/fst_data_all.csv")
unique_inversions <- unique(fst_data$inversion) # get unique inversions

pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/figures/fst_plots.pdf", height=4.5, width=4.5)
for (inv in unique_inversions) {
  k = 1

  current_inv <- subset(fst_data, inversion == inv)
  current_inv <- na.omit(current_inv)
  current_inv$index <- 1:nrow(current_inv)
  current_inv$WEIGHTED_FST[current_inv$WEIGHTED_FST < 0] <- 0
  current_inv$MEAN_FST[current_inv$MEAN_FST < 0] <- 0
  
  # add window midpoint value to data for plotting
  for (i in 1:nrow(current_inv)) {
    current_inv$mid[i]=floor((current_inv$BIN_START[i]+current_inv$BIN_END[i])/2)
  }
  
  # determine optimal loess smoothing values
  fst.lo <- loess.as(current_inv$index, current_inv$WEIGHTED_FST, degree = 2, criterion ="aicc", user.span = NULL, plot = F)
  fst.lo.pred <- predict(fst.lo)
  
  #add loess smoothing values to the dataset 
  fst.dat <- cbind(current_inv,fst.lo.pred)
  
  #get current inversion start and end, to add to plot
  current_inversion <- inv
  
  split_string <- unlist(strsplit(current_inversion, "[:-]"))
  inv_start <- as.numeric(split_string[2])
  inv_end <- as.numeric(split_string[3])
  
  current_pop <- unique(current_inv$population)
  
  k = k + 1 
  
  #plot
  fst.plot <- ggplot(current_inv, aes(x=mid/1000000, y=WEIGHTED_FST)) + theme_bw() +
    geom_point(col="grey") +
    geom_line(data=fst.dat, aes(x=mid/1000000, y=fst.lo.pred),col="navyblue") +
    annotate("rect", fill = "lightblue", alpha = 0.5, 
             xmin = inv_start/1000000, xmax = inv_end/1000000,
             ymin = 0, ymax = 1) +
    xlab("Positions (Mbp)") +
    ylab("Fst (weighted)") +
    (ggtitle(paste(current_pop, " ", inv)) ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"),
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12),
          legend.position = "none")
  print(fst.plot)
}
dev.off()