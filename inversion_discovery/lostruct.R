#localPCA w/ lostruct

library(lostruct)
library(doParallel)
library(tools)
library(stringr)
library(tidyverse)

# 1000 SNP window size 
window_size <- 1000 

# list of population pairs 
pop_pairs <- list("D00 H00", "D01 H01", "D03 H02", "D04 H05", "D05 H06", "D32 H12")

foreach(p=1:length(pop_pairs)) %do% {
  ### subset the data to the current population pair 
  p1 <- str_split(pop_pairs[p], " ", n = 2, simplify=TRUE)[1]
  p2 <- str_split(pop_pairs[p], " ", n = 2, simplify=TRUE)[2]
  
  # subset the samples
  samples <- read.table("/g/data/ht96/McLay_UQ/.../samples.tsv")
  colnames(samples) <- "name"
  sample.list <- samples[substr(samples$name, 1, 3) == p1 | substr(samples$name, 1, 3) == p2, ]
  file_path <- paste("/g/data/ht96/McLay_UQ/.../", p1, p2, ".txt", sep = "")
  write.table(sample.list, file_path, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # subset the bcf file
  x <- paste("bcftools view --threads 12 -O b -S /g/data/ht96/McLay_UQ/.../", p1,p2, ".txt -o /g/data/ht96/McLay_UQ/.../",p1,p2,".bcf /g/data/ht96/McLay_UQ/.../VO_merged_final_NS.recode.bcf", sep = "")
  system(x)
  index <- paste("bcftools index --threads 24 /g/data/ht96/McLay_UQ/.../",p1,p2,".bcf", sep = "")
  system(index)
  
  ### run lostruct, removing NA values in any matrix 
  bcf_file <- paste("/g/data/ht96/McLay_UQ/.../",p1,p2, ".bcf", sep = "")
  sites <- vcf_positions(bcf_file)
  windows <- vcf_windower(bcf_file, size=window_size, type="snp", sites=sites)
  eigen_win <- eigen_windows(windows, k=2, mc.cores=24) 
  pc_dist <- pc_dist(eigen_win, mc.cores=24) 
  na_win <- is.na(eigen_win[,1])
  pc_dist <- pc_dist[!na.wins, !na.wins]
  nan.wins <- pc_dist[,1]=="NaN"
  pc_dist <- pc_dist[!nan.wins, !nan.wins]
  
  # generate MDS coordinates 
  mds <- cmdscale(pc_dist, eig=TRUE, k=2) 
  mds_coords <- mds$points
  colnames(mds_coords) <- paste("MDS coordinate", 1:ncol(mds_coords))
 
 # create a dataset of windows with the chromosome, start, end and midpoint of the window as the first four columns
  win_data <- region(windows)()
  win_data <- win_data[!na.wins,][!nan.wins,]
  win_data$mid <- (win_data$start + win_data$end)/2

  # add the mds values to win_data data, skipping the first four columns 
  for (i in 1:k_kept){
    j = i + 4
    win_data[,j] <- mds_coords[,i]
  }
  # add a window count column to the win_data dataset 
  win_data$n <- 1:nrow(win_data)
  
  # write the data to file
  write.table(win_data, paste("/g/data/ht96/McLay_UQ/.../", p1,p2, "_win_regions.txt", sep = ""), sep="\t")
}