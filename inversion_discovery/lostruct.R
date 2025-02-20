#Local PCA + MDS outlier analysis from Todesco et al. 2020 

library(lostruct)
library(doParallel)
library(tools)
library(stringr)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(grid)
select=dplyr::select

window_size <- 100 # window size 
k_kept <- 40 # Number of MDS to use
min_windows <- 4 # Minimum number of windows to call an outlier regions
max_distance_between_outliers <- 10 # Max distance between oulier windows to be included in the same region 
n_permutations <- 1000 # Number of permutation for chromosomal clustering test
min_cor <- 0.8 # Correlation theshold for collapsing MDS

pop_pairs <- list("D00 H00", "D01 H01", "D03 H02", "D04 H05", "D05 H06", "D32 H12")

foreach(p=1:length(pop_pairs)) %do% {
  ### subset the data to the current population pair 
  p1 <- str_split(pop_pairs[p], " ", n = 2, simplify=TRUE)[1]
  p2 <- str_split(pop_pairs[p], " ", n = 2, simplify=TRUE)[2]
  
  samples <- read.table("/home/564/km6006/Scripts/inversion_paper/samples.tsv")
  colnames(samples) <- "name"
  sample.list <- samples[substr(samples$name, 1, 3) == p1 | substr(samples$name, 1, 3) == p2, ]
  file_path <- paste("/home/564/km6006/Scripts/inversion_paper/", p1, p2, ".txt", sep = "")
  write.table(sample.list, file_path, col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  x <- paste("bcftools view --threads 12 -O b -S /home/564/km6006/Scripts/inversion_paper/", p1,p2, ".txt -o /g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/",p1,p2,".bcf /g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/VO_merged_final_NS.recode.bcf", sep = "")
  system(x)
  index <- paste("bcftools index --threads 24 /g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/",p1,p2,".bcf", sep = "")
  system(index)
  
  ### run lostruct, removing NA values in any matrix 
  bcf.file <- paste("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/",p1,p2, ".bcf", sep = "")
  sites <- vcf_positions(bcf.file)
  win.fn.snp <- vcf_windower(bcf.file, size=window_size, type="snp", sites=sites)
  snp.pca <- eigen_windows(win.fn.snp, k=2, mc.cores=24) 
  pcdist <- pc_dist(snp.pca, mc.cores=24) 
  na.wins <- is.na(snp.pca[,1])
  pcdist <- pcdist[!na.wins, !na.wins]
  nan.wins <- pcdist[,1]=="NaN"
  pcdist <- pcdist[!nan.wins, !nan.wins]
  # generate MDS coordinates 
  mds <- cmdscale(pcdist, eig=TRUE, k=k_kept) 
  mds.coords <- mds$points
  colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
  # create a dataset of windows with the chromosome, start, end and midpoint of the window as the first four columns
  win.regions <- region(win.fn.snp)()
  win.regions <- win.regions[!na.wins,][!nan.wins,]
  win.regions %>% mutate(mid=(start+end)/2) -> win.regions
  # rename the MDS coordinate columns in win.regions with format "mds01"
  for (k in 1:k_kept){
    name = paste("mds", str_pad(k, 2, pad = "0"), sep="")
    win.regions$tmp <- "NA"
    win.regions <- win.regions %>% rename(!!name := tmp)
  }
  #add the mds values to win.regions data, skipping the first four columns 
  for (i in 1:k_kept){
    j = i + 4
    win.regions[,j] <- mds.coords[,i]
  }
  # add a window count column to the win.regions dataset 
  win.regions$n <- 1:nrow(win.regions)
  
  #write the data to file
  win_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_1_win_regions.txt", sep = "")
  write.table(win.regions, win_tn, sep="\t")
}