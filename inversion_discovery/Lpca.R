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

  ### Look for "chromosomally clustered" windows along each MDS coordinate (in both positive and negetive directions)
  
  # extract the mds column names (columns 5-40) from win.regions
  mds_pcs <- colnames(win.regions)[5:(ncol(win.regions)-1)]
  # create an empty data frame to hold the clustering data 
  mds_clustering <- tibble(mds_coord=character(), direction=character(), clust_pvalue=numeric(), outliers=numeric(), n1_outliers=numeric(), 
                          high_cutoff=numeric(), lower_cutoff=numeric(), chr=character())
  # for all 40 mds calculate SD and set the MDS outlier threshold to 4x the SD 
  for (mds_chosen in mds_pcs){
    print(paste("Processing",mds_chosen))
    win.regions %>%
      mutate_(the_mds=mds_chosen) %>%
      summarize(sd_mds=sd(the_mds)) %>% pull() -> sd_mds
    mds_high_cutoff <- sd_mds*4
    mds_low_cutoff <- sd_mds*3
    
    # positive windows 
    win.regions %>%
      mutate_(the_mds=mds_chosen) %>%
      mutate(sd_mds=sd(the_mds)) %>%
      filter(the_mds > (sd_mds*4)) -> pos_windows
    if (nrow(pos_windows) >= min_windows){
      permutations <- matrix(nrow=n_permutations, ncol=1)
      for (i in 1:n_permutations){
        max_1 <- win.regions %>% sample_n(nrow(pos_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
          ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
        permutations[i,1] <- max_1
      }
      pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
      pos_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
        pull(chrom) -> clustered_chr
      x <- as.tibble(permutations) %>% filter(V1 >= sampled_max_1) %>% nrow() 
      pvalue <- (x+1)/(n_permutations+1)
      tmp <- tibble(mds_coord=as.character(mds_chosen), direction=as.character("pos"), clust_pvalue=as.numeric(pvalue), outliers=as.numeric(nrow(pos_windows)),
                    n1_outliers=as.numeric(sampled_max_1), high_cutoff=as.numeric(mds_high_cutoff), lower_cutoff=as.numeric(mds_low_cutoff), chr=as.character(clustered_chr))
      mds_clustering <- rbind(mds_clustering, tmp)
    }else{
      tmp <- tibble(mds_coord=as.character(mds_chosen), direction=as.character("pos"), clust_pvalue=as.numeric(NA), outliers=as.numeric(nrow(pos_windows)),
                    n1_outliers=as.numeric(NA), high_cutoff=as.numeric(NA), lower_cutoff=as.numeric(NA), chr=as.numeric(NA))
      mds_clustering <- rbind(mds_clustering, tmp)
    }
    
    #negative windows
    win.regions %>%
      mutate_(the_mds=mds_chosen) %>%
      mutate(sd_mds=sd(the_mds)) %>%
      filter(the_mds < -(sd_mds*4)) -> neg_windows
    if (nrow(neg_windows) >= min_windows){
      permutations <- matrix( nrow=n_permutations, ncol=1)
      for (i in 1:n_permutations){
        max_1 <- win.regions %>% sample_n(nrow(neg_windows)) %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
          ungroup() %>% summarize(sum=sum(count)) %>% pull(sum)
        permutations[i,1] <- max_1
      }
      neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>% 
        ungroup() %>% summarize(sum=sum(count)) %>% pull(sum) -> sampled_max_1
      neg_windows %>% group_by(chrom) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head(n=1) %>%
        pull(chrom) -> clustered_chr
      x <- as.tibble(permutations) %>%  filter(V1 >= sampled_max_1) %>% nrow()
      pvalue <- (x+1)/(n_permutations+1)
      tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(pvalue), outliers=as.numeric(nrow(neg_windows)),
                    n1_outliers=as.numeric(sampled_max_1), high_cutoff = as.numeric(-mds_high_cutoff),lower_cutoff=as.numeric(-mds_low_cutoff),
                    chr=as.character(clustered_chr))
      mds_clustering <- rbind(mds_clustering, tmp)
    }else{
      tmp <- tibble(mds_coord = as.character(mds_chosen),direction = as.character("neg"),clust_pvalue = as.numeric(NA), outliers=as.numeric(nrow(neg_windows)),
                    n1_outliers=as.numeric(NA), high_cutoff = as.numeric(NA),lower_cutoff=as.numeric(NA),
                    chr=as.numeric(NA))
      mds_clustering <- rbind(mds_clustering, tmp)
    }
  }
  mds_clustering %>% filter(clust_pvalue < 0.01) -> sig_mds_clusters

  sig_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_2_sig_mds_clusters.txt", sep= "")
  write.table(sig_mds_clusters, sig_tn, sep="\t")

  #make a list of the outlier windows, and kmeans clustering for each "chromosomally clustered" region of windows 
  outlier_windows <- tibble(chrom=character(), start=numeric(), end=numeric(), mid=numeric(), 
                            the_mds=numeric(), mds_coord=character(), outlier=character(), 
                            n=numeric())
  cluster_genotypes <- tibble(mds_coord=character(), name=character(), PC1=numeric(), genotype=character())

  for (i in 1:nrow(sig_mds_clusters)) {
    coord <- pull(sig_mds_clusters[i,1]) #the first column is the MDS no. 
    direction <- pull(sig_mds_clusters[i,2])
    high_cutoff <- pull(sig_mds_clusters[i,6])
    low_cutoff <- pull(sig_mds_clusters[i,7])
    cluster_chr <- pull(sig_mds_clusters[i,8])
    coord_direction <- paste(coord, "-", direction, sep="")
    print(paste("Testing",coord_direction))
    
    if (direction == "pos"){
      current_windows <- win.regions %>%
        mutate_(the_mds=coord ) %>% #it's looking at the mds column for the first MDS listed in sig_mds_clusters
        mutate(outlier=case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                                TRUE ~ "Non-outlier")) %>%
        filter(outlier != "Non-outlier") %>%
        select(chrom,start,end,mid,the_mds,outlier,n) %>% 
        mutate(mds_coord=coord_direction) %>%
        mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
        mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
        filter(min_dist < max_distance_between_outliers) %>%
        select(-ahead_n,-behind_n,-min_dist)
      windows <- current_windows %>% pull(n)
      outlier_windows <- rbind(outlier_windows, current_windows)
    }else{
      current_windows <- win.regions %>%
        mutate_(the_mds=coord) %>% 
        mutate(outlier=case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                                TRUE ~ "Non-outlier")) %>%
        filter(outlier != "Non-outlier") %>%
        select(chrom,start,end,mid,the_mds,outlier,n) %>%
        mutate(mds_coord=coord_direction) %>%
        mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
        mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
        filter(min_dist < max_distance_between_outliers) %>%
        select(-ahead_n,-behind_n,-min_dist)
      windows <- current_windows %>% pull(n)
      outlier_windows <- rbind(outlier_windows, current_windows)
    }
    out <- cov_pca(win.fn.snp(windows), k=2)
    matrix.out <- t(matrix(out[4:length(out)], ncol=nrow(samples), byrow=T))
    out <- as_tibble(cbind(samples, matrix.out)) %>% 
      rename(name=name, PC1="1", PC2="2") %>% 
      mutate(PC1=as.double(PC1), PC2=as.double(PC2))
    try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
    if("try-error" %in% class(try_3_clusters)){
      kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
    }else{
      kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
    }
    out$cluster <- kmeans_cluster$cluster - 1
    out$cluster <- as.character(out$cluster)
    out$mds_coord <- paste(coord, direction, sep="-")
    genotype.out <- out %>% select(mds_coord,name,PC1,cluster) %>% rename(genotype=cluster)
    cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
  }

  outlier_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_3_outlier_windows.txt", sep= "")
  write.table(outlier_windows, outlier_tn, sep="\t")
  cluster_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_4_cluster_genotypes.txt", sep= "")
  write.table(cluster_genotypes, cluster_tn, sep="\t")

  #Count outlier windows

  mds_counts <- tibble(mds_coord=character(), n_outliers=numeric())
  for (mds in unique(outlier_windows$mds_coord)){
    count_outliers <- outlier_windows %>%
      filter(mds_coord == mds) %>% nrow()
    tmp_tibble <- tibble(mds_coord=as.character(mds), n_outliers=as.numeric(count_outliers))
    mds_counts <- rbind(mds_counts, tmp_tibble)
  }

  mds_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_5_mds_counts.txt", sep= "")
  write.table(mds_counts, mds_tn, sep="\t")

  # Collapsing correlated MDS
  correlated_mds <- tibble(mds1=character(), mds2=character(), correlation=numeric())
  total_mds_coords <- unique(outlier_windows$mds_coord)
  for (i in 1:(length(total_mds_coords)-1)){
    for (j in (i+1):length(total_mds_coords)){
      mds1 <- total_mds_coords[i]
      mds2 <- total_mds_coords[j]
      chr1 <- outlier_windows %>% filter(mds_coord == mds1) %>% select(chrom) %>% unique() %>% pull()
      chr2 <- outlier_windows %>% filter(mds_coord == mds2) %>% select(chrom) %>% unique() %>% pull()
      if (chr1 != chr2){next;}
      cluster_genotypes %>% mutate(mds_coord=gsub("_", "-", mds_coord)) %>% filter(mds_coord == mds1 | mds_coord == mds2) %>%
        select(-PC1) %>%
        spread(mds_coord, genotype) %>% select(-name) -> tmp
      x <- tmp %>% pull(1) %>% as.numeric()
      y <- tmp %>% pull(2) %>% as.numeric()
      test_result <- cor.test(x, y, na.rm=T) #calculate Pearson's product moment correlation coefficient
      tmp_tibble <- tibble(mds1=as.character(mds1), mds2=as.character(mds2), correlation=as.numeric(abs(test_result$estimate)))
      correlated_mds <- rbind(correlated_mds, tmp_tibble)
    }
  }

  cor_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_6_correlated_mds.txt", sep= "")
  write.table(correlated_mds, cor_tn, sep="\t")

  #Check pairs with high correlation and remove the mds that has fewer outlier windows.
  for (i in 1:nrow(correlated_mds)){
    if (correlated_mds[i,3] >= min_cor){ 
      count1 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,1])) %>% pull(n_outliers)
      count2 <- mds_counts %>% filter(mds_coord == as.character(correlated_mds[i,2])) %>% pull(n_outliers)
      if (count1 < count2){
        total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,1]))] -> total_mds_coords
      }else{
        total_mds_coords[which(total_mds_coords != as.character(correlated_mds[i,2]))] -> total_mds_coords
      }
    }
  }

  totmds_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_7_total_mds_coords.txt", sep= "")
  write.table(total_mds_coords, totmds_tn, sep="\t")


  #output final data much of this analysis is the same as above, but uses total_mds_coords

  #cm <- read_tsv("/home/564/km6006/Scripts/outliers/centromeres.tsv",col_names=c("chrom", "centromere"))

  cluster_genotypes <- tibble(mds_coord=character(), name=character(), 
                              PC1=numeric(), PC2=numeric(), genotype=character())
  mds_info <- tibble(mds_coord=character(), mds=character(), chromosome=character(), 
                    start=numeric(), end=numeric(), n_outliers=numeric(),
                    PC1_perc=numeric(), PC2_perc=numeric(), 
                    betweenSS_perc=numeric())
  het_data <- tibble(inversion=character(), name=character(), 
                              PC1=numeric(), PC2=numeric(), genotype=character(), het=numeric())

  for (i in 1:nrow(sig_mds_clusters)){
    coord <- pull(sig_mds_clusters[i,1])
    direction <- pull(sig_mds_clusters[i,2])
    if(! paste(coord,"-", direction, sep="") %in% total_mds_coords){
      next;
    }
    high_cutoff <- pull(sig_mds_clusters[i,6])
    low_cutoff <- pull(sig_mds_clusters[i,7])
    cluster_chr <- pull(sig_mds_clusters[i,8])
    coord_direction <- paste(coord, "-", direction, sep="")
    print(paste("Testing",coord_direction))
    
    if (direction == "pos"){
      current_windows <- win.regions %>% #all MDS 
        mutate_(the_mds=coord ) %>% 
        mutate(outlier=case_when((the_mds > high_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                                TRUE ~ "Non-outlier")) %>%
        filter(outlier != "Non-outlier") %>%
        select(chrom,start,end,mid,the_mds,outlier,n) %>%
        mutate(mds_coord=coord_direction) %>%
        mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
        mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
        filter(min_dist < max_distance_between_outliers ) %>%
        select(-ahead_n,-behind_n,-min_dist)
      windows <- current_windows %>% pull(n)
    }else{
      current_windows <- win.regions %>%
        mutate_(the_mds=coord) %>% 
        mutate(outlier=case_when((the_mds < high_cutoff) & (chrom == cluster_chr) ~ "Outlier", 
                                TRUE ~ "Non-outlier")) %>%
        filter(outlier != "Non-outlier") %>%
        select(chrom,start,end,mid,the_mds,outlier,n) %>%
        mutate(mds_coord=coord_direction) %>%
        mutate(ahead_n=n-lag(n), behind_n=abs(n-lead(n))) %>%
        mutate(min_dist=pmin(ahead_n, behind_n, na.rm=T)) %>%
        filter(min_dist < max_distance_between_outliers) %>%
        select(-ahead_n,-behind_n,-min_dist)
      windows <- current_windows %>% pull(n)
    }
    chromosome <- current_windows %>% head(1) %>% pull(chrom)
    start <- current_windows %>% summarize(s=min(start)) %>% pull()
    end <- current_windows %>% summarize(e=max(end)) %>% pull()
    
    # cm %>% 
    #   mutate(centromere=as.double(centromere)) %>%
    #   filter(chrom == cluster_chr) -> cent
    
    # output data tables 
    out <- cov_pca(win.fn.snp(windows), k=2)
    PC1_perc <- out[2]/out[1]
    PC2_perc <- out[3]/out[1]
    matrix.out <- t(matrix(out[4:length(out)], ncol=nrow(samples), byrow=T))
    out <- as_tibble(cbind(samples, matrix.out)) %>% 
      rename(name=name, PC1="1", PC2="2") %>% 
      mutate(PC1=as.double(PC1), PC2=as.double(PC2))
    try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
    if("try-error" %in% class(try_3_clusters)){
      kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
    }else{
      kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
    }
    out$cluster <- kmeans_cluster$cluster - 1
    out$cluster <- as.character(out$cluster)
    out$mds_coord <- paste(coord, direction, sep="_")
    betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss

    genotype.out <- out %>% select(mds_coord,name,PC1,PC2,cluster) %>% rename(genotype=cluster) 
    cluster_genotypes <- rbind(cluster_genotypes, genotype.out)
    tmp_info <- tibble(mds_coord=as.character(coord_direction), mds=as.character(toupper(coord)), chromosome=as.character(chromosome), 
                      start=as.numeric(start), end=as.numeric(end), n_outliers=as.numeric(length(windows)),
                      PC1_perc=as.numeric(round(PC1_perc*100,2)), PC2_perc=as.numeric(round(PC2_perc*100,2)),
                      betweenSS_perc=as.numeric(round(betweenSS_perc,4)))
    mds_info <- rbind(mds_info, tmp_info)

  # heterozygosity
    win.fn.snp(windows) %>% as.tibble() -> snps 
    colnames(snps) <- pull(samples) 
    snps %>% gather("name","genotype",1:ncol(snps)) %>% group_by(name, genotype) %>% 
      summarize(count=n()) %>% 
      spread(genotype, count) %>% 
      summarize(het=`1`/(`0` + `1` + `2`)) -> heterozygosity 
    
    inversion <- paste0(chromosome, ":", start, "-", end)
    het_out <- out %>% select(inversion,name,PC1,PC2,cluster) %>% rename(genotype=cluster) 
    heterozygosity <- inner_join(het_out, heterozygosity)
    het_data <- rbind(het_data, heterozygosity)
  }


  mdsinfo_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_8_mds_info.txt", sep= "")
  write.table(mds_info, mdsinfo_tn, sep="\t")
  het_tn <- paste("/g/data/ht96/McLay_UQ/inversion_paper/local_pca/", p1,p2, "_9_het_gen.txt", sep= "")
  write.table(het_data, het_tn, sep="\t")

}