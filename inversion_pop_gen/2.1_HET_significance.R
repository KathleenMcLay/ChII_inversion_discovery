## Mann-W U test 
het_data <- read.csv("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/data/Het/het_data_all_inversions.csv")
het_data$genotype <- as.character(het_data$genotype)
het_data$genotype <- trimws(het_data$genotype, whitespace = '"')

man_u01 <- data.frame(pop=character(), inversion=character(), group=character(), p_val=numeric()) 
man_u12 <- data.frame(pop=character(), inversion=character(), group=character(), p_val=numeric()) 

unique_inversions <- unique(het_data$inversion) # get unique inversions

#0-1 genotype comparison 
for (inv in unique_inversions) {
  current_inv <- subset(het_data, inversion == inv) # subset the data to current inversion
  current_pop <- unique(current_inv$pop)
  het01 <- subset(current_inv, genotype %in% c('0','1')) # subset the data again to only include rows with 0 or 1 in the genotype column 
  if (sum(het01$genotype == '0') == 0 || sum(het01$genotype == '1') == 0) {
    next # skip to next iteration if either '0' or '1' is missing in the subset
  }
  test <- wilcox.test(het01$het ~ het01$genotype) # perform wilcox.test
  p_val <- test$p.value # get p-value for wilcox 
  new_01<- data.frame(pop = current_pop, inversion = inv, group = "0-1", p_val = p_val)
  man_u01 <- rbind(man_u01, new_01) #rbind current result 
}

#1-2 genotype comparison 
for (inv in unique_inversions) {
  current_inv <- subset(het_data, inversion == inv) # subset the data to current inversion
  current_pop <- unique(current_inv$pop)
  het12 <- subset(current_inv, genotype %in% c('1','2')) # subset the data again to only include rows with 0 or 1 in the genotype column 
  if (sum(het12$genotype == '1') == 0 || sum(het12$genotype == '2') == 0) {
    next # skip to next iteration if either '0' or '1' is missing in the subset
  }
  
  test <- wilcox.test(het12$het ~ het12$genotype) # perform wilcox.test
  p_val <- test$p.value # get p-value for wilcox 
  new_12 <- data.frame(pop = current_pop, inversion = inv, group = "1-2", p_val = p_val)
  man_u12 <- rbind(man_u12, new_12) #rbind current result 
}

het_sig <- rbind(man_u01, man_u12)
write.csv(het_sig, "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /3_population_level_inversions/data/Het/het_significance.csv", sep = "\t")
