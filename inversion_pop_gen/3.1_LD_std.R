ld_data <- read.csv("new_LD.csv")
ld_data <- drop_na(ld_data) 
ld_cg_data <- read.csv("LD_cmn_gen_data.csv")
ld_cg_data <- drop_na(ld_cg_data)

unique_inversions <- unique(ld_cg_data$inversion) # get unique inversions

### compare mean LD for the inversion between cmn genotype and all samples 
mean_r2 <- mean(ld_cg_data[,5])
sd_r2 <- sd(ld_cg_data[,5])
  
# Normalize the data
for (i in 1:nrow(ld_cg_data)) {
  ld_cg_data$std_mean_r2[i] <- (ld_cg_data[i,5] - mean_r2) / sd_r2
  print(i)
}
 
result_df <- data.frame(population = character(), inv = character(), mean_all_inv = numeric(), mean_cmn_inv = numeric(), check = character())


for (inv in unique_inversions) {
  ld_dat <- subset(ld_data, inversion == inv) # subset the data to current inversion
  ld_cg_dat <- subset(ld_cg_data, inversion == inv) 
  current_pop <- unique(ld_dat$pop)
  split_string <- unlist(strsplit(inv, "[:-]"))
  inv_start <- as.numeric(split_string[2])
  inv_end <- as.numeric(split_string[3])
  
  # subset the data to the inversion
  win1_cg_inv <- which(ld_cg_dat$win1 >= inv_start & lag(ld_cg_dat$win1) <= inv_end) 
  sub_dat_cg <- subset(ld_cg_dat[win1_cg_inv,]) 
  sub_dat_cg <- sub_dat_cg[order(sub_dat_cg$win2, decreasing = FALSE),] 
  win2_cg_inv <- which(sub_dat_cg$win2 <= inv_end) 
  sub_dat_cg <- subset(sub_dat_cg[win2_cg_inv,])
  
  win1_inv <- which(ld_dat$win1 >= inv_start & lag(ld_dat$win1) <= inv_end) 
  sub_dat <- subset(ld_dat[win1_inv,]) 
  sub_dat <- sub_dat[order(sub_dat$win2, decreasing = FALSE),] 
  win2_inv <- which(sub_dat$win2 <= inv_end) 
  sub_dat <- subset(sub_dat[win2_inv,])
  
  if (nrow(sub_dat) == 0 || nrow(sub_dat_cg) == 0) {
    print("no windows")
  } else {
    # Calculate means
    mean_all_inv <- mean(sub_dat[, 5])
    mean_cmn_inv <- mean(sub_dat_cg[, 5])
    
    if (mean_cmn_inv <= mean_inv) {
      check <- "pass"
    } else {
      check <- "fail"
    }
    
    result_df <- rbind(result_df, data.frame(population = current_pop, inv = inv, mean_all_inv = mean_all_inv, mean_cmn_inv = mean_cmn_inv, check = check))
  }
  
}