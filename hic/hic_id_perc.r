library(data.table)
library(tidyverse)
library(grid)
library(gridExtra)
library(tools)


# read in the inversion data
inversions <- read.csv("/Users/.../inversions_all.csv")
#create a column for the HiC 1Mb start and end windows
inversions$start_window <- floor(inversions$start / 1000000) * 1000000
inversions$end_window <- floor(inversions$end / 1000000) * 1000000
data_directory <- "/Users/.../int_tables_1Mb"

# Initialize data frames
inv_genosets <- data.frame(
  inversion = integer(),
  scaffold = character(),
  start_window = integer(),
  end_window = integer(),
  `2` = character(),   
  `0` = character(),   
  stringsAsFactors = FALSE
)

perc_peaks <- tibble(
  chr2 = character(),
  bin2 = integer(),
  chr1 = character(),
  bin1 = integer(),
  inv_links = numeric(),
  genotype = character(),
  noninv_links = numeric(),
  genotype2 = character(),
  link_comparison = character(),
  abs_link_comparison = character(),
  distance = numeric(),
  percent_ranks = numeric(),
  inversion = character(),
)  

# Loop through each row of the data frame
for (inv in 1:nrow(inversions)) {
  row_values <- inversions[inv, 12:19]
  columns_2 <- colnames(inversions)[which(row_values == 2) + 11]
  columns_0 <- colnames(inversions)[which(row_values == 0) + 11]
  
  # Add a new row to the result data frame
  inv_genosets <- rbind(inv_genosets, data.frame(
    inversion = inversions[inv, 1],
    scaffold = inversions[inv, 3],
    start_window = inversions[inv, 20],
    end_window = inversions[inv, 21],
    `gen2` = ifelse(length(columns_2) > 0, paste(columns_2, collapse = " "), ""),
    `gen0` = ifelse(length(columns_0) > 0, paste(columns_0, collapse = " "), ""),
    stringsAsFactors = FALSE
  ))
}

#remove any rows that do not contain both individuals with 0 and 2 genotypes 
inv_genosets <- inv_genosets[inv_genosets$gen2 != "" & inv_genosets$gen0 != "", ]

# Loop through each row in inv_genosets
for (i in 1:nrow(inv_genosets)) {
  # Get the current scaffold, gen2, and gen0 values
  scaffold <- inv_genosets$scaffold[i]
  inversion <- inv_genosets$inversion[i]
  gen2_values <- unlist(strsplit(inv_genosets$gen2[i], " "))  # Split gen2 into individual values
  gen0_values <- unlist(strsplit(inv_genosets$gen0[i], " "))  # Split gen0 into individual values
  
  # Generate all possible combinations of gen2 and gen0 values
  combinations <- expand.grid(gen2 = gen2_values, gen0 = gen0_values, stringsAsFactors = FALSE)
  
  # Process each combination
  for (j in 1:nrow(combinations)) {
    # get HiC data for current genotype 2 individual
    gen2 <- combinations$gen2[j] 
    file2 <- paste0(gen2, "_", scaffold, "_hic_table.txt")
    data_pop2 <- setDT(read.table(file.path(data_directory, file2), header = TRUE))
    colnames(data_pop2) <- c("chr2", "bin2", "chr1", "bin1", "inv_links", "genotype")
    # get HiC data for current genotype 0 individual
    gen0 <- combinations$gen0[j]
    file0 <- paste0(gen0, "_", scaffold, "_hic_table.txt")
    data_pop0 <- setDT(read.table(file.path(data_directory, file0), header = TRUE))
    colnames(data_pop0) <- c("chr2", "bin2", "chr1", "bin1", "noninv_links", "genotype")
    # Combine the data from the two tables, and calculate link comparisons (difference in links between pop1 and pop2)
    data_pop2[,noninv_links:=data_pop0[,noninv_links]]
    data_pop2[,genotype2:=data_pop0[,genotype]] 
    data_pop2[,link_comparison:=inv_links - noninv_links] 
    data_pop2$abs_link_comparison <- abs(data_pop2$link_comparison) #absolute difference in interactions between the two individuals for ranking with percent_rank() which ranks highest to lowest 
    data_pop2$inversion <- inversion
    write.table(data_pop2, paste("/Users/.../pair_tables/", inversion, "_", gen2, "_", gen0, ".txt", sep = ""), sep = "\t", row.names = FALSE)
  
    data_pop2 %>%  
      filter(bin1 >= bin2) %>%
      mutate(distance = bin1 - bin2) %>%
      group_by(distance) %>%
      mutate(percent_ranks = percent_rank(abs_link_comparison)) %>%
      ungroup() -> pop_pair
    
    peak <- subset(pop_pair, bin2 == inv_genosets$start_window[i] & bin1==inv_genosets$end_window[i]) 
    peak$perc <- floor(peak$percent_ranks*100)
    perc_peaks <- rbind(perc_peaks, peak)
    }
}
write.csv(perc_peaks, "/Users/.../hi_interaction_perc_rank.csv")


