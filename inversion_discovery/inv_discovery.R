library(cluster) 
library(dplyr)    

# Load data into a dataframe `df` with columns: chrom, position, mds01, mds02, window_index
df <- read.csv("MDS_win_regions_data_all.csv", header = TRUE)
df <- df[, c("chrom", "start", "end", "mid", "mds01", "mds02", "n", "population")]

# Get unique populations
pops <- unique(df$population)

# Initialize an empty dataframe to store results
selected_regions <- data.frame()

# Loop over each population
for (pop in pops) {
  df_pop <- df[df$population == pop, ]
  
  # Get unique chromosomes
  chromosomes <- unique(df_pop$chrom)
  
  for (chrom in chromosomes) {
    df_chrom <- df_pop[df_pop$chrom == chrom, ] 
    
    # Ensure there are enough points for clustering
    if (nrow(df_chrom) < 10) {
      message(paste("Skipping", chrom, "in population", pop, "- not enough windows for clustering."))
      next
    }
    
    # Determine the best k using silhouette score
    set.seed(123)  # For reproducibility
    sil_scores <- numeric()
    
    for (k in 2:10) {
      km <- kmeans(df_chrom[, c("mds01", "mds02")], centers = k, nstart = 25)
      sil <- silhouette(km$cluster, dist(df_chrom[, c("mds01", "mds02")]))
      sil_scores[k] <- mean(sil[, 3])  # Store mean silhouette width
    }
    
    best_k <- which.max(sil_scores)  # Get k with highest silhouette score
    
    # Perform k-means clustering with optimal k
    final_kmeans <- kmeans(df_chrom[, c("mds01", "mds02")], centers = best_k, nstart = 25)
    df_chrom$cluster <- final_kmeans$cluster  # Assign clusters
    
    # Compute z-scores for MDS1
    df_chrom$MDS1_z <- scale(df_chrom$mds01)
    
    # Identify high z-score windows
    df_chrom$high_z <- df_chrom$MDS1_z > 1.5
    
    # Detect consecutive windows in the same cluster with high z-score
    df_chrom_TRUE <- df_chrom %>%
      arrange(n) %>%  # Ensure sorted order
      group_by(cluster) %>%
      mutate(run_id = cumsum(c(1, diff(n) != 1))) %>%  # Identify consecutive runs
      group_by(cluster, run_id) %>%
      filter(high_z) %>%
      filter(n() >= 10)  %>%  # Require at least 10 consecutive windows
      ungroup()
    
    # Store results
    if (nrow(df_chrom_TRUE) > 0) {
      df_chrom_TRUE$population <- pop  # Retain population info
      selected_regions <- bind_rows(selected_regions, df_chrom_TRUE)
    }
  }
}

# Output results
print(selected_regions)
write.csv(selected_regions, file = "selected_regions.csv", row.names = FALSE)

inversion_summary <- selected_regions %>%
  group_by(population, chrom, run_id) %>%  # Group by population, chromosome, and unique inversion (run_id)
  summarise(
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

inversion_summary <- inversion_summary %>%
  mutate(
    length = (end - start + 1)/1000000 # Calculate the length of each inversion
  )
write.csv(inversion_summary, file = "inversion_summary.csv", row.names = FALSE)


has_consecutive_10 <- function(n_values) {
  diffs <- diff(n_values)  # Compute differences between consecutive n values
  consec_streaks <- rle(diffs == 1)  # Identify runs where n increases by exactly 1
  any(consec_streaks$lengths[consec_streaks$values] >= 9)  # Check if any run has at least 9 differences (10 consecutive values)
}

# Filtering dataset to keep all rows for a run_id if it has at least 10 consecutive n values
filtered_data <- selected_regions %>%
  group_by(run_id) %>%
  filter(has_consecutive_10(n)) %>%
  ungroup()

filter_summary <- filtered_data %>%
  group_by(population, chrom, run_id) %>%  # Group by population, chromosome, and unique inversion (run_id)
  summarise(
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

filter_summary <- filter_summary %>%
  mutate(
    length = (end - start + 1)/1000000  # Calculate the length of each inversion
  )
write.csv(filter_summary, file = "region_summary_allwindows.csv", row.names = FALSE)

