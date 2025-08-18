library(data.table)
library(tidyverse)
library(dplyr)

# Read in the file containing the 'inversion' details, create inversion id 
inv_data <- read.csv("4_inversion_discovery/final_inversion_list.csv", header=TRUE)

# Directory for hic paired interaction tables 
hic_files <- list.files(path = "6_inversion_hic_interactions/HiC_int_tables_1Mb/pair_tables", pattern=".txt", all.files=FALSE)

for (i in 1:length(hic_files)) {
  pop_pair <- read.table(paste("6_inversion_hic_interactions/HiC_int_tables_1Mb/pair_tables/", hic_files[i], sep = ""), header=TRUE)
  
  window_size <- 1000000 # 1 MBP WINDOWS - each diamond on the graph is 1Mbp wide
  chosen_inversion <- unique(pop_pair$inversion)
  chosen_chr <-  unique(pop_pair$chr1)
  
  #filter inv-dat # to only include the chosen inversion
  inv_d <- inv_data %>%
    filter(inversion == chosen_inversion)
  
  chosen_start <- inv_d %>% pull(inversion_start_bp) %>% min() 
  chosen_end <- inv_d %>% pull(inversion_end_bp) %>% max() 
  chosen_middle = ((floor(chosen_start/1000000) + floor(chosen_end/1000000))/2)*1000000 #FIND INVERSION MIDPOINT 
  chosen_height = ((floor(chosen_end/1000000) - floor(chosen_start/1000000))/2) #FIND INVERSION HEIGHT (above midpoint)
  inv_code <- inv_d %>% pull(inversion_code)
  ind1 <- pop_pair %>% pull(genotype) %>% unique()
  ind2 <- pop_pair %>% pull(genotype2) %>% unique()
  ylab <- paste(inv_code, " comparision - ",ind1, ":", ind2, sep = "")
  scaff <- inv_d %>% pull(scaffold) #replace _ with space
  scaff <- gsub("_", " 0", scaff) #replace _ with space
  
  #Set the beginning and end of the plot as the start and end of the scaffold
  view_start <- 0 #inv_d %>% pull(scaffold_start_bp)
  view_end <- inv_d %>% pull(scaffold_end_bp)
  
  #CALCULATE PERCENTAGE RANK, PREPARE DATA FOR PLOTTING 
  plot_dat <- as_tibble(pop_pair) %>%
    filter(chr2 == chosen_chr & chr1 == chosen_chr & 
             bin2 >= view_start & bin2 <= view_end &
             bin1 >= view_start & bin1 <= view_end) %>% 
    filter(bin1 >= bin2) %>%
    mutate(distance = bin1 - bin2) %>%
    group_by(distance) %>%
    mutate(percent_ranks = percent_rank(link_comparison)) %>%
    ungroup() %>%
    mutate(link_comparison = case_when(bin1 == bin2 ~ 99, #all bins that are the same are marked as 99
                                       link_comparison > 0.3 ~ 0.3, #all links > 0.3 are converted to 0.3
                                       link_comparison < -0.3 ~ -0.3, #all links < -0.3 are converted to 0.3
                                       link_comparison == 0 ~ 0.00001, #all links that are 0, are changed to 0.00001, presumably to allow for plotting
                                       TRUE ~ link_comparison)) %>% 
    mutate(link_comparison = na_if(link_comparison, link_comparison == 99))  %>% #bins that are the same are excluded (shown as grey triangles on plot)
    mutate(x = (bin1 + bin2)/2,y=(bin1-bin2)/2000000) %>% #CREATE X and Y variables 
    group_by(y,x,link_comparison,percent_ranks, distance) %>%
    expand(count = seq(1:4)) %>% #THIS STUFF IS TO USE GEOM_POLYGON! - it has four sides! 
    mutate(bin_id = paste0(x,"_",y)) %>%
    ungroup() %>%
    mutate(y = case_when(count == 2 ~ y+0.5,
                         count == 4 ~ y - 0.5,
                         TRUE ~ y),
           x = case_when(count == 1 ~ x -(window_size/2),
                         count == 3 ~ x + (window_size/2),
                         TRUE ~ x)) %>%
    mutate(y = case_when(y < 0 ~ 0,
                         TRUE ~ y)) 
  
  # Make the plot 
  plot_name <- paste("6_inversion_hic_interactions/figures/", tools::file_path_sans_ext(hic_files[i]), "_HiCleg.tiff", sep = "")
  tiff(file=plot_name, units="in", width=8, height=5, res=300)
  print(ggplot(plot_dat, aes()) + geom_polygon(aes(x=x/1000000,y=y,fill=link_comparison,group = bin_id),size=0) +
          scale_fill_gradient2(low = "red", mid = "white",
                               high = "#39568CFF", midpoint = 0,name="HiC interaction\ndifference",limits=c(-0.3,0.3)) +
          annotate("segment",x=floor(chosen_start/1000000),y=0,
                   yend=chosen_height,xend=chosen_middle/1000000,
                   alpha=1,color="black",size=0.5) +
          annotate("segment",x=floor(chosen_end/1000000),y=0,
                   yend=chosen_height,xend=chosen_middle/1000000,
                   alpha=1,color="black",size=0.5) +
          theme_linedraw() + ylab(ylab) + xlab(paste(scaff, "position (Mb)")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_blank(),
                panel.border = element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()
                ))
  dev.off()
}
