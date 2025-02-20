library(SNPRelate)
library(data.table)
library(stringr)
#library(ggplot2)

#read in the file containing the inversions(chr:start-end) [,1] and the population in which they were discovered (P01P02) [,2]
inversions <- read.csv("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/4_regions_results_summary.csv", header = FALSE)

#function to perform kmeans clustering on new samples
perform_kmeans <- function(data, centers_list) {
  # Predefined identifiers for each center
  center_ids <- list(
    c(0, 1, 2),  # Identifiers for 3 clusters (min, midpoint, max)
    c(0, 2),     # Identifiers for 2 clusters (min, max)
    c(0, 1),     # Identifiers for 2 clusters (min, midpoint)
    c(1, 2)      # Identifiers for 2 clusters (midpoint, max)
  )
  
  for (i in seq_along(centers_list)) {
    centers <- centers_list[[i]]
    result <- try(kmeans(data, length(centers), centers = centers), silent = TRUE)
    if (!inherits(result, "try-error") && all(result$size > 0)) {
      # Map cluster assignments to predefined identifiers
      mapping <- order(centers)  # Order centers to determine their identifiers
      result$cluster <- sapply(result$cluster, function(cluster) center_ids[[i]][mapping[cluster]])
      return(result)
    }
  }
  
  # If all attempts fail, return NULL and exit the function gracefully
  return(NULL)
}

#open a pdf to save the plots 
#pdf("/scratch/user/uqkmcla4/PCA_projection/PCA_plots.pdf", height=7, width=7) 

for (i in 1:nrow(inversions)) {
    ### Create input file for current inversion 
    #isolate the current inversions coordinates chr:start-end
    current_inv <- inversions[i,1]
    print(current_inv)
    #subset the vcf.gz file to the current inversion with all populations
    # x <- paste("bcftools view --threads 12 -O z -r ",current_inv, " -o /scratch/user/uqkmcla4/PCA_projection/",current_inv,".vcf.gz /scratch/user/uqkmcla4/sf8_SNPs_reheader_nsng.vcf.gz", sep="")
    # system(x)
    # index <- paste("bcftools index --threads 12 /scratch/user/uqkmcla4/PCA_projection/",current_inv,".vcf.gz", sep = "")
    # system(index)
    
    ### Principal component analysis - using the samples from the inversion discovery population
    #create a vector of the discovery populations  
    dis_pop <- unlist(str_extract_all(inversions[i, 2], "[A-Z][0-9]{2}"))
    # Create a regular expression to filter the samples to those from the discovery populations
    dis_pop_fil <- paste0("^(", paste(dis_pop, collapse = "|"), ")")
    print(paste("Discovery population: ", dis_pop_fil))

    #current inversion vcf.gz file
    vcf_fn <- paste("/scratch/user/uqkmcla4/PCA_projection/",current_inv,".vcf.gz", sep = "")
    #convert VCF file to GDS 
    gds_file <- paste("/scratch/user/uqkmcla4/PCA_projection/",current_inv,".gds", sep = "")
    GDS <- snpgdsVCF2GDS(vcf_fn, gds_file) 
    #assign the new file to an open object 
    genofile <- snpgdsOpen(GDS)
    #create a vector for all samples in the vcf file  
    all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
    print("start of all_samples")
    print(all_samples)
    #create a vector of just samples from the discovery population
    PCA_samples <- all_samples[grep(dis_pop_fil, all_samples)]
    print("start of PCA_samples")
    print(PCA_samples)
    #run PCA 
    senecio_pca <- snpgdsPCA(genofile, autosome.only=FALSE, sample.id=PCA_samples)
    #create matrix of eigenvectors for original samples - to use as centre for calculating new 
    OG_PCA <- senecio_pca$eigenvect[, 1:2]
    #extract SNP loadings 
    SnpLoadings <- snpgdsPCASNPLoading(senecio_pca, genofile)
    
    ### Calculate the eigenvectors for each of the remaining samples
    # specify samples to be projected 
    proj_samples <- all_samples
    #make a dataframe of the samples
    proj_samples_df <- data.frame(Samples = proj_samples)
    print("start of proj_samples_df")
    print(head(proj_samples_df))
    # calculate the eigenvectors 
    proj_EV <-  snpgdsPCASampLoading(SnpLoadings, genofile, sample.id=proj_samples)
    print("start of proj_EV")
    print(proj_EV)
    
    ### Genotype remaining samples 
    #create matrix of eigenvectors to compute kmean clusters
    PCA <- proj_EV$eigenvect[, 1:2]
    colnames(PCA) <- c("PC1", "PC2")
    print("start of PCA")
    print(head(PCA))
    matrixPCA <- as.matrix(PCA) 
    print("start of matrixPCA")
    print(head(matrixPCA))
    
    ###kmeans clustering to identify genotype groups - get the three cluster centres from the original PCA sample eigenvectors 
    # Define the fallback centers
    centers_list <- list(
      c(min(OG_PCA[, 1]), (min(OG_PCA[, 1]) + max(OG_PCA[, 1])) / 2, max(OG_PCA[, 1])), # 3 clusters
      c(min(OG_PCA[, 1]), max(OG_PCA[, 1])),                                            # 2 clusters (min, max)
      c(min(OG_PCA[, 1]), (min(OG_PCA[, 1]) + max(OG_PCA[, 1])) / 2),                   # 2 clusters (min, midpoint)
      c((min(OG_PCA[, 1]) + max(OG_PCA[, 1])) / 2, max(OG_PCA[, 1]))                    # 2 clusters (midpoint, max)
    )
    # Perform k-means clustering
    kmeans_cluster <- perform_kmeans(matrixPCA[, 1], centers_list)
    
    #output genotype and eigenvector values for new samples 
    PCA <- cbind(proj_samples_df, PCA)
    #check if the kmeans_cluster has assigned genotypes
    if (!is.null(kmeans_cluster)) {
      PCA$genotype <- kmeans_cluster$cluster
    }
    write.table(PCA, paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/pca_projection/new_genotypes_", current_inv, ".txt", sep = ""))

    ### Plot PCA of projected samples 
    #Set up variable for plot theme settings 
    #p_theme <- theme(axis.title.x = element_text(margin = margin(t = 20), size = 12), 
    #             axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 12),
    #             plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
    #             axis.ticks = element_line(colour="black"), panel.background=element_rect(colour="black"), panel.border=element_rect(colour="black", size=0.75),
    #             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #             legend.position = "none") 

    #make plot! 
    #PCA_plot <- ggplot(PCA, aes(x=PC1, y=PC2, col=genotype)) + geom_point(size = 2) + theme_bw() +
    #scale_color_manual(name="Genotype", values=c("red","purple","blue")) +
    #ggtitle(current_inv, " projected samples") + xlab("PC1") + ylab("PC2") +
    #p_theme   
    #print(PCA_plot)
    
    #close the file 
    snpgdsClose(genofile)
}
#dev.off()