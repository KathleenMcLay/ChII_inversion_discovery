# kmeans clustering of all SNPs from putative inversion region to identify genotype clustering 
library(cluster)  # For silhouette analysis
library(dplyr)  
library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggplot2)

inversion_list <- read.csv("/QRISdata/Q6656/chapter_II/new_inv_discovery/2_regions_summary_final.csv", header=TRUE)

# Create a loop to run PCA on each inversion region
for (row in 1:nrow(inversion_list)) {
    # create a variable to filter the vcf
    region <- inversion_list[row, "vcf_filter_alt"]
    # create a variable for the population name/s
    pop <- inversion_list[row, "population"]

    # create the correct population prefix to create a sample list from the vcf for populations with two or more subpopulations
    if (nchar(pop) > 3) { 
        pop_filter <- paste0(substr(pop, 1, 3), "|", substr(pop, 4, 6))
        } else {
        pop_filter <- pop  # Keep original if <= 2 characters
    }
    print(paste("CURRENT POP FILTER VALUE IS:", pop_filter, " for inversion", region))
    # extract a list of samples for the current population

    extract <- paste(
    "bcftools query -l /QRISdata/Q6656/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz | grep -E \"^",
    pop_filter,
    "\" > /QRISdata/Q6656/chapter_II/new_inv_discovery/alt/",
    pop, "_sample_list.txt",
    sep = ""
    )
    system(extract)
    # filter the vcf 
    filter <- paste("bcftools view --threads 12 -O z -S /QRISdata/Q6656/chapter_II/new_inv_discovery/alt/", pop, "_sample_list.txt -r ", region, " -o /QRISdata/Q6656/chapter_II/new_inv_discovery/alt/inv_vcfs/", pop, "_", region, ".vcf.gz /QRISdata/Q6656/sf8_no_variants_noD1_reheader_final_nsng.vcf.gz", sep = "")
    system(filter)
    # index the vcf
    index <- paste("bcftools index --threads 24 /QRISdata/Q6656/chapter_II/new_inv_discovery/alt/inv_vcfs/", pop, "_", region, ".vcf.gz", sep = "")
    system(index)

    # read in list of samples for the population 
    samples <- read.table(paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/", pop, "_sample_list.txt", sep=""), header=FALSE)

    # read in the vcf for the current inversion
    vcf_fn <- paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/inv_vcfs/", pop, "_", region, ".vcf.gz", sep="")
    # Convert VCF file
    snpgdsVCF2GDS(vcf_fn, paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/pca/", pop, "_", region, "_senecio.gds", sep=""))
    # Assign the new file to an object 
    genofile <- snpgdsOpen(paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/pca/", pop, "_", region, "_senecio.gds", sep=""))
    # Run PCA 
    senecio_pca <- snpgdsPCA(genofile, autosome.only=FALSE)

    # Calcualte PVE% for PC1 and PC2
    PC1_perc <- head(round(senecio_pca$varprop[1]*100, 2))
    PC2_perc <- head(round(senecio_pca$varprop[2]*100, 2))

    # Create matrix to compute kmean clusters 
    PC1 <- senecio_pca$eigenvect[,1]
    PC2 <- senecio_pca$eigenvect[,2]
    PCA <- cbind(PC1, PC2)
    colnames(PCA) <- c("PC1", "PC2")
    PCA <- as.matrix(PCA) 

    # Compute kmeans - makes either 3 or 2 genotype clusters
    try_3_clusters <-try(kmeans(PCA[,1], 3, centers=c(min(PCA[,1]), (min(PCA[,1])+max(PCA[,1]))/2, max(PCA[,1]))))
    if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(PCA[,1], 2, centers=c(min(PCA[,1]), max(PCA[,1])))
    }else{
    kmeans_cluster <- kmeans(PCA[,1], 3, centers=c(min(PCA[,1]), (min(PCA[,1])+max(PCA[,1]))/2, max(PCA[,1])))
    }

    # Calculate silhouette score
    sil <- silhouette(kmeans_cluster$cluster, dist(PCA[, c("PC1", "PC2")]))
    mean_sil_score <- mean(sil[, 3])

    # Make genotype cluster dataset   
    cluster_genotypes <- as_tibble(cbind(samples, PCA)) %>% 
    mutate(PC1=as.double(PC1), PC2=as.double(PC2))
    cluster_genotypes$cluster <- kmeans_cluster$cluster - 1
    cluster_genotypes$cluster <- as.character(cluster_genotypes$cluster)
    cluster_name <- paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/pca/", pop, "_", region, "_PCA_genotypes.txt", sep= "")
    write.table(cluster_genotypes, cluster_name, sep="\t")

    # Make percentage stats dataset 
    betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss
    stat_table <- tibble(PC1_perc=as.numeric(PC1_perc), PC2_perc=as.numeric(PC2_perc),
                        betweenSS_perc=as.numeric(round(betweenSS_perc,4)), mean_sil_score=as.numeric(round(mean_sil_score,4)))
    stat_name <- paste("/QRISdata/Q6656/chapter_II/new_inv_discovery/alt/pca/", pop, "_", region, "_PCA_stats.txt", sep= "")
    write.table(stat_table, stat_name, sep="\t")

    snpgdsClose(genofile)
}
dev.off()