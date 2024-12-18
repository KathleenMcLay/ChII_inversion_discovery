# kmeans clustering of all SNPs from putative inversion region to identify genotype clustering 

library(tidyverse)
library(SNPRelate)
library(Matrix)
library(ggrepel)

# list of inversion regions vcfs 
pca.files <- list.files(path = "/QRISdata/Q6656/new_lpca/inversion_vcfs", pattern=".vcf.gz", all.files=FALSE)
sample.files <- list.files(path = "/QRISdata/Q6656/new_lpca/sample_lists", pattern=".txt", all.files=FALSE) #samples must be in same order as samples in vcf 

pdf("/QRISdata/Q6656/new_lpca/pca_out/PCA_plots.pdf", height=7, width=7) 

for (j in 1:length(pca.files)) {
    pop_base <- sub("\\.\\w+$", "", pca.files[j])
    pop_base <- sub("_.+$", "", pop_base)
    print(pop_base)
    for (i in 1:length(sample.files)) {
        sample_base <- sub("\\.\\w+$", "", sample.files[i])  
        if (pop_base == sample_base) {
        print(sample_base)
        samples <- read.table(paste("/QRISdata/Q6656/new_lpca/sample_lists/", sample.files[i], sep=""), header=FALSE) 
        colnames(samples) <- "sample"
        }
    }
    vcf_fn <- paste("/QRISdata/Q6656/new_lpca/inversion_vcfs/", pca.files[j], sep="")

    # Convert VCF file 
    snpgdsVCF2GDS(vcf_fn, "/QRISdata/Q6656/new_lpca/pca_out/senecio.gds") 
    # Assign the new file to an object 
    genofile <- snpgdsOpen("/QRISdata/Q6656/new_lpca/pca_out/senecio.gds")
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

    # Compute kmeans - makes either 3 or 2 genotype clusters, 3 clusters consistent with inversion
    try_3_clusters <-try(kmeans(PCA[,1], 3, centers=c(min(PCA[,1]), (min(PCA[,1])+max(PCA[,1]))/2, max(PCA[,1]))))
    if("try-error" %in% class(try_3_clusters)){
    kmeans_cluster <-kmeans(PCA[,1], 2, centers=c(min(PCA[,1]), max(PCA[,1])))
    }else{
    kmeans_cluster <- kmeans(PCA[,1], 3, centers=c(min(PCA[,1]), (min(PCA[,1])+max(PCA[,1]))/2, max(PCA[,1])))
    }

    # Make genotype cluster dataset   
    cluster_genotypes <- as_tibble(cbind(samples, PCA)) %>% 
    mutate(PC1=as.double(PC1), PC2=as.double(PC2))
    cluster_genotypes$cluster <- kmeans_cluster$cluster - 1
    cluster_genotypes$cluster <- as.character(cluster_genotypes$cluster)
    cluster_name <- paste("/QRISdata/Q6656/new_lpca/pca_out/", tools::file_path_sans_ext(pca.files[j]), "_PCA_genotypes.txt", sep= "")
    write.table(cluster_genotypes, cluster_name, sep="\t")

    # Make percentage stats dataset 
    betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss
    stat_table <- tibble(PC1_perc=as.numeric(PC1_perc), PC2_perc=as.numeric(PC2_perc),
                        betweenSS_perc=as.numeric(round(betweenSS_perc,4)))
    stat_name <- paste("/QRISdata/Q6656/new_lpca/pca_out/", tools::file_path_sans_ext(pca.files[j]), "_PCA_stats.txt", sep= "")
    write.table(stat_table, stat_name, sep="\t")

    # PCA Plot 
    title <- tools::file_path_sans_ext(pca.files[j])
    PCA_plot <- ggplot(cluster_genotypes, aes(x=PC1, y=PC2, col=cluster)) + geom_point(size = 2) + theme_bw() +
    scale_color_manual(name="Cluster", values=c("red","purple","blue")) +
    ggtitle(title) +
    xlab(paste("PC1 (",PC1_perc, "%)")) +
    ylab(paste("PC2 (",PC2_perc, "%)")) +
    theme(legend.title = element_text(size=18, face="bold"), 
          legend.text=element_text(size=15),
          plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"), 
          axis.title.x = element_text(margin = margin(t = 20), face="bold", size = 12), 
          axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), face="bold", size = 12),
          plot.title = element_text(vjust = 4, size = 10), axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          axis.text = element_text(size = 12),
          legend.position = "none")    
    print(PCA_plot)
    snpgdsClose(genofile)
    file.remove("/QRISdata/Q6656/new_lpca/pca_out/senecio.gds")
}
dev.off()