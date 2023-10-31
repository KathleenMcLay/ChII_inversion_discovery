### Density plots for VariantFiltration annotations

library(data.table)
library(ggplot2)
library(doParallel)

annfiles <- list.files(path = "/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only/ann_tables", pattern="merged", all.files=FALSE)

registerDoParallel(4)

# Loop through the .txt files and read data into variables
foreach (i=1:length(annfiles)) %dopar% {
  variable_name <- gsub(".txt", "", sub(".+_([^_]+)\\.txt", "\\1", annfiles[i]))
  data <- fread(paste("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only/ann_tables/", annfiles[i], sep = ""))
  plot_name=paste("/g/data/ht96/McLay_UQ/inversion_paper/6_SNP_filtering/variant_only/ann_tables/", variable_name, "_LC_VO.jpeg", sep = "")
  jpeg(file=plot_name)
  ggplot(data, aes(x=V1)) + geom_density()
  dev.off()
}