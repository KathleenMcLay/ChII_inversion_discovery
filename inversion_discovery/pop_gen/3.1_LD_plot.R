LD_files <- list.files(path = "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /wholegen_LD",pattern=".txt", all.files=FALSE)
chromosomes <- read.csv("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /wholegen_LD/chromosomes.csv")

pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /wholegen_LD/wholegen_LD_plots.pdf", height=7, width=7)
for (i in 1:nrow(chromosomes)) {
  chosen_chr <- chromosomes$chr[i]
  print(chosen_chr)
  for (j in 1:length(LD_files)){
    if (tools::file_path_sans_ext(LD_files[j]) == chosen_chr) {
      print(paste0("reading file for: ", LD_files[j]))
      ld_chrom <- read.table(paste("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /wholegen_LD/",LD_files[j], sep=""), sep="\t",header=T,stringsAsFactors=FALSE)
    }
  }
ld_plot <- ggplot(ld_chrom,aes(x=win1/1000000,y=win2/1000000)) + theme_classic() + 
  geom_tile(aes(fill=max_2_r2)) +
  #geom_tile(data=ld_mds, aes(x=win2/1000000,y=win1/1000000,fill=max_2_r2)) +
  scale_fill_gradientn(colours=c("grey95","blue","red"), values=c(0,0.5,1), name="LD") +
  #geom_segment(mapping=aes(x=start/1000000,xend=end/1000000,y=-4,yend=-4), col="purple",size = 2.5) +
  #geom_segment(mapping=aes(x=-4,xend=-4,y=start/1000000,yend=end/1000000), col="purple",size = 2.5) +
  scale_x_continuous(expand=c(0.02,0)) +
  scale_y_continuous(expand=c(0.02,0)) +
  coord_fixed(ratio = 1) +
  xlab("Mbp") + ylab("Mbp") +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches")) + ggtitle(chosen_chr)
print(ld_plot)
}
dev.off()