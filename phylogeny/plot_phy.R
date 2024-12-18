### Visualise the phylogeny! 

library(devtools)
install_github("https://github.com/YuLab-SMU/ggtree")

library(tidyverse)
library(ggtree)
library(treeio)
library(ggrepel)
library(phangorn)


raxml_file <- "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /2_phylogeny_pop_structure/data/WAout/RAxML_bipartitionsBranchLabels.ML"
raxml <- read.raxml(raxml_file)

# create annotation details for populations of each ecotype ie. Dunes, Headlands, Alpines 
D_dat <- data.frame(
  node = c(333, 385, 394, 412, 475, 533, 524, 2, 1, 306),
  name = c("D01 Lennox Head (NSW)", "D03 Cabarita Beach (QLD)", "D00 Stradbroke Island (QLD)", "D04 Coffs Harbour (NSW)", "D05 Hat Head (NSW)", "D35 Flinders Peninsula (WA)", "D09 Leeuwin-Naturaliste National Park (WA)", "D32-254 Discovery Bay (VIC)", "D32-288 Discovery Bay (VIC)", "D12 Bermagui Dune (NSW)")
)

H_dat <- data.frame(
  node = c(366, 403, 430, 320, 442, 484, 493, 502, 276, 521, 511),
  name = c("H02 Cabarita Beach (QLD)", "H00 Stradbroke Island (QLD)", "H05 Coffs Harbour (NSW)","H04 Byron Bay (NSW)", "H01 Lennox Head (NSW)", "H06 Hat Head (NSW)", "H07 Port Macquarie (NSW)", "H03 Kiama Blowhole (NSW)", "H12 Portland, Cape Bridgewater (VIC)", "H15 Port Aurthur (TAS)", "H14 Green Cape (NSW)")
)

AP_dat <- data.frame(
  node = c(250, 301, 297, 291),
  name = c("A14 Ben Lomond NP (TAS)", "A10 Prussian Creek - Kosciuszko NP (NSW)", "A07 Mount Kosciuszko NP (NSW)", "A11 Falls Creek (VIC)")
) 

AE_dat <- data.frame(
  node = c(251, 294),
  name = c("A04 Mount Wellington (TAS)", "A03 Falls Creek (VIC)")
) 

O_dat <- data.frame(
  node = c(196, 197),
  name = c("W02 Upper Brookfield (QLD)", "T01 Lamington NP (QLD)")
)

# create the ggtree tree
tree <- ggtree(raxml) +
  geom_rootpoint(size=4) +
  theme(legend.position="none") +
  hexpand(.3, direction = 1)

#save tree data as table
data <- tree$data
write.table(data, "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /2_phylogeny_pop_structure/data/WAout/phy_sf7_pruned_PCA_10kb_noD1_tree_data_WAout.txt", sep = "\t")

# add the clade annotations to the tree (all samples tree only)
tree <- tree + geom_cladelab(data = D_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#FFA500", barcolour="#FFA500", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = H_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#32CD32", barcolour="#32CD32", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = AP_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#333399", barcolour="#333399", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = AE_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#333399", barcolour="#333399", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = O_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour=c("#663300", "#6699CC"), barcolour=c("#663300", "#6699CC"), fontsize=5, fontface=2) #add offset=.8 to move clade labels to the right



### List of nodes to display bootstrap values
bootlist <- c(523, 274, 275, 273, 520, 522, 285, 286, 287, 288, 289, 290, 296, 314, 473, 474, 315, 316, 317, 318, 319, 329, 411, 330, 331, 365, 332)

### list of nodes to collapse
nodelist <- c(333, 385, 394, 412, 475, 533, 524, 306, 366, 403, 430, 320, 442, 484, 493, 502, 276, 521, 511, 251, 250, 301, 297, 291, 294, 196, 197)

### collapse tree
for (cn in nodelist) {
  tree <- collapse(tree, node=cn) 
}

### add labels to population nodes 
tree <- tree + geom_cladelab(data = D_dat, mapping = aes(node = node, label = name), offset=.02, fontsize=16)
tree <- tree + geom_cladelab(data = H_dat, mapping = aes(node = node, label = name), offset=.02, fontsize=16) 
tree <- tree + geom_cladelab(data = AE_dat, mapping = aes(node = node, label = name), offset=.02, fontsize=16) 
tree <- tree + geom_cladelab(data = AP_dat, mapping = aes(node = node, label = name), offset=.02, fontsize=16) 
tree <- tree + geom_cladelab(data = O_dat, mapping = aes(node = node, label = name), offset=.02, fontsize=16) #add offset=.8 to move clade labels to the right

### update colours of node points 
tree <- tree + geom_point2(aes(subset=(node==333)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==385)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==394)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==412)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==475)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==533)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==524)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==306)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==1)), shape=19, size=14, colour='#FFA500')
tree <- tree + geom_point2(aes(subset=(node==2)), shape=19, size=14, colour='#FFA500')

tree <- tree + geom_point2(aes(subset=(node==366)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==403)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==430)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==320)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==442)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==484)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==493)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==502)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==276)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==521)), shape=19, size=14, colour='#32CD32')
tree <- tree + geom_point2(aes(subset=(node==511)), shape=19, size=14, colour='#32CD32')

tree <- tree + geom_point2(aes(subset=(node==251)), shape=19, size=14, colour='#66CCFF')
tree <- tree + geom_point2(aes(subset=(node==250)), shape=19, size=14, colour='#333399')
tree <- tree + geom_point2(aes(subset=(node==301)), shape=19, size=14, colour='#333399')
tree <- tree + geom_point2(aes(subset=(node==297)), shape=19, size=14, colour='#333399')
tree <- tree + geom_point2(aes(subset=(node==291)), shape=19, size=14, colour='#333399')
tree <- tree + geom_point2(aes(subset=(node==294)), shape=19, size=14, colour='#66CCFF')

tree <- tree + geom_point2(aes(subset=(node==196)), shape=19, size=14, colour='#663300')
tree <- tree + geom_point2(aes(subset=(node==197)), shape=19, size=14, colour='#6699CC')

### add bootstrap values to select nodes 
tree <- tree + geom_label(data = subset(data, node %in% bootlist), aes(label = bootstrap), size=8)

### print
tiff("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/3_results /2_phylogeny_pop_structure/figures/test2.tiff", units="in", width=60, height=35, res=200)
print(tree)
dev.off()


# # find 'population' nodes for annotations 
# MRCA(tree, "D01415", "D01243") #D01 315
# MRCA(tree, "H02286", "H02206") #H02 348
# MRCA(tree, "D03294", "D03270") #D03 367
# MRCA(tree, "D00217", "D00210") #D00 376
# MRCA(tree, "H00252", "H00215") #H00 385
# MRCA(tree, "D04201", "D04364") #D04 394
# MRCA(tree, "H05331", "D04193") #H05 412
# MRCA(tree, "H04253", "H04222") #H04 302
# MRCA(tree, "H01276", "H01281") #H01 424
# MRCA(tree, "H06641", "H06220") #H06 466
# MRCA(tree, "D05268", "D05305") #D05 457
# MRCA(tree, "H07233", "H07270") #H07 475
# MRCA(tree, "H03219", "H03273") #H03 484
# MRCA(tree, "D35097", "D35061") #D35 518
# MRCA(tree, "D09309", "D09225") #D09 509
# MRCA(tree, "H12205", "H12264") #H12 529
# MRCA(tree, "H15058", "H15051") #H15 505
# MRCA(tree, "H14270", "H14301") #H14 494
# MRCA(tree, "D12212", "D12221") #D12 287
# MRCA(tree, "A1005", "A1003") #A10 281
# MRCA(tree, "A0702", "A07p7") #A7 277
# MRCA(tree, "A1103", "A1102") #A11 538
# MRCA(tree, "A0315", "A3") #A03 538