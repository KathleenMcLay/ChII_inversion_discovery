### Visualise the phylogeny! 

library(devtools)
install_github("https://github.com/YuLab-SMU/ggtree")
install.packages("phangorn")

library(tidyverse)
library(ggtree)
library(treeio)
library(ggrepel)
library(phangorn)

### 1: Create tree with midpoint root 

# read in raxml file 
raxml_file <- "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/data/RAxML_bipartitionsBranchLabels.ML"
raxml <- read.raxml(raxml_file)
pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/figures/Mid_R_BL_phy_sf7_pruned_PCA_10kb_noD1.pdf", width=25, height=60)

# get midpoint root for tree with phangorn package 
raxml <- read.tree(raxml_file)
mid <- midpoint(raxml)
getRoot(mid)

# create the ggtree tree using the midpoint rooted phylogeny create with phangorn
tree <- ggtree(mid) +
  geom_rootpoint(size=4) +
  theme(legend.position="none") +
  hexpand(.3, direction = 1)

# save tree data as table
data <- tree$data
write.table(data, "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/data/phy_sf7_pruned_PCA_10kb_noD1_tree_data_MID_R.txt", sep = "\t")

# create annotation details for populations of each ecotype ie. Dunes, Headlands, Alpines 
D_dat <- data.frame(
  node = c(316, 367, 376, 394, 457, 518, 509, 256, 255, 288),
  name = c("D01 Lennox Head (NSW)", "D03 Cabarita Beach (QLD)", "D00 Stradbroke Island (QLD)", "D04 Coffs Harbour (NSW)", "D05 Hat Head (NSW)", "D35 Flinders Peninsula (WA)", "D09 Leeuwin-Naturaliste National Park (WA)", "D32-254 Discovery Bay (VIC)", "D32-288 Discovery Bay (VIC)", "D12 Bermagui Dune (NSW)")
)

H_dat <- data.frame(
  node = c(348, 385, 412, 303, 424, 466, 475, 484, 529, 505, 494),
  name = c("H02 Cabarita Beach (QLD)", "H00 Stradbroke Island (QLD)", "H05 Coffs Harbour (NSW)","H04 Byron Bay (NSW)", "H01 Lennox Head (NSW)", "H06 Hat Head (NSW)", "H07 Port Macquarie (NSW)", "H03 Kiama Blowhole (NSW)", "H12 Portland, Cape Bridgewater (VIC)", "H15 Port Aurthur (TAS)", "H14 Green Cape (NSW)")
)

A_dat <- data.frame(
  node = c(234, 233, 282, 278, 538, 274),
  name = c("A04 Mount Wellington (TAS)", "A14 Ben Lomond NP (TAS)", "A10 Prussian Creek - Kosciuszko NP (NSW)", "A07 Kosciuszko NP (NSW)", "A11 Falls Creek (VIC", "A03 Falls Creek (VIC)")
) #A03, A04 = upright

O_dat <- data.frame(
  node = c(179, 180),
  name = c("W02 Upper Brookfield (QLD)", "T01 Lamington NP (QLD)")
)

# add the annotations to the tree 
tree <- tree + geom_tiplab() +
  geom_cladelab(data = D_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour="#FF9933", barcolour="#FF9933", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = H_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour="#009933", barcolour="#009933", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = A_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour="#333399", barcolour="#333399", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = O_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour=c("#663300", "#6699CC"), barcolour=c("#663300", "#6699CC"), fontsize=5, fontface=2) 

print(tree)
dev.off()

### 2: Create tree just from raxml file (WA OUTGROUP)

# annotate tree with populations 
raxml_file <- "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/data/WAout/RAxML_bipartitionsBranchLabels.ML"
raxml <- read.raxml(raxml_file)
pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/figures/WA_out_BL_phy_sf7_pruned_PCA_10kb_noD1.pdf", width=25, height=60)

# create the ggtree tree
tree <- ggtree(raxml) +
  geom_rootpoint(size=4) +
  theme(legend.position="none") +
  hexpand(.3, direction = 1)

#save tree data as table
data <- tree$data
write.table(data, "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/data/WAout/phy_sf7_pruned_PCA_10kb_noD1_tree_data_WAout.txt", sep = "\t")

# create annotation details for populations of each ecotype ie. Dunes, Headlands, Alpines 
D_dat <- data.frame(
  node = c(333, 385, 394, 412, 475, 533, 524, 2, 1, 306),
  name = c("D01 Lennox Head (NSW)", "D03 Cabarita Beach (QLD)", "D00 Stradbroke Island (QLD)", "D04 Coffs Harbour (NSW)", "D05 Hat Head (NSW)", "D35 Flinders Peninsula (WA)", "D09 Leeuwin-Naturaliste National Park (WA)", "D32-254 Discovery Bay (VIC)", "D32-288 Discovery Bay (VIC)", "D12 Bermagui Dune (NSW)")
)

H_dat <- data.frame(
  node = c(366, 403, 430, 320, 442, 484, 493, 502, 276, 521, 511),
  name = c("H02 Cabarita Beach (QLD)", "H00 Stradbroke Island (QLD)", "H05 Coffs Harbour (NSW)","H04 Byron Bay (NSW)", "H01 Lennox Head (NSW)", "H06 Hat Head (NSW)", "H07 Port Macquarie (NSW)", "H03 Kiama Blowhole (NSW)", "H12 Portland, Cape Bridgewater (VIC)", "H15 Port Aurthur (TAS)", "H14 Green Cape (NSW)")
)

A_dat <- data.frame(
  node = c(251, 250, 301, 297, 291, 294),
  name = c("A04 Mount Wellington (TAS)", "A14 Ben Lomond NP (TAS)", "A10 Prussian Creek - Kosciuszko NP (NSW)", "A07 Mount Kosciuszko NP (NSW)", "A11 Falls Creek (VIC)", "A03 Falls Creek (VIC)")
) #A03, A04 = upright

O_dat <- data.frame(
  node = c(196, 197),
  name = c("W02 Upper Brookfield (QLD)", "T01 Lamington NP (QLD)")
)

# add the annotations to the tree 
tree <- tree + geom_tiplab() +
  geom_cladelab(data = D_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#FF9933", barcolour="#FF9933", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = H_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#009933", barcolour="#009933", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = A_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour="#333399", barcolour="#333399", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = O_dat, mapping = aes(node = node, label = name), offset=.05, barsize=0.5, textcolour=c("#663300", "#6699CC"), barcolour=c("#663300", "#6699CC"), fontsize=5, fontface=2) #add offset=.8 to move clade labels to the right

print(tree)
dev.off()

# find 'population' nodes for annotations 
MRCA(tree, "D01215", "D01243") #D01 315
MRCA(tree, "H02286", "H02206") #H02 348
MRCA(tree, "D03294", "D03270") #D03 367
MRCA(tree, "D00217", "D00210") #D00 376
MRCA(tree, "H00252", "H00215") #H00 385
MRCA(tree, "D04201", "D04364") #D04 394
MRCA(tree, "H05331", "D04193") #H05 412
MRCA(tree, "H04253", "H04222") #H04 302
MRCA(tree, "H01276", "H01281") #H01 424
MRCA(tree, "H06641", "H06220") #H06 466
MRCA(tree, "D05268", "D05305") #D05 457
MRCA(tree, "H07233", "H07270") #H07 475
MRCA(tree, "H03219", "H03273") #H03 484
MRCA(tree, "D35097", "D35061") #D35 518
MRCA(tree, "D09309", "D09225") #D09 509
MRCA(tree, "H12205", "H12264") #H12 529
MRCA(tree, "H15058", "H15051") #H15 505
MRCA(tree, "H14270", "H14301") #H14 494
MRCA(tree, "D12212", "D12221") #D12 287
MRCA(tree, "A1005", "A1003") #A10 281
MRCA(tree, "A0702", "A07p7") #A7 277
MRCA(tree, "A1103", "A1102") #A11 538
MRCA(tree, "A0315", "A3") #A03 538
