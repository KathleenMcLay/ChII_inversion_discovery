raxml_file <- "/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/data/RAxML_bipartitionsBranchLabels.ML"
raxml <- read.raxml(raxml_file)
pdf("/Users/kathleenmclay/Google Drive/PhD/Chapter_1_inversions/results /phylogeny_pop_structure/figures/phy.pdf", width=25, height=60)


D_dat <- data.frame(
  node = c(315, 367, 376, 394, 457, 518, 509, 256, 255, 287),
  name = c("D01 Lennox Head (NSW)", "D03 Cabarita Beach (QLD)", "D00 Stradbroke Island (QLD)", "D04 Coffs Harbour (NSW)", "D05 Hat Head (NSW)", "D35 Flinders Peninsula (WA)", "D09 Leeuwin-Naturaliste National Park (WA)", "D32-254 Discovery Bay (VIC)", "D32-288 Discovery Bay (VIC)", "D12 Bermagui Dune (NSW)")
)

H_dat <- data.frame(
  node = c(348, 385, 412, 302, 424, 466, 475, 484, 529, 505, 495),
  name = c("H02 Cabarita Beach (QLD)", "H00 Stradbroke Island (QLD)", "H05 Coffs Harbour (NSW)","H04 Byron Bay (NSW)", "H01 Lennox Head (NSW)", "H06 Hat Head (NSW)", "H07 Port Macquarie (NSW)", "H03 Kiama Blowhole (NSW)", "H12 Portland, Cape Bridgewater (VIC)", "H15 Port Aurthur (TAS)", "H14 Green Cape (NSW)")
)

#A03, A04 = upright
A_dat <- data.frame(
  node = c(234, 233, 281, 277, 538, 2, 271, 1),
  name = c("A04 Mount Wellington (TAS)", "A14", "A10", "A07 Mount Kosciuszko NP (SA)", "A11", "A03 Falls Creek (SA)", "A03-23 Falls Creek (SA)", "A03-15 Falls Creek (SA)")
)

O_dat <- data.frame(
  node = c(179, 180),
  name = c("W02 Upper Brookfield (QLD)", "T01 Lamington NP (QLD)")
)

tree <- ggtree(raxml, branch.length = "none") +
  geom_rootpoint(size=4) +
  theme(legend.position="none") +
  hexpand(.3, direction = 1)
tree <- tree + geom_tiplab() + 
  geom_cladelab(data = D_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour="#FF9933", barcolour="#FF9933", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = H_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour="#009933", barcolour="#009933", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = A_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour="#333399", barcolour="#333399", fontsize=5, fontface=2) 
tree <- tree + geom_cladelab(data = O_dat, mapping = aes(node = node, label = name), offset=.8, align=TRUE, barsize=2, textcolour=c("#663300", "#6699CC"), barcolour=c("#663300", "#6699CC"), fontsize=5, fontface=2) 

print(tree)
dev.off()