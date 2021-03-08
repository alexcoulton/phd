library(ape)
library(phylobase)
library(adegenet)
library(adephylo)

tree = read.tree("original_data/phylogeny/snphylo.output.ml.tree")
tree.dist = distTips(tree)

tree.dist2 = as.data.frame(as.matrix(tree.dist))
tree.dist3 = tree.dist2[grep("Paragon", rownames(tree.dist2))[1], ]
