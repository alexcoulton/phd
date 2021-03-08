library(treeio)
library(ape)

tree1 = read.beast("rotation1scripts_v4/original_data/watkins.phylogenies/full_sequences_w_tauschii_NO_indels_mcc.trees")
get.fields(tree1)
tree1
tree1.data = get.data(tree1)



tree.nexus = read.nexus("rotation1scripts_v4/original_data/watkins.phylogenies/full_sequences_w_tauschii_NO_indels_mcc.trees")
distances1 = as.data.frame(cophenetic.phylo(tree.nexus))
#grab the distance between tips and then divide by 2 for the TMRCA
distances2 = (distances1$ref / 2) * 1000000

distances3 = data.frame(rownames(distances1), distances2)

colnames(distances3) = c("sample", "tmrca")

sample.names1 = gsub(".*/", "", output.paths1)
sample.names1 = gsub("\\..*", "", sample.names1)
sample.names1 = gsub(".*_", "", sample.names1)

distances3$watkin.id = gsub(".*_", "", distances3$sample)

match(distances3$watkin.id, sample.names1)

distances3 = distances3[-c(105, 106), ]

distances3$n.missense = sapply(sift.scores.all, length)

mean(distances3$n.missense / distances3$tmrca)
sd(distances3$n.missense / distances3$tmrca)

plot(distances3$tmrca, distances3$n.missense)

sapply(sift.scores.all, mean)
boxplot(sift.scores.all)
