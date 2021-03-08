setwd("E:/phd.project.main/")

exome.tree = readLines("rotation1scripts_v4/original_data/watkins.phylogenies/exome.sequences.fa.treefile")
array.tree = readLines("rotation1scripts_v4/original_data/watkins.phylogenies/watkins.array.exome.subset.fa.treefile")

beast.tree = readLines("rotation1scripts_v4/original_data/watkins.phylogenies/sequences_w_tauschii_v2_output.trees")




library(readxl)
wilkins.data = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 2)

wilkins.data$`Watkins Number` = gsub("Watkins_", "", wilkins.data$`Watkins Number`)
wilkins.data$`Watkins Number` = paste0("1190", wilkins.data$`Watkins Number`)
wilkins.data$`Country of Origin` = gsub(" ", "_", wilkins.data$`Country of Origin`)


for(i in 1:nrow(wilkins.data)){
    exome.tree = gsub(wilkins.data$`Watkins Number`[i], paste0(wilkins.data$`Country of Origin`[i], "_", wilkins.data$`Watkins Number`[i]), exome.tree)
}


for(i in 1:nrow(wilkins.data)){
    array.tree = gsub(wilkins.data$`Watkins Number`[i], paste0(wilkins.data$`Country of Origin`[i], "_", wilkins.data$`Watkins Number`[i]), array.tree)
}

for(i in 1:nrow(wilkins.data)){
    beast.tree = gsub(wilkins.data$`Watkins Number`[i], paste0(wilkins.data$`Country of Origin`[i], "_", wilkins.data$`Watkins Number`[i]), beast.tree)
}

beast.tree = gsub("Sample_.*?_", "", beast.tree)
beast.tree = gsub("/$", "", beast.tree)
beast.tree = gsub("/,$", ",", beast.tree)
writeLines(beast.tree, "rotation1scripts_v4/original_data/watkins.phylogenies/beast.tree.w.geo.trees")
