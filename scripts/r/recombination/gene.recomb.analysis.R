functional.anno = read.csv("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.csv", header = T, fill = T, stringsAsFactors = F)

dmc1 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("dmc1", x, ignore.case = T)
})), ]

rad51 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("rad51", x, ignore.case = T)
})), ]

zyp1 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("zyp1", x, ignore.case = T)
})), ]

asy1 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("asy1", x, ignore.case = T)
})), ]

msh5 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("msh5", x, ignore.case = T)
})), ]

msh4 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("msh4", x, ignore.case = T)
})), ]

meiosis = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("meiosis", x, ignore.case = T)
})), ]

spo11 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("spo11", x, ignore.case = T)
})), ]

fancm = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("fancm", x, ignore.case = T)
})), ]

h2ax = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("h2ax", x, ignore.case = T)
})), ]

mlh1 = functional.anno[unlist(lapply(functional.anno, function(x){
    grep("mlh1", x, ignore.case = T)
})), ]

allmeiosis = list(dmc1, rad51, zyp1, asy1, msh5, msh4, meiosis, spo11, fancm, h2ax, mlh1)
allmeiosis2 = allmeiosis[which(unlist(lapply(allmeiosis, function(x) nrow(x) > 0)))]

library(tibble)
library(dplyr)
allmeiosis3 = bind_rows(allmeiosis2)
allm3 = allmeiosis3 


write.csv(allmeiosis3, "rotation1scripts_v4/original_data/recombination.gene.analysis/recombination.gene.analysis.csv", row.names = F)

meiosis.genes = data.frame(unique(multi.str.split(allm3$Gene.ID, "\\.", 1)))
write.csv(meiosis.genes, "rotation1scripts_v4/original_data/recombination.gene.analysis/meiosis.genes.csv", row.names = F)

g1 = allmeiosis2$Gene.ID

g = functional.anno[sort(unique(unname(unlist(lapply(functional.anno, function(x) grep("meiosis", x, ignore.case = T)))))), ]
write.csv(g, "rotation1scripts_v4/original_data/recombination.gene.analysis/g.csv", row.names = F)






##### misc code 


mgenes = read.csv("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/meiosis.genes.csv")
mgenes2 = unlist(lapply(as.character(mgenes[, 1]), function(x) grep(x, names(genes))))

lapply(mgenes2, function(x){
	writeXStringSet(genes[x], paste0("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/genes/", x, ".fa"))
	})




 # for i in ~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/genes/* ; do fng=$(echo $i | sed 's/.*\///') echo $fng; done
