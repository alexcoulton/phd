list.of.nonmonotonic.markers = read.table("processed_data/nonmonotonic.markers.txt", sep = "\n", stringsAsFactors = F)
list.of.nonmonotonic.markers = as.character(list.of.nonmonotonic.markers$V1)

path = function(quadrant) paste("original_data/genotype.data/alexq", quadrant, "_man_cur_sk_map_genogeno.csv", sep = "")
genotype.data.paths = list(path(1), path(2), path(3), path(4))

list.of.genotypes = lapply(genotype.data.paths, function(x) read.csv(x, header = T, stringsAsFactors = F))

list.of.genotypes = lapply(list.of.genotypes, function(x){ 
    x = x[, -match(list.of.nonmonotonic.markers, colnames(x))]
    return(x)
    })

#remove NA value from coord 1,1
list.of.genotypes = lapply(list.of.genotypes, function(x){
    x[1,1] = ""
    return(x)
})

for(i in seq_along(list.of.genotypes)){
    write.csv(list.of.genotypes[[i]], paste("processed_data/genotypes/alexq", i, "_man_cur_nonmonotonic_removed.csv", 
                                                                                    sep = ""), row.names = F)
}
