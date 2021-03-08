#AutoCloner alignment comparisons


library(Biostrings)


gene1 = readDNAStringSet("rotation1scripts_v4/original_data/autocloner/TraesCS5A01G531300_set1.fa")
gene2 = readDNAStringSet("rotation1scripts_v4/original_data/autocloner/TraesCS5A01G531700_set1.fa")
gene3 = readDNAStringSet("rotation1scripts_v4/original_data/autocloner/TraesCS5A01G530800_set1.fa")


analyse.genes = function(gene){
    g1.m = as.matrix(gene)
    gene.range = min(which(g1.m[1, ] != "-")):max(which(g1.m[1, ] != "-"))
    
    g1.m2 = g1.m[, gene.range]
    
    
    lapply(3:nrow(g1.m2), function(i){
        gene.length = ncol(g1.m2)
        non.hyphen.coords = which(g1.m2[i, ] != '-' & g1.m2[1, ] != '-')
        
        identity1 = length(which(g1.m2[1, ] == g1.m2[i, ])) / gene.length
        identity2 = length(which(g1.m2[1, non.hyphen.coords] == g1.m2[i, non.hyphen.coords])) / length(non.hyphen.coords)
        
        gc.content = length(which(g1.m2[i, non.hyphen.coords] == "G" | g1.m2[i, non.hyphen.coords] == "C")) / length(non.hyphen.coords)
        
        coverage = length(non.hyphen.coords) / ncol(g1.m2)
        c(identity1, identity2, coverage, gc.content)
    })
    
    
}

g = strsplit(as.character(gene1[[1]]), "")[[1]]
min(which(g != "-"))

g = strsplit(as.character(gene1[[5]]), "")[[1]]
min(which(g != "-"))

analyse.genes(gene1)
analyse.genes(gene2)
analyse.genes(gene3)
