#coming up with test of the sequences in the multiple sequence alignment... remove any duds

multalign1 = readDNAMultipleAlignment("rotation1scripts_v4/auto.primer.picker/Traes700.full.test/seq/extended/alignments/all.align.rev.fa.trimmed.fa")
g = as.matrix(multalign1)
start.base.coord = min(grep("A|G|T|C", mult.align1.mat[1, ]))
matrix.before.start = mult.align1.mat[, 1:start.base.coord]

#get number of insertions (hyphens) in each sequence
num.hyphens = apply(matrix.before.start, 1, function(x) length(which(x == "-")))[2:nrow(matrix.before.start)]
num.hyphens[which(num.hyphens == max(num.hyphens))]


#exclude first sequence from comparison as it only contains hyphens
combinations.mat = as.data.frame(t(combn(2:nrow(matrix.before.start), 2)))
combinations.mat$n.same = ""

for(i in 1:nrow(combinations.mat)){
    matrix.before.start = as.data.frame(matrix.before.start)
    matrix.only.bases = matrix.before.start[, -which(matrix.before.start[combinations.mat[i, 1], ] == "-")]
    browser()
    
    combinations.mat$n.same[i] = length(which(matrix.before.start[combinations.mat[i, 1], ] == matrix.before.start[combinations.mat[i, 2]]))
}

g2 = expand.grid(1:nrow(matrix.before.start), 1:nrow(matrix.before.start))
g2$n.same = ""

for(i in 1:nrow(g2)){
    
    print(g2[i, 1])
    print(g2[i, 2])
    browser()
    g2$n.same[i] = length(which(matrix.before.start[g2[i, 1], ] == matrix.before.start[g2[i, 2]]))
    browser()
}

for(i in unique(g2[1, ])){
    

}

lapply(unique(g2[, 1]), function(x){
    temp.df = filter(g2, Var1 == x)
    mean(as.numeric(temp.df$n.same))
})



g = strsplit(primer.files, "")
num.coords = lapply(g, function(x){
    q = grep("p", x)
    q[1]
})

primer.files.order = unlist(Map(function(x, y){
    as.numeric(paste(x[1:(y - 1)], collapse = ""))
}, g, num.coords))

primer.files[match(1:length(primer.files.order), primer.files.order)]




