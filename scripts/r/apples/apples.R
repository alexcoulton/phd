library(readxl)
apples = read_excel("rotation1scripts_v4/original_data/apple.project/The girls.xlsx")

which(apples$JAMESGREVE == apples$WORCESTPEAR)


data1 = newdf("variety", "num.markers.same", no.rows = T)
for(i in 10:ncol(apples)){
    g = lapply(apples[, 10:ncol(apples)], function(x){
        length(which(apples[, i] == x))
        
        
    })
    # browser()
    
    
    g1 = colnames(apples)[i]
    g2 = sort(unlist(g), decreasing = T)[2]
    g3 = data.frame(g1, g2, stringsAsFactors = F)    
    
    data1 = rbind(data1, g3)
}

colnames(apples)[1:40]
girls = apples[, 1:40]

parents = apples[, c(10, 11, 13, 14)]
girls = apples[, 15:40]
girls = girls[, -grep("X_", colnames(girls))]

v(parents)

uninformative.coords = which(parents$JAMESGREVE == parents$WORCESTPEAR & parents$JAMESGREVE == parents$MICHELIN & parents$JAMESGREVE == parents$DABINET)
parents = parents[-uninformative.coords, ]
girls = girls[-uninformative.coords, ]


parents.split = lapply(parents, function(x){
    allele1 = multi.str.split(x, "/", 1)
    allele2 = multi.str.split(x, "/", 2)
    data.frame(allele1, allele2)
    
})

parents.split2 = do.call(cbind, parents.split)

girls.split = lapply(girls, function(x){
    allele1 = multi.str.split(x, "/", 1)
    allele2 = multi.str.split(x, "/", 2)
    data.frame(allele1, allele2)
    
})


girls.split2 = do.call(cbind, girls.split)

combgirls = cbind(parents, girls)
all.possible = unique(unname(unlist(lapply(combgirls, unique))))
to.rm1 = all.possible[-grep("[A-Z]?/[A-Z]?$|\\.", all.possible)]
to.keep1 = c(all.possible[grep("[A-Z]?/[A-Z]?$|\\.", all.possible)], NA)

rows.containing = apply(combgirls, 1, function(x){
    any(x %in% to.rm1    )
})

combgirls2 = combgirls[-which(rows.containing), ]

comb.girls3 = lapply(combgirls2, function(x){
    g = gsub("G/A", "R", x)
    g = gsub("T/C", "Y", g)
    g = gsub("A/G", "R", g)
    g = gsub("C/A", "M", g)
    g = gsub("C/T", "Y", g)
    g = gsub("T/G", "K", g)
    g = gsub("G/T", "K", g)
    g = gsub("A/C", "M", g)
    g = gsub("A/T", "W", g)
    g = gsub("T/A", "W", g)
    g = gsub("A/A", "A", g)
    g = gsub("G/G", "G", g)
    g = gsub("C/C", "C", g)
    g = gsub("T/T", "T", g)
    g
})

comb.girls4 = do.call(cbind, comb.girls3)
comb.girls4 = as.data.frame(comb.girls4)


girls.sequences = lapply(comb.girls4, function(x){
    
    paste(x, collapse = "")
})


girls.sequences = lapply(girls.sequences, function(x){
    DNAString(x)
})

do.call(DNAStringSet, girls.sequences)

girls.sequences = DNAStringSet(girls.sequences)
writeXStringSet(girls.sequences, "rotation1scripts_v4/original_data/apple.project/appleseq.fa")


