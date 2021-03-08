sacha.all.snps = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.sacha/aas.all.snps.output.txt", c("P04_Apogee 2.CEL_call_code", "P18_Paragon 4.CEL_call_code"))
grep("Paragon", sacha.all.snps$probeset_id)
grep("Apogee", sacha.all.snps$probeset_id)

sacha.all.snps = sacha.all.snps[-c(350, 353, 344, 346), ]

sacha.all.snps2 = sacha.all.snps[, c(1, 2, match(allmf2, colnames(sacha.all.snps)))]

source("rotation1scripts_v4/scripts/r/axiom.data.processing/generic.axiom.data.preparation.functions.R")

sacha.all.snps3 = cleanup(sacha.all.snps2, c(2, 3))
sacha.all.snps3 = convert.to.character.data.frame(sacha.all.snps3)

sacha.all.snps3 = assign.parental.genotypes(sacha.all.snps3, 2:3)

write.csv(sacha.all.snps3, "rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.sacha/sacha.map.reconstructed.csv", row.names = F)


axiomdata_map = mstmap(axiomdata, trace = T, dist.fun = "kosambi", id = "probeset_id", p.value = 1e-45, anchor = F)
sachapullmap = pull.map(axiomdata_map, as.table = T)

sachapullmap$notshared = ""
sachapullmap$marker = rownames(sachapullmap)
sachapullmap$marker = switch.affy.format(sachapullmap$marker)

sachapullmap[which(sachapullmap$marker %in% markers.not.shared), ]$notshared = T

sp2 = split(sachapullmap, sachapullmap$chr)


sp2 = lapply(sp2, function(x){
    g = grab.consensus.chromosome(x$marker, T)
    x$chromo = names(g[1])
    x
})

sp2 = sp2[sort(sapply(sp2, nrow), decreasing = T, index.return = T)$ix]

map.comparison = data.frame(sapply(sp2, nrow),

sapply(sp2, function(x){
    (length(which(x$notshared == T)) / nrow(x)) * 100
}),

unlist(lapply(sp2, function(x){
    q = unique(x$chromo)[1]
    if(is.null(q)) q = "none"
    q
}))

)

colnames(map.comparison) = c("num.markers", "percent.ns.marker", "chr")
sum(map.comparison[which(map.comparison$percent.ns.marker == 100), ]$num.markers)


sp3 = bind_rows(sp2)

allen.split1 = split(allen.axp.map4, allen.axp.map4$chromo)

lapply(allen.split1, function(x){
    g = sort(table(na.omit(sachapullmap[match(x$marker, rownames(sachapullmap)), ]$chr)), decreasing = T)
    g[-which(g == 0)]
})




c5.split = split(c5, c5$chr1)
mapping1 = lapply(c5.split, function(x){
    g = sort(table(sachapullmap[match(x$marker, rownames(sachapullmap)), ]$chr), decreasing = T)
    g[-which(g == 0)][-1]
})


sort(sapply(sp2, nrow), decreasing = T)






sort(table(sachapullmap[na.omit(match(markers.not.shared, rownames(sachapullmap))), ]$chr), decreasing = T)

sort(unlist(unname(mapping1)), decreasing = T)

#### RECONSTRUCT AXP F2 GENETIC MAP ####

