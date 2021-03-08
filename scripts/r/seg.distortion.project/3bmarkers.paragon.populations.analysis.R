#other paragon datasets 
p.x.baj = read_csv("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/p.x.baj.cleaned.csv")
p.x.baj = read_csv("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/p.x.becard.cleaned.csv")

all.p.crosses = lapply(list.files("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/", pattern = ".csv"), function(x){
    read_csv(p("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/", x))
    })

names(all.p.crosses) = list.files("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/", pattern = ".csv")

files1 = list.files("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/", pattern = ".csv")


load("rotation1scripts_v4/saved.objects/all.markers.distort")

#all.markers.distort from avni.reanalysis.R



chi.stats1 = lapply(all.p.crosses, function(x){
    g = x[, which(colnames(x) %in% all.markers.distort)]
    q = lapply(g, geno.count.chi)
    mean(unlist(lapply(q, function(x) x[3])))
    # browser()
    # parent1 = x[2, 1]
    # parent2 = x[3, 1]
    # haplotype1 = g[2, ]
    # haplotype2 = g[3, ]
    # haps = list(haplotype1, haplotype2)
    # names(haps) = c(parent1, parent2)
    # 
    # lapply(g, geno.count.chi)
    
    # haps
    
})

lapply(chi.stats1, length)

chi.stats1

do.call(rbind, chi.stats1)

chi.stats1[[1]]
chi.stats1[[2]]


#analysis of haplotypes

all.pops = lapply(list.files("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/", pattern = ".txt"), function(x){
    convert.aas.to.rqtl(p("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/", x))
})



g = lapply(all.pops, function(x){
    cellnames1 = as.character(x[, 2])
    # browser()
    cellnames1.char = nchar(cellnames1)
    coords1 = sort(cellnames1.char, index.return = T)$ix
    cellnames1[coords1]
    
})

g[6:7]

all.pops = all.pops[-(6:7)]
all.pops = lapply(all.pops, convert.to.character.data.frame)

parents1 = list(c("P04_Apogee", "P08_Apogee"),
                                c("M23_CIM47", "M24_CIM47", "N23_CIM47", "P24_PARAGON10", "P23_PARAGON14", "P23_PARAGON9"),
                                c("P23_Paragon18", "C22_PAR", "N23_CIM49", "N24_CIM49"),
                                c("N23_MISR1"),
                                c("M01_Starke", "O05_Starke"),
                                c("N23_BAJ"),
                                c("N23_Bacard_kachu"),
                                c("N24_Pfau"),
                                c("M23_Super152"),
                                c("N24_Waxwing"),
                                c("M23_Wyalkatchem"))

parent.haplotypes1 = Map(function(x, y){
    g = sapply(y, function(z) x[grep(z, x[, 2]), all.markers.distort])
    as.data.frame(t(g))
}, all.pops, parents1)

parent.haplotypes1 = do.call(rbind, parent.haplotypes1)

parent.haplotypes1[] = lapply(parent.haplotypes1, as.numeric)

#perform sorting to find haplotypes most similar to apogee haplotyp
vec1 = as.numeric()
for(i in 2:nrow(parent.haplotypes1)){
    q = length(which(as.numeric(parent.haplotypes1[1, ]) == as.numeric(parent.haplotypes1[i, ])))
    vec1 = c(vec1, q)
}

parent.haplotypes1.sorted = parent.haplotypes1[c(1, sort(vec1, decreasing = T, index.return = T)$ix + 1), ]


names(all.pops) = list.files("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/", pattern = ".txt")[-(6:7)]

#grab ratios of raw AAS export data for all other paragon crosses
chi.stats2 = lapply(all.pops, function(x){
    g = x[, which(colnames(x) %in% all.markers.distort)]
    g = g[-1, ]
    q = lapply(g, function(y){
        a = length(which(y == 0))
        b = length(which(y == 2))
        ab = list(a, b)
        names(ab) = c("A", "B")
        ab
    })
    do.call(rbind, q)
    # browser()
})

#remove monomorphic SNPs
chi.stats2 = lapply(chi.stats2, function(x){
    x = as.data.frame(x)
    coords1 = unlist(lapply(x, function(y){
        
        which(y == 0 | y == 1 | y == 2 | y == 3)
    }))
    
    coords1 = unique(unname(coords1))
    if(length(coords1 > 0)) x = x[-coords1, ]
    
    # browser()
    x
})

lapply(chi.stats2, nrow)

seg.ratios.para.crosses = lapply(chi.stats2, function(x){
    x = as.data.frame(t(x))
    mean(unlist(lapply(x, function(y) chisq.test(c(y[[1]], y[[2]]))[3])))
})

seg.ratios.para.crosses[which(is.na(seg.ratios.para.crosses))] = 1000

seg.ratios.para.crosses.sorted = sort(unlist(seg.ratios.para.crosses))

seg.ratios.para.crosses.sorted2 = as.data.frame(seg.ratios.para.crosses.sorted)

write.csv(parent.haplotypes1.sorted, "rotation1scripts_v4/processed_data/3b.segregation.distortion.apogee.x.paragon/parent.haplotypes1.sorted.csv")
write.csv(seg.ratios.para.crosses.sorted2, "rotation1scripts_v4/processed_data/3b.segregation.distortion.apogee.x.paragon/seg.ratios.para.crosses.sorted2.csv")
