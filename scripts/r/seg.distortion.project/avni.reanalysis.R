library(readxl)
setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
avni.data = read_excel("rotation1scripts_v4/original_data/genotype.data/avni.etal.2014/avni.genotype.data.xls")

avni.data2 = split.df.into.list.of.dfs.by.column(avni.data, "chr")
avni.data2 =lapply(avni.data2, function(x){
    as.data.frame(t(x))
})

source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

lapply(avni.data2, function(x){
    g = convert.recomb.to.seg2(x)
    
    g2 = unlist(lapply(g, function(y) y[[3]]))
    
    g3 = p.adjust(g2, "BH")
    g4 = p.adjust(g2, "bonferroni")
    g2.stat = unlist(length(which(g2 < 0.05)))
    g3.stat = unlist(length(which(g3 < 0.05)))
    g4.stat = unlist(length(which(g4 < 0.05)))
    c(g2.stat, g3.stat, g4.stat)
})


avni.data.flip = as.data.frame(t(avni.data))
avni.data.flip = convert.to.character.data.frame(avni.data.flip)

calc1 = function(x){
    g = convert.recomb.to.seg2(x)
    g2 = unlist(lapply(g, function(y) y[[3]]))
    g3 = p.adjust(g2, "BH")
    g4 = p.adjust(g2, "bonferroni")
    g2.stat = unlist(length(which(g2 < 0.05)))
    g3.stat = unlist(length(which(g3 < 0.05)))
    g4.stat = unlist(length(which(g4 < 0.05)))
    c(g2.stat, g3.stat, g4.stat)
}

calc2 = function(x, fdr){
    if(missing(fdr)) fdr = "bonf"
    
    g = convert.recomb.to.seg2(x)
    g2 = unlist(lapply(g, function(y) y[[3]]))
    g3 = p.adjust(g2, "BH")
    g4 = p.adjust(g2, "bonferroni")
    g2.stat = unlist(length(which(g2 < 0.05)))
    g3.stat = unlist(length(which(g3 < 0.05)))
    g4.stat = unlist(length(which(g4 < 0.05)))
    if(fdr == "bonf"){
        which(g4 < 0.05)    
    } else {
        which(g3 < 0.05)    
    }
    
}

calc1(avni.data.flip)
unique(as.character(avni.data.flip[, calc2(avni.data.flip)][3, ]))
avni.bonf.snps = avni.data.flip[, calc2(avni.data.flip)]
avni.bonf.3b = avni.bonf.snps[, which(avni.bonf.snps[3, ] == "3B")]

iselect.durum.probes = read_excel("E:/ac14037.backup/Google Drive/University/PhD/Seg dist simulation/avni et al iselect 90k array/pbi12183-sup-0010-tables5.xlsx")
iselect.durum.probes = reset.colnames(iselect.durum.probes)
iselect.durum.probes$V1

avni.3b.probe.seqs = iselect.durum.probes[match(as.character(avni.bonf.3b[1, ]), iselect.durum.probes$V1), ]



avni.3b.probe.seqs2 = avni.3b.probe.seqs$V2
#replace snp ambiguity codes with other format
avni.3b.probe.seqs2 = gsub("\\[A\\/G\\]", "R", avni.3b.probe.seqs2)
avni.3b.probe.seqs2 = gsub("\\[T\\/C\\]", "Y", avni.3b.probe.seqs2)
avni.3b.probe.seqs2 = gsub("\\[T\\/G\\]", "K", avni.3b.probe.seqs2)
avni.3b.probe.seqs2 = gsub("\\[A\\/C\\]", "M", avni.3b.probe.seqs2)

avni.3b.probe.seqs3 = na.omit(avni.3b.probe.seqs2)

library(Biostrings)

avni.3b.seqs = DNAStringSet(avni.3b.probe.seqs3)
names(avni.3b.seqs) = na.omit(avni.3b.probe.seqs$V1)

#writeXStringSet(avni.3b.seqs, "rotation1scripts_v4/original_data/fasta/avni.probe.sequences.fa")

avni.blast = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/avni.probes.vs.genome.blast")
avni.blast2 = grab.best.hits(avni.blast)
avni.blast3 = avni.blast2[which(avni.blast2$sseqid == "chr3B"), ]

#allen reanalysis

allen.axc = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", 2)
allen.sxr = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", 3)
allen.oxs = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", 4)
allen.axp = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", 5)
allen.axp[allen.axp == "AA"] = "A"
allen.axp[allen.axp == "BB"] = "B"
allen.axp[allen.axp == "AB"] = "H"
allen.csxp = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", 6)

allen.maps = list(allen.axc, allen.sxr, allen.oxs, allen.axp, allen.csxp)

allen.maps2 = lapply(allen.maps, function(x){
    g = as.data.frame(t(x))
    convert.to.character.data.frame(g)
})

lapply(allen.maps2, calc1)
calc2.results = lapply(allen.maps2, calc2)
calc2.results.fdr = lapply(allen.maps2, calc2, fdr = "fdr")





as.data.frame(allen.axc[calc2.results[[1]], ])[, 1:3]

as.data.frame(allen.axp[calc2.results[[4]], ][, 1:3])

big.distort.markers = allen.axp[calc2.results[[4]], ]
big.distort.markers = as.data.frame(big.distort.markers)



g = as.data.frame(t(big.distort.markers[big.distort.markers$chr == "3B", ]))
v(g)
g = reset.colnames(g)
lapply(g, function(x){
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    abs((a / (a + b)) - 0.5)
})



big.distort.markers3b = g

v(big.distort.markers3b)
big.distort.markers3b = convert.to.character.data.frame(big.distort.markers3b)

#grab coordinates of lines that contain no recombination events between these markers 
coords = as.numeric()
for(i in 4:nrow(big.distort.markers3b)){
    
    q = unique(as.character(big.distort.markers3b[i, ]))
    # browser()
    if(length(q) == 1){
        coords = c(coords, i)
    }
    if(length(q) == 2 & "-" %in% q){
        coords = c(coords, i)
    }
}

#remove those coordinates from the data frame
big.distort.markers3b.2 = big.distort.markers3b[-coords, ]
v(big.distort.markers3b.2)


allen.axp3b = allen.axp[which(allen.axp$chr == "3B"), ]

allen.axp3b = as.data.frame(t(allen.axp3b))

allen.axp3b.2 = allen.axp3b[-coords, ]



allen.axp3b.2 = convert.to.character.data.frame(allen.axp3b.2)

allen.axp3b.3 = rbind(1:ncol(allen.axp3b.2), allen.axp3b.2)

allen.axp3b.3[1, ] = unlist(lapply(allen.axp3b, function(x){
    q = geno.count(x)
    q[1] / q[2]
}))




v(allen.axp3b.3)

lapply(big.distort.markers3b, geno.count)

# write.csv(allen.axp3b.3, "rotation1scripts_v4/processed_data/genotypes/apogee.x.paragon.sacha/selected.lines.all.3b.csv", row.names = T)

q = attributes(genetictophysical("3B", as.character(convert.to.character.data.frame(big.distort.markers3b)[1, ]), "-"))
q$bp



v(big.distort.markers3b)
b.d.3b.v2 = as.data.frame(t(big.distort.markers3b))
v(b.d.3b.v2)

b.d.3 = b.d.3b.v2[9:14, ]

b.d.4 = lapply(b.d.3, function(x){
    g = unique(as.character(x))
    length(g)
})

b.d.3[, which(unlist(b.d.4) > 1)]


g2 = as.data.frame(t(big.distort.markers3b))
g3 = lapply(g2, function(x){
    "H" %in% x
})

g4 = which(unlist(g3))

selected.lines = g2[, g4]





# write.csv(selected.lines, "rotation1scripts_v4/temp/selected.lines.csv", row.names = T)

all.markers.distort = big.distort.markers[big.distort.markers$chr == "3B", ]$marker
# s(all.markers.distort, "rotation1scripts_v4/saved.objects/all.markers.distort", "avni.reanalysis.R")

markernames3b = big.distort.markers[big.distort.markers$chr == "3B", ]$marker[9:14]


markerpositions = genetictophysical("3B", markernames3b, "-")
markerpositionsbp = attributes(markerpositions)$bp

markerpositionsbp = sort(as.integer(na.omit(markerpositionsbp)))

genes3b = grab.genes.and.functions("3B", markerpositionsbp[1], markerpositionsbp[length(markerpositionsbp)])


genetictophysical("3B", c("AX-94457592", "AX-94531806", "AX-94452872", "AX-94452872"), "-")


genetictophysical("3B", c("AX-94785822", "AX-94716868", "AX-94452872"))
g2 = grab.genes.and.functions("3B", 707441696, 723032202)


g1 = grab.genes.and.functions("3B", 707441696, 731264570)

g1 = grab.genes.and.functions("3B", 707441696, 731264570)

write(g1$Gene.ID, "rotation1scripts_v4/processed_data/3b.segregation.distortion.apogee.x.paragon/list.of.genes.txt")



processing1 = function(geno.df1){
    geno.df2 = split.df.into.list.of.dfs.by.column(geno.df1, "chr")
    lapply(geno.df2, function(x){
        g = as.data.frame(t(x))
        g = convert.to.character.data.frame(g)
        
        g1 = convert.recomb.to.seg2(g)
        g2 = unlist(lapply(g, function(y) y[[3]]))
        g3 = p.adjust(g2, "BH")
        g4 = p.adjust(g2, "bonferroni")
        browser()
        g2.stat = unlist(length(which(g2 < 0.05)))
        g3.stat = unlist(length(which(g3 < 0.05)))
        g4.stat = unlist(length(which(g4 < 0.05)))
        c(g2.stat, g3.stat, g4.stat) 
    })
}


lapply(allen.maps, processing1)





allmarkers3b = allen.axp$marker[which(allen.axp$chr == "3B")]

range1 = match(all.markers.distort, allmarkers3b)
range1 = seq(160, 195)


phys2 = genetictophysical("3B", allmarkers3b)
na.omit(phys2)
# write.csv(as.data.frame(t(as.data.frame(list(round(phys2, digits = 2), attributes(phys2)$bp)))), "rotation1scripts_v4/processed_data/genotypes/apogee.x.paragon.sacha/3bphys.csv")


g = grab.genes.and.functions("3B", 707441696, 772406808)

# s(probe.sequences, "rotation1scripts_v4/saved.objects/probe.sequences", "NA")
load("rotation1scripts_v4/saved.objects/probe.sequences")
probesubset3b = allmarkers3b[range1]
probe.seq.subset = probe.sequences[na.omit(match(switch.affy.format(probesubset3b), names(probe.sequences)))]
# writeXStringSet(probe.seq.subset, "rotation1scripts_v4/original_data/fasta/3bprobesubset.seg.dist.markers.fa")


blast3bprobes = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/3bprobes.subset.paragon.genome.blast")
blast3bprobes2 = grab.best.hits(blast3bprobes)
paragon.3b.scaffolds = blast3bprobes2$sseqid

s(paragon.3b.scaffolds, "rotation1scripts_v4/saved.objects/paragon.3b.scaffolds", "avni.reanalysis.R")

blastpara3b.vs.iwgsc.3b.truncated = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/paragon3b.vs.iwgsc3b.blast")






