#generic axiom data preparation

setwd("E:/phd.project.main")
source("rotation1scripts_v4/scripts/r/functions.R")
library(dplyr)
library(tibble)


#     ____________________________________________________________________________
#     INITIAL PROCESSING                                                                                                            ####

p.x.baj = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/alex.aas.export/p.x.baj2.txt",
                                                            c("X23_Paragon", "N23_BAJ"))

p.x.becard = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/alex.aas.export/PxBacard.txt",
                                                                 c("X23_Paragon", "N23_Bacard_kachu"))

p.x.super152 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/alex.aas.export/PxSuper152.txt",
                                                                 c("X23_Paragon", "M23_Super"))

p.x.waxwing = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/alex.aas.export/PxWaxwing.txt",
                                                                     c("Paragon", "N24_Waxwing"))

p.x.wyalkatchem = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/alex.aas.export/PxWyalkatchem.txt",
                                                                            c("Paragon", "M23_Wyalkatchem"))

p.x.pfau = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/other.paragon.mapping.populations/alex.aas.export/P.Pfau2.txt",
                                                                            c("Paragon", "N24_Pfau"))

apogee.x.paragon = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.sacha/rec.snps.txt", c("Paragon", "Apogee"))

apogee.x.paragon = apogee.x.paragon[-(3:6), ]

m.data = convert.aas.to.rqtl("rotation1scripts_v4/temp/Re-exportedphighres.txt", c("E09_PBW", "O13_PBW"), call.code.numeric = T)


c.x.a.f2 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/c.x.a.f2.txt", c("Cadenza", "Avalon"))

    

axc51_2 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/axc51-2.all.snps.txt")
axc61_1 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/axc61-1.all.snps.txt")

cxa65_1 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/c.x.a65-1.all.snps.txt")
cxa65_3 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/c.x.a65-3.all.snps.txt")


#grab only the probes that were selected in the c.x.a.f2 cross
axc51_2 = axc51_2[, c(match((colnames(c.x.a.f2)), colnames(axc51_2)))]
axc61_1 = axc61_1[, c(match((colnames(c.x.a.f2)), colnames(axc61_1)))]

cxa65_1 = cxa65_1[, c(match((colnames(c.x.a.f2)), colnames(cxa65_1)))]
cxa65_3 = cxa65_3[, c(match((colnames(c.x.a.f2)), colnames(cxa65_3)))]


#add cadenza and avalon parents to dataframe
axc51_2 = rbind(c.x.a.f2[2:3, ], axc51_2)
axc51_2 = axc51_2[c(3, 1, 2, 5:nrow(axc51_2)), ]

axc61_1 = rbind(c.x.a.f2[2:3, ], axc61_1)
axc61_1 = axc61_1[c(3, 1, 2, 4:nrow(axc61_1)), ]

cxa65_1 = rbind(c.x.a.f2[2:3, ], cxa65_1)
cxa65_1 = cxa65_1[c(3, 1, 2, 4:nrow(cxa65_1)), ]

cxa65_3 = rbind(c.x.a.f2[2:3, ], cxa65_3)
cxa65_3 = cxa65_3[c(3, 1, 2, 4:nrow(cxa65_3)), ]

#prepare double haploid genotyping data
a.x.cdh.arraya = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza/avalon.x.cad.all.snps.arraya.txt")
a.x.cdh.arrayb = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza/avalon.x.cad.all.snps.arrayb.txt")

a.x.cdh.arraya$probeset_id = as.character(a.x.cdh.arraya$probeset_id)
a.x.cdh.arrayb$probeset_id = as.character(a.x.cdh.arrayb$probeset_id)

#the individuals on array a and b are not in the same row order - here i am fixing this
gb = gsub("_B", "", multi.str.split(multi.str.split(a.x.cdh.arrayb$probeset_id, "\\.", 1), "^[A-Z0-9]*_", 2))
ga = gsub("_A", "", multi.str.split(multi.str.split(a.x.cdh.arraya$probeset_id, "\\.", 1), "^[A-Z0-9]*_", 2))

gb = gsub("_", "", multi.str.split(g, "AxC", 2))
ga = gsub("_", "", multi.str.split(g1, "AxC", 2))

coords = match(ga, gb)

#this is the fixed arrayb that matches array a in terms of individuals
a.x.cdh.arraybv2 = a.x.cdh.arrayb[coords, ]

colnames(a.x.cdh.arraybv2)[1:5]

combined.a.x.c = cbind(a.x.cdh.arraya, a.x.cdh.arraybv2[, 3:ncol(a.x.cdh.arraybv2)])

na.omit(match(colnames(axc61_1), colnames(combined.a.x.c)))

#refine probe set to match f2 pop (note not all probes are present in the DH pop)
combined.a.x.c2 = combined.a.x.c[, na.omit(match(colnames(axc61_1), colnames(combined.a.x.c)))]

#add parents from f2 pop
combined.a.x.c3 = rbind(axc61_1[2:3, match(colnames(combined.a.x.c2), colnames(axc61_1))], combined.a.x.c2)
combined.a.x.c4 = combined.a.x.c3[c(3, 1, 2, 4:nrow(combined.a.x.c3)), ]





#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####


#remove columns in which parents are the same

rm.interparental.homozygosities = function(df, parent.rows){
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    
    if(length((which(df[parent.rows[[1]], 3:ncol(df)] == df[parent.rows[[2]], 3:ncol(df)]) + 2)) != 0) {
        df = df[, -(which(df[parent.rows[[1]], 3:ncol(df)] == df[parent.rows[[2]], 3:ncol(df)]) + 2)]    
    }
    return(df)
}

#remove parental heterozygosities
rm.parental.het = function(df, parent.rows){
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    if(length(which(df[parent.rows[1], ] == "H")) != 0) df = df[, -which(df[parent.rows[1], ] == "H")]
    if(length(which(df[parent.rows[2], ] == "H")) != 0) df = df[, -which(df[parent.rows[2], ] == "H")]
    return(df)
}


#remove parental nocalls
rm.parental.nocall = function(df, parent.rows){
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    if(length(which(df[parent.rows[1], ] == "-")) != 0) df = df[, -which(df[parent.rows[1], ] == "-")]
    if(length(which(df[parent.rows[2], ] == "-")) != 0) df = df[, -which(df[parent.rows[2], ] == "-")]
    return(df)
}

cleanup = function(df, parent.rows){
    df = rm.interparental.homozygosities(df, parent.rows)
    df = rm.parental.het(df, parent.rows)
    df = rm.parental.nocall(df, parent.rows)
    return(df)
}

assign.parental.genotypes = function(df, parent.rows){
    #performs "flipping algorithm" to assign parental genotypes to data. 
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    
    #make new df w/ only parents (2 rows total)
    pc2 = df[c(parent.rows[[1]], parent.rows[[2]]), ]
    pc2 = as.data.frame(pc2)
    pc2 = convert.to.character.data.frame(pc2)
    
    df = as.data.frame(df)
    
    #make list of columns to flip
    n = which(lapply(pc2, function(x){
        if(x[1] == "B" & x[2] == "A"){
            return(T)
        } else return(F)
    }) == T)
    
    #flip columns in main dataframe
    df[, n] = lapply(df[, n], function(x){
        g = x
        
        g[g == "A"] = "X"
        g[g == "B"] = "A"
        g[g == "X"] = "B"
        
        return(g)
    })
    
    
    return(df)
}

remove.excess.parents = function(df, parent.name){
    #removes all but one of the parent individuals for one of the parents in the cross (only operates on one parent at a time)
    #df - df in rqtl format
    #parent.name - character string of the parent to perform operation on
    to.remove = grep(parent.name, df$probeset_id)
    df = df[-to.remove[2:length(to.remove)], ]
    return(df)
}


#     ____________________________________________________________________________
#     PROCESS DATA                                                                                                                        ####




combined.geno.data3 = cleanup(cxa65_3, 2:3)

combined.geno.data3 = convert.to.character.data.frame(combined.geno.data3)

combined.geno.data4 = assign.parental.genotypes(combined.geno.data3, 2:3)

# combined.geno.data4 = remove.excess.parents(combined.geno.data4, "Paragon")


write.csv(combined.geno.data4, "rotation1scripts_v4/processed_data/genotypes/avalon.x.cadenza/a.x.c.dh.flipped.probeset.matching.axcf2.csv", row.names = F)
write.csv(combined.geno.data4, "rotation1scripts_v4/processed_data/genotypes/a.x.c.f2/axc51_2.csv", row.names = F)
write.csv(combined.geno.data4, "rotation1scripts_v4/processed_data/genotypes/a.x.c.f2/axc61_1.csv", row.names = F)
write.csv(combined.geno.data4, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa653_redone.csv", row.names = F)



#     ____________________________________________________________________________
#     CHECK MARKERS OF INTEREST                                                                                             ####

markers.seg.5a = c("AX-94405899", "AX-94472591", "AX-94742567", "AX-94519690", "AX-94453668", "AX-94386526", "AX-94729059")

p.x.baj.clean = combined.geno.data4
s(p.x.baj.clean, "rotation1scripts_v4/saved.objects/p.x.baj.clean", "generic.axiom.data.preparation.R")

check.5a.markers = function(geno.data){
    g = geno.data[, unlist(lapply(markers.seg.5a, function(x) which(colnames(geno.data) == x)))]
    
    # browser()
    g1 = lapply(g, function(x){
        g2 = table(x)
        g2[["A"]] / g2[["B"]]
        })    
    
    abs(1 - mean(unlist(g1)))
}

p.x.becard.clean = combined.geno.data4
s(p.x.becard.clean, "rotation1scripts_v4/saved.objects/p.x.becard.clean", "generic.axiom.data.preparation.R")

p.x.super152.clean = combined.geno.data4
s(p.x.super152.clean, "rotation1scripts_v4/saved.objects/p.x.super152.clean", "generic.axiom.data.preparation.R")

p.x.waxwing.clean = combined.geno.data4
s(p.x.waxwing.clean, "rotation1scripts_v4/saved.objects/p.x.waxwing.clean", "generic.axiom.data.preparation.R")

apogee.x.paragon.clean = combined.geno.data4
s(apogee.x.paragon.clean, "rotation1scripts_v4/saved.objects/apogee.x.paragon.clean", "generic.axiom.data.preparation.R")

p.x.wyalkatchem.clean = combined.geno.data4
s(p.x.wyalkatchem.clean, "rotation1scripts_v4/saved.objects/p.x.wyalkatchem.clean", "generic.axiom.data.preparation.R")

p.x.pfau.clean = combined.geno.data4
s(p.x.pfau.clean, "rotation1scripts_v4/saved.objects/p.x.pfau.clean", "generic.axiom.data.preparation.R")

check.5a.markers(p.x.becard.clean)
check.5a.markers(p.x.baj.clean)
check.5a.markers(p.x.super152.clean)
check.5a.markers(p.x.waxwing.clean)
check.5a.markers(p.x.wyalkatchem.clean)
check.5a.markers(p.x.pfau.clean)
check.5a.markers(apogee.x.paragon.clean)

all.crosses = list(p.x.baj.clean, p.x.super152.clean, p.x.waxwing.clean, p.x.wyalkatchem.clean, p.x.pfau.clean, apogee.x.paragon.clean)
names(all.crosses) = c("p.x.baj", "p.x.super152", "p.x.waxwing", "p.x.wyalkatchem", "p.x.pfau", "apogee.x.paragon")

sum.crosses = sort(sapply(all.crosses, function(x){
    check.5a.markers(x)
}, USE.NAMES = T))




#     ____________________________________________________________________________
#     CHECK IF PROBE SEQUENCE MAPS TO GENE                                                                        ####

probe.sequences = readDNAStringSet("rotation1scripts_v4/original_data/fasta/35kprobes.fa")

probe.sequences2 = probe.sequences[unlist(lapply(markers.seg.5a, function(x) which(switch.affy.format(names(probe.sequences)) == x)))]

rev.probe.sequences2 = reverseComplement(probe.sequences2)


probe.align = pairwiseAlignment(rev.probe.sequences2, ms2.gene.cs.cds)
writePairwiseAlignments(probe.align, "rotation1scripts_v4/processed_data/TILLING/rev.probe.align.vs.ms2.txt")



#     ____________________________________________________________________________
#     MISC FUNCTIONS                                                                                                                    ####

replace.aas.genotype.format.w.rqtl.format = function(dataframe){
    dataframe[dataframe == "AA"] = "A"
    dataframe[dataframe == "BB"] = "B"
    dataframe[dataframe == "AB"] = "H"
    dataframe[dataframe == "BA"] = "H"
    dataframe[dataframe == "NoCall"] = "-"
    return(dataframe)
}

grab.columns.in.which.parents.are.concordant = function(dataframe, parent){
    #supply dataframe in rqtl format and name of parent variety (only processes one parent at a time)
    parent.row.coords = grep(parent, dataframe[, 2])
    parent.df = dataframe[parent.row.coords, ]
    parent.df = as.data.frame(parent.df)
    parent.df[] = lapply(parent.df, as.character)
    
    g = lapply(parent.df, function(x){
        b = x[[1]] == x[[2]]
        c = x[[2]] == x[[3]]
        g = all(c(c, b))
        return(g)
    })
    
    which(unlist(g))
    
}

convert.cerealsdb.format.to.rqtl = function(df){
    g = df[, 1:3]
    g = t(g)
    colnames(g) = g[2, ]
    
    g = as.data.frame(g)
    g = add_column(g, probeset_id = g[1, 1], .before = colnames(g)[1])
    g = add_column(g, V1 = 1, .before = colnames(g)[1])
    g = g[3, ]
    g = convert.to.character.data.frame(g)
    g = replace.aas.genotype.format.w.rqtl.format(g)
    return(g)
}
