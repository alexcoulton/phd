


#read files for marker designation for each cxa / axc subpopulation (originally outputted from Axiom Analysis Suite)
marker.conversion.types = lapply(paste0("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/", list.files("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/", pattern = "marker.conversion")), read.delim, sep = "\t", stringsAsFactors = F)

#select only polyhighres markers in each sub population
polyhighres = lapply(marker.conversion.types, function(x){
    x[which(x$ConversionType == "PolyHighResolution"), ]$probeset_id
})

#get markers that are polyhighres in all subpopulations
highly.conservative.marker.list = Reduce(intersect, polyhighres)

#remove some markers that I've manually identified as erroneous (see OneNote 09/07/2019 "Genotyping errors effecting recombination distribution")
highly.conservative.marker.list = highly.conservative.marker.list[-which(highly.conservative.marker.list %in% c("AX-94450697", "AX-94656347", "AX-94756234", "AX-94501087", "AX-94691753", "AX-94529552", "AX-94535421", "AX-94559013", "AX-94623475", "AX-95091073", "AX-95173034", "AX-95244086", "AX-94474254"))]

#write highly conservative marker list to file to be reread back into axiom analysis suite
# snp.list.high.cons = data.frame(highly.conservative.marker.list)
# colnames(snp.list.high.cons) = "probeset_id"
# write.table(snp.list.high.cons, "rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/snp.list.high.cons.txt", sep = "\t", row.names = F, quote = F)

#load genetic map that was made using only cxa651 CEL files as input to axiom analysis suite
cxaf2.c651.asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa.axc.all.replicates.mapv3.newprobeset.csv")
# cxaf2.c651.asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa651-redone-map.csv")

# g = split(cxaf2.c651.asmap.old1, cxaf2.c651.asmap.old1$chromo)
# lapply(g, nrow)
# 
# g2 = split(cxaf2.c651.asmap, cxaf2.c651.asmap$chromo)
# lapply(g2, nrow)


cxaf2.c651.asmap.highly.conserved = cxaf2.c651.asmap[which(cxaf2.c651.asmap$marker %in% highly.conservative.marker.list), ]
cxaf2.c651.asmap.not.hc = cxaf2.c651.asmap[which(!cxaf2.c651.asmap$marker %in% switch.affy.format(highly.conservative.marker.list)), ]



# cxaf2.c651.asmap = cxaf2.c651.asmap.highly.conserved
# cxaf2.c651.asmap$X = switch.affy.format(cxaf2.c651.asmap$X)
cxaf2.c651.asmap$marker = switch.affy.format(cxaf2.c651.asmap$marker)

cxaf2.c651.asmap.list = split.df.into.list.of.dfs.by.column(cxaf2.c651.asmap, "chromo", do.sort = F)

c.x.af2 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.csv", stringsAsFactors = F)

to.rm1 = which(!cxaf2.c651.asmap$marker %in% colnames(c.x.af2))

if(length(to.rm1) > 0) cxaf2.c651.asmap = cxaf2.c651.asmap[-to.rm1, ]




#genotyping data with asmap genetic map order
c.x.af2.asmap.ord = c.x.af2[, c(1, 2, na.omit(match(cxaf2.c651.asmap$marker, colnames(c.x.af2))))]
c.x.af2.asmap.ord[1, 3:ncol(c.x.af2.asmap.ord)] = cxaf2.c651.asmap$chromo

c.x.af2v2 = c.x.af2.asmap.ord

a.x.c512 = read.csv("rotation1scripts_v4/processed_data/genotypes/a.x.c.f2/axc51_2.csv", stringsAsFactors = F, header = T)
a.x.c611 = read.csv("rotation1scripts_v4/processed_data/genotypes/a.x.c.f2/axc61_1.csv", stringsAsFactors = F, header = T)

a.x.c512 = a.x.c512[, match(colnames(c.x.af2v2), colnames(a.x.c512))]
a.x.c512[1, ] = c.x.af2v2[1, ]

a.x.c611 = a.x.c611[, match(colnames(c.x.af2v2), colnames(a.x.c611))]
a.x.c611[1, ] = c.x.af2v2[1, ]

c.x.a651 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa651_redone.csv", stringsAsFactors = F, header = T)
c.x.a653 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa653_redone.csv", stringsAsFactors = F, header = T)

#add chromosome assignments to genotype data
c.x.a651 = c.x.a651[, match(colnames(c.x.af2v2), colnames(c.x.a651))]
c.x.a651[1, ] = c.x.af2v2[1, ]

#add chromosome assignments again
c.x.a653 = c.x.a653[, match(colnames(c.x.af2v2), colnames(c.x.a653))]
c.x.a653[1, ] = c.x.af2v2[1, ]

# write.csv(c.x.a651, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa651_redone_w_chromos.csv", row.names = F)
# write.csv(c.x.a653, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa653_redone_w_chromos.csv", row.names = F)


c.x.acomb = rbind(c.x.a651, c.x.a653[3:nrow(c.x.a653), ])
asmaplgs = unique(cxaf2.c651.asmap$chromo)



#check heterozygosity in genotyping dataframe (any hidden parents?)
# which(sapply(as.data.frame(t(a.x.c611)), function(x){
#     length(which(x == "H"))
# }) < 100)

#remove parents from genotyping data
parent.coords = c(2, 3, grep("M23_CxA 65-1 G12.CEL", c.x.a651$probeset_id), grep("O23_CxA 65-1 H12.CEL", c.x.a651$probeset_id))
c.x.a651 = c.x.a651[-parent.coords, ]

parent.coords2 = c(2, 3, grep("O23_AxC 61-1 H12.CEL", a.x.c611$probeset_id), grep("M23_AxC 61-1 G12.CEL", a.x.c611$probeset_id))
a.x.c611 = a.x.c611[-parent.coords2, ]

c.x.a653 = c.x.a653[-c(2, 3), ]

a.x.c512 = a.x.c512[-c(2, 3), ]

list.of.geno.dfs.each.pop = list(c.x.a651, c.x.a653, a.x.c512, a.x.c611)

#randomly sample individuals from each population so they have the same number of individuals
list.of.geno.dfs.each.pop2 = lapply(list.of.geno.dfs.each.pop, function(x) x[c(1, sample(2:(nrow(x)), 94)), ])

#modify genotyping df to only include highly conserved markers
list.of.geno.dfs.each.pop3 = lapply(list.of.geno.dfs.each.pop2, function(x){
    g = x[, c(1, 2, which(colnames(x) %in% switch.affy.format(highly.conservative.marker.list)))]
    g
})

list.of.geno.dfs.each.pop3 = list.of.geno.dfs.each.pop2

recomb2 = lapply(list.of.geno.dfs.each.pop3, detect.recombination)

#get number of recombination events for each individual
individual.events = lapply(recomb2, function(x){
    sort(sapply(unique(x$individual), function(y){
        # browser()
        sum(x[which(x$individual == y), ]$num.events)
    }))
})

indivi.lar1 = names(which(unlist(individual.events) > 90))


#remove individuals with excessive (erroneous) number of recombination events
list.of.geno.dfs.each.pop = lapply(list.of.geno.dfs.each.pop, function(x){
    to.rm1 = which(x$probeset_id %in% indivi.lar1)
    
    if(length(to.rm1) > 0){
        x = x[-to.rm1, ]
    }
    
    x
})

# lapply(list.of.geno.dfs.each.pop, nrow)

#randomly sample individuals from each population so they have the same number of individuals
# list.of.geno.dfs.each.pop2 = lapply(list.of.geno.dfs.each.pop, function(x){
#     browser()
#     x[c(1, sample(2:(nrow(x)), length(2:(nrow(x) - 1)))), ]
# })

#modify genotyping df to only include highly conserved markers
# list.of.geno.dfs.each.pop3 = lapply(list.of.geno.dfs.each.pop, function(x){
#     g = x[, c(1, 2, which(colnames(x) %in% switch.affy.format(highly.conservative.marker.list)))]
#     g
# })

#make initial recombination dataframes for cluster plot examination
recombination.dfs1 = lapply(list.of.geno.dfs.each.pop3, detect.recombination)
histlists1 = lapply(recombination.dfs1, makehistlist.new)

# setwd("C:/Users/Public/Documents/AxiomAnalysisSuite/Output/AxC 61-1/")



axiom.folders1 = c("E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/c.x.a65-1/", 
                                     "E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/c.x.a65-3/", 
                                     "E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/AxC 51-2/", 
                                     "E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/AxC 61-1/")


asmaplgs = unique(as.character(list.of.geno.dfs.each.pop3[[1]][1, ]))[-1]

axc.lg.lengths = table(as.character(list.of.geno.dfs.each.pop3[[1]][1, ]))
axc.lg.lengths = axc.lg.lengths[-c(1, length(axc.lg.lengths))]