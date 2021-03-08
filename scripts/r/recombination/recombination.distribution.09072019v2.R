library(rmarkdown)
library(gridExtra)
library(readr)
library(ASMap)
library(bookdown)
library(cowplot)
library(data.table)
library(BatchMap)
# library(SNPolisher)
library(zoo) #for na.approx()
library(dplyr)
library(sp)
library(ggplot2)

#setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R")

#### SETUP REPLICATE AxC ####

#read files for marker designation for each cxa / axc subpopulation (originally outputted from Axiom Analysis Suite)
marker.conversion.types = lapply(paste0("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/", list.files("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/", pattern = "marker.conversion")), read.delim, sep = "\t")

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
cxaf2.c651.asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa651-redone-map.csv")

cxaf2.c651.asmap.highly.conserved = cxaf2.c651.asmap[which(cxaf2.c651.asmap$marker %in% switch.affy.format(highly.conservative.marker.list)), ]

cxaf2.c651.asmap$X = switch.affy.format(cxaf2.c651.asmap$X)
cxaf2.c651.asmap$marker = switch.affy.format(cxaf2.c651.asmap$marker)

cxaf2.c651.asmap.list = split.df.into.list.of.dfs.by.column(cxaf2.c651.asmap, "chromo", do.sort = F)

c.x.af2 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.csv", stringsAsFactors = F)

#genotyping data with asmap genetic map order
c.x.af2.asmap.ord = c.x.af2[, c(1, 2, match(cxaf2.c651.asmap$marker, colnames(c.x.af2)))]
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
list.of.geno.dfs.each.pop3 = lapply(list.of.geno.dfs.each.pop, function(x){
    g = x[, c(1, 2, which(colnames(x) %in% switch.affy.format(highly.conservative.marker.list)))]
    g
})

#make initial recombination dataframes for cluster plot examination
recombination.dfs1 = lapply(list.of.geno.dfs.each.pop3, detect.recombination)
histlists1 = lapply(recombination.dfs1, makehistlist.new)

# setwd("C:/Users/Public/Documents/AxiomAnalysisSuite/Output/AxC 61-1/")

axiom.folders1 = c("E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/c.x.a65-1/", 
                                     "E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/c.x.a65-3/", 
                                     "E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/AxC 51-2/", 
                                     "E:/windows.backup/Users/Public/Documents/AxiomAnalysisSuite/Output/AxC 61-1/")


geno.to.change = Map(function(x, y){
    # browser()
    evaluate.cluster.plots(x, y)    
}, recombination.dfs1, axiom.folders1)

# .orig = without change in chromosome order to be numerical/alphabetical
l.geno.dfs3.w.change.orig = Map(function(x, y){
    for(z in 1:nrow(y)){
        x[which(x$probeset_id == y[z, 1]), which(colnames(x) == y[z, 2])] = "-"
    }
    x
}, list.of.geno.dfs.each.pop3, geno.to.change)

l.geno.dfs3.w.change = Map(function(x, y){
    for(z in 1:nrow(y)){
        x[which(x$probeset_id == y[z, 1]), which(colnames(x) == y[z, 2])] = "-"
    }
    x
}, list.of.geno.dfs.each.pop3, geno.to.change)
l.geno.dfs3.w.change = lapply(l.geno.dfs3.w.change, put.geno.in.order.of.chromo)

#remove erroneous-looking chromosomes
l.geno.dfs3.w.change = lapply(l.geno.dfs3.w.change, function(x){
    x[, -which(as.character(x[1, ]) %in% c("4B", "5A", "5B"))]
})



g = lapply(l.geno.dfs3.w.change, function(x) as.character(x[1, ]))

asmaplgs = unique(as.character(l.geno.dfs3.w.change[[1]][1, ]))[-1]

axc.lg.lengths = table(as.character(l.geno.dfs3.w.change[[1]][1, ]))
axc.lg.lengths = axc.lg.lengths[-c(1, length(axc.lg.lengths))]


#test whether geno.to.change worked
# q = which(l.geno.dfs3.w.change[[1]]$probeset_id == geno.to.change[[1]][2, 1])
# q2 = which(colnames(l.geno.dfs3.w.change[[1]]) == geno.to.change[[1]][2, 2])
# 
# 
# l.geno.dfs3.w.change[[1]][q, (q2 - 3):(q2 + 3)]
# list.of.geno.dfs.each.pop3[[1]][q, (q2 - 3):(q2 + 3)]

#### PREPARE HIGHER QC GENOTYPING DATA FOR AxC / CxA ####

#examine lower qc sample metrics


low.qc.metrics.files = paste0("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/lower.qc.metrics/", 
                                                    list.files("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/lower.qc.metrics/", pattern = ".txt"))

low.qc.metrics.files = low.qc.metrics.files[c(3, 4, 1, 2)]
library(readr)
lowqc.met2 = lapply(low.qc.metrics.files, read_delim, delim = "\t")

lowqc.met3 = bind_rows(lowqc.met2)

source("rotation1scripts_v4/scripts/r/recombination/axc.analysis.R")

#average mrd across chromosomes for all individuals
mrd.avg = aggregate(axc.avg2$phys.dist, list(factor(axc.avg2$individual, levels = unique(axc.avg2$individual))), function(x) mean(na.omit(x)))
mrd.avg$Group.1 = multi.str.split(as.character(mrd.avg$Group.1), "_call_code", 1)

axc.avg.chromo = arrange(axc.avg2, chromo)
axc.avg.chromo2 = split.df.into.list.of.dfs.by.column(axc.avg.chromo, "chromo")

#match up phys.dist and dqc/qc information
axc.avg.chromo3 = lapply(axc.avg.chromo2, function(x){
    x$individual = multi.str.split(as.character(x$individual), "_call_code", 1)
    lowqc.temp = lowqc.met3[-which(!lowqc.met3$`Sample Filename` %in% x$individual), ]
    lowqc.temp2 = bind_cols(lowqc.temp, x[match(lowqc.temp$`Sample Filename`, x$individual), ])
    lowqc.temp2[c("Sample Filename", "DQC", "QC call_rate", "phys.dist", "marker.dist", "chromo")]
})

#linear models of MRD vs DQC
lm1 = lapply(axc.avg.chromo3, function(x){
    # browser()
    summary(lm(x$phys.dist ~ x$DQC))
})

lm2 = lapply(axc.avg.chromo3, function(x){
    # browser()
    summary(lm(x$phys.dist ~ x$`QC call_rate`))
})


dqcp1 = p.adjust(sapply(lm1, function(x){
    x[[4]][2, 4]
}), method = "fdr")


dqcr1 = round(unlist(sapply(lm1, function(x){
    x[9]
})), digits = 4)

qcp2 = p.adjust(sapply(lm2, function(x){
    x[[4]][2, 4]
}), method = "fdr")

qcr2 = round(unlist(sapply(lm2, function(x){
    x[9]
})), digits = 4)

data.frame(dqcp1, dqcr1, qcp2, qcr2)

count1 = 1
dqc.mrd.plot1 = lapply(axc.avg.chromo3, function(x){
    g = ggplot(x, aes(x = phys.dist, y = DQC)) + geom_point() + theme_bw() + ggtitle(asmaplgs[count1])
    count1 <<- count1 + 1
    g
})

count1 = 1
qc.mrd.plot2 = lapply(axc.avg.chromo3, function(x){
    g = ggplot(x, aes(x = phys.dist, y = `QC call_rate`)) + geom_point() + theme_bw() + ggtitle(asmaplgs[count1])
    count1 <<- count1 + 1
    g
})



# axcc2 = lapply(axc.avg.chromo2, compare.qc.phys.dist, qc.type = 2)
# 
# do.call(grid.arrange, axcc1)
# 
# do.call(grid.arrange, axcc2)



lowqc.met3$mrd = ""

lowqc.met4 = lowqc.met3[-which(!lowqc.met3$`Sample Filename` %in% mrd.avg$Group.1), ]

lowqc.met5 = bind_cols(lowqc.met4, mrd.avg[match(lowqc.met4$`Sample Filename`, mrd.avg$Group.1), ])


hist(lowqc.met5$x)

lm1 = lm(lowqc.met5$x ~ lowqc.met5$DQC)
summary(lm1)

lm2 = lm(lowqc.met5$x ~ lowqc.met5$`QC call_rate`)
summary(lm2)





lm2 = lm(lowqc.met5$x ~ lowqc.met5$`QC call_rate`)
summary(lm2)


geno.order




#examine sample QC metrics

qc.metrics.files = paste0("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.sample.metrics/", 
                                                    list.files("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.sample.metrics/", pattern = ".txt"))

qc.metrics.files = qc.metrics.files[c(3, 4, 1, 2)]

qc.met2 = lapply(qc.metrics.files, read_delim, delim = "\t")

geno.order = c("c.x.a651", "c.x.a653", "a.x.c512", "a.x.c611")

result.qc.met3 = sapply(qc.met2, function(x){
    dqc = c(mean(x$DQC), sd(x$DQC))
    qc.cr = c(mean(x$`QC call_rate`), sd(x$`QC call_rate`))
    
    qc.pass = mean(x[which(x$`Pass/Fail` == "Pass"), ]$`QC call_rate`)
    
    data.frame(dqc[1], dqc[2], qc.cr[1], qc.cr[2], qc.pass)
    
})

colnames(result.qc.met3) = geno.order
rownames(result.qc.met3) = c("dqc.mean", "dqc.sd", "qc.cr.mean", "qc.cr.sd", "avg.qc.for.passing.samples")




qc.met4 = comb.treatments(qc.met2)
qc.met4$treatment = as.factor(qc.met4$treatment)



dqc.hist = ggplot(qc.met4, aes(x = DQC)) + geom_histogram(bins = 100) + facet_grid(vars(treatment)) + theme_bw()

printplot("dqc.hist2", dqc.hist)

exportresults(result.qc.met3, "result.qc.met3")

dqc.anova = aov(DQC ~ treatment, data = qc.met4)
summary(dqc.anova)

TukeyHSD(dqc.anova)



hist(qc.met2[[1]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
hist(qc.met2[[2]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
hist(qc.met2[[3]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
hist(qc.met2[[4]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)






#are dqc values significantly different between treatments?
kruskal.test(sapply(qc.met2, function(x) x$DQC))

wilcox.test(sapply(qc.met2, function(x) x$DQC)[[1]], sapply(qc.met2, function(x) x$DQC)[[3]])


#examine SNP QC metrics

snp.qc.metrics.files = paste0("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.probe.metrics/", 
                                                    list.files("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.probe.metrics//", pattern = ".txt"))

snp.qc.metrics.files = snp.qc.metrics.files[c(3, 4, 1, 2)]

snp.met2 = lapply(snp.qc.metrics.files, read.table, sep = "\t", header = T)

g = read.table("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.probe.metrics/axc512.txt", sep = "\t")

g[which(g$V1 == "AX-94796264"), ]

snp.met3 = lapply(snp.met2, function(x){
    x[which(x$probeset_id %in% highly.conservative.marker.list), ]
})

lapply(snp.met3, function(x){
    mean(na.exclude(x$HomFLD))
})

lapply(snp.met3, function(x){
    mean(na.exclude(x$FLD))
})

lapply(snp.met3, function(x){
    mean(na.exclude(x$HetSO))
})

lapply(snp.met3, function(x){
    mean(na.exclude(x$HomRO))
})


#prepare genotyping data for higher qc samples

higher.qc.genotypes = paste0("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/", 
             list.files("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/", pattern = ".txt"))

higher.qc.genotypes = higher.qc.genotypes[c(3, 4, 1, 2)]

higher.qc2 = lapply(higher.qc.genotypes, convert.aas.to.rqtl)

higher.qc3 = lapply(higher.qc2, function(x){
    x = x[, match(switch.affy.format(colnames(l.geno.dfs3.w.change[[1]])), colnames(x))]
    x
})

#get parent coordinates
p1 = grep("M23", higher.qc3[[1]]$probeset_id)
p2 = grep("O23", higher.qc3[[1]]$probeset_id)
higher.qc3[[1]] = higher.qc3[[1]][c(1, p1, p2, 2:76, 78:88), ]

#add parents to all genotyping dataframes
higher.qc4 = lapply(2:4, function(x){
    x = higher.qc3[[x]]
    g1 = convert.to.character.data.frame(as.data.frame(rbind(higher.qc3[[1]][2:3, ], x)))
    g1[c(3, 1:2, 4:nrow(g1)), ]
})

higher.qc4 = c(list(higher.qc3[[1]]), higher.qc4)

which(higher.qc4[[1]][2, ] == higher.qc4[[1]][3, ])
which(higher.qc4[[1]][2, ] == "-")
which(higher.qc4[[1]][2, ] == "H")
which(higher.qc4[[1]][3, ] == "-")
which(higher.qc4[[1]][3, ] == "H")

#assign.parental.genotypes from generic.axiom.data.preparation.R

higher.qc5 = lapply(higher.qc4, cleanup, parent.rows = 2:3)
higher.qc5 = lapply(higher.qc5, convert.to.character.data.frame)
higher.qc6 = lapply(higher.qc5, assign.parental.genotypes, parent.rows = 2:3)


chromos1 = na.exclude(cxaf2.c651.asmap$chromo[match(switch.affy.format(colnames(higher.qc6[[1]])), cxaf2.c651.asmap$marker)])

#chromosome assignments
higher.qc7 = lapply(higher.qc6, function(x){
    x[1, 3:ncol(x)] = chromos1
    colnames(x)[3:ncol(x)] = switch.affy.format(colnames(x)[3:ncol(x)])
    x
})





#make random mixtures of populations

c.x.a.comb = as.data.frame(rbind(c.x.a651[2:nrow(c.x.a651), ], c.x.a653[2:nrow(c.x.a653), ]))

vec1 = 1:190

samp1 = sample(190, 95)
samp2 = vec1[-match(samp1, vec1)]

c.x.a.random1 = as.data.frame(rbind(c.x.a651[1, ], c.x.a.comb[samp1, ]))
c.x.a.random2 = as.data.frame(rbind(c.x.a651[1, ], c.x.a.comb[samp2, ]))

a.x.c.comb = as.data.frame(rbind(a.x.c512[2:nrow(a.x.c512), ], a.x.c611[2:nrow(a.x.c611), ]))

vec1 = 1:190

samp1 = sample(190, 95)
samp2 = vec1[-match(samp1, vec1)]

a.x.c.random1 = as.data.frame(rbind(a.x.c512[1, ], a.x.c.comb[samp1, ]))
a.x.c.random2 = as.data.frame(rbind(a.x.c512[1, ], a.x.c.comb[samp2, ]))

list.of.geno.dfs.each.pop.random = list(c.x.a.random1, c.x.a.random2, a.x.c.random1, a.x.c.random2)

#remove individuals with excessive (erroneous) number of recombination events
list.of.geno.dfs.each.pop.random = lapply(list.of.geno.dfs.each.pop.random, function(x){
    to.rm1 = which(x$probeset_id %in% indivi.lar1)
    
    if(length(to.rm1) > 0){
        x = x[-to.rm1, ]
    }
    
    x
})


#setting all linkage groups to the same to see if any individuals are skewing recombination distribution by a large amount
list.of.geno.dfs.each.pop.random.all.one.lg = lapply(list.of.geno.dfs.each.pop.random, function(x){
    x[1, 3:ncol(x)] = "1"
    x
})

#### PREPARE COMPARISON OF PARENTS ON AxC ARRAY TO OTHER AxC GENOTYPES ####

#first do comparison for the axc611 array

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')
avalon = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Avalon'")
avalon2 = fetch(avalon, n=-1)

cadenza = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Cadenza'")
cadenza2 = fetch(cadenza, n=-1)


avalon651 = read.delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/cxa651.txt", sep = "\t", comment.char = "#")
avalon651.2 = avalon651[, c(1, grep("O23", colnames(avalon651)))]

cadenza651.2 = avalon651[, c(1, grep("M23", colnames(avalon651)))]


avalon611 = read.delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/axc611.low.qc.txt", sep = "\t", comment.char = "#")
avalon611.2 = avalon611[, c(1, grep("O23", colnames(avalon611)))]

cadenza611.2 = avalon611[, c(1, grep("M23", colnames(avalon611)))]

# match(avalon2$Probe_row, avalon611.2$probeset_id)

# length(which(avalon2$Matrix_value == avalon611.2$O23_AxC.61.1.H12.CEL_call_code)) / 35143

# length(which(cadenza2$Matrix_value == cadenza611.2$M23_AxC.61.1.G12.CEL_call_code)) / 35143

# length(which(avalon2$Matrix_value == avalon651.2$O23_CxA.65.1.H12.CEL_call_code)) / 35143

# length(which(cadenza2$Matrix_value == cadenza651.2$M23_CxA.65.1.G12.CEL_call_code)) / 35143

t11 = as.data.frame(table(cadenza2$Matrix_value))
t11$var = "cadenza_cereals"
t12 = as.data.frame(table(avalon2$Matrix_value))
t12$var = "avalon_cereals"

t13 = as.data.frame(table(cadenza611.2$M23_AxC.61.1.G12.CEL_call_code))
t13$var = "cadenza_axc_array"
t14 = as.data.frame(table(avalon611.2$O23_AxC.61.1.H12.CEL_call_code))
t14$var = "avalon_axc_array"

t15 = as.data.frame(table(cadenza651.2$M23_CxA.65.1.G12.CEL_call_code))
t15$var = "cadenza_cxa_array"

t16 = as.data.frame(table(avalon651.2$O23_CxA.65.1.H12.CEL_call_code))
t16$var = "avalon_cxa_array"

geno.freq.comp = bind_rows(t11, t12, t13, t14, t15, t16)
colnames(geno.freq.comp) = c("Genotype", "Frequency", "Sample")
geno.freq.comp = dcast(geno.freq.comp, Genotype ~ Sample, value.var = "Frequency")
geno.freq.comp = geno.freq.comp[, c(1, 3, 6, 2, 4, 5, 7)]
colnames(geno.freq.comp) = c("Genotype", "Avalon Cereals", "Cadenza Cereals", "Avalon AxC", "Avalon CxA", "Cadenza AxC", "Cadenza CxA")


exportresults(bind_rows(t11, t12, t13, t14, t15, t16), "genotype.frequency.comparison")


#do comparison only for highly.conservative.marker.list

avalon651.3 = avalon651.2[which(avalon651.2$probeset_id %in% highly.conservative.marker.list), ]
cadenza651.3 = cadenza651.2[which(cadenza651.2$probeset_id %in% highly.conservative.marker.list), ]

avalon3 = avalon2[which(avalon2$Probe_row %in% highly.conservative.marker.list), ]
avalon611.3 = avalon611.2[which(avalon611.2$probeset_id %in% highly.conservative.marker.list), ]





qq1 = length(which(avalon3$Matrix_value == avalon611.3$O23_AxC.61.1.H12.CEL_call_code))
qq2 = length(highly.conservative.marker.list)

qq1 / qq2

cadenza3 = cadenza2[which(cadenza2$Probe_row %in% highly.conservative.marker.list), ]
cadenza611.3 = cadenza611.2[which(cadenza611.2$probeset_id %in% highly.conservative.marker.list), ]

qq3 = length(which(cadenza3$Matrix_value == cadenza611.3$M23_AxC.61.1.G12.CEL_call_code))
qq4 = length(highly.conservative.marker.list)

qq3 / qq4


mismatch.coords1 = which(cadenza3$Matrix_value != cadenza611.3$M23_AxC.61.1.G12.CEL_call_code)

cad3.1 = cadenza3[mismatch.coords1, ]
cad3.2 = cadenza611.3[mismatch.coords1, ]
cad3.3 = bind_cols(cad3.1, cad3.2)
t22 = as.data.frame(table(cad3.3$Matrix_value))
t22$var = "cad_cerealsdb"
t23 = as.data.frame(table(cad3.3$M23_AxC.61.1.G12.CEL_call_code))
t23$var = "cad_axc_array"

mismatch.coords2 = which(avalon3$Matrix_value != avalon611.3$O23_AxC.61.1.H12.CEL_call_code)

ava3.1 = avalon3[mismatch.coords2, ]
ava3.2 = avalon611.3[mismatch.coords2, ]
ava3.3 = bind_cols(ava3.1, ava3.2)
t24 = as.data.frame(table(ava3.3$Matrix_value))
t24$var = "ava_cerealsdb"
t25 = as.data.frame(table(ava3.3$O23_AxC.61.1.H12.CEL_call_code))
t25$var = "ava_axc_array"




qq5 = length(which(avalon3$Matrix_value == avalon651.3$O23_CxA.65.1.H12.CEL_call_code))

qq5 / qq4

qq6 = length(which(cadenza3$Matrix_value == cadenza651.3$M23_CxA.65.1.G12.CEL_call_code))

qq6 / qq4

length(which(cadenza651.3$M23_CxA.65.1.G12.CEL_call_code == cadenza611.3$M23_AxC.61.1.G12.CEL_call_code)) / nrow(avalon651.3)

length(which(avalon651.3$O23_CxA.65.1.H12.CEL_call_code == avalon611.3$O23_AxC.61.1.H12.CEL_call_code)) / nrow(cadenza651.3)



mismatch.coords3 = which(avalon3$Matrix_value != avalon651.3$O23_CxA.65.1.H12.CEL_call_code)

ava3.11 = avalon3[mismatch.coords3, ]
ava3.12 = avalon651.3[mismatch.coords3, ]
ava3.13 = bind_cols(ava3.11, ava3.12)

t26 = as.data.frame(table(ava3.13$O23_CxA.65.1.H12.CEL_call_code))
t26$var = "ava_cxa_array"


mismatch.coords4 = which(cadenza3$Matrix_value != cadenza651.3$M23_CxA.65.1.G12.CEL_call_code)

cad3.11 = cadenza3[mismatch.coords4, ]
cad3.12 = cadenza651.3[mismatch.coords4, ]
cad3.13 = bind_cols(cad3.11, cad3.12)

t28 = as.data.frame(table(cad3.13$M23_CxA.65.1.G12.CEL_call_code))
t28$var = "cad_cxa_array"

geno.freq.comp2 = bind_rows(t22, t24, t23, t25, t26, t28)
colnames(geno.freq.comp2) = c("Genotype", "Frequency", "Sample")
library(reshape2)
geno.freq.comp2 = dcast(geno.freq.comp2, Genotype ~ Sample, value.var = "Frequency")
geno.freq.comp2 = geno.freq.comp2[, c(1, 3, 6, 2, 4, 5, 7)]
colnames(geno.freq.comp2) = c("Genotype", "Avalon Cereals", "Cadenza Cereals", "Avalon AxC", "Avalon CxA", "Cadenza AxC", "Cadenza CxA")

exportresults(geno.freq.comp2, "comparison.genotypes.mismatch.markers.highly.cons")




#compare chinese spring samples for reference

cs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Chinese%'")
cs2 = fetch(cs, n=-1)

unique(cs2$Var_col)

cs3 = cs2[which(cs2$Var_col == "Chinese_Spring"), ]

cs4 = cs2[which(cs2$Var_col == "192_Chinese_Spring"), ]


length(which(cs3$Matrix_value == cs4$Matrix_value)) / 35143





#### EXAMINE POSTERIOR AREAS ####

posteriors1 = lapply(axiom.folders1, examine.posteriors)
# s(posteriors1, "rotation1scripts_v4/saved.objects/posteriors1", "recombination.distribution.09072019v2.R")

post2 = lapply(posteriors1, function(x){
    colnames(x)[1] = "marker"
    x[which(x$marker %in% highly.conservative.marker.list), ]
})

treatments1 = c("c.x.a651", "c.x.a653", "a.x.c512", "a.x.c611")

for(i in 1:4){
    post2[[i]]$treatment = treatments1[i]
}

post3 = bind_rows(post2)
post3$treatment = as.factor(post3$treatment)

#run some statistics on the areas of the ovals
post3.lm = aov(log(post3$area) ~ post3$treatment)

TukeyHSD(post3.lm)

anova(post3.lm)





#are posterior areas bigger in the third treatment compared to other treatments - indicating more genotyping error?
# answer: no, they are smaller
lapply(post2, function(x){
    mean(x$area)
})


kruskal.test(post2[[1]]$area, post2[[2]]$area, post2[[3]]$area, post2[[4]]$area)

t.test(log(post2[[1]]$area), log(post2[[2]]$area))

hist(log(post2[[1]]$area))
hist(log(post2[[2]]$area))
hist(log(post2[[3]]$area))
hist(log(post2[[4]]$area))

#### SAVANNAH X RIALTO #### 

library(data.table)




s.x.r.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/sxr.allen.csv", stringsAsFactors = F)

sxr2 = prepare.allen.map(s.x.r.allen)

#get physical positions for markers without interpolating between NAs
sxr.phys1 = get.phys.for.chromos(sxr2, F, T, "-")

#extract genetic map from allen genotyping dataframe
sxr.gen.map = get.allen.genetic.map(s.x.r.allen, sxr2)

sxr.qtl1 = qtlanalysis.all(sxr2, sxr.gen.map, "RIL10")
sxr.qtl1.step0 = qtlanalysis.all(sxr2, sxr.gen.map, "RIL10", geno.prob.step = 0)

#analyse recombination frequency without taking LISS (longest increasing subsequence)
sxr.qtl1.step0.noliss = qtlanalysis.all(sxr2, sxr.gen.map, "RIL10", geno.prob.step = 0, skip.liss = T)

lapply(sxr.qtl1.step0.noliss[[1]], function(x) x[[1]])


sxr.qtl1[[1]][[4]]



sxr.qtl1[[1]][[4]]



sxr.qtl1[[1]][[4]][[2]]


sxr.qtl.anal12 = qtlanalysis.all(sxr2, sxr.gen.map, "RIL5")

lapply(sxr.qtl.anal12[[1]], function(x) x[[1]])

sxr.qtl.anal12[[1]][[4]][[2]]




#add physical positions to genetic map dataframe
sxr.gen.map$phys = sxr.phys1$phys.pos[match(sxr.gen.map$marker, sxr.phys1$marker)]


sxr.gen.map2 = order.by.physical.dist.within.bins(sxr.gen.map, "chr", "cM", "phys")

#get longest increasing subsequence of physical positions
sxr.gen.map3 = split.df.into.list.of.dfs.by.column(sxr.gen.map2, "chr", F)
sxr.gen.map4 = lapply(sxr.gen.map3, function(x){
    # browser()
    x = x[-which(is.na(x$phys)), ]
    x = x[longest_subseq.R(na.omit(x$phys)), ]
    x
})

sxr.gen.map5 = bind_rows(sxr.gen.map4)


#grab skeleton markers from longest increasing subsequence (take the central marker in each bin)
#need to take only skeleton markers as QTL analysis downstream doesn't like markers in the same position
sxr.gen.map.ske = only.skeleton.markers(sxr.gen.map5)


#make new genotyping dataframe containing only skeleton markers in longest increasing subsequence of physical positions
sxr.ls = sxr2[, c(1, 2, match(sxr.gen.map.ske$marker, colnames(sxr2)))]

#recalculate centimorgan distances between markers
sxr.cm2 = extract.centimorgan.from.genotypes(sxr.ls)

#are any of the markers in the savannah x rialto map monomorphic?
t.len2 = sapply(sxr.ls[2:nrow(sxr.ls), 3:ncol(sxr.ls)], function(x){
    if(length(which(x == "B")) == 0 | length(which(x == "A")) == 0){
        return(0)
    } else {
        return(1)
    }
})

which(t.len2 == 0) #no, they are not


sxr.phys = get.phys.for.chromos(sxr.ls, F, T, "-")

sxr.recomb = detect.recombination(sxr.ls)

sxr.recomb$phys.pos = sxr.phys[match(sxr.recomb$probe.before.transition, sxr.phys$marker), ]$phys.pos

sxr.ind = ind.avg(sxr.recomb)

sxr.freq = get.recomb.freq(sxr.recomb)

sxr.qtl = prepareqtl(sxr.ls, sxr.ind, sxr.freq, 6, sxr.gen.map.ske)

sxr.qtl[[7]][[2]]
sxr.qtl[[10]][[2]]
sxr.qtl[[1]][[2]]

sxr.qtl[[6]][[4]]

#analyse qtls

#check if there are any significant qtls
lapply(sxr.qtl, function(x) x[[1]])
#get phenotypes
sapply(sxr.qtl, function(x) x[[3]])

#get markers that are significant qtls
qtl1 = sxr.qtl[[6]][[4]]
qtl1.real.markers = qtl1[grep("AX-[0-9]*", rownames(qtl1)), ]

#get rqtl object
qtl1c = sxr.qtl[[6]][[5]]

#make multi-qtl model with qtl-qtl interactions
qtl2 = makeqtl(qtl1c, chr = c("3B", "3B", "3B"), pos = c(11.15, 12.75, 81.54), what = "prob")
qtl2 = makeqtl(qtl1c, chr = qtl1.real.markers$chr, pos = qtl1.real.markers$pos, what = "prob")
fitqtl(qtl1c, qtl = qtl2, formula = y~Q1*Q2, method = "hk", pheno.col = 6)


#get additive effects
summary(fitqtl(qtl1c, qtl = qtl2, formula = y~Q1+Q2+Q3, method = "hk", pheno.col = 6, get.ests = T, dropone = F))

#get effect strengths of qtls
eff2 = effectplot(sxr.qtl[[6]][[5]], mname1 = "AX-94754126", pheno.col = 6, draw = F)

effectplot(sxr.qtl[[6]][[5]], mname1 = "AX-94754126", pheno.col = 6)

effectplot(sxr.qtl[[1]][[5]], mname1 = "AX-95206556", pheno.col = 1)


eff = effectplot(sxr.qtl[[1]][[5]], mname1 = "AX-95206556", pheno.col = 1, draw = F)

pdf("rotation1scripts_v4/plots/recombination.distribution11072019/sxr.mdr.qtl3b.pdf", 10, 7)
sxr.qtl[[6]][[2]]
dev.off()





#### OPATA X SYNTHETIC #### 

oxs.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/o.x.s.allen.csv", stringsAsFactors = F)

oxs2 = prepare.allen.map(oxs.allen)

oxs.map = get.allen.genetic.map(oxs.allen, oxs2)

# oxs.qtl1 = qtlanalysis.all(oxs2, oxs.map, "RIL7")
oxs.qtl1.step0 = qtlanalysis.all(oxs2, oxs.map, "RIL7", geno.prob.step = 0)
#analyse recombination frequency
oxs.qtl1.step0.noliss = qtlanalysis.all(oxs2, oxs.map, "RIL7", geno.prob.step = 0, skip.liss = T)

oxs.qtl1.step0.noliss[[1]][[1]]

v(oxs.qtl1[[2]])


lapply(oxs.qtl1[[1]], function(x) x[[1]])



#### CHINESE SPRING X PARAGON ####


csxp.map = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.cs.x.p.map.csv", stringsAsFactors = F)
csxp.map.orig = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.cs.x.p.map.csv", stringsAsFactors = F)

csxp.geno = prepare.allen.map(csxp.map.orig)
csxp.geno.nofilter = prepare.allen.map(csxp.map.orig, F)

csxp.map1 = get.allen.genetic.map(csxp.map.orig, csxp.geno)
csxp.map1.nofilter = get.allen.genetic.map(csxp.map.orig, csxp.geno.nofilter)

# csxp.qtl1.1 = qtlanalysis.all(csxp.geno, csxp.map1, "RIL7", geno.prob.step = 1)
# csxp.qtl1.1v2 = qtlanalysis.all(csxp.geno, csxp.map1, "RIL7", geno.prob.step = 1)

csxp.qtl1.0 = qtlanalysis.all(csxp.geno, csxp.map1, geno.prob.step = 0)
csxp.qtl1.0.noliss = qtlanalysis.all(csxp.geno, csxp.map1, geno.prob.step = 0, skip.liss = T)

csxp.qtl1.0.noliss.nofilter = qtlanalysis.all(csxp.geno.nofilter, csxp.map1.nofilter, geno.prob.step = 0, skip.liss = T)

lapply(csxp.qtl1.0[[1]], function(x) x[[1]])
lapply(csxp.qtl1.0.noliss[[1]], function(x) x[[1]])
lapply(csxp.qtl1.0.noliss.nofilter[[1]], function(x) x[[1]])

csxp.qtl1.0.noliss.nofilter[[1]][[1]][[2]]


lapply(csxp.qtl1.1[[1]], function(x) x[[1]])
lapply(csxp.qtl1.1v2[[1]], function(x) x[[1]])



csxp.liss.map = qtlanalysis.all(csxp.geno, csxp.map1, geno.prob.step = 0, skip.qtl.anal = T)


csxp.qtl1.0[[1]][[1]][[2]]

csxp.qtl1.0[[1]][[1]]

csxp.qtl1.0[[1]][[12]][[2]]


png("rotation1scripts_v4/plots/recombination.distribution11072019/lod.calcgeno.prob.steps.png", 1500, 1000)
grid.arrange(csxp.qtl1.0[[1]][[16]][[2]], csxp.qtl1.1[[1]][[16]][[2]])
dev.off()

csxp.qtl1.1[[1]][[5]][[2]]


csxp.qtl1[[1]][[12]][[2]]

csxp.sig.plots1 = list(csxp.qtl1[[1]][[5]][[2]], csxp.qtl1[[1]][[6]][[2]], csxp.qtl1[[1]][[7]][[2]], csxp.qtl1[[1]][[9]][[2]], csxp.qtl1[[1]][[13]][[2]], csxp.qtl1[[1]][[15]][[2]], csxp.qtl1[[1]][[16]][[2]])

do.call(grid.arrange, csxp.sig.plots1)


v(csxp.qtl1[[2]])


#### AVALON X CADENZA #### 

axc.map = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.axc.csv")
axc.geno = prepare.allen.map(axc.map, F)

axc.geno[1, 1] = ""

write.csv(axc.geno, "rotation1scripts_v4/temp/dhrqtltest/axcgeno.csv", row.names = F)
g1 = rqtl.read("rotation1scripts_v4/temp/dhrqtltest/axcgeno.csv", crosstype = "dh")

read.cross("csv", file = "rotation1scripts_v4/temp/dhrqtltest/axcgeno.csv", crosstype = "dh")



write.cross(g1, "csv", "rotation1scripts_v4/temp/dhrqtltest/test2")


axc.map1 = get.allen.genetic.map(axc.map, axc.geno)

# axc.qtl1 = qtlanalysis.all(axc.geno, axc.map1, "RIL10")


# axc.qtl1.step0 = qtlanalysis.all(axc.geno, axc.map1, "RIL10", geno.prob.step = 0, fgen = "dh")

lapply(axc.qtl1.step0[[1]], function(x) x[[1]])

axc.qtl1.step0[[1]][[6]][[2]]
axc.qtl1.step0[[1]][[17]][[2]]

v(axc.qtl1.step0[[2]])

# axc.qtl1.step0.noliss = qtlanalysis.all(axc.geno, axc.map1, "RIL10", geno.prob.step = 0, skip.liss = T, fgen = "dh")

lapply(axc.qtl1[[1]], function(x) x[[1]])
lapply(axc.qtl1.step0[[1]], function(x) x[[1]])

lapply(axc.qtl1.step0.noliss[[1]], function(x) x[[1]])
axc.qtl1.step0.noliss[[1]][[1]][[2]]

axc.qtl1[[1]][[17]][[2]]
axc.qtl1[[1]][[17]]

# 
# colnames(csxp.map.orig)[1] = "marker"
# csxp.map.orig2 = split.df.into.list.of.dfs.by.column(csxp.map.orig, "chr")
# 
# csxp.phys = lapply(csxp.map.orig2, function(x){
#     genetictophysical(unique(x$chr), x$marker, "-")
# })
# 
# lapply(csxp.phys, function(x){
#     na.omit(x)
# })
# 
# 
# find.max.gap = function(phys.positions){
#     #examines how well distributed the SNPs are across the chromosome by physical position
#     #the smaller the number the better (0 = perfect, 100 = worst)
#     max(unlist(lapply(1:100, function(x){
#         min(abs(phys.positions - x))
#     })))
# }
# 
# concordance = sapply(csxp.phys, function(y){
#     find.max.gap(na.omit(y))
# })
# 
# chr.to.keep = names(which(concordance <= 18))
# 
# 
# 
# 
# 
# 
# 
# csxp.map = as.data.frame(t(csxp.map))
# csxp.map = convert.to.character.data.frame(csxp.map)
# colnames(csxp.map) = csxp.map[1, ]
# csxp.map = csxp.map[-1, ]
# csxp.map = csxp.map[-2, ]
# csxp.map = cbind(csxp.map[, 1:2], csxp.map)
# csxp.map[, 2] = rownames(csxp.map)
# csxp.map[, 1] = ""
# csxp.map[2:nrow(csxp.map), 1] = 2:nrow(csxp.map)
# csxp.map[1, 2] = ""
# csxp.map[1, 1] = NA
# coln1 = colnames(csxp.map)[1:2]
# colnames(csxp.map)[1:2] = ""
# colnames(csxp.map)[3:4] = coln1
# csxp.map = csxp.map[, c(1, 2, which(csxp.map[1, ] %in% chr.to.keep))]
# colnames(csxp.map)[2] = "probeset_id"
# 
# csxp.map.phys1 = get.phys.for.chromos(csxp.map, T, T, ".")
# 
# 
# csxp.recomb = detect.recombination(csxp.map)
# 
# csxp.recomb$phys.pos = csxp.map.phys1[match(csxp.recomb$probe.before.transition, csxp.map.phys1$marker), ]$phys.pos
# 
# csxp.recomb.ind = ind.avg(csxp.recomb)
# 
# get.recomb.freq = function(x){
#     x = arrange(x, individual)
#     aggregate(x$num.events, list(factor(x$individual, levels = unique(x$individual))), sum)
# }
# 
# csxp.recomb.freq2 = get.recomb.freq(csxp.recomb)
# 
# csxp.qtl = prepareqtl(csxp.map, csxp.recomb.ind, csxp.recomb.freq2)
# 
# csxp.qtl[[1]][[2]]
# 
# lapply(csxp.qtl, function(x) x[[1]])
# 
# csxp.qtl[[1]][[2]]
# 
# 
# 
# 
# 
# 
# samp.coord1 = sample(2:nrow(csxp.map), 94)
# samp.temp = 2:nrow(csxp.map)
# samp.temp2 = samp.temp[-which(samp.temp %in% samp.coord1)]
# samp.coord2 = sample(samp.temp2, 94)
# 
# csxp.map.samp1 = csxp.map[c(1, samp.coord1), ]
# csxp.map.samp2 = csxp.map[c(1, samp.coord2), ]
# 
# csxp.samps = list(csxp.map.samp1, csxp.map.samp2)
# 
# csxp.recomb = lapply(csxp.samps, detect.recombination)
# csxp.hist = lapply(csxp.recomb, makehistlist.new)
# 
# p.adjust(round(perform.test(csxp.hist, c(1, 2)), digits = 5), "bonf")
# 
# hist(csxp.hist[[1]][[3]], breaks = 100)
# hist(csxp.hist[[2]][[3]], breaks = 100)
# 


#### APOGEE X PARAGON F5 ####

a.x.p.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic map AxP.csv", stringsAsFactors = F)

axp2 = prepare.allen.map(a.x.p.allen)



#extract genetic map from allen genotyping dataframe
axp.gen.map = a.x.p.allen[match(colnames(axp2), a.x.p.allen$marker), 1:3][-c(1, 2), ]


# axpqtl12 = qtlanalysis.all(axp2, axp.gen.map) #does specifying RIL5 effect the outcome? No.

axpqtl13 = qtlanalysis.all(axp2, axp.gen.map, "RIL5")


axpqtl13.step0.riself = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk")

grab.sig.phenotype.qtl(axpqtl13.step0.bcfst)
grab.sig.phenotype.qtl(axpqtl13.step0)

axpqtl13.step0.bcfst[[1]][[2]][[2]]
axpqtl13.step0.bcfst[[1]][[4]][[2]]
axpqtl13.step0.bcfst[[1]][[4]][[2]] + coord_cartesian(ylim = c(0, 10))

axpqtl13.step0.bcfst[[1]][[4]][[5]]$pheno



axpqtl13.step0[[1]][[4]][[2]]


axpqtl13.step0.bcfst[[1]][[11]][[2]]
axpqtl13.step0[[1]][[4]][[2]]



axpqtl13.step0 = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0)
axpqtl13.step0

#analyse recombination frequency
axpqtl13.step0.noliss = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0, skip.liss = T)

axpqtl13.step0[[1]][[3]][[2]]

grid.arrange(axpqtl13.step0[[1]][[3]][[2]], axpf2.qtl.analysis.step0[[1]][[3]][[2]])


lapply(axpqtl13.step0[[1]], function(x) x[[1]])
lapply(axpqtl13.step0.noliss[[1]], function(x) x[[1]])

lapply(axpqtl13.step0.noliss[[1]], function(x) x[[1]])

axpqtl13.step0.noliss[[1]][[9]][[2]]


axpqtl13.step0.noliss
axpqtl13.step0.noliss[[1]][[9]]

axpqtl13.step0.noliss[[1]]

lapply(axpqtl13[[1]], function(x) x[[3]])

axpqtl13[[1]][[3]][[2]]

axpqtl13[[1]][[2]][[2]]

axpqtl13[[1]][[7]][[2]]

do.call(grid.arrange, lapply(axpqtl13[[1]], function(x) x[[2]]))





# #### ALL POPULATION QTL ANALYSIS #### 
# # sxr.qtl1.step0
# # axpqtl13.step0
# # oxs.qtl1.step0
# # csxp.qtl1.0
# # axc.qtl1.step0
# # axcf2.qtl.analysis.step0
# # axpf2.qtl.analysis.step0
# # 
# # all.pops.qtl = list(sxr.qtl1.step0, axpqtl13.step0, oxs.qtl1.step0, csxp.qtl1.0, axc.qtl1.step0, axcf2.qtl.analysis.step0, axpf2.qtl.analysis.step0)
# # names(all.pops.qtl) = c("sxr.qtl1.step0", "axpqtl13.step0", "oxs.qtl1.step0", "csxp.qtl1.0", "axc.qtl1.step0", "axcf2.qtl.analysis.step0", "axpf2.qtl.analysis.step0")
# 
# 
# 
# 
# # lapply(all.pops.qtl, get.sig.qtl)
# 
# #significant QTL for step value of 0 in calc.genoprob()
# 
# 
# 
# 
# # sig.qtl1 = list(all.pops.qtl$axpqtl13.step0[[1]][[2]], all.pops.qtl$csxp.qtl1.0[[1]][[1]], all.pops.qtl$csxp.qtl1.0[[1]][[12]], all.pops.qtl$axcf2.qtl.analysis.step0[[1]][[8]], all.pops.qtl$axpf2.qtl.analysis.step0[[1]][[1]], all.pops.qtl$axpf2.qtl.analysis.step0[[1]][[3]])
# 
# 
# pop.list = c("A x P", "CS x P", "CS x P", "A x C F2", "A x P F2", "A x P F2")
# # all.qtl.dat1 = lapply(sig.qtl1, parse.qtl.data)
# 
# # all.qtl.dat2 = Map(function(x, y){
# #     x$pop = y
# #     x
# # }, all.qtl.dat1, pop.list)
# 
# 
# 
# 
# # all.qtl.dat3 = bind_rows(all.qtl.dat2)
# 
# # allq4 = all.qtl.dat3[, -which(colnames(all.qtl.dat3) %in% c("p-value Chi", "p-value F", "add.t", "dom.t", "pos", "pos.after", "pos.before", "lod.before", "lod.after", "dom.est", "dom.se"))]
# 
# # allq4
# 
# # allq4 = add.column.at.position(allq4, 0)
# # allq4 = add.column.at.position(allq4, 0)
# # 
# # allq4[, 1] = allq4$pop
# # allq4[, 2] = allq4$phenotype
# # 
# # allq4 = allq4[, -which(colnames(allq4) %in% c("phenotype", "pop"))]
# # colnames(allq4)[1:2] = c("pop", "phenotype")
# 
# #grab only the numeric columns (from all being character)
# numeric.cols1 = which(!sapply(lapply(allq4, as.numeric), function(x) all(is.na(x))))
# for(i in numeric.cols1){
#     allq4[, i] = round(as.numeric(allq4[, i]), digits = 3)
# }
# 
# 
# 
# 
# multi.qtl1.var = as.data.frame(unclass(sig.qtl1[[1]]$multi.qtl.models[[3]])$result.full)$`%var`[1]
# 
# 
# multi.qtl.drop = as.data.frame(unclass(sig.qtl1[[1]]$multi.qtl.models[[3]])$result.drop)
# 
# multi.qtl.drop[which(multi.qtl.drop$`Pvalue(F)` < 0.05), ]
# 
# 
# 
# # ----- PLOTS
# 
# all.pops.qtl$csxp.qtl1.0[[1]][[1]][[5]]
# 
# qtlplot1 = all.pops.qtl$axpqtl13.step0[[1]][[2]][[2]] + ggtitle("Apogee x Paragon - MRD 1B")
# qtlplot2 = all.pops.qtl$csxp.qtl1.0[[1]][[1]][[2]] + ggtitle("Chinese Spring x Paragon - Recombination Frequency")
# qtlplot3 = all.pops.qtl$csxp.qtl1.0[[1]][[12]][[2]] + ggtitle("Chinese Spring x Paragon - MRD 5B")
# qtlplot4 = all.pops.qtl$axcf2.qtl.analysis.step0[[1]][[8]][[2]] + ggtitle("Avalon x Cadenza - MRD 3A")
# qtlplot5 = all.pops.qtl$axpf2.qtl.analysis.step0[[1]][[1]][[2]] + ggtitle("Apogee x Paragon - Recombination Frequency")
# qtlplot6 = all.pops.qtl$axpf2.qtl.analysis.step0[[1]][[3]][[2]] + ggtitle("Apogee x Paragon - MRD 2A")
# 
# qtlplots = list(qtlplot1, qtlplot2, qtlplot3, qtlplot4, qtlplot5, qtlplot6)
# 
# do.call(grid.arrange, qtlplots)
# 
# 
# grid.arrange(all.pops.qtl$axpqtl13.step0[[1]][[2]][[2]], all.pops.qtl$csxp.qtl1.0[[1]][[1]][[2]], all.pops.qtl$csxp.qtl1.0[[1]][[12]][[2]], all.pops.qtl$axcf2.qtl.analysis.step0[[1]][[8]][[2]], all.pops.qtl$axpf2.qtl.analysis.step0[[1]][[1]][[2]], all.pops.qtl$axpf2.qtl.analysis.step0[[1]][[3]][[2]])
# 


#### QTL ANALYSIS W/O HETS ####

axp2[axp2 == "H"] = "-"
pxw492[pxw492 == "H"] = "-"
pxw942[pxw942 == "H"] = "-"
initial.data.comb[initial.data.comb == "H"] = "-"
lgeno3comb[lgeno3comb == "H"] = "-"


axpqtl13.step0.riself.het = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk")
pxw49.qtl1.step0.het = qtlanalysis.all(pxw492, pxw49.map, "RIL7", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk")
pxw94.qtl1.step0.het = qtlanalysis.all(pxw942, pxw94.map, "RIL7", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk")
axpf2.qtl.analysis.step0.het = qtlanalysis.all(initial.data.comb, initial.data.comb.map, "RIL2", F, ".", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk")
axcf2.qtl.analysis.step0.het = qtlanalysis.all(lgeno3comb, lgeno3comb.map, "RIL2", F, ".", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk")


grab.sig.phenotype.qtl(axpqtl13.step0.riself.het)
grab.sig.phenotype.qtl(pxw49.qtl1.step0.het)
grab.sig.phenotype.qtl(pxw94.qtl1.step0.het)
pxw94.qtl1.step0.het[[1]][[7]][[2]]

grab.sig.phenotype.qtl(axpf2.qtl.analysis.step0.het)
grab.sig.phenotype.qtl(axcf2.qtl.analysis.step0.het)


axpqtl13.step0.riself.het.noliss = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk", skip.liss = T)
pxw49.qtl1.step0.het.noliss = qtlanalysis.all(pxw492, pxw49.map, "RIL7", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk", skip.liss = T)
pxw94.qtl1.step0.het.noliss = qtlanalysis.all(pxw942, pxw94.map, "RIL7", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk", skip.liss = T)
axpf2.qtl.analysis.step0.het.noliss = qtlanalysis.all(initial.data.comb, initial.data.comb.map, "RIL2", F, ".", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk", skip.liss = T)
axcf2.qtl.analysis.step0.het.noliss = qtlanalysis.all(lgeno3comb, lgeno3comb.map, "RIL2", F, ".", geno.prob.step = 0, fgen = "riself", qtl.algorithm = "ehk", skip.liss = T)



#### MISC ####


# 
# #### testing resizing by scaling points
# plot.new()
# plot.window(xlim = c(-5, 5), ylim = c(-5, 5))
# 
# #generate circle vertices
# csamp <- function(n,rad=1,centre=c(0,0)){ 
#     x0 <- centre[1] ; y0 <- centre[2] 
#     u <- 2*pi*seq(0, 1, 0.01)
#     rad*cbind(x=cos(u)+x0, y=sin(u)+y0) 
# } 
# 
# q2 = as.data.frame(csamp(100))
# 
# x1 = c(0, 1, 2, 2, 0)
# y1 = c(0, 1, 2, 5, 0)
# lines(x = q2$x*1.5, y = q2$y*1.5)
# 
# 
# dev.off()
# 
# 
# # trace(Ps_Visualization, edit = T)
# 
# 
# render("rotation1scripts_v4/rmarkdown/recombination_distribution.Rmd")
# 
# 





#### qtl function ####

# s.x.r.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/sxr.allen.csv", stringsAsFactors = F)
# 
# sxr2 = prepare.allen.map(s.x.r.allen)
# 
# #extract genetic map from allen genotyping dataframe
# sxrmap1 = s.x.r.allen[match(colnames(sxr2), s.x.r.allen$marker), 1:3][-c(1, 2), ]
# 
# g123 = qtlanalysis.all(sxr2, sxrmap1)








##### FUNCTIONAL ANNOTATION ANALYSIS ####


functional.anno = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", sep = "\t", header = T, fill = T)

v(functional.anno[grep("SPO", functional.anno$Human.Readable.Description), ])




#### PHYLOGENETIC ANALYSIS ####

library(readxl)
alab1 = read_excel("rotation1scripts_v4/original_data/recombination.gene.analysis/alabdullah.genes1.xlsx")
alab1 = alab1[-1, ]
colnames(alab1) = alab1[1, ]
alab1 = alab1[-1, ]
alab2 = alab1[, which(colnames(alab1) %in% c("Gene ID (v1.0)", "Gene name", "Expected function in wheat"))]
alab3 = alab2[match(unique(alab2$`Gene name`), alab2$`Gene name`), ]

alab.subset = alab2[match(unique(alab2$`Expected function in wheat`), alab2$`Expected function in wheat`), ]


alab.tree = readLines("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/alab1/all.alignments/cds_nogaps/tree.lengths.txt")
alab.tree2 = multi.str.split(alab.tree, ": ", 2)
alab.tree2 = as.numeric(alab.tree2)

random.tree = readLines("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/cds_nogaps/tree.lengths.txt")
random.tree2 = multi.str.split(random.tree, ": ", 2)
random.tree2 = as.numeric(random.tree2)
report_t(t.test(alab.tree2, random.tree2))

boxplot(alab.tree2, random.tree2, names = c("Meiotic genes", "Random gene set"))

# ------ WITHOUT BARLEY

alab.tree.nobar = readLines("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/alab1/all.alignments/cds_nogaps/nobarley/tree.lengths.txt")
alab.tree.nobar2 = multi.str.split(alab.tree.nobar, ": ", 2)
alab.tree.nobar2 = as.numeric(alab.tree.nobar2)

random.tree.nobar = readLines("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/cds_nogaps/nobarley/tree.lengths.txt")
random.tree.nobar2 = multi.str.split(random.tree.nobar, ": ", 2)
random.tree.nobar2 = as.numeric(random.tree.nobar2)

report_t(t.test(alab.tree.nobar2, random.tree.nobar2))
report_mann(wilcox.test(alab.tree.nobar2, random.tree.nobar2))

boxplot(alab.tree.nobar2, random.tree.nobar2, names = c("Meiotic genes", "Random gene set"))

# ------ MORE PHYLOGENETIC ANALYSIS


library(ape)
library(adephylo)

#read in trees calculated from superalignments
alab.tree = read.tree("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/alab1/all.alignments/cds_nogaps/superalignment.fa.treefile")
random.tree = read.tree("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/cds_nogaps/superalignment.fa.treefile")

alab.tree.nobar = read.tree("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/alab1/all.alignments/cds_nogaps/nobarley/nobar.superalign.fa.treefile")
random.tree.nobar = read.tree("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/cds_nogaps/nobarley/nobar.superalign.fa.treefile")


alab.tree$tip.label = c("Chinese Spring", "Cadenza", "Robigus", "Claire", "Paragon", "Barley")
alab.tree.nobar$tip.label = c("Chinese Spring", "Cadenza", "Robigus", "Claire", "Paragon")
random.tree$tip.label = c("Chinese Spring", "Cadenza", "Paragon", "Robigus", "Claire", "Barley")
random.tree.nobar$tip.label = c("Chinese Spring", "Cadenza", "Paragon", "Robigus", "Claire")

alabdist = as.matrix(distTips(alab.tree))
alabdist[upper.tri(alabdist)] <- 0
alabdist = as.data.frame(alabdist)
alabdist[alabdist == "0"] = ""



randomdist = as.matrix(distTips(random.tree))[c(1, 2, 4, 5, 3, 6), c(1, 2, 4, 5, 3, 6)]
randomdist[upper.tri(randomdist)] <- 0
randomdist = as.data.frame(randomdist)
randomdist[randomdist == "0"] = ""



par(mfrow = c(1, 2))
plot(alab.tree) + title("Meiosis Gene Set")
plot(random.tree) + title("Random Gene Set")

plot(unroot(alab.tree.nobar), type = "unrooted", cex = 0.8) + title("Meiosis Gene Set")
add.scale.bar(length = 0.0005)
plot(unroot(random.tree.nobar), type = "unrooted", cex = 0.8) + title("Random Gene Set")
add.scale.bar(length = 0.0005)


#### EM QTL ALGORITHM ####
sxr.qtl1.step0.em = qtlanalysis.all(sxr2, sxr.gen.map, "RIL10", geno.prob.step = 0, qtl.algorithm = "em")

#analyse recombination frequency without taking LISS (longest increasing subsequence, qtl.algorithm = "em")
sxr.qtl1.step0.noliss.em = qtlanalysis.all(sxr2, sxr.gen.map, "RIL10", geno.prob.step = 0, skip.liss = T, qtl.algorithm = "em")

axpqtl13.step0.em = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0, qtl.algorithm = "em")
axpqtl13.step0.noliss.em = qtlanalysis.all(axp2, axp.gen.map, "RIL5", geno.prob.step = 0, skip.liss = T, qtl.algorithm = "em")

oxs.qtl1.step0.em = qtlanalysis.all(oxs2, oxs.map, "RIL7", geno.prob.step = 0, qtl.algorithm = "em")
#analyse recombination frequency
oxs.qtl1.step0.noliss.em = qtlanalysis.all(oxs2, oxs.map, "RIL7", geno.prob.step = 0, skip.liss = T, qtl.algorithm = "em")

csxp.qtl1.0.em = qtlanalysis.all(csxp.geno, csxp.map1, geno.prob.step = 0, qtl.algorithm = "em")
csxp.qtl1.0.noliss.em = qtlanalysis.all(csxp.geno, csxp.map1, geno.prob.step = 0, skip.liss = T, qtl.algorithm = "em")

axc.qtl1.step0.em = qtlanalysis.all(axc.geno, axc.map1, "RIL10", geno.prob.step = 0, fgen = "dh", qtl.algorithm = "em")
axc.qtl1.step0.noliss.em = qtlanalysis.all(axc.geno, axc.map1, "RIL10", geno.prob.step = 0, skip.liss = T, fgen = "dh", qtl.algorithm = "em")


axcf2.qtl.analysis.step0.em = qtlanalysis.all(lgeno3comb, lgeno3comb.map, "RIL2", F, ".", geno.prob.step = 0, fgen = 2, qtl.algorithm = "em")


axpf2.qtl.analysis.step0.em = qtlanalysis.all(initial.data.comb, initial.data.comb.map, "RIL2", F, ".", geno.prob.step = 0, fgen = 2, qtl.algorithm = "em")



grab.sig.phenotype.qtl(oxs.qtl1.step0.em)
oxs.qtl1.step0.em[[1]][[11]][[2]]


grab.sig.phenotype.qtl(axpf2.qtl.analysis.step0.em)


axc.qtl1.step0.em[[1]][[17]][[2]]
axc.qtl1.step0.noliss.em[[1]][[1]][[2]]

axpf2.qtl.analysis.step0.em[[1]][[7]][[2]]

#### GENE DISTRIBUTION / RECOMBINATION ANALYSIS #### 

gff1 = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", sep = "\t", stringsAsFactors = F, header = T)
gff1 = gff1[, -1]
gff1 = reset.colnames(gff1)
gff1$V4 = gff1$V4 / 1000000

gff1$V1 = multi.str.split(gff1$V1, "chr", 2)
colnames(gff1)[1] = "chr"
gff1 = gff1[-which(gff1$chr == "Un"), ]
gff2 = split(gff1, gff1$chr)

centromeres = read.table("rotation1scripts_v4/original_data/IWGSC/centromere.positions.txt", sep = " ", stringsAsFactors = F)
#taken from S11 of supplementary data from 2018 IWGSC paper in Science. Where more than one region has been defined,
#the region that overlaps the highest centromeric TE family (either Cereba or Quinta; figure S7) has been taken.
centromere.positions = sapply(centromeres$V2, function(x){
    g = strsplit(x, "-")
    g = as.numeric(g[[1]])
    mean(g)
})

centromere.positions.bp = centromere.positions * 1000000

centromeres$pos = centromere.positions
centromeres$pos.bp = centromere.positions.bp
colnames(centromeres)[1] = "chr"

gene.dist.hists = Map(function(x, y){
    
    g = make.ggplot.histogram(x$V4, num.bins = 1000, breaks = seq(1, max(x$V4), (max(x$V4) / 1000)), xlabel = "Physical Position", ylabel = "Frequency", plot.title = "Histogram") +
        coord_cartesian(ylim = c(0, 50))
    g + geom_vline(xintercept = y) + ggtitle(x$chr[[1]]) + theme_bw()
    
}, gff2[1:21], centromeres$pos)

do.call(grid.arrange, gene.dist.hists)







#chinese spring recombination plots


cslm1 = csxp.liss.map[[4]]
cslm2 = split(cslm1, cslm1$chr)
cslm2 = lapply(cslm2, function(x){
    x$phys.mb = x$phys.bp / 1000000
    x
})


cslm3 = lapply(cslm2, function(x){
    x$cmdiff = c(diff(x$cM), 0)
    x$physdiff = c(diff(x$phys.mb), 0)
    x$cmmb = x$cmdiff / x$physdiff
    x
})

cslm4 = unsplit(cslm3, cslm1$chr)
cslm4 = cslm4[-which(cslm4$chr == "1B" | cslm4$chr == "5A" | cslm4$chr == "7B"), ]
centromeres.cs = centromeres[which(centromeres$chr %in% unique(cslm4$chr)), ]

chinese.spring.recombination = ggplot(cslm4, aes(x = phys.mb, y = cmmb)) + geom_line() +
    facet_wrap(. ~ chr, scales = "free_y", ncol = 1, strip.position = 'right') + 
    geom_vline(aes(xintercept = pos), data = centromeres.cs) + geom_point(shape = 16, size = 1) +
    ylab("Recombination (cM / Mb)") + xlab("Physical position (Mb)") + theme_bw(base_size = 20) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2))



plots111 = lapply(cslm3, function(x){
    ggplot(x, aes(x = phys.mb, y = cmmb)) + geom_line()
})
do.call(grid.arrange, plots111)

#### COMPARISON OF ARRAY TO GENOME ####

probes35k = readDNAStringSet("rotation1scripts_v4/original_data/fasta/35kprobes.fa")
probes820k = readDNAStringSet("rotation1scripts_v4/original_data/fasta/820kprobes.fa")

luzie.probes = read.csv("rotation1scripts_v4/original_data/IWGSC/Luzie.probes.JIC_RefSeq1_axiom_order_Nov18.csv", stringsAsFactors = F)

unique(luzie.probes$WGA.Chr)


luzie.probes = luzie.probes[-which(luzie.probes$WGA.Chr == "unknown"), ]
luzie.probes$mb = luzie.probes$WGA.bp / 1000000
colnames(luzie.probes)[2] = "chr"
luzie.probes2 = split(luzie.probes, luzie.probes$WGA.Chr)


probe.dist.hists = Map(function(x, y){
    g = make.ggplot.histogram(x$mb, num.bins = 1000, breaks = seq(1, max(na.omit(x$mb)), (max(na.omit(x$mb)) / 1000)), xlabel = "Physical Position", ylabel = "Frequency", plot.title = "Histogram") + 
        coord_cartesian(ylim = c(0, 50))
    g + geom_vline(xintercept = y) + ggtitle(x$WGA.Chr[[1]]) + theme_bw()
}, luzie.probes2, centromeres$pos)

gene.dist.hists[[1]]
probe.dist.hists[[1]]


all.genes.plot = ggplot(gff1, aes(x = V4)) + geom_histogram(binwidth = 1) + facet_grid(rows = vars(chr)) + geom_vline(aes(xintercept = pos), data = centromeres) + coord_cartesian(ylim = c(0, 60))
all.probes.plot = ggplot(luzie.probes, aes(x = mb)) + geom_histogram(binwidth = 1) + facet_grid(rows = vars(chr)) + geom_vline(aes(xintercept = pos), data = centromeres) + coord_cartesian(ylim = c(0, 60))

grid.arrange(all.genes.plot, all.probes.plot, ncol = 2)


#### CORRELATION CHROMOSOME LENGTH AND MRD ####


#axc analysis
chromosomecounts2 = chromosomecounts
chromosomecounts2$V2 = multi.str.split(chromosomecounts2$V2, "chr", 2)

chromomrd2 = cbind(chromosomecounts2[match(means.per.chromo2$Group.2, chromosomecounts2$V2), ], means.per.chromo2)

chromomrd2$V1 = chromomrd2$V1 / 1000000

cor.test(chromomrd2$V1, chromomrd2$`4`)

#axp analysis
chromosomecounts2 = chromosomecounts
chromosomecounts2$V2 = multi.str.split(chromosomecounts2$V2, "chr", 2)

chromomrd2 = cbind(chromosomecounts2[match(axpmeans.per.chromo2$Group.2, chromosomecounts2$V2), ], axpmeans.per.chromo2)

chromomrd2$V1 = chromomrd2$V1 / 1000000

cor.test(chromomrd2$V1, chromomrd2$`4`)


