#simulated genotype data

#install.packages("ASMap")
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/")
library(ASMap)
library(qtl)
library(dplyr)
library(tibble)

axiom.data.csv = read.csv("original_data/simulated.genotype.data/very.simple.h.csv", header = T, stringsAsFactors = F)

axiom.data = read.cross(format = "csv", dir = "./", file = "original_data/simulated.genotype.data/very.simple.csv", estimate.map = F)
axiom.data = read.cross(format = "csv", dir = "./", file = "original_data/simulated.genotype.data/very.simple.even.more.markers.csv", estimate.map = F)
axiom.data = read.cross(format = "csv", dir = "./", file = "original_data/simulated.genotype.data/very.simple.h.csv", estimate.map = F)

axiom.data = read.cross(format = "csv", dir = "./", file = "original_data/simulated.genotype.data/very.simple.v2.csv", estimate.map = F)
axiom.data = read.cross(format = "csv", dir = "./", file = "original_data/simulated.genotype.data/very.simple.v2.dist.csv", estimate.map = F)

#     ____________________________________________________________________________
#     ESTIMATE MAPS (PRESERVE ORDER)                                                                                    ####

#estimate the maps
b = quickEst(axiom.data, map.function = "kosambi")
b = est.map(axiom.data)
pull.map(b, as.table = T)

b = est.rf(axiom.data)
pull.rf(b)


test.genolist = makegenotypelist(axiom.data.csv)

test.histlist = makehistlist(1, test.genolist, return_recombinationcounts = T)

events = lapply(unique(test.histlist$individual), function(x){
    fil.events = filter(test.histlist, individual == x)
    num.events = nrow(fil.events)
    return(num.events)
})

num.missing = length(unique(test.genolist$individual)) - nrow(test.histlist)

events = c(events, rep(0, num.missing))



View(initial.data[[1]])
View(axiom.data.csv)


#     ____________________________________________________________________________
#     testing a subset of the axiom data                                                                            ####


load.data = function(x) read.csv(paste("processed_data/genotypes/alexq", x, "_man_cur_sk_map_genogeno.csv", sep = ""),
                                                                 header = T, stringsAsFactors = F)

initial.data = lapply(1:4, load.data)

test.data = initial.data[[1]][1:5,]
test.data[1, 1] = ""
test.data.noparents = test.data[-c(2,3),]

write.csv(test.data.noparents, "processed_data/genotypes/test.reduced.csv", row.names = F)

geno.test = makegenotypelist(test.data)

geno.hist.test = makehistlist(linkagegroup = 1, genolist2 = geno.test, T)

events = lapply(unique(geno.hist.test$individual), function(x){
    fil.events = filter(geno.hist.test, individual == x)
    num.events = nrow(fil.events)
    return(num.events)
})

#removed NA value in excel before this next line
axiom.test = read.cross(format = "csv", dir = "./", file = "processed_data/genotypes/test.reduced.csv", estimate.map = F)
axiom.test = read.cross(format = "csv", dir = "./", file = "processed_data/genotypes/alexq1_man_cur_sk_map_genogeno.csv", estimate.map = F)

g = quickEst(axiom.test, map.function = "kosambi")
g = est.map(axiom.test, map.function = "kosambi")
g.table = pull.map(g, as.table = T)

calculate.map.length = function(genotypedata){
    
}

