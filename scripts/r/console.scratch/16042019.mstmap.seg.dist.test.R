library(ASMap)
library(parallel)

load("rotation1scripts_v4/saved.objects/cxa.6b.cm")
load("rotation1scripts_v4/saved.objects/axc.chr1a.cm")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/pedigree.sim.functions.R")
#### FUNCTIONS ####

make.cluster.comparison = function(selec.str, skip.run, selec.pos1, selec.pos2){
    
    if(missing(skip.run)) skip.run = F
    if(missing(selec.pos1)) selec.pos1 = 30
    if(missing(selec.pos2)) selec.pos2 = 200
    
    if(skip.run == F){
    
    run.1000.sims(p("axcf2.6b.selec", selec.str, "pop300pos", selec.pos1), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = selec.pos1, pop.size = 300, cm.pos = cxa.6b.cm)    
    
    
    run.1000.sims(p("axcf2selec", selec.str, "pop300pos", selec.pos2), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = selec.pos2, pop.size = 300, cm.pos = axc.chr1a.cm)
    
    }
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec", selec.str, "pop300pos", selec.pos2))
    axc1a = all.sim.files
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec", selec.str, "pop300pos", selec.pos1))
    axc6b = all.sim.files
    
    
    axc6b = lapply(axc6b, function(y){
        colnames(y) = paste0("6b", colnames(y))
        y
    })
    
    
    sim1.2comb = Map(function(sim6b, sim1a){
        g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)
        
    }, axc1a, axc6b)
    
    sim1.2comb
    
}

make.cluster.comparison.one.selec = function(selec.str, skip.run){
    
    if(missing(skip.run)) skip.run = F
    
    if(skip.run == F){
        
        run.1000.sims(p("axcf2.6b.selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm, F)    
        
        
        run.1000.sims(p("axcf2noselecpop300"), 2, num.sims = 1, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)
    }
    
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2noselecpop300"))
    axc1a = all.sim.files
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec", selec.str, "pop300"))
    axc6b = all.sim.files
    
    
    axc6b = lapply(axc6b, function(y){
        colnames(y) = paste0("6b", colnames(y))
        y
    })
    
    
    sim1.2comb = Map(function(sim6b, sim1a){
        g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)
        
    }, axc1a, axc6b)
    
    sim1.2comb
    
}

make.cluster.comparison.no.selec = function(skip.run){
    
    if(missing(skip.run)) skip.run = F
    
    if(skip.run == F){
        
        run.1000.sims(p("axcf2.6b.noselecpop300"), 2, num.sims = 1, w.selec = F, pop.size = 300, cm.pos = cxa.6b.cm, F)    
        
        
        run.1000.sims(p("axcf2noselecpop300"), 2, num.sims = 1, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm, F)
    }
    
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2noselecpop300"))
    axc1a = all.sim.files
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.noselecpop300"))
    axc6b = all.sim.files
    
    
    axc6b = lapply(axc6b, function(y){
        colnames(y) = paste0("6b", colnames(y))
        y
    })
    
    
    sim1.2comb = Map(function(sim6b, sim1a){
        g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)
        
    }, axc1a, axc6b)
    
    sim1.2comb
    
}

gen.asmap.for.sim = function(sim.dat1, p.value1, RIL){
    if(missing(p.value1)) p.value1 = 4
    if(missing(RIL)) RIL = 2
    sim1.rqtl2 = sim.dat1
    sim1.rqtl2[sim1.rqtl2 == "H"] = "X"
    sim1.rqtl2 = as.data.frame(t(sim1.rqtl2))
    sim1.rqtl2 = convert.to.character.data.frame(sim1.rqtl2)
    sim1.rqtl3 = mstmap.data.frame(sim1.rqtl2, p("RIL", RIL), "kosambi", p.value = p.value1)
    sim1.rqtl3
}

do.order.comp = function(x){
    g = gen.asmap.for.sim(convert.to.character.data.frame(x))
    g2 = pull.map(g, as.table = T)
    
    g3 = rownames(g2)
    
    g4 = as.numeric(multi.str.split(g3, "marker", 2))
    
    224:1
    
    mean((match(g4, 224:1) - 1:224)^2)	
}

f1 = function(x, p){
    #get number of linkage groups from genotyping data
    g = gen.asmap.for.sim(convert.to.character.data.frame(x), p.value1 = p)
    g1 = pull.map(g, as.table = T)
    q = length(unique(g1$chr))
    attr(q, "map") = g1
    q
} 


f2 = function(x, p){
    #get genetic map from genotyping data
    gen.asmap.for.sim(convert.to.character.data.frame(x), p.value1 = p)
} 

order.within.bins = function(cluster.output, decreasing, mapmap){
    if(missing(mapmap)) mapmap = F
    if(mapmap == F){
        map1 = attributes(cluster.output)$map
    } else {
        map1 = cluster.output
    }
    
    
    #order markers within bins (asmap flips the order for some reason)
    map2 = lapply(unique(map1$chr), function(x){
        # browser()
        g = map1[which(map1$chr == x), ]
        
        g1 = lapply(unique(g$pos), function(y){
            # browser()
            g1 = g[which(g$pos == y),]
            g2 = as.numeric(multi.str.split(rownames(g1), "marker", 2))
            g1[sort(g2, decreasing = decreasing, index.return = T)$ix, ]
        })
        
        do.call(rbind, g1)
    })
    
    map3 = do.call(rbind, map2)
    map3
}


examine.order1 = function(cluster.output){    
    g = rownames(order.within.bins(cluster.output, T))
    g1 = rownames(order.within.bins(cluster.output, F))

    g1.1 = g1[grep("^marker", g1)]
    g2 = g[grep("^6bmarker", g)]

    print(g1.1)
    print(g2)


    g2.nums = as.numeric(multi.str.split(g2, "marker", 2))
    res1 = all(g2.nums == 222:1)
    bad.coord1 = which(!g2.nums == 222:1)


    g1.nums = as.numeric(multi.str.split(g1.1, "marker", 2))
    res2 = all(g1.nums == 1:224)
    bad.coord2 = which(!g1.nums == 1:224)
    list(res1, res2, bad.coord1, bad.coord2)
}



do.t.test = function(map1, map2){
    t1 = sapply(map1, function(x){
        g = attributes(x)$map
        max(g$pos)
    })
    
    t2 = sapply(map2, function(x){
        g = attributes(x)$map
        max(g$pos)
    })
    
    t.test(t1, t2)
    
}


#### MAIN CODE #### 



# all.sims1 = mclapply(all.sims1, function(x){
# 	lapply(x, function(y){
# 			gen.asmap.for.sim(y)
# 		})
# 	}, mc.cores = 5)



# g = gen.asmap.for.sim(convert.to.character.data.frame(axcf2selec2pop300[[1]]))
# pull.map(g, as.table = T)
# 
# g = gen.asmap.for.sim(convert.to.character.data.frame(axcf2selec20pop300[[1]]))
# pull.map(g, as.table = T)
# 


# load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop300")
# 
# axcf2.noselec.pop300 = all.sim.files
# 






# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec2pop300")
# 
# axcf2.6b.selec2pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec4pop300")
# 
# axcf2.6b.selec4pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec6pop300")
# 
# axcf2.6b.selec6pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec8pop300")
# 
# axcf2.6b.selec8pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec10pop300")
# 
# axcf2.6b.selec10pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec12pop300")
# 
# axcf2.6b.selec12pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec14pop300")
# 
# axcf2.6b.selec14pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec16pop300")
# 
# axcf2.6b.selec16pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec18pop300")
# 
# axcf2.6b.selec18pop300 = all.sim.files
# 
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/6b.sims/axcf2.6b.selec20pop300")
# 
# axcf2.6b.selec20pop300 = all.sim.files

# all.sim6b = list(axcf2.6b.selec2pop300, axcf2.6b.selec4pop300, axcf2.6b.selec6pop300, axcf2.6b.selec8pop300, axcf2.6b.selec10pop300, axcf2.6b.selec12pop300, axcf2.6b.selec14pop300, axcf2.6b.selec16pop300, axcf2.6b.selec18pop300, axcf2.6b.selec20pop300)

# all.sim6b = lapply(all.sim6b, function(x){
# 	lapply(x, function(y){
# 		colnames(y) = paste0("6b", colnames(y))
# 		y
# 		})
# 	})


# all.sim.comb = Map(function(z, q){
# 	Map(function(sim6b, sim1a){
# 	g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)
# 
# 	}, z, q)
# 
# 	}, all.sim6b, all.sims1)




# g = gen.asmap.for.sim(convert.to.character.data.frame(all.sim.comb[[1]][[1]]), p.value1 = 1e-20)
# convert.recomb.to.seg.f2(all.sim.comb[[1]][[1]], T) < 0.000000001


# g1 = gen.asmap.for.sim(convert.to.character.data.frame(all.sim.comb[[10]][[1]]), p.value1 = 1e-20)

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec1.2pop300")
axcf2selec1.2pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec1.2pop300")
axcf2.6b.selec1.2pop300 = all.sim.files

#rename columns for later identification
axcf2.6b.selec1.2pop300 = lapply(axcf2.6b.selec1.2pop300, function(y){
		colnames(y) = paste0("6b", colnames(y))
		y
		})

#combine both genotyping dataframes
sim1.2comb = Map(function(sim6b, sim1a){
	g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)

	}, axcf2selec1.2pop300, axcf2.6b.selec1.2pop300)

#generate genetic map from combined genotyping data
g1 = gen.asmap.for.sim(convert.to.character.data.frame(sim1.2comb[[1]]), p.value1 = 1e-40)


####    old code?    ####

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec1.01pop300")
axcf2selec1.01pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec1.01pop300")
axcf2.6b.selec1.01pop300 = all.sim.files


axcf2.6b.selec1.01pop300 = lapply(axcf2.6b.selec1.01pop300, function(y){
		colnames(y) = paste0("6b", colnames(y))
		y
		})


sim1.2comb = Map(function(sim6b, sim1a){
	g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)

	}, axcf2selec1.01pop300, axcf2.6b.selec1.01pop300)

g1 = gen.asmap.for.sim(convert.to.character.data.frame(sim1.2comb[[1]]), p.value1 = 1e-45)

#### old code?    ####

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec1.2pop300")
axcf2selec1.2pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec1.2pop300")
axcf2.6b.selec1.2pop300 = all.sim.files


axcf2.6b.selec1.2pop300 = lapply(axcf2.6b.selec1.2pop300, function(y){
		colnames(y) = paste0("6b", colnames(y))
		y
		})


sim1.2comb = Map(function(sim6b, sim1a){
	g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)

	}, axcf2selec1.2pop300, axcf2.6b.selec1.2pop300)

g1 = gen.asmap.for.sim(convert.to.character.data.frame(sim1.2comb[[1]]), p.value1 = 1e-45)

g1 = gen.asmap.for.sim(convert.to.character.data.frame(cluster1.15[[1]]), p.value1 = 1e-45)

#### cluster comparisons ####

cluster1.15 = make.cluster.comparison(1.15)

cluster1.12 = make.cluster.comparison(3, F, 30, 100)



geno.count(cluster1.12[[1]][[1]]["marker100"]) 

f1(cluster1.12)

cluster1.1 = make.cluster.comparison(1.1)

f1(cluster1.1, 1.5e-48)

f1(cluster1.12[[1]], 1.5e-48)

selec.pressures = seq(1.1, 1.2, 0.005)

selec.pressures.debug = rep(1, 10)

selec.pressures.lr = seq(1.1, 2, 0.1)

#move one selection position to a region of low recombination
clusters = lapply(selec.pressures.lr, function(x){
	make.cluster.comparison(x, F, 30, 100)
	})

#move both selection pressures to positions of low recombination
clusters = lapply(selec.pressures.lr, function(x){
    make.cluster.comparison(x, F, 110, 100)
    })

#scramble order of genotypes
clusters2 = lapply(clusters, function(x){
	x[[1]][, sample(ncol(x[[1]]), ncol(x[[1]]))]
	})





clus3 = f1(clusters2[[3]], 1e-45)
clus3.map = order.within.bins(clus3, F)
multi.str.split(rownames(clus3.map), "marker", 2)

geno.counts1 = lapply(clusters[[3]][[1]], geno.count)
geno.counts1[c("6bmarker30", "marker200")]

lapply(clusters[1:2], function(z){
    geno.counts1 = lapply(z[[1]], geno.count)
    geno.counts1[c("6bmarker30", "marker100")]    
    })

lapply(clusters[1:2], function(z){
    geno.counts1 = lapply(z[[1]], geno.count)
    geno.counts1[c("6bmarker110", "marker100")]    
    })



seq2 = 10^-seq(40, 50, 1)
#adjust the value of X in clusters2[[X]] to change the selection pressure
lgs1.1 = lapply(seq2, function(x){
	f1(clusters2[[1]], x)
	})

lgs1.105 = lapply(seq2, function(x){
    f1(clusters2[[2]], x)
    })


examine.order1(lgs1.1[[1]])
examine.order1(lgs1.105[[1]])

lgs1.105 = lapply(seq2, function(x){
	f2(clusters2[[2]], x)
	})

newseq = min(which(lgs1.105 >= 2))
newseq2 = seq(1, 10, 1)*seq2[newseq]

lgs2 = lapply(newseq2, function(x){
	f1(clusters[[2]], x)
	})


lgs1.11 = lapply(seq2, function(x){
	f1(clusters2[[3]], x)
	})


lapply(seq2, function(x){
    f1(clusters2[[2]], x)
    })



pull.rf(attributes(lgs1.11[[3]]))


lgs1.11 = lapply(seq2, function(x){
	f2(clusters2[[3]], x)
	})


pdf("rotation1scripts_v4/temp/lg1.11.pdf", 10, 10)
heatMap(lgs1.11[[3]], what = "rf")
dev.off()

pdf("rotation1scripts_v4/temp/lg1.105.pdf", 10, 10)
heatMap(lgs1.105[[6]], what = "rf")
dev.off()



newseq = min(which(lgs1.105 >= 2))
newseq2 = seq(1, 10, 1)*seq2[newseq]

lgs2 = lapply(newseq2, function(x){
	f1(clusters[[2]], x)
	})






clusters[[2]]

cluster1.2 = make.cluster.comparison(1.2)

clus1.11 = lapply(1:100, function(x){
    clus1 = make.cluster.comparison(1.11)
    clus1
})


clus1.11.maps = mclapply(clus1.11, function(x){
    clus1.scram2 = x[[1]][, sample(ncol(x[[1]]), ncol(x[[1]]))]
    q2 = 1
    
    seq2 = 10^-seq(45, 40, -1)
    counter = 1
    
    while(q2 != 2 & counter < 7){
        q2 = f1(list(clus1.scram2), seq2[counter])    
        counter = counter + 1
    }
    q2
}, mc.cores = 50)



geno.counts3 = lapply(clus1.11[[2]][[1]], geno.count)
geno.counts3[c("6bmarker30", "marker200")]

#### ONE SELECTION PRESSURE #### 

clus1s = lapply(1:100, function(x){
    clus1 = make.cluster.comparison.one.selec(1)
    clus1
})


clus1s.maps = mclapply(clus1s, function(x){
    clus1.scram2 = x[[1]][, sample(ncol(x[[1]]), ncol(x[[1]]))]
    q2 = 1
    
    seq2 = 10^-seq(43, 38, -1)
    counter = 1
    
    while(q2 != 2 & counter < 7){
        q2 = f1(list(clus1.scram2), seq2[counter])    
        counter = counter + 1
    }
    q2
}, mc.cores = 50)

## lets try removing the distorted markers 

clus1s[[1]][[1]]

lapply(clus1s[[1]][[1]], geno.count)

seg.coords = which(sapply(clus1s[[1]][[1]], function(x) geno.count.chi.f2(x)['p-value']) < 0.0001)

temp.map = clus1s[[1]][[1]][, -seg.coords]

f1(temp.map, 1e-32)

##

map.lengths1 = as.data.frame(t(sapply(clus1s.maps, function(x){
    x = attributes(x)$map
    
    g = rownames(x[which(x$chr == "L1"), ])
    g2 = grep("6b", g)
    if(length(g2) > 0){
        sixb.len = max(x[which(x$chr == "L1"), ]$pos)
        onea.len = max(x[which(x$chr == "L2"), ]$pos)
    } else {
        sixb.len = max(x[which(x$chr == "L2"), ]$pos)
        onea.len = max(x[which(x$chr == "L1"), ]$pos)
    }
    
    
    c(onea.len, sixb.len)
})))


clus1 = make.cluster.comparison.one.selec(1, T)

geno.counts3 = lapply(clus1[[1]], geno.count)
geno.counts3[c("6bmarker30", "marker200")]

clus1.scram = clus1[[1]][, sample(ncol(clus1[[1]]), ncol(clus1[[1]]))]
q = f1(list(clus1.scram), 1e-43)

q = lgs[[3]]


clus2s = lapply(1:100, function(x){
    clus2 = make.cluster.comparison.no.selec()    
    clus2
})

geno.counts4 = lapply(clus2[[1]], geno.count)
geno.counts4[c("6bmarker30", "marker200")]

clus2s.maps = mclapply(clus2s, function(x){
    clus2.scram2 = x[[1]][, sample(ncol(x[[1]]), ncol(x[[1]]))]    
    q2 = f1(list(clus2.scram2), 1e-31)
}, mc.cores = 50)


map.lengths = as.data.frame(t(sapply(clus2s.maps, function(x){
    x = attributes(x)$map
    
    g = rownames(x[which(x$chr == "L1"), ])
    g2 = grep("6b", g)
    if(length(g2) > 0){
        sixb.len = max(x[which(x$chr == "L1"), ]$pos)
        onea.len = max(x[which(x$chr == "L2"), ]$pos)
    } else {
        sixb.len = max(x[which(x$chr == "L2"), ]$pos)
        onea.len = max(x[which(x$chr == "L1"), ]$pos)
    }
    
    
    c(onea.len, sixb.len)
})))


map.lengths = as.data.frame(t(map.lengths))

q2 = lgs[[3]]


#### one selec. w/ smaller str ####

clus3s = lapply(1:100, function(x){
    clus1 = make.cluster.comparison.one.selec(4)
    clus1
})


clus1s.maps = mclapply(clus1s, function(x){
    clus1.scram2 = x[[1]][, sample(ncol(x[[1]]), ncol(x[[1]]))]
    q2 = 1
    
    seq2 = 10^-seq(43, 38, -1)
    counter = 1
    
    while(q2 != 2 & counter < 7){
        q2 = f1(list(clus1.scram2), seq2[counter])    
        counter = counter + 1
    }
    q2
}, mc.cores = 50)



#### EFFECT ON MAP LENGTH #### 

#CHR 6B
g = seq(1, 10, 2)
lapply(g, function(x){
    run.1000.sims(p("axcf2.6b.selec", x, "pop300"), 2, num.sims = 100, w.selec = T, selec.str = x, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)        
})

run.1000.sims(p("axcf2.6b.noselecpop300"), 2, num.sims = 100, w.selec = F, pop.size = 300, cm.pos = cxa.6b.cm)     

make.maps1 = function(genotypes){
    mclapply(genotypes, function(x){
        clus1.scram2 = x[, sample(ncol(x), ncol(x))]
        q2 = 2
        
        seq2 = 10^-seq(43, 34, -1)
        counter = 1
        
        while(q2 != 1 & counter <= length(seq2)){
            q2 = f1(list(clus1.scram2), seq2[counter])    
            counter = counter + 1
        }
        q2
    }, mc.cores = 50)
}

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec1pop300")

axcf2.6b.selec1pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec3pop300")

axcf2.6b.selec3pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec5pop300")

axcf2.6b.selec5pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec7pop300")

axcf2.6b.selec7pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.selec9pop300")

axcf2.6b.selec9pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.noselecpop300")

axcf2.6b.noselecpop300 = all.sim.files

sixb.geno = list(axcf2.6b.selec1pop300, axcf2.6b.selec3pop300, axcf2.6b.selec5pop300, axcf2.6b.selec7pop300, axcf2.6b.selec9pop300, axcf2.6b.noselecpop300)

lapply(sixb.geno, function(x) geno.count(x[[1]][, 30]))

geno.count(axcf2.1a.selec1pop300[[1]][, 30])
geno.count(axcf2.1a.selec9pop300[[1]][, 30])


s1.maps = make.maps1(axcf2.6b.selec1pop300)
s3.maps = make.maps1(axcf2.6b.selec3pop300)
s5.maps = make.maps1(axcf2.6b.selec5pop300)
s7.maps = make.maps1(axcf2.6b.selec7pop300)
s9.maps = make.maps1(axcf2.6b.selec9pop300)
sno.maps = make.maps1(axcf2.6b.noselecpop300)

smaps = list(s1.maps, s3.maps, s5.maps, s7.maps, s9.maps, sno.maps)
s(smaps, "rotation1scripts_v4/saved.objects/smaps", "16042019.mstmap.seg.dist.test.R")


lapply(smaps, function(y){
    mean(sapply(y, function(x) max(attributes(x)$map$pos)))    
})

#1A

g = seq(1, 10, 2)
lapply(g, function(x){
    run.1000.sims(p("axcf2.1a.selec", x, "pop300"), 2, num.sims = 100, w.selec = T, selec.str = x, selec.pos = 200, pop.size = 300, cm.pos = axc.chr1a.cm)        
})

run.1000.sims(p("axcf2.1a.noselecpop300"), 2, num.sims = 100, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)        

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a.selec1pop300")

axcf2.1a.selec1pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a.selec3pop300")

axcf2.1a.selec3pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a.selec5pop300")

axcf2.1a.selec5pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a.selec7pop300")

axcf2.1a.selec7pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a.selec9pop300")

axcf2.1a.selec9pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a.noselecpop300")

axcf2.1a.noselecpop300 = all.sim.files

onea.geno = list(axcf2.1a.selec1pop300, axcf2.1a.selec3pop300, axcf2.1a.selec5pop300, axcf2.1a.selec7pop300, axcf2.1a.selec9pop300, axcf2.1a.noselecpop300)

a1.maps = make.maps1(axcf2.1a.selec1pop300)
a3.maps = make.maps1(axcf2.1a.selec3pop300)
a5.maps = make.maps1(axcf2.1a.selec5pop300)
a7.maps = make.maps1(axcf2.1a.selec7pop300)
a9.maps = make.maps1(axcf2.1a.selec9pop300)
ano.maps = make.maps1(axcf2.1a.noselecpop300)

amaps = list(a1.maps, a3.maps, a5.maps, a7.maps, a9.maps, ano.maps)
s(amaps, "rotation1scripts_v4/saved.objects/amaps", "16042019.mstmap.seg.dist.test.R")

lapply(amaps, function(y){
    mean(sapply(y, function(x) max(attributes(x)$map$pos)))    
})


lapply(onea.geno, function(x) geno.count(x[[1]][, 200]))






run.1000.sims(p("axcf2selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)


# 1A DIFF SELEC POS 

g = seq(1, 10, 2)
lapply(g, function(x){
    run.1000.sims(p("axcf2.1a2.selec", x, "pop300"), 2, num.sims = 100, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 300, cm.pos = axc.chr1a.cm)        
})

run.1000.sims(p("axcf2.1a2.noselecpop300"), 2, num.sims = 100, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)        

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a2.selec1pop300")

axcf2.1a2.selec1pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a2.selec3pop300")

axcf2.1a2.selec3pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a2.selec5pop300")

axcf2.1a2.selec5pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a2.selec7pop300")

axcf2.1a2.selec7pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a2.selec9pop300")

axcf2.1a2.selec9pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a2.noselecpop300")

axcf2.1a2.noselecpop300 = all.sim.files

oneaa.geno = list(axcf2.1a2.selec1pop300, axcf2.1a2.selec3pop300, axcf2.1a2.selec5pop300, axcf2.1a2.selec7pop300, axcf2.1a2.selec9pop300, axcf2.1a2.noselecpop300)



aa1.maps = make.maps1(axcf2.1a2.selec1pop300)
aa3.maps = make.maps1(axcf2.1a2.selec3pop300)
aa5.maps = make.maps1(axcf2.1a2.selec5pop300)
aa7.maps = make.maps1(axcf2.1a2.selec7pop300)
aa9.maps = make.maps1(axcf2.1a2.selec9pop300)
aano.maps = make.maps1(axcf2.1a2.noselecpop300)

aamaps = list(aa1.maps, aa3.maps, aa5.maps, aa7.maps, aa9.maps, aano.maps)

get.mean.map.length = function(z){
        mean(sapply(z, function(x) max(attributes(x)$map$pos)))    
}

get.sd.map.length = function(z){
    sd(sapply(z, function(x) max(attributes(x)$map$pos)))    
}

lapply(aamaps, function(y){
    mean(sapply(y, function(x) max(attributes(x)$map$pos)))    
})

lapply(oneaa.geno, function(x) geno.count(x[[3]][, 100]))

par(mfrow = c(2, 1))
plot(unique(axc.chr1a.cm))
plot(unique(cxa.6b.cm))

# selec 1 vs no selec # 
run.1000.sims(p("axcf2.1a3.selec1pop300"), 2, num.sims = 50, w.selec = T, selec.str = 1, selec.pos = 100, pop.size = 300, cm.pos = axc.chr1a.cm) 

run.1000.sims(p("axcf2.1a3.noselecpop300"), 2, num.sims = 50, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)    

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a3.selec1pop300")

axcf2.1a3.selec1pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.1a3.noselecpop300")

axcf2.1a3.noselecpop300 = all.sim.files

qmaps1 = make.maps1(axcf2.1a3.selec1pop300)
qmaps2 = make.maps1(axcf2.1a3.noselecpop300)



do.t.test(qmaps1, qmaps2)

get.mean.map.length(qmaps1)
get.mean.map.length(qmaps2)



q = convert.to.character.data.frame(axcf2.1a3.selec1pop300[[2]])
q[q == "H"] = "X"
pull.map(gen.asmap.for.sim(q), as.table = T)

q2 = convert.to.character.data.frame(axcf2.1a3.noselecpop300[[2]])
q2[q2 == "H"] = "X"
pull.map(gen.asmap.for.sim(q2), as.table = T)





#normal cm test # 

g = seq(1, 10, 2)
lapply(g, function(x){
    run.1000.sims(p("normalcm.selec", x, "pop300"), 2, num.sims = 100, w.selec = T, selec.str = x, selec.pos = 20, pop.size = 300, cm.pos = seq(1, 100, 2))        
})

run.1000.sims(p("normalcm.noselecpop300"), 2, num.sims = 100, w.selec = F, pop.size = 300, cm.pos = seq(1, 100, 2))        

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/normalcm.selec1pop300")

normalcm.selec1pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/normalcm.selec3pop300")

normalcm.selec3pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/normalcm.selec5pop300")

normalcm.selec5pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/normalcm.selec7pop300")

normalcm.selec7pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/normalcm.selec9pop300")

normalcm.selec9pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/normalcm.noselecpop300")

normalcm.noselecpop300 = all.sim.files

onea.geno = list(normalcm.selec1pop300, normalcm.selec3pop300, normalcm.selec5pop300, normalcm.selec7pop300, normalcm.selec9pop300, normalcm.noselecpop300)



n1.maps = make.maps1(normalcm.selec1pop300)
n3.maps = make.maps1(normalcm.selec3pop300)
n5.maps = make.maps1(normalcm.selec5pop300)
n7.maps = make.maps1(normalcm.selec7pop300)
n9.maps = make.maps1(normalcm.selec9pop300)
nno.maps = make.maps1(normalcm.noselecpop300)

nmaps = list(n1.maps, n3.maps, n5.maps, n7.maps, n9.maps, nno.maps)

lapply(nmaps, function(y){
    mean(sapply(y, function(x) max(attributes(x)$map$pos)))    
})


#### make heatmap for no selection ####

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.6b.noselecpop300")
axcf2.6b.noselecpop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop300")
axcf2.noselec.pop300 = all.sim.files

noselec1 = cbind(axcf2.6b.noselecpop300[[1]], axcf2.noselec.pop300[[1]])
noselec1 = as.data.frame(noselec1)
noselec1 = convert.to.character.data.frame(noselec1)
noselec2 = gen.asmap.for.sim(noselec1, p.value1 = 1e-39)

noselec2.rf = pull.rf(noselec2, what = "rf")

c4 = melt(noselec2.rf)
c4[is.na(c4)] = 0
c4$Var1 = unlist(lapply(1:ncol(noselec2.rf), function(x) rep(x, ncol(noselec2.rf))))
c4$Var2 = rep(1:ncol(noselec2.rf), ncol(noselec2.rf))
colnames(c4)[3] = "RF"

#### plot7 ####

plot7 = ggplot(c4, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = RF)) + 
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() + ggtitle("(a)") +
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.05),
                legend.position = "none") 

plot7 = plot7 + geom_segment(aes(x = 222.5, xend = 222.5, y = 1, yend = 447), color = "#FFFFFF", size = 0.7) +
    geom_segment(aes(x = 1, xend = 447, y = 222.5, yend = 222.5), color = "#FFFFFF", size = 0.7) +
    scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) 
plot7


#### MAKE PUBLICATION HEATMAPS ####

cluster1.11 = make.cluster.comparison(1.11, T)
s(cluster1.11, "rotation1scripts_v4/saved.objects/cluster1.11", "16042019.mstmap.seg.dist.test.R")
load("rotation1scripts_v4/saved.objects/cluster1.11")

seq2 = 10^-seq(40, 50, 1)

lgs1.11 = lapply(seq2, function(x){
    f2(cluster1.11[[1]], x)
})


c1.11.rf = pull.rf(lgs1.11[[3]], what = "rf")

q = pull.map(lgs1.11[[3]], as.table = T)
q = as.data.frame(rbind(q[which(q$chr == "L2"), ], q[which(q$chr == "L1"), ]), stringsAsFactors = F)
q = convert.to.character.data.frame(q)
q[q == "L1"] = "qq1"
q[q == "L2"] = "L1"
q[q == "qq1"] = "L2"

c1.11.rf = c1.11.rf[match(rownames(q), colnames(c1.11.rf)), match(rownames(q), colnames(c1.11.rf))]


c2 = melt(c1.11.rf)
c2[is.na(c2)] = 0
c2$Var1 = unlist(lapply(1:ncol(c1.11.rf), function(x) rep(x, ncol(c1.11.rf))))
c2$Var2 = rep(1:ncol(c1.11.rf), ncol(c1.11.rf))
colnames(c2)[3] = "RF"

# s(c2, "rotation1scripts_v4/saved.objects/c2", "16042019.mstmap.seg.dist.test.R")
load("rotation1scripts_v4/saved.objects/c2")
#### plot5 ####
plot5 = ggplot(c2, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = RF)) + 
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() + ggtitle("(b)") +
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.05),
                legend.position = "none") 

plot5 = plot5 + geom_segment(aes(x = 222.5, xend = 222.5, y = 1, yend = 447), color = "#FFFFFF", size = 0.7) +
    geom_segment(aes(x = 1, xend = 447, y = 222.5, yend = 222.5), color = "#FFFFFF", size = 0.7) +
    scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) 
plot5

heatmap.rf.plot.1.11 = plot5

#second plot

cluster1.105 = make.cluster.comparison(1.105, T)
s(cluster1.105, "rotation1scripts_v4/saved.objects/cluster1.105", "16042019.mstmap.seg.dist.test.R")
load("rotation1scripts_v4/saved.objects/cluster1.105")

seq2 = 10^-seq(40, 50, 1)

lgs1.105 = lapply(seq2, function(x){
    f2(cluster1.105[[1]], x)
})

cluster1.105[[1]] = cluster1.105[[1]][, match(rownames(attributes(lgs1.105[[6]])$map), colnames(cluster1.105[[1]]))]

c1.105.rf = pull.rf(lgs1.105[[6]], what = "rf")



q = pull.map(lgs1.105[[6]], as.table = T)
q = as.data.frame(rbind(q[which(q$chr == "L2"), ], q[which(q$chr == "L1"), ]), stringsAsFactors = F)
q = convert.to.character.data.frame(q)
q[q == "L1"] = "qq1"
q[q == "L2"] = "L1"
q[q == "qq1"] = "L2"

c1.105.rf = c1.105.rf[match(rownames(q), colnames(c1.105.rf)), match(rownames(q), colnames(c1.105.rf))]




c3 = melt(c1.105.rf)
c3[is.na(c3)] = 0
c3$Var1 = unlist(lapply(1:ncol(c1.105.rf), function(x) rep(x, ncol(c1.105.rf))))
c3$Var2 = rep(1:ncol(c1.105.rf), ncol(c1.105.rf))
colnames(c3)[3] = "RF"

# s(c3, "rotation1scripts_v4/saved.objects/c3", "16042019.mstmap.seg.dist.test.R")
load("rotation1scripts_v4/saved.objects/c3")

# c3$Var2 = rep(1:446, 446)
#### plot 6 ####
plot6 = ggplot(c3, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = RF)) + 
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() + ggtitle("(c)") +
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.05)) 

plot6 = plot6 + geom_segment(aes(x = 279.5, xend = 279.5, y = 1, yend = 447), color = "#FFFFFF", size = 0.7) +
    geom_segment(aes(x = 1, xend = 447, y = 279.5, yend = 279.5), color = "#FFFFFF", size = 0.7) +
    scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) 
plot6


library(cowplot)




genetic.rf.plot1 = plot_grid(plot7, plot5, plot6, ncol = 3, rel_widths = c(0.2925, 0.2925, 0.415))


y.grob = textGrob("Marker number", gp=gpar(fontsize=11), rot = 90) #https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
x.grob = textGrob("Marker number", gp=gpar(fontsize=11))

genetic.rf.plot1 = arrangeGrob(genetic.rf.plot1, left = y.grob, bottom = x.grob)

grid.arrange(arrangeGrob(genetic.rf.plot1, left = y.grob, bottom = x.grob))

tiff("rotation1scripts_v4/plots/simulation/pedsim/effect.on.genetic.mapv7_plos.tiff", units = "cm", width = 17.4, height = 6, res = 450)
grid.arrange(arrangeGrob(genetic.rf.plot1, left = y.grob, bottom = x.grob))
dev.off()



plot_grid(plot7, plot5, plot6, ncol = 3, rel_widths = c(0.318, 0.318, 0.364))

#### f8 testing ####

f8.test = make.cluster.comparison.f8(2)


do.processing = function(x){
	colnames(x[[2]][[1]]) = paste0("6b", colnames(x[[2]][[1]]))
	f8.comb = cbind(x[[1]][[1]], x[[2]][[1]])
	f8.comb = as.data.frame(f8.comb)
	f8.comb = convert.to.character.data.frame(f8.comb)
	f8.comb	
}

f8.comb2 = gen.asmap.for.sim(f8.comb, 1e-42, 8)	


f8.1.11 = make.cluster.comparison.f8(1.11)

f8.1.11.2 = do.processing(f8.1.11)

f81.11.map = gen.asmap.for.sim(f8.1.11.2, 1e-46, 8)

do.geno.counts = function(x){
	geno.counts3 = lapply(x, geno.count)
	geno.counts3[c("6bmarker30", "marker200")]	
}





f8.2 = make.cluster.comparison.f8(1.3)

f8.2.2 = do.processing(f8.2)

f82.map = gen.asmap.for.sim(f8.2.2, 1e-46, 8)

order.within.bins(pull.map(f82.map, as.table = T), T, mapmap = T)

do.geno.counts(f8.2.2)


geno.counts3 = lapply(f8.2.2, geno.count)
geno.counts3[c("6bmarker30", "marker200")]


#### PLOTTING UPDATES 24/02/2021 ####


plot7 = plot7 + theme_classic(base_size = 22) + xlab('') + ylab('') +
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.05),
                legend.position = "none") 
plot6 = plot6 + theme_classic(base_size = 22) + xlab('') + ylab('') + 
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.05)) 
plot5 = plot5 + theme_classic(base_size = 22) + xlab('') + ylab('') + 
    theme(axis.title = element_blank(), plot.title = element_text(hjust = 0.05),
                legend.position = "none") 

sgchap.plot4.11a = plot7
sgchap.plot4.11b = plot5
sgchap.plot4.11c = plot6

save(sgchap.plot4.11a, file = '~/project.phd.main/rotation1scripts_v4/saved.objects/sgchap.plot4.11a')
save(sgchap.plot4.11b, file = '~/project.phd.main/rotation1scripts_v4/saved.objects/sgchap.plot4.11b')
save(sgchap.plot4.11c, file = '~/project.phd.main/rotation1scripts_v4/saved.objects/sgchap.plot4.11c')

