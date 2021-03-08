#### ANALYSIS ####

# axc.chr1a.cm = cxaf2.c651.asmap.list[[15]]$pos
# s(axc.chr1a.cm, "rotation1scripts_v4/saved.objects/axc.chr1a.cm", "05042019.pedigree.sim.R")

mclapply(1:2, function(x){
    complete.f5.selection.procedure(length(axc.chr1a.cm), 300, 6, axc.chr1a.cm, "rotation1scripts_v4/original_data/simulation/pedigreesim/", paste0("axc1", x), "KOSAMBI", T, 5, 200)    
}, mc.cores = 2)

library(parallel)

ptm = proc.time()


lapply(c(96, 300, 1000, 10000), function(x){
    run.1000.sims(p("axcf2.noselec.pop", x), 2, num.sims = 1000, w.selec = F, pop.size = x)    
})

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop96"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 200, pop.size = 96)    
})

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop300"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 200, pop.size = 300)    
})

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop1000"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 200, pop.size = 1000)    
})

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop10000"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 200, pop.size = 10000)    
})

lapply(c(96, 300, 1000, 10000), function(x){
    run.1000.sims(p("axcf6selec20.pop", x), f.gen = 6, num.sims = 1100, w.selec = T, selec.str = 20, pop.size = x)    
})

lapply(c(96, 300, 1000, 10000), function(x){
    run.1000.sims(p("axcf6.noselec.pop", x), f.gen = 6, num.sims = 1000, w.selec = F, pop.size = x)    
})

# axc.chr1a.cm

lapply(c(96, 300, 1000, 10000), function(x){
    run.1000.sims(p("axcf2.sparse.noselec.pop", x), 2, num.sims = 1000, w.selec = F, pop.size = x, cm.pos = 60)    
})

lapply(c(96, 300, 1000), function(x){
    run.1000.sims(p("axcf2.sparse.2markers.noselec.pop", x), 2, num.sims = 1000, w.selec = F, pop.size = x, cm.pos = c(20, 60))
})

lapply(c(96, 300, 1000), function(x){
    run.1000.sims(p("axcf2.sparse.2markers_2.noselec.pop", x), 2, num.sims = 1000, w.selec = F, pop.size = x, cm.pos = c(59, 60))
})

lapply(c(96), function(x){
    run.1000.sims(p("axcf2.ultradense.noselec.pop", x), 2, num.sims = 1000, w.selec = F, pop.size = x, cm.pos = seq(0, 130, 0.01))    
})






lapply(c(96, 300, 1000, 10000), function(x){
    run.1000.sims(p("axcf2.dense.noselec.pop", x), 2, num.sims = 1000, w.selec = F, pop.size = x, cm.pos = seq(0, 130, 0.1))    
})



# run more sims to make up to 1000 in the sims that had some errors before

run.1000.sims(p("axcf6selec20.pop", 96, ".batch2"), f.gen = 6, num.sims = 500, w.selec = T, selec.str = 20, pop.size = 96)    
run.1000.sims(p("axcf6selec20.pop", 1000, ".batch2"), f.gen = 6, num.sims = 20, w.selec = T, selec.str = 20, pop.size = 1000)    


lapply(seq(12, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop96.batch2"), 2, num.sims = 600, w.selec = T, selec.str = x, selec.pos = 200, pop.size = 96)    
})


run.1000.sims(p("axcf2.dense.noselec.pop", 1000, ".batch2"), 2, num.sims = 50, w.selec = F, pop.size = 1000, cm.pos = seq(0, 130, 0.1))    



lapply(1:1000, function(x){
    complete.f5.selection.procedure(length(axc.chr1a.cm), 300, 6, axc.chr1a.cm, "rotation1scripts_v4/original_data/simulation/pedigreesim/", paste0("test", x), "KOSAMBI", T, 5, 200)        
})





run.1000.sims(p("axcf2noselecpop300.het.test"), 2, num.sims = 1000, pop.size = 300)    
run.1000.sims(p("axcf3noselecpop300.het.test"), 3, num.sims = 1000, pop.size = 300)    
run.1000.sims(p("axcf4noselecpop300.het.test"), 4, num.sims = 1000, pop.size = 300)    
run.1000.sims(p("axcf5noselecpop300.het.test"), 5, num.sims = 1000, pop.size = 300)    
run.1000.sims(p("axcf6noselecpop300.het.test"), 6, num.sims = 1000, pop.size = 300)    





run.1000.sims("axcf6noselec", 6, num.sims = 100, w.selec = F)

run.1000.sims("axcf6selec", 6, num.sims = 100, w.selec = T, selec.str = 10, selec.pos = 200)


load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec")
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf6selec")



#### make additional chromosomes to check segregation distortion effect on clustering ####

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2.6b.selec", x, "pop300"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    
})

run.1000.sims(p("axcf2.6b.noselecpop300"), 2, num.sims = 1000, w.selec = F, pop.size = 300, cm.pos = cxa.6b.cm)    


run.1000.sims(p("axcf2.6b.selec1.2pop300"), 2, num.sims = 1000, w.selec = T, selec.str = 1.2, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    


run.1000.sims(p("axcf2selec1.2pop300"), 2, num.sims = 1000, w.selec = T, selec.str = 1.2, selec.pos = 200, pop.size = 300, cm.pos = axc.chr1a.cm)    


run.1000.sims(p("axcf2.6b.selec1.01pop300"), 2, num.sims = 1000, w.selec = T, selec.str = 1.01, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    


run.1000.sims(p("axcf2selec1.01pop300"), 2, num.sims = 1000, w.selec = T, selec.str = 1.01, selec.pos = 200, pop.size = 300, cm.pos = axc.chr1a.cm)    




pop.sizes = seq(10, 1000, 5)
selec.strengths = seq(2, 50, 2)

Map(function(x, y){
    run.1000.sims(p("axcf2.heatmap.selec, ", y, "pop", x), 2, num.sims = 10, w.selec = T, selec.str = x, selec.pos = 10, pop.size = 300, cm.pos = seq(1, 100, 5))
}, pop.sizes, selec.strengths)





#### two selection pressures on different chromosome functions ####




make.cluster.comparison = function(selec.str){
    run.1000.sims(p("axcf2.6b.selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    
    
    
    run.1000.sims(p("axcf2selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = 200, pop.size = 300, cm.pos = axc.chr1a.cm)
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec", selec.str, "pop300"))
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

make.cluster.comparison.one.selec = function(selec.str){
    run.1000.sims(p("axcf2.6b.selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    
    
    
    run.1000.sims(p("axcf2selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec", selec.str, "pop300"))
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


run.1000.sims(p("axcf2.6b.noselecpop300"), 2, num.sims = 5, w.selec = F, pop.size = 300, cm.pos = cxa.6b.cm)




make.cluster.comparison.one.selec = function(selec.str){
    run.1000.sims(p("axcf2.6b.selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = T, selec.str = selec.str, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    
    
    
    run.1000.sims(p("axcf2selec", selec.str, "pop300"), 2, num.sims = 1, w.selec = F, pop.size = 300, cm.pos = axc.chr1a.cm)
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec", selec.str, "pop300"))
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

make.cluster.comparison.f8 = function(selec.str){
    run.1000.sims(p("axcf8.6b.selec", selec.str, "pop300"), 8, num.sims = 2, w.selec = T, selec.str = selec.str, selec.pos = 30, pop.size = 300, cm.pos = cxa.6b.cm)    
    
    
    run.1000.sims(p("axcf8selec", selec.str, "pop300"), 8, num.sims = 2, w.selec = T, selec.str = selec.str, selec.pos = 200, pop.size = 300, cm.pos = axc.chr1a.cm)
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf8selec", selec.str, "pop300"))
    axc1a = all.sim.files
    
    load(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf8.6b.selec", selec.str, "pop300"))
    axc6b = all.sim.files
    
    list(axc1a, axc6b)
    
    # axc6b = lapply(axc6b, function(y){
    #     colnames(y) = paste0("6b", colnames(y))
    #     y
    # })
    # 
    # 
    # sim1.2comb = Map(function(sim6b, sim1a){
    #     g = as.data.frame(cbind(sim6b, sim1a), stringsAsFactors = F)
    #     
    # }, axc1a, axc6b)
    # 
    # sim1.2comb
    
}

10^seq(40, 50, 1)


unique()



g = lapply(all.sim.files, convert.recomb.to.seg)

par(mfrow = c(2, 2))
lapply(sample(100, 4), function(x){
    plot(g[[x]], ylim = c(0.3, 0.7))
    abline(h = 0.5)
})


g2 = sapply(all.sim.files, function(x){
    min(convert.recomb.to.seg2(x, p.value1 = T))
})



g = read.ped("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axc.pop300.w.selec5.f6.1.output_genotypes.dat_wselec.f6.dat")
g = convert.ped.sim.to.zygote(g)
g = reset.rownames(g)
rownames(g) = paste0("marker", rownames(g))
g = as.data.frame(t(g))




seg.ratios1 = sapply(g, function(x){
    q = geno.count(x)[c("A", "B")]
    q[1] / (q[1] + q[2])
})

plot(seg.ratios1, ylim = c(0.3, 0.7))

#
proc.time() - ptm




#selection procedure for    F2
make.pedsim.config.files.and.run()

#make 10k gametes
make.pedsim.config.files.and.run(length(cxaf2.c651.asmap.list[[15]]$pos), 10000, 2, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "axc1a.pop10k", "KOSAMBI", T)





# todo: check that selection.test.f3 contains enough gametes with a "B" at locus 200 to replace the gametes that contain an 
# "A" in the previous simulation for those lines that are to undergo selection





################


ped.file.1 = strsplit(colnames(ped.selec), "F")
ped.file.1 = sapply(ped.file.1, function(x) x[[2]])
ped.file.2 = strsplit(ped.file.1, "_")
as.numeric(sapply(ped.file.2, function(x) x[[1]])) + 1



lapply(as.data.frame(t(ped.zygote[, grep("F3", colnames(ped.zygote))])), geno.count)
lapply(as.data.frame(t(ped.zygote[, grep("F2", colnames(ped.zygote))])), geno.count)




ped2 = ped1[, grep("F3", colnames(ped1))]

names(selection.ind)









#testing runtimes 
ptm <- proc.time()
make.pedsim.config.files.and.run(length(cxaf2.c651.asmap.list[[15]]$pos), 300, 3, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "testf3", map.function = "KOSAMBI", T)
proc.time() - ptm

ptm <- proc.time()
read.delim("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axc1a.pop300.output_genotypes.dat", sep = "\t", stringsAsFactors = F)
proc.time() - ptm



ptm <- proc.time()
g = gen.data(225, 300, prob.vec = prob.vec3)
proc.time() - ptm

g2 = convert.sim.data.to.rqtl.format(g)
write.csv(g2, "rotation1scripts_v4/temp/g2.csv", row.names = F)
g3 = rqtl.read("rotation1scripts_v4/temp/g2.csv")
heatMap(g3, what = "rf")
mean(countXO(g3))

make.pedsim.config.files.and.run(length(cxaf2.c651.asmap.list[[15]]$pos), 10000, 6, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "axc1a.pop10k", T)

cxa1a.cm = cxaf2.c651.asmap.list[[15]]$pos
recomb.probs = diff(cxa1a.cm) / sum(diff(cxa1a.cm))
c(0, recomb.probs)

prob.vec3



library(readr)
big.ped.sim.pop = read_delim("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axc1a.pop10k.output_genotypes.dat", 
                                                         delim = "\t")



ped.sim = read.delim("rotation1scripts_v4/original_data/simulation/PedigreeSim v2.0/axc1a.pop300_results1_out_genotypes.dat", sep = "\t", stringsAsFactors = F)


rownames(ped.sim) = ped.sim[, 1]
ped.sim = ped.sim[, -1]
ped.sim = convert.to.character.data.frame(ped.sim)

ped.sim2 = lapply(seq(1, ncol(ped.sim), 2), function(q){
    combined.geno1 = unlist(Map(function(x, y){
        if(x == "A" & y == "A") z = "A"
        if(x == "A" & y == "B") z = "H"
        if(x == "B" & y == "A") z = "H"
        if(x == "B" & y == "B") z = "B"
        z
    }, ped.sim[, q], ped.sim[, (q + 1)]))    
    combined.geno1
})

ped.sim3 = do.call(cbind, ped.sim2)
ped.sim3 = as.data.frame(t(ped.sim3))

ped.sim4 = convert.sim.data.to.rqtl.format(ped.sim3)
write.csv(ped.sim4, "rotation1scripts_v4/temp/ped.sim.csv", row.names = F)


ped.sim.rqtl = read.cross("csv", "rotation1scripts_v4/temp/", "ped.sim.csv", genotypes=c("A", "H", "B"),estimate.map=F)

mean(countXO(ped.sim.rqtl))

par(mfrow = c(1, 1))
heatMap(ped.sim.rqtl, what = "rf")

gen1.rf = as.data.frame(pull.rf(ped.sim.rqtl))
max(sapply(gen1.rf, function(x) max(na.omit(x))))
# pop10k.max.rf = max(sapply(gen1.rf, function(x) max(na.omit(x))))


rqtl.read = function(rqtl.csv.path, ...){
    g = strsplit(rqtl.csv.path, "/")
    filename = g[[1]][length(g[[1]])]
    path1 = paste(g[[1]][-length(g[[1]])], collapse = "/")
    path1 = paste0(path1, "/")
    read.cross("csv", path1, filename, genotypes = c("A", "H", "B"), estimate.map = F, ...)
}

rqtl1 = rqtl.read("rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/cs.x.para.flipped.csv", F.gen = 7)

rqtl1.rf = pull.rf(rqtl1)
rqtl2.rf = as.data.frame(rqtl1.rf)
hist(rqtl1.rf)
ncol(rqtl1.rf)
nrow(rqtl1.rf)
# dir.create("rotation1scripts_v4/plots/recombination.frequency")
png("rotation1scripts_v4/plots/recombination.frequency/cs.x.paragon.rf.png", 500, 500)
hist(as.vector(rqtl1.rf), breaks = 200, main = "Chinese Spring X Paragon F7, 300 Individuals", xlab = "Recombination Frequency")
dev.off()
mean(na.omit(as.vector(rqtl1.rf)))
sd(na.omit(as.vector(rqtl1.rf)))



rqtl2 = rqtl.read("rotation1scripts_v4/processed_data/genotypes/avalon.x.cadenza/a.x.c.flipped.csv", crosstype = "dh")
rqtl2.rf = pull.rf(rqtl2)

rqtl2.lod = pull.rf(rqtl2, what = "lod")

png("rotation1scripts_v4/plots/recombination.frequency/a.x.c.rf.png", 500, 500)
hist(rqtl2.rf, breaks = 200, main = "Avalon X Cadenza Double Haploid, 130 Individuals", xlab = "Recombination Frequency")
dev.off()


mean(na.omit(as.vector(rqtl2.rf)))
sd(na.omit(as.vector(rqtl2.rf)))


#### MISC ####

col.w.b = sapply(geno10k, function(x){
    if(x[200] == "B") return(T)
    if(x[200] == "A") return(F)
})

col.w.b = unlist(col.w.b)

geno10k.w.b = geno10k[, which(col.w.b == T)]
geno10k.w.b.f2 = geno10k.w.b[, grep("F2", colnames(geno10k.w.b))]







#### test f3 selection procedure ####

make.pedsim.config.files.and.run(length(cxaf2.c651.asmap.list[[15]]$pos), 300, 3, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "testf3", map.function = "KOSAMBI", T)

f3.selection.procedure("testf3", "fj843", 3)


q = read.ped("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/testf3.output_genotypes.dat_wselec.dat")
q2 = read.ped("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/testf3.output_genotypes.dat")

length(which(sapply(q, function(x) x[[200]] == "B")))
length(which(sapply(q2, function(x) x[[200]] == "B")))

length(which(sapply(q, function(x) x[[200]] == "A")))
length(which(sapply(q2, function(x) x[[200]] == "A")))


make.pedsim.config.files.and.run(length(cxaf2.c651.asmap.list[[15]]$pos), 300, 2, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "testf2", map.function = "KOSAMBI", T)

f2.selection.procedure("testf2", 10, 200)

q = read.ped("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/testf2.output_genotypes.dat")
q2 = read.ped("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/testf2.output_genotypes.dat_wselec.dat")




make.pedsim.config.files.and.run(length(cxaf2.c651.asmap.list[[15]]$pos), 300, 2, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "axc1a.f2.300pop", "KOSAMBI", T)




ptm <- proc.time()
perform.1000.simulations.fgen(300, 20, 200, 225, recombination.profiles1a, 1, 6, number.simulations = 1)
proc.time() - ptm


ptm <- proc.time()
complete.f5.selection.procedure(length(cxaf2.c651.asmap.list[[15]]$pos), 300, 4, cxaf2.c651.asmap.list[[15]]$pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", "new.testf5v2", "KOSAMBI", T, 5, 200)
proc.time() - ptm



