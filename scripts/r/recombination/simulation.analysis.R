#to be used in conjunction with recombination.distribution.09072019v2.R

#### RECOMBINATION SIMULATION ####
source('~/project.phd.main/rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R')
source('~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R')
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/recombination.distribution.test.6b")

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/recombination.distribution.test.sparse")

recomb.test.sparse = all.sim.files
# recombination.distribution.test = all.sim.files
# recombination.distribution.test6b = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/recombination.distribution.test.6b.500.pop")
recombination.distribution.test6b.500.pop = all.sim.files


load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/recombination.distribution.test.axp3b.500.pop.500.sims")
recombination.distribution.test.axp3b.500.pop = all.sim.files


sim1 = convert.sim.to.rqtl.ske.format(recomb.test.sparse[[1]])


#plot simulated genotyping datasets as genetic map comparisons
sim.plots = lapply(1:10, function(x){
    sim2 = lapply(recomb.test.sparse[sample(99, 4)], convert.sim.to.rqtl.ske.format)
    
    sim.centi = lapply(sim2, extract.centimorgan.from.genotypes)
    
    sim.plot = prepare.genetic.map.comparison.plots(sim2, sim.centi)
    sim.plot[[1]]
})


do.call(grid.arrange, sim.plots)




#### sparse marker set ####
# 
# r.sp1 = lapply(recomb.test.sparse, process.sim.recomb, sample.size = 50)
# r.sp1.2 = lapply(r.sp1, makehistlist.new)
# 
# r.sp2 = lapply(recomb.test.sparse, process.sim.recomb, sample.size = 100)
# r.sp2.2 = lapply(r.sp2, makehistlist.new)
# 
# r.sp3 = lapply(recomb.test.sparse, process.sim.recomb, sample.size = 1000)
# r.sp3.2 = lapply(r.sp3, makehistlist.new)
# 
# 
# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(r.sp1.2[[sample(99, 1)]])[[1]], dist.from.cent(r.sp1.2[[sample(99, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 
# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(r.sp2.2[[sample(99, 1)]])[[1]], dist.from.cent(r.sp2.2[[sample(99, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 
# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(r.sp3.2[[sample(99, 1)]])[[1]], dist.from.cent(r.sp3.2[[sample(99, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 


# rd1 = process.sim.recomb(recombination.distribution.test6b[[1]])
# rd2 = process.sim.recomb(recombination.distribution.test6b[[2]])
# rd3 = process.sim.recomb(recombination.distribution.test6b[[3]])
# rd4 = process.sim.recomb(recombination.distribution.test6b[[4]])

#### AxP 3B marker distribution #####

# rd.all = list(rd1, rd2, rd3, rd4)


axp.3b.map = axp.liss[[2]][which(axp.liss[[2]]$chr == "3B"), ]

axp.3b.map$marker = paste0("marker", 1:nrow(axp.3b.map))

rd.all = mclapply(recombination.distribution.test.axp3b.500.pop, process.sim.recomb, mc.cores = 50)

rd.all2 = lapply(rd.all, function(x){
	 x$phys.pos = axp.3b.map[match(x$probe.before.transition, axp.3b.map$marker), ]$phys
	 x
 })

rd.all250 = mclapply(recombination.distribution.test.axp3b.500.pop, process.sim.recomb, sample.size = 250, 
mc.cores = 50)
rd.all100 = mclapply(recombination.distribution.test.axp3b.500.pop, process.sim.recomb, sample.size = 100)
rd.all30 = mclapply(recombination.distribution.test.axp3b.500.pop, process.sim.recomb, sample.size = 30)

hs.all = lapply(rd.all, makehistlist.new)
hs.all250 = lapply(rd.all250, makehistlist.new)
hs.all100 = lapply(rd.all100, makehistlist.new)
hs.all30 = lapply(rd.all30, makehistlist.new)

par(mfrow = c(2, 2))
hist(hs.all[[sample(100, 1)]][[1]], breaks = 100)
hist(hs.all250[[sample(100, 1)]][[1]], breaks = 100)
hist(hs.all100[[sample(100, 1)]][[1]], breaks = 100)
hist(hs.all30[[sample(100, 1)]][[1]], breaks = 100)


hs.all.hist1 = make.ggplot.histogram(hs.all[[1]], num.bins = 100, breaks = seq(0, 223, 1), ylabel = "Frequency") + xlab("Recombination Position") + theme(axis.title = element_text(size = 8)) + ggtitle("a")
hs.all.hist2 = make.ggplot.histogram(hs.all250[[1]], num.bins = 100, breaks = seq(0, 223, 1), ylabel = "Frequency") + xlab("Recombination Position") + theme(axis.title = element_text(size = 8)) + ggtitle("b")
hs.all.hist3 = make.ggplot.histogram(hs.all100[[1]], num.bins = 100, breaks = seq(0, 223, 1), ylabel = "Frequency") + xlab("Recombination Position") + theme(axis.title = element_text(size = 8)) + ggtitle("c")
hs.all.hist4 = make.ggplot.histogram(hs.all30[[1]], num.bins = 100, breaks = seq(0, 223, 1), ylabel = "Frequency") + xlab("Recombination Position") + theme(axis.title = element_text(size = 8)) + ggtitle("d")


grid.arrange(hs.all.hist1, hs.all.hist2, hs.all.hist3, hs.all.hist4)


rd.all.ind = mclapply(rd.all, ind.avg, use.centromere = T, scale.by.arm = T, mc.cores = 50)
rd.all250.ind = mclapply(rd.all250, ind.avg, use.centromere = T, scale.by.arm = T, mc.cores = 50)
rd.all100.ind = mclapply(rd.all100, ind.avg, use.centromere = T, scale.by.arm = T, mc.cores = 50)
rd.all30.ind = mclapply(rd.all30, ind.avg, use.centromere = T, scale.by.arm = T, mc.cores = 50)


rd.all.means = sapply(rd.all.ind, function(x) mean(x$marker.dist))
rd.all250.means = sapply(rd.all250.ind, function(x) mean(x$marker.dist))
rd.all100.means = sapply(rd.all100.ind, function(x) mean(x$marker.dist))
rd.all30.means = sapply(rd.all30.ind, function(x) mean(x$marker.dist))

rd.means = data.frame(rd.all.means, rd.all250.means, rd.all100.means, rd.all30.means)

rd.means2 = melt(rd.means)
rd.means2$variable = as.factor(rd.means2$variable)

levels(rd.means2$variable) = c("500", "250", "100", "30")


rd.anova = aov(rd.means2$value ~ rd.means2$variable)
summary(rd.anova)




hist(rd.all.means, breaks = 30)
hist(rd.all250.means, breaks = 30)
hist(rd.all100.means, breaks = 30)
hist(rd.all30.means, breaks = 30)

histogram.sim.means1 = ggplot(rd.means2, aes(x = value)) + geom_histogram(binwidth = 0.2) + facet_grid(rows = vars(variable)) + theme_bw() + ylab("Frequency") +
    xlab("Mean MRD")





histpop500 = make.ggplot.histogram(rd.all.means, num.bins = 100, breaks = seq(65, 85, 0.2), ylabel = "Frequency") + coord_cartesian(ylim = c(0, 10))
histpop250 = make.ggplot.histogram(rd.all250.means, num.bins = 100, breaks = seq(65, 85, 0.2), ylabel = "Frequency") + coord_cartesian(ylim = c(0, 10))
histpop100 = make.ggplot.histogram(rd.all100.means, num.bins = 100, breaks = seq(65, 85, 0.2), ylabel = "Frequency") + coord_cartesian(ylim = c(0, 10))
histpop30 = make.ggplot.histogram(rd.all30.means, num.bins = 100, breaks = seq(65, 85, 0.2), ylabel = "Frequency") + coord_cartesian(ylim = c(0, 10))



all.rd.ind = list(rd.all.ind, rd.all250.ind, rd.all100.ind, rd.all30.ind)

s(all.rd.ind, "rotation1scripts_v4/saved.objects/all.rd.ind", "simulation.analysis.R")

combinations1 = as.data.frame(combn(1:200, 4))

combinations1 = combinations1[, sample(ncol(combinations1), 10000)]

#TO RUN ON WILKINS #SEE recomb.simulation.R
num.sim.sig2 = mclapply(all.rd.ind, function(x){
    length(which(sapply(combinations1, function(y){
        
        
        m1 = x[[y[1]]]$phys.dist.short
        m2 = x[[y[2]]]$phys.dist.short
        m3 = x[[y[3]]]$phys.dist.short
        m4 = x[[y[4]]]$phys.dist.short
        
        m1 = data.frame(m1, "1")
        m2 = data.frame(m2, "2")
        m3 = data.frame(m3, "3")
        m4 = data.frame(m4, "4")
        
        colnames(m1) = c("phys.dist", "treatment")
        colnames(m2) = c("phys.dist", "treatment")
        colnames(m3) = c("phys.dist", "treatment")
        colnames(m4) = c("phys.dist", "treatment")
        
        m5 = bind_rows(m1, m2, m3, m4)
        # browser()
        # abs(mean(m1) - mean(m2))
        
        # browser()
        # hist(m1)
        # hist(m2)
        # wilcox.test(m1, m2)$p.value
        
        kruskal.test(m5$phys.dist ~ m5$treatment)$p.value
        
    }) < 0.05))
}, mc.cores = 50)






num.sim.sig2 = round(((unlist(num.sim.sig1) / ncol(combinations1)) * 100), digits = 2)




rdstat500 = length(which(round(unlist(sapply(1:100, function(x){
    to.compare = sample(100, 2)
    # browser()
    
    g = wilcox.test(rd.all.ind[[to.compare[1]]][[2]], rd.all.ind[[to.compare[2]]][[2]])
    g[3]
})), digits = 5) < 0.05))


rdstat250 = length(which(round(unlist(sapply(1:100, function(x){
    to.compare = sample(100, 2)
    # browser()
    
    g = wilcox.test(rd.all250.ind[[to.compare[1]]][[2]], rd.all250.ind[[to.compare[2]]][[2]])
    g[3]
})), digits = 5) < 0.05))

rdstat100 = length(which(round(unlist(sapply(1:100, function(x){
    to.compare = sample(100, 2)
    # browser()
    
    g = wilcox.test(rd.all100.ind[[to.compare[1]]][[2]], rd.all100.ind[[to.compare[2]]][[2]])
    g[3]
})), digits = 5) < 0.05))

rdstat30 = length(which(round(unlist(sapply(1:100, function(x){
    to.compare = sample(100, 2)
    # browser()
    
    g = wilcox.test(rd.all30.ind[[to.compare[1]]][[2]], rd.all30.ind[[to.compare[2]]][[2]])
    g[3]
    
    list(mean(rd.all30[[to.compare[1]]]$marker.dist), mean(rd.all30[[to.compare[2]]]$marker.dist))
    
    
})), digits = 5) < 0.05))






ind.avg(rd.all[[1]], manual.center = 111)




# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(rd.all[[sample(999, 1)]])[[1]], dist.from.cent(rd.all[[sample(999, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 
# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(rd.all40[[sample(999, 1)]])[[1]], dist.from.cent(rd.all40[[sample(999, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 
# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(rd.all20[[sample(999, 1)]])[[1]], dist.from.cent(rd.all20[[sample(999, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 
# length(which(round(unlist(sapply(1:1000, function(x){
#     g = wilcox.test(dist.from.cent(rd.all10[[sample(999, 1)]])[[1]], dist.from.cent(rd.all10[[sample(999, 1)]])[[1]])    
#     g[3]
# })), digits = 5) < 0.05))
# 
# hist(dist.from.cent(rd.all10[[10]])[[1]], breaks = 100)
# hist(dist.from.cent(rd.all10[[200]])[[1]], breaks = 100)
# 
# wilcox.test(dist.from.cent(rd.all10[[10]])[[1]], dist.from.cent(rd.all10[[200]])[[1]])



#### examine Sacha apogee X Paragon map #### 

allenaxp = read_csv("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic map AxP.csv")

allenaxp2 = t(allenaxp)
allenaxp2 = allenaxp2[-3, ]
allenaxp2 = convert.to.character.data.frame(allenaxp2)
colnames(allenaxp2) = allenaxp2[1, ]
allenaxp2 = allenaxp2[-1, ]
allenaxp2[allenaxp2 == "BB"] = "B"
allenaxp2[allenaxp2 == "AA"] = "A"
allenaxp2[allenaxp2 == "AB"] = "H"
allenaxp2[allenaxp2 == "BB"] = "B"
allenaxp3 = cbind(allenaxp2[, 1:2], allenaxp2)
allenaxp3[, 2] = rownames(allenaxp3)
allenaxp3[1, 2] = ""
allenaxp3[, 1] = c(NA, 2:nrow(allenaxp3))
tempcolname1 = colnames(allenaxp3)[1:2]
colnames(allenaxp3)[1:2] = ""
colnames(allenaxp3)[3:4] = tempcolname1


allenaxp4 = allenaxp3[, c(1, 2, which(as.character(unique(allenaxp3[1, ])) %in% listofwheatchromosomes))]




allenaxpphys = get.phys.for.chromos(allenaxp4, T, T, "-")

allenaxp.recomb = detect.recombination(allenaxp4)


allenaxp.recomb$phys.pos = allenaxpphys[match(allenaxp.recomb$probe.before.transition, allenaxpphys$marker), ]$phys.pos

allenaxpind.avg = ind.avg(allenaxp.recomb)

hist(allenaxpind.avg$phys.dist)


run.comparison.samples = function(x){
    individuals1 = sample(unique(allenaxpind.avg$individual), x)
    
    allenaxpind.avg2 = allenaxpind.avg[which(allenaxpind.avg$individual %in% individuals1[1:(length(individuals1) / 2)]), ]
    allenaxpind.avg3 = allenaxpind.avg[which(allenaxpind.avg$individual %in% individuals1[((length(individuals1) / 2) + 1):length(individuals1)]), ]
    
    wilcox.test(allenaxpind.avg2$phys.dist, allenaxpind.avg3$phys.dist)
    
}


p.vals30 = sapply(1:1000, function(x){
    run.comparison.samples(30)$p.value    
})

length(which(p.vals30 < 0.05))

p.vals60 = sapply(1:1000, function(x){
    run.comparison.samples(60)$p.value    
})

length(which(p.vals60 < 0.05))


p.vals100 = sapply(1:1000, function(x){
    run.comparison.samples(100)$p.value    
})

length(which(p.vals100 < 0.05))


p.vals300 = sapply(1:1000, function(x){
    run.comparison.samples(300)$p.value    
})

length(which(p.vals300 < 0.05))




# hist(allenaxpind.avg2$phys.dist)








