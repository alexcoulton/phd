load("rotation1scripts_v4/saved.objects/axcf6selec20.pop96.comb")

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6selec20.pop300")

axcf6selec20.pop300 = all.sim.files

load("rotation1scripts_v4/saved.objects/axcf6selec20.pop1000.comb")

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6selec20.pop10000")

axcf6selec20.pop10000 = all.sim.files

axcf6selec20.pop96.seg = lapply(axcf6selec20.pop96.comb, convert.recomb.to.seg)

axcf6selec20.pop96.seg2 = do.call(rbind, axcf6selec20.pop96.seg)

axcf6selec20.pop300.seg = lapply(axcf6selec20.pop300, convert.recomb.to.seg)

axcf6selec20.pop300.seg2 = do.call(rbind, axcf6selec20.pop300.seg)

axcf6selec20.pop1000.seg = lapply(axcf6selec20.pop1000.comb, convert.recomb.to.seg)

axcf6selec20.pop1000.seg2 = do.call(rbind, axcf6selec20.pop1000.seg)

axcf6selec20.pop10000.seg = lapply(axcf6selec20.pop10000, convert.recomb.to.seg)

axcf6selec20.pop10000.seg2 = do.call(rbind, axcf6selec20.pop10000.seg)

all.pop.seg.ratios = list(axcf6selec20.pop96.seg2, axcf6selec20.pop300.seg2, axcf6selec20.pop1000.seg2, axcf6selec20.pop10000.seg2)

axcf6selec20.pop96.seg2.mean = sapply(as.data.frame(axcf6selec20.pop96.seg2), mean)
axcf6selec20.pop300.seg2.mean = sapply(as.data.frame(axcf6selec20.pop300.seg2), mean)
axcf6selec20.pop1000.seg2.mean = sapply(as.data.frame(axcf6selec20.pop1000.seg2), mean)
axcf6selec20.pop10000.seg2.mean = sapply(as.data.frame(axcf6selec20.pop10000.seg2), mean)

axcf6selec20.pop96.seg2.sd = sapply(as.data.frame(axcf6selec20.pop96.seg2), sd)
axcf6selec20.pop300.seg2.sd = sapply(as.data.frame(axcf6selec20.pop300.seg2), sd)
axcf6selec20.pop1000.seg2.sd = sapply(as.data.frame(axcf6selec20.pop1000.seg2), sd)
axcf6selec20.pop10000.seg2.sd = sapply(as.data.frame(axcf6selec20.pop10000.seg2), sd)



all.stats1 = list(axcf6selec20.pop96.seg2.mean, axcf6selec20.pop300.seg2.mean, axcf6selec20.pop1000.seg2.mean, axcf6selec20.pop10000.seg2.mean, axcf6selec20.pop96.seg2.sd, axcf6selec20.pop300.seg2.sd, axcf6selec20.pop1000.seg2.sd, axcf6selec20.pop10000.seg2.sd)

all.stats2 = as.data.frame(do.call(cbind, all.stats1))
colnames(all.stats2) = c("axcf6selec20.pop96.seg2.mean", "axcf6selec20.pop300.seg2.mean", "axcf6selec20.pop1000.seg2.mean", "axcf6selec20.pop10000.seg2.mean", "axcf6selec20.pop96.seg2.sd", "axcf6selec20.pop300.seg2.sd", "axcf6selec20.pop1000.seg2.sd", "axcf6selec20.pop10000.seg2.sd")

write.csv(all.stats2, "rotation1scripts_v4/temp/all.stats2.csv", row.names = F)


#histogram of peak of distortion

pop96.anal3 = as.data.frame(t(axcf6selec20.pop96.seg2))
nrow(pop96.anal3)

pop96.peak = sapply(pop96.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})

pop300.anal3 = as.data.frame(t(axcf6selec20.pop300.seg2))
nrow(pop300.anal3)

pop300.peak = sapply(pop300.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})


pop1000.anal3 = as.data.frame(t(axcf6selec20.pop1000.seg2))
nrow(pop1000.anal3)

pop1000.peak = sapply(pop1000.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})

pop10000.anal3 = as.data.frame(t(axcf6selec20.pop10000.seg2))
nrow(pop10000.anal3)

pop10000.peak = sapply(pop10000.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})

all.pop.peak = list(pop96.peak, pop300.peak, pop1000.peak, pop10000.peak)

all.pop.peak.pedsim = all.pop.peak
s(all.pop.peak.pedsim, "rotation1scripts_v4/saved.objects/all.pop.peak.pedsim", "15042019.make.seg.dist.pedsim.plot.R")
















#make plot data for no selec


load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6.noselec.pop96")

axcf6.noselec.pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6.noselec.pop300")

axcf6.noselec.pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6.noselec.pop1000")

axcf6.noselec.pop1000 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6.noselec.pop10000")

axcf6.noselec.pop10000 = all.sim.files



axcf6.noselec.pop96.seg = lapply(axcf6.noselec.pop96, convert.recomb.to.seg)
axcf6.noselec.pop300.seg = lapply(axcf6.noselec.pop300, convert.recomb.to.seg)
axcf6.noselec.pop1000.seg = lapply(axcf6.noselec.pop1000, convert.recomb.to.seg)
axcf6.noselec.pop10000.seg = lapply(axcf6.noselec.pop10000, convert.recomb.to.seg)

axcf6.noselec.pop96.seg2 = do.call(rbind, axcf6.noselec.pop96.seg)
axcf6.noselec.pop300.seg2 = do.call(rbind, axcf6.noselec.pop300.seg)
axcf6.noselec.pop1000.seg2 = do.call(rbind, axcf6.noselec.pop1000.seg)
axcf6.noselec.pop10000.seg2 = do.call(rbind, axcf6.noselec.pop10000.seg)







axcf6.noselec.pop96.seg2.mean = sapply(as.data.frame(axcf6.noselec.pop96.seg2), mean)
axcf6.noselec.pop300.seg2.mean = sapply(as.data.frame(axcf6.noselec.pop300.seg2), mean)
axcf6.noselec.pop1000.seg2.mean = sapply(as.data.frame(axcf6.noselec.pop1000.seg2), mean)
axcf6.noselec.pop10000.seg2.mean = sapply(as.data.frame(axcf6.noselec.pop10000.seg2), mean)

axcf6.noselec.pop96.seg2.sd = sapply(as.data.frame(axcf6.noselec.pop96.seg2), sd)
axcf6.noselec.pop300.seg2.sd = sapply(as.data.frame(axcf6.noselec.pop300.seg2), sd)
axcf6.noselec.pop1000.seg2.sd = sapply(as.data.frame(axcf6.noselec.pop1000.seg2), sd)
axcf6.noselec.pop10000.seg2.sd = sapply(as.data.frame(axcf6.noselec.pop10000.seg2), sd)

noselec.list = list(axcf6.noselec.pop96.seg2.mean, axcf6.noselec.pop300.seg2.mean, axcf6.noselec.pop1000.seg2.mean, axcf6.noselec.pop10000.seg2.mean, axcf6.noselec.pop96.seg2.sd, axcf6.noselec.pop300.seg2.sd, axcf6.noselec.pop1000.seg2.sd, axcf6.noselec.pop10000.seg2.sd)

noselec.data = as.data.frame(do.call(cbind, noselec.list))

colnames(noselec.data) = c("axcf6.noselec.pop96.seg2.mean", "axcf6.noselec.pop300.seg2.mean", "axcf6.noselec.pop1000.seg2.mean", "axcf6.noselec.pop10000.seg2.mean", "axcf6.noselec.pop96.seg2.sd", "axcf6.noselec.pop300.seg2.sd", "axcf6.noselec.pop1000.seg2.sd", "axcf6.noselec.pop10000.seg2.sd")


write.csv(noselec.data, "rotation1scripts_v4/temp/noselec.data.csv", row.names = F)