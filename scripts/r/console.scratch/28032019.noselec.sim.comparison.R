source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")
load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb1a.pop96")
load("rotation1scripts_v4/saved.objects/recomb.sim/recomb.profile1a.no.selec.sims")
load("rotation1scripts_v4/saved.objects/all.recomb.prob2")
load("rotation1scripts_v4/saved.objects/all.recomb.prob.pop10000")

list.of.no.selec = list(list.of.first.gen.recomb1a.pop96, recomb.profile1a.no.selec.sims, all.recomb.prob2, all.recomb.prob.pop10000)

list.of.no.selec2 = lapply(list.of.no.selec, function(x){
	x[1:10]
	})

library(parallel)

pop96.seg = mclapply(list.of.first.gen.recomb1a.pop96, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop300.seg = mclapply(recomb.profile1a.no.selec.sims, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop1000.seg = mclapply(all.recomb.prob2, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop10000.seg = mclapply(all.recomb.prob.pop10000, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop96.seg2 = as.data.frame(do.call(rbind, pop96.seg))
pop300.seg2 = as.data.frame(do.call(rbind, pop300.seg))
pop1000.seg2 = as.data.frame(do.call(rbind, pop1000.seg))
pop10000.seg2 = as.data.frame(do.call(rbind, pop10000.seg))


pop96.means = unlist(lapply(pop96.seg2, mean))
pop96.sds = unlist(lapply(pop96.seg2, sd))


pop300.means = lapply(pop300.seg2, mean)
pop300.sds = lapply(pop300.seg2, sd)

pop1000.means = lapply(pop1000.seg2, mean)
pop1000.sds = lapply(pop1000.seg2, sd)

pop10000.means = lapply(pop10000.seg2, mean)
pop10000.sds = lapply(pop10000.seg2, sd)

pop.stats = cbind(pop96.means, pop96.sds, pop300.means, pop300.sds, pop1000.means, pop1000.sds, pop10000.means, pop10000.sds)

s(pop.stats, "rotation1scripts_v4/saved.objects/pop.stats", "28032019.noselec.sim.comparison.R")
