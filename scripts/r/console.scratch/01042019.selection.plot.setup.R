source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

selec.96.pop.sims = lapply(seq(12, 20, 2), function(x){
	perform.1000.simulations(96, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = x, selection.pos = 200, corestouse = 60)	
	})


selec.96.pop.sims.analyses = Map(function(x, y){
	perform.analysis(x, T, 200, y, "recombination.profiles1a", T, 2)
	}, selec.96.pop.sims, seq(12, 20, 2))



load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop300.selec20")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec20")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop10000.selec20")



pop96.seg.ratios = lapply(selec.96.pop.sims[[5]], function(x){
		convert.recomb.to.seg(x)
	})

pop96.seg.ratios2 = do.call(rbind, pop96.seg.ratios)

pop300.seg.ratios = lapply(pop300.selec20, function(x){
		convert.recomb.to.seg(x)
	})

pop300.seg.ratios2 = do.call(rbind, pop300.seg.ratios)

pop1000.seg.ratios = lapply(pop1000.selec20, function(x){
		convert.recomb.to.seg(x)
	})

pop1000.seg.ratios2 = do.call(rbind, pop1000.seg.ratios)

pop10000.seg.ratios = lapply(pop10000.selec20, function(x){
		convert.recomb.to.seg(x)
	})

pop10000.seg.ratios2 = do.call(rbind, pop10000.seg.ratios)



	
all.pop.seg.ratios = list(pop96.seg.ratios2, pop300.seg.ratios2, pop1000.seg.ratios2, pop10000.seg.ratios2)


p1 = plotter2(g2, g2$`96mean`, g2$`96sd`, title = "Population Size: 96", ylimits = y2, hline = T)
p2 = plotter2(g2, g2$`300mean`, g2$`300sd`, title = "Population Size: 300", ylimits = y2, hline = T)
p3 = plotter2(g2, g2$`1000mean`, g2$`1000sd`, title = "Population Size: 1000", ylimits = y2, hline = T)
p4 = plotter2(g2, g2$`10000mean`, g2$`10000sd`, title = "Population Size: 10000", ylimits = y2, hline = T)


pop96.f6.sim = perform.1000.simulations.fgen(96, 20, 200, 225, recombination.profiles1a, 60, 6)
pop300.f6.sim = perform.1000.simulations.fgen(300, 20, 200, 225, recombination.profiles1a, 60, 6)
pop1000.f6.sim = perform.1000.simulations.fgen(1000, 20, 200, 225, recombination.profiles1a, 60, 6)
pop10000.f6.sim = perform.1000.simulations.fgen(10000, 20, 200, 225, recombination.profiles1a, 60, 6)

pop96.f6.sim.w.selection = pop96.f6.sim
pop300.f6.sim.w.selection = pop300.f6.sim
pop1000.f6.sim.w.selection = pop1000.f6.sim
pop10000.f6.sim.w.selection = pop10000.f6.sim

s(pop96.f6.sim.w.selection, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/pop96.f6.sim.w.selection", "01042019.selection.plot.setup.R")
s(pop300.f6.sim.w.selection, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/pop300.f6.sim.w.selection", "01042019.selection.plot.setup.R")
s(pop1000.f6.sim.w.selection, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/pop1000.f6.sim.w.selection", "01042019.selection.plot.setup.R")
s(pop10000.f6.sim.w.selection, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/pop10000.f6.sim.w.selection", "01042019.selection.plot.setup.R")



all.f6.pop.w.selection = list(pop96.f6.sim, pop300.f6.sim, pop1000.f6.sim, pop10000.f6.sim)

pop96.anal = mclapply(pop96.f6.sim.w.selection, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop96.anal2 = combine.list.of.data.frames(pop96.anal)

s(pop96.anal2, "rotation1scripts_v4/saved.objects/pop96.anal2", "01042019.selection.plot.setup.R")

pop96.mean = sapply(pop96.anal2, mean)
pop96.sd = sapply(pop96.anal2, sd)

pop300.anal = mclapply(pop300.f6.sim.w.selection, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop300.anal2 = combine.list.of.data.frames(pop300.anal)

pop300.mean = sapply(pop300.anal2, mean)
pop300.sd = sapply(pop300.anal2, sd)

pop1000.anal = mclapply(pop1000.f6.sim.w.selection, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop1000.anal2 = combine.list.of.data.frames(pop1000.anal)

pop1000.mean = sapply(pop1000.anal2, mean)
pop1000.sd = sapply(pop1000.anal2, sd)


pop10000.anal = mclapply(pop10000.f6.sim.w.selection, function(x){
	convert.recomb.to.seg(x)
	}, mc.cores = 60)

pop10000.anal2 = combine.list.of.data.frames(pop10000.anal)

pop10000.mean = sapply(pop10000.anal2, mean)
pop10000.sd = sapply(pop10000.anal2, sd)


all.pop.dat.f6.w.selec = data.frame(pop96.mean, pop300.mean, pop1000.mean, pop10000.mean, pop96.sd, pop300.sd, pop1000.sd, pop10000.sd)




pop96.f6.sim.w.selection.analysis = perform.analysis(pop96.f6.sim.w.selection, T, 200, 20, "recombination.profiles1a", F, 6)
pop300.f6.sim.w.selection.analysis = perform.analysis(pop300.f6.sim.w.selection, T, 200, 20, "recombination.profiles1a", F, 6)
pop1000.f6.sim.w.selection.analysis = perform.analysis(pop1000.f6.sim.w.selection, T, 200, 20, "recombination.profiles1a", F, 6)
pop10000.f6.sim.w.selection.analysis = perform.analysis(pop10000.f6.sim.w.selection, T, 200, 20, "recombination.profiles1a", F, 6)

all.f6.analyses = combine.list.of.data.frames(list(pop96.f6.sim.w.selection.analysis, pop300.f6.sim.w.selection.analysis, pop1000.f6.sim.w.selection.analysis, pop10000.f6.sim.w.selection.analysis))

write.csv(all.f6.analyses, "rotation1scripts_v4/temp/all.f6.analyses", row.names = F)

#check the levels of heterozygosity - this is actually an F5
evaluate.hets(pop96.f6.sim.w.selection[[1]])
evaluate.hets(pop300.f6.sim.w.selection[[1]])
evaluate.hets(pop1000.f6.sim.w.selection[[1]])
evaluate.hets(pop10000.f6.sim.w.selection[[1]])

#grab the average number of A and B genotypes to calculate chi-square significance limits
sum(sapply(as.data.frame(t(sapply(pop96.f6.sim.w.selection[[1]], geno.count))), mean)[c("A", "B")])
sum(sapply(as.data.frame(t(sapply(pop300.f6.sim.w.selection[[1]], geno.count))), mean)[c("A", "B")])
sum(sapply(as.data.frame(t(sapply(pop1000.f6.sim.w.selection[[1]], geno.count))), mean)[c("A", "B")])
sum(sapply(as.data.frame(t(sapply(pop10000.f6.sim.w.selection[[1]], geno.count))), mean)[c("A", "B")])


#histogram of peak of distortion

pop96.anal3 = as.data.frame(t(pop96.anal2))
nrow(pop96.anal3)

pop96.peak = sapply(pop96.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})

pop300.anal3 = as.data.frame(t(pop300.anal2))
nrow(pop300.anal3)

pop300.peak = sapply(pop300.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})


pop1000.anal3 = as.data.frame(t(pop1000.anal2))
nrow(pop1000.anal3)

pop1000.peak = sapply(pop1000.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})

pop10000.anal3 = as.data.frame(t(pop10000.anal2))
nrow(pop10000.anal3)

pop10000.peak = sapply(pop10000.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
})

all.pop.peak = list(pop96.peak, pop300.peak, pop1000.peak, pop10000.peak)
s(all.pop.peak, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/all.pop.peak", "01042019.selection.plot.setup.R")



lapply(all.f6)

s(all.f6.pop.w.selection, "rotation1scripts_v4/saved.objects/")



f6.no.sele.1000 = perform.1000.simulations.fgen(96, recombination.profile = recombination.profiles1a, corestouse = 60, fgen = 6)
f6.no.sele.1000 = perform.1000.simulations.fgen(300, recombination.profile = recombination.profiles1a, corestouse = 60, fgen = 6)
f6.no.sele.1000 = perform.1000.simulations.fgen(1000, recombination.profile = recombination.profiles1a, corestouse = 60, fgen = 6)
f6.no.sele.1000 = perform.1000.simulations.fgen(10000, recombination.profile = recombination.profiles1a, corestouse = 60, fgen = 6)

