#!/home/ac14037/bin/Rscript

#execution of recombination / segregation distortion simulation
setwd("/home/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")
load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

library(parallel)

#### ANALYSIS ####

corestouse = 50


pop.sizes = seq(50, 1000, 50)
pop.sizes2 = seq(1050, 4000, 50)
pop.sizes3 = seq(2050, 3000, 50)



list.of.pop.size.sims = list()
for(i in pop.sizes){
	print(p("doing pop size: ", i))
	list.of.first.gen.recomb.anywhere = mclapply(1:100, function(x){
  		g = gen.data(225, i)
  		g
	}, mc.cores = corestouse)	


	list.of.pop.size.sims = c(list.of.pop.size.sims, list.of.first.gen.recomb.anywhere)
	print(p("done: ", i))
}

count1 = make.counter()
ptm = proc.time()
proc.time() - ptm
simulations.recomb.any.pop.50.to.10000 = lapply(pop.sizes, function(y){
	print(p("doing pop size ", y))

	list.of.first.gen.recomb.anywhere = mclapply(1:100, function(x){
  		g = gen.data(225, y)
  		g
	}, mc.cores = corestouse)		
})

simulations.recomb.any.pop.1050.to.4000 = lapply(pop.sizes2, function(y){
	print(p("doing pop size ", y))

	list.of.first.gen.recomb.anywhere = mclapply(1:100, function(x){
  		g = gen.data(225, y)
  		g
	}, mc.cores = corestouse)		
})

s(simulations.recomb.any.pop.50.to.10000, "rotation1scripts_v4/saved.objects/recomb.sim/simulations.recomb.any.pop.50.to.10000", "recombination.sim.execution5.R")


simulations.recomb.any.pop.50.to.1000.mag = lapply(simulations.recomb.any.pop.50.to.1000, check.magnitude.of.distortion)

simulations.recomb.any.pop.1050.to.4000.mag = lapply(simulations.recomb.any.pop.1050.to.4000, check.magnitude.of.distortion)