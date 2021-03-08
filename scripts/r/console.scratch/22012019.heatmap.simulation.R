setwd("/home/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")
load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

pop.sizes = seq(10, 1000, 5)
selec.strengths = seq(2, 50, 2)

selec.strengths.smaller = seq(2, 20, 2)

all.param2 = expand.grid(pop.sizes, selec.strengths.smaller)
colnames(all.param2) = c("pop.sizes", "selec.strengths") 

all.param = expand.grid(pop.sizes, selec.strengths) #9975 rows
colnames(all.param) = c("pop.sizes", "selec.strengths")

library(parallel)

x = 10
y = 10

count1 = make.counter()

# all.simulations = mcMap(function(x, y){
# 	gen.data(225, x, recombination.profiles = recombination.profiles1a, selection = y, selection.pos = 200)
# 	print(p("done ", count1()))
# 	}, all.param$pop.sizes, all.param$selec.strengths, mc.cores = 40)


all.sims2 = mcMap(function(x, y){
	gen.data(50, x, selection = y, selection.pos = 25)
	}, all.param$pop.sizes[1:2000], all.param$selec.strengths[1:2000], mc.cores = 50)


all.sims3 = mcMap(function(x, y){
	gen.data(50, x, selection = y, selection.pos = 25)
	}, all.param$pop.sizes[2001:5000], all.param$selec.strengths[2001:5000], mc.cores = 60)

all.sims4 = mcMap(function(x, y){
	gen.data(50, x, selection = y, selection.pos = 25)
	}, all.param$pop.sizes[5001:9975], all.param$selec.strengths[5001:9975], mc.cores = 60)

all.sims.comb = c(all.sims2, all.sims3, all.sims4)


all.seg1 = lapply(all.sims.comb, function(x){
	g = convert.recomb.to.seg(x)
	sum((g - 0.5)^2)

	})


all.sims.replicates = lapply(1:4, function(x){
	all.sims.replicates = mcMap(function(x, y){
		gen.data(50, x, selection = y, selection.pos = 25)
		}, all.param2$pop.sizes, all.param2$selec.strengths, mc.cores = 50)
	all.sims.replicates		
	})



# all.sims.rep.calc = lapply(all.sims.replicates, function(x){
# 	lapply(x, function(y){
# 		g = convert.recomb.to.seg(x)
# 		sum((g - 0.5)^2)
# 		})
# 	})


all.seg2 = unlist(lapply(all.sims.replicates[[1]], function(x){
	g = convert.recomb.to.seg(x)
	sum((g - 0.5)^2)

	}))


all.seg3 = unlist(lapply(all.sims.replicates[[2]], function(x){
	g = convert.recomb.to.seg(x)
	sum((g - 0.5)^2)

	}))

all.seg4 = unlist(lapply(all.sims.replicates[[3]], function(x){
	g = convert.recomb.to.seg(x)
	sum((g - 0.5)^2)

	}))

all.seg5 = unlist(lapply(all.sims.replicates[[4]], function(x){
	g = convert.recomb.to.seg(x)
	sum((g - 0.5)^2)

	}))






all.seg.average = (all.seg2 + all.seg3 + all.seg4 + all.seg5) / 4

all.param2$heat.value = as.numeric(all.seg.average)

heatmap.sim.df2 = all.param2
s(heatmap.sim.df2, "rotation1scripts_v4/saved.objects/heatmap.sim.df2", "22012019.heatmap.simulation.R")

all.param$heat.value = as.numeric(unlist(all.seg1))


library(ggplot2)

