source("rotation1scripts_v4/scripts/r/modelling/pedigree.sim.functions.R")

pop.sizes = seq(10, 1000, 5)
selec.strengths.smaller = seq(4, 20, 2)

all.param2 = expand.grid(pop.sizes, selec.strengths.smaller)
colnames(all.param2) = c("pop.sizes", "selec.strengths") 


Map(function(x, y){
    run.1000.sims(p("axcf2.heatmap.selec", y, "pop", x), 2, num.sims = 40, w.selec = T, selec.str = y, selec.pos = 10, pop.size = x, cm.pos = seq(1, 100, 5), F)
}, all.param2$pop.sizes, all.param2$selec.strengths)

files1 = list.files("~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/simresults", pattern = "heatmap")

all.sims = lapply(files1, function(x){
	load(p("~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", x))
	all.sim.files
	})


g = min(sapply(all.sims, length))

all.sims2 = lapply(all.sims, function(x){
	x[1:g]
	})


all.sims3 = lapply(all.sims2, function(x){
	lapply(x, function(y){
		g2 = convert.recomb.to.seg(y)
		sum((g2 - 0.5)^2)
		})
	})

all.sims4 = sapply(all.sims3, function(x){
	mean(unlist(x))
	})

selec.strengths1 = as.numeric(multi.str.split(multi.str.split(files1, "selec", 2), "pop", 1))

pop1 = as.numeric(multi.str.split(files1, "pop", 2))

dat1 = data.frame(all.sims4, selec.strengths1, pop1)
colnames(dat1) = c("dev", "s.str", "pop")

pedsim.heatmap.data = dat1

s(pedsim.heatmap.data, "rotation1scripts_v4/saved.objects/pedsim.heatmap.data", "24042019.pedsim.heatmap.R")