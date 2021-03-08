source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")
load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop300")
load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop1000")
load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop10000")


chi300 = convert.recomb.to.seg2(list.of.first.gen.recomb.anywherepop300)
chi1000 = convert.recomb.to.seg2(list.of.first.gen.recomb.anywherepop1000)
chi10000 = convert.recomb.to.seg2(list.of.first.gen.recomb.anywherepop10000)

chi300 = lapply(list.of.first.gen.recomb.anywherepop300, function(x){
	g = convert.recomb.to.seg2(x)
	max(unlist(lapply(g, function(y) y[[1]])))
})

chi1000 = unlist(lapply(list.of.first.gen.recomb.anywherepop1000, function(x){
	g = convert.recomb.to.seg2(x)
	max(unlist(lapply(g, function(y) y[[1]])))
}))

chi10000 = unlist(lapply(list.of.first.gen.recomb.anywherepop10000, function(x){
	g = convert.recomb.to.seg2(x)
	max(unlist(lapply(g, function(y) y[[1]])))
}))

chi300.1 = unlist(lapply(chi300, function(x){
	x[[1]]
}))

chi1000.1 = unlist(lapply(chi1000, function(x){
	x[[1]]
}))

chi10000.1 = unlist(lapply(chi10000, function(x){
	x[[1]]
}))

check.magnitude.of.distortion = function(x){
	#what is the magnitude of distortion in the simulations at the peak
	#args:
	#x - a list of genotyping dataframes 
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		max(g1)
	}))

	g
}

mag300 = check.magnitude.of.distortion(list.of.first.gen.recomb.anywherepop300)
mag1000 = check.magnitude.of.distortion(list.of.first.gen.recomb.anywherepop1000)
mag10000 = check.magnitude.of.distortion(list.of.first.gen.recomb.anywherepop10000)

g = convert.recomb.to.seg2(list.of.first.gen.recomb.anywherepop300[[986]])
g1 = convert.recomb.to.seg2(list.of.first.gen.recomb.anywherepop1000[[986]])

gg1 = (unlist(lapply(g, function(x) x[[1]])))
which(gg1 == max(gg1))

gg1 = (unlist(lapply(g1, function(x) x[[1]])))
which(gg1 == max(gg1))

g = list.of.first.gen.recomb.anywherepop300[[986]][, 95]
g1 = list.of.first.gen.recomb.anywherepop1000[[986]][, 109]

numa.numb = function(x){
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    g = c(a, b)
    names(g) = c("a", "b")
    g
}

numa.numb(g1)

chisq.test(c(77, 56))
chisq.test(c(274, 233))

library(gridExtra)
pdf("rotation1scripts_v4/plots/simulation/selec300pop.pdf", 20, 20)
do.call(grid.arrange, plots1[1:10])
dev.off()

	pdf("rotation1scripts_v4/plots/simulation/selec1000pop.pdf", 20, 20)
	do.call(grid.arrange, plots2[1:10])
	dev.off()

	pdf("rotation1scripts_v4/plots/simulation/selec10000pop.pdf", 20, 20)
	do.call(grid.arrange, plots3[1:10])
	dev.off()


##### PERFORM SIMULATION WITH 96 INDIVIDUALS ####
pdf("rotation1scripts_v4/plots/simulation/plot.w.selec.300.pdf", 20, 20)pdf
corestouse = 30

list.of.first.gen.recomb1a.pop96 = mclapply(1:1000, function(x){
    g = gen.data(225, 96, recombination.profiles = recombination.profiles1a)
    g
}, mc.cores = corestouse)


##### do some plotting #####

pdf("rotation1scripts_v4/plots/simulation/plot.w.selec.300.pdf", 20, 20)
do.call(grid.arrange, plot300)
dev.off()

pdf("rotation1scripts_v4/plots/simulation/plot.w.selec.10000.pdf", 20, 20)
do.call(grid.arrange, plot10000)
dev.off()

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop8")

pdf("rotation1scripts_v4/plots/simulation/plot.w.selec2.10000.pdf", 20, 20)
do.call(grid.arrange, plotss)
dev.off()

smallerpop1 = perform.1000.simulations(96, 10, 200, 225)


pdf("rotation1scripts_v4/plots/simulation/plot.w.selec10.96.pdf", 20, 20)
do.call(grid.arrange, plots96popwselec)
dev.off()


popsizes1 = c(50, 100, 500, 1000)

selectionpressures = c(10, 8, 7, 6, 5, 4, 3, 2)


sims1 = lapply(selectionpressures, function(x) perform.1000.simulations(96, x, 200, 225))

count1 = make.counter()
lapply(sims1, function(x){
	plots1 = make.plots(x[1:10])
	pdf(p("rotation1scripts_v4/plots/simulation/selection96pop", selectionpressures[count1()], ".pdf"), 20, 20)
	do.call(grid.arrange, plots1)
	dev.off()
	})
