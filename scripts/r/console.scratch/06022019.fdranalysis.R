source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")


# all.recomb.prob2 = list()



# for(i in 1:50){
# 	load(p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob", i))
# 	all.recomb.prob2 = c(all.recomb.prob2, list.of.first.gen.recomb.prob)
# }

# s(all.recomb.prob2, "rotation1scripts_v4/saved.objects/recomb.sim/all.recomb.prob2")

# all.recomb.prob.pop10000 = list()



# for(i in 1:50){
# 	load(p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob10000pop", i))
# 	all.recomb.prob.pop10000 = c(all.recomb.prob.pop10000, list.of.first.gen.recomb.prob.10000pop)
# }

# s(all.recomb.prob.pop10000, "rotation1scripts_v4/saved.objects/recomb.sim/all.recomb.prob.pop10000")

# recomb.profile1a.no.selec.sims = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a)

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec2")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec3")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec4")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec5")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec6")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec7")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec8")





all.selec = list(list.of.first.gen.w.prob.w.selec, list.of.first.gen.w.prob.w.selec2, list.of.first.gen.w.prob.w.selec3, list.of.first.gen.w.prob.w.selec4, list.of.first.gen.w.prob.w.selec5, list.of.first.gen.w.prob.w.selec6, list.of.first.gen.w.prob.w.selec7, list.of.first.gen.w.prob.w.selec8)

count1 = make.counter()
selection.strengths = c(10, 8, 7, 6, 5, 4, 3, 2)
selec.analysis = lapply(all.selec, function(x){
	i = count1()	
	perform.analysis(x, T, 200, selection.strengths[i], "recombination.profiles1a", 1000, F2 = T)
	})




load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop2")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop3")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop4")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop5")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop6")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop7")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.1000pop8")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop2")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop3")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop4")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop5")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop6")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop7")
load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop8")

all.selec1000pop = list(list.of.first.gen.w.prob.w.selec.1000pop, list.of.first.gen.w.prob.w.selec.1000pop2, list.of.first.gen.w.prob.w.selec.1000pop3, list.of.first.gen.w.prob.w.selec.1000pop4, list.of.first.gen.w.prob.w.selec.1000pop5, list.of.first.gen.w.prob.w.selec.1000pop6, list.of.first.gen.w.prob.w.selec.1000pop7, list.of.first.gen.w.prob.w.selec.1000pop8)

all.selec10000pop = list(list.of.first.gen.w.prob.w.selec.10000pop, list.of.first.gen.w.prob.w.selec.10000pop2, list.of.first.gen.w.prob.w.selec.10000pop3, list.of.first.gen.w.prob.w.selec.10000pop4, list.of.first.gen.w.prob.w.selec.10000pop5, list.of.first.gen.w.prob.w.selec.10000pop6, list.of.first.gen.w.prob.w.selec.10000pop7, list.of.first.gen.w.prob.w.selec.10000pop8)

count1 = make.counter()
selection.strengths = c(10, 8, 7, 6, 5, 4, 3, 2)
selec.analysis1000pop = lapply(all.selec1000pop, function(x){
	i = count1()	
	perform.analysis(x, T, 200, selection.strengths[i], "recombination.profiles1a", 1000, F2 = T)
	})


count1 = make.counter()
selection.strengths = c(10, 8, 7, 6, 5, 4, 3, 2)
selec.analysis10000pop = lapply(all.selec10000pop, function(x){
	i = count1()	
	perform.analysis(x, T, 200, selection.strengths[i], "recombination.profiles1a", 1000, F2 = T)
	})

all.1000pop = combine.list.of.data.frames(selec.analysis1000pop)
all.10000pop = combine.list.of.data.frames(selec.analysis10000pop)

write.csv(all.1000pop, "rotation1scripts_v4/temp/all.1000pop.csv", row.names = F)
write.csv(all.10000pop, "rotation1scripts_v4/temp/all.10000pop.csv", row.names = F)




selection.strengths2 = c(20, 18, 16, 14, 12)

lapply(selection.strengths2, function(x){
	 = 

})


pop300.selec20 = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 20, selection.pos = 200, corestouse = 60)	
s(pop300.selec20, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop300.selec20", "06022019fdranalysis.R")
pop300.selec18 = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 18, selection.pos = 200, corestouse = 60)	
s(pop300.selec18, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop300.selec18", "06022019fdranalysis.R")
pop300.selec16 = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 16, selection.pos = 200, corestouse = 60)	
s(pop300.selec16, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop300.selec16", "06022019fdranalysis.R")
pop300.selec14 = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 14, selection.pos = 200, corestouse = 60)	
s(pop300.selec14, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop300.selec14", "06022019fdranalysis.R")
pop300.selec12 = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 12, selection.pos = 200, corestouse = 60)	
s(pop300.selec12, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop300.selec12", "06022019fdranalysis.R")

pop1000.selec20 = perform.1000.simulations(1000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 20, selection.pos = 200, corestouse = 60)	
s(pop1000.selec20, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec20", "06022019fdranalysis.R")
pop1000.selec18 = perform.1000.simulations(1000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 18, selection.pos = 200, corestouse = 60)	
s(pop1000.selec18, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec18", "06022019fdranalysis.R")
pop1000.selec16 = perform.1000.simulations(1000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 16, selection.pos = 200, corestouse = 60)	
s(pop1000.selec16, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec16", "06022019fdranalysis.R")
pop1000.selec14 = perform.1000.simulations(1000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 14, selection.pos = 200, corestouse = 60)	
s(pop1000.selec14, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec14", "06022019fdranalysis.R")
pop1000.selec12 = perform.1000.simulations(1000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 12, selection.pos = 200, corestouse = 60)	
s(pop1000.selec12, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec12", "06022019fdranalysis.R")


list.of.pop300 = list(pop300.selec20, pop300.selec18, pop300.selec16, pop300.selec14, pop300.selec12)
list.of.pop1000 = list(pop1000.selec20, pop1000.selec18, pop1000.selec16, pop1000.selec14, pop1000.selec12)

i = 1
pop300.analyses = lapply(list.of.pop300, function(x){
	g = perform.analysis(x, T, 200, selection.strengths2[i], "recombination.profiles1a", 1000)
	i <<- i + 1
	g
	})

i = 1
pop1000.analyses = lapply(list.of.pop1000, function(x){
	g = perform.analysis(x, T, 200, selection.strengths2[i], "recombination.profiles1a", 1000)
	i <<- i + 1
	g
	})


pop10000.selec20 = perform.1000.simulations(10000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 20, selection.pos = 200, corestouse = 60)	
s(pop10000.selec20, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop10000.selec20", "06022019fdranalysis.R")
pop10000.selec20.analysis = perform.analysis(pop10000.selec20, T, 200, 20, "recombination.profiles1a", 1000)
rm(pop10000.selec20)


pop10000.selec18 = perform.1000.simulations(10000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 18, selection.pos = 200, corestouse = 60)	
s(pop10000.selec18, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop10000.selec18", "06022019fdranalysis.R")
pop10000.selec18.analysis = perform.analysis(pop10000.selec18, T, 200, 18, "recombination.profiles1a", 1000)
rm(pop10000.selec18)
gc()


pop10000.selec16 = perform.1000.simulations(10000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 16, selection.pos = 200, corestouse = 60)	
s(pop10000.selec16, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop10000.selec16", "06022019fdranalysis.R")
pop10000.selec16.analysis = perform.analysis(pop10000.selec16, T, 200, 16, "recombination.profiles1a", 1000)
rm(pop10000.selec16)
gc()


pop10000.selec14 = perform.1000.simulations(10000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 14, selection.pos = 200, corestouse = 60)	
s(pop10000.selec14, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop10000.selec14", "06022019fdranalysis.R")
pop10000.selec14.analysis = perform.analysis(pop10000.selec14, T, 200, 14, "recombination.profiles1a", 1000)
rm(pop10000.selec14)
gc()

pop10000.selec12 = perform.1000.simulations(10000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 12, selection.pos = 200, corestouse = 60)	
s(pop10000.selec12, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop10000.selec12", "06022019fdranalysis.R")
pop10000.selec12.analysis = perform.analysis(pop10000.selec12, T, 200, 12, "recombination.profiles1a", 1000)
rm(pop10000.selec12)
gc()



analyses.list = list(pop10000.selec20.analysis, pop10000.selec18.analysis, pop10000.selec16.analysis, pop10000.selec14.analysis, pop10000.selec12.analysis)

analyses1 = combine.list.of.data.frames(analyses.list)
write.csv("rotation1scripts_v4/temp/analyses1.csv", row.names = F)





pop1000.selec12 = perform.1000.simulations(1000, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = 12, selection.pos = 200, corestouse = 60)	
s(pop1000.selec12, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/pop1000.selec12", "06022019fdranalysis.R")

selec.96.pop.sims = lapply(seq(12, 20, 2), function(x){
	perform.1000.simulations(96, number.of.markers = 225, recombination.profile = recombination.profiles1a, selection.strength = x, selection.pos = 200, corestouse = 60)	
	})