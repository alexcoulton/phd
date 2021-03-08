setwd("~/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

all.recomb.prob2 = list()
for(i in 1:50){
  
  load(p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob", i))
  all.recomb.prob2 = c(all.recomb.prob2, list.of.first.gen.recomb.prob)
  
}

all.recomb.prob.pop10000 = list()
for(i in 1:50){
  
  load(p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob10000pop", i))
  all.recomb.prob.pop10000 = c(all.recomb.prob.pop10000, list.of.first.gen.recomb.prob.10000pop)
  
}

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec2")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec3")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec4")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec5")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec6")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec7")

load("rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec8")

load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop300")

load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop1000")

load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop10000")

load("rotation1scripts_v4/saved.objects/list.of.first.gen.recomb.5b.225.300")