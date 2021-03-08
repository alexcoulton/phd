source("rotation1scripts_v4/scripts/r/functions.R")

load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop300")
load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop1000")
load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop10000")

source("rotation1scripts_v4/scripts/r/modelling/analysisv2.R")


checklist1 = function(y){
	lapply(c(0.05, 0.01, 0.001, 0.0001), function(x){
		number.dist.markers(y, x)
	})	
}

markers1 = checklist1(list.of.first.gen.recomb.anywherepop300)
markers2 = checklist1(list.of.first.gen.recomb.anywherepop1000)
markers3 = checklist1(list.of.first.gen.recomb.anywherepop10000)


