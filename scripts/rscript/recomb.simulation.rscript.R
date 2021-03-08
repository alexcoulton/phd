#!/home/ac14037/R/bin/Rscript
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

library(parallel)

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

corestouse = 16


listfirstgen.rec.anywhere.225.300 = mclapply(1:1000, function(x){
  g = gen.data(225, 300)
  g
}, mc.cores = corestouse)

list.of.first.gen.recomb.1a.225.300 = mclapply(1:1000, function(x){
  g = gen.data(225, 300, recombination.profiles = recombination.profiles1a)
  g
}, mc.cores = corestouse)

s(listfirstgen.rec.anywhere.225.300, "rotation1scripts_v4/saved.objects/listfirstgen.rec.anywhere.225.300", "recomb.simulation.rscript.R")

s(list.of.first.gen.recomb.1a.225.300, "rotation1scripts_v4/saved.objects/list.of.first.gen.recomb.1a.225.300", "recomb.simulation.rscript.R")



