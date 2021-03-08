#!/home/ac14037/R/bin/Rscript
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

library(parallel)

load("rotation1scripts_v4/saved.objects/recombination.profiles5b")

corestouse = 16

list.of.first.gen.recomb.5b.225.300 = mclapply(1:1000, function(x){
  g = gen.data(399, 300, recombination.profiles = recombination.profiles5b)
  g
}, mc.cores = corestouse)



s(list.of.first.gen.recomb.5b.225.300, "rotation1scripts_v4/saved.objects/list.of.first.gen.recomb.5b.225.300", "recomb.sim.5b.R")



