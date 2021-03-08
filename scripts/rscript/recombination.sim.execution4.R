#!/home/ac14037/bin/Rscript

#execution of recombination / segregation distortion simulation
setwd("/home/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")
load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

library(parallel)

#### ANALYSIS ####

corestouse = 30

list.of.first.gen.w.prob.w.selec.10000pop5 = mclapply(1:1000, function(x){
  g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a, selection = 5, selection.pos = 200)
  g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec.10000pop5, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop5", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec.10000pop6 = mclapply(1:1000, function(x){
  g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a, selection = 4, selection.pos = 200)
  g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec.10000pop6, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop6", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec.10000pop7 = mclapply(1:1000, function(x){
  g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a, selection = 3, selection.pos = 200)
  g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec.10000pop7, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop7", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec.10000pop8 = mclapply(1:1000, function(x){
  g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a, selection = 2, selection.pos = 200)
  g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec.10000pop8, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec.10000pop8", "recombination.sim.execution2.R")
