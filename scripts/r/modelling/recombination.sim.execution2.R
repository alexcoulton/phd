#!/home/ac14037/bin/Rscript

#execution of recombination / segregation distortion simulation
setwd("/home/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")
source("rotation1scripts_v4/scripts/r/functions.R")
load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

#### ANALYSIS ####

corestouse = 30

list.of.first.gen.w.prob.w.selec = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 10, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec2 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 8, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec2, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec2", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec3 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 7, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec3, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec3", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec4 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 6, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec4, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec4", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec5 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 5, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec5, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec5", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec6 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 4, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec6, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec6", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec7 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 3, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec7, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec7", "recombination.sim.execution2.R")

list.of.first.gen.w.prob.w.selec8 = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = 2, selection.pos = 200)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.w.prob.w.selec8, "rotation1scripts_v4/saved.objects/recomb.sim/wselection/list.of.first.gen.w.prob.w.selec8", "recombination.sim.execution2.R")
