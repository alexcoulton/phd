#!/home/ac14037/bin/Rscript

setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

library(parallel)

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

args = commandArgs(trailingOnly = T)

# list.of.first.gen.recomb.prob = lapply(1:20, function(x){
#   g = gen.data(225, 1000, recombination.profiles = recombination.profiles1a)
#   g
# })

# s(list.of.first.gen.recomb.prob, p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob", args[1]), "recomb.simulation2.rscript.R")


list.of.first.gen.recomb.prob.10000pop = lapply(1:20, function(x){
  g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a)
  g
})

s(list.of.first.gen.recomb.prob.10000pop, p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob10000pop", args[1]), "recomb.simulation2.rscript.R")