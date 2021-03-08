#!/home/ac14037/R/bin/Rscript
#Approximate Bayesian Computation

setwd("/home/ac14037/project.phd.main")

source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

#COMMENTS ARE FOR 1A SIMULATION - SAVED TO sim1.txt
# load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

# load("rotation1scripts_v4/saved.objects/segratios1a")

load("rotation1scripts_v4/saved.objects/recombination.profiles5b")

load("rotation1scripts_v4/saved.objects/segratios5b")

seg.ratios = unname(segratios5b[[1]]$g651)
library(parallel)


write(c("simstat1", "selection.position", "selection.strength"), "rotation1scripts_v4/original_data/abc.simulation/sim5b.txt", append = T)

abc.run1 = mclapply(1:20000, function(x){
  # selection.strength = round(rexp(1, 0.2)) + 4
  # selection.position = round(runif(1, 0, 225))
  
  selection.strength = round(rexp(1, 0.1)) + 4
  selection.position = round(runif(1, 0, 399))


  # g = gen.data(225, 94, recombination.profiles = recombination.profiles1a, selection = selection.strength, selection.pos = selection.position)

  g = gen.data(399, 94, recombination.profiles = recombination.profiles5b, selection = selection.strength, selection.pos = selection.position)
  seg1 = unname(convert.recomb.to.seg(g))
  
  sumstat1 = sum((seg.ratios - seg1)^2)
  
  final = c(sumstat1, selection.position, selection.strength)
  names(final) = c("sumstat1", "selection.position", "selection.strength")
  write(final, "rotation1scripts_v4/original_data/abc.simulation/sim5b.txt", append = T)
  final
  
}, mc.cores = 16)



#lapply(abc.run1, function(x){
#	write(x, "rotation1scripts_v4/original_data/abc.simulation/sim1.txt", append = T)
#})



