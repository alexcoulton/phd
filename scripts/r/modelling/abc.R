#Approximate Bayesian Computation

source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")


load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

load("rotation1scripts_v4/saved.objects/segratios1a")

#get data to compare simulation to with summary statistic 
seg.ratios = unname(segratios1a$g651[2:226])

abc.run1 = mclapply(1:20000, function(x){
    selection.strength = round(rexp(1, 0.2)) + 4
    selection.position = round(runif(1, 0, 225))
    
    g = gen.data(225, 94, recombination.profiles = recombination.profiles1a, selection = selection.strength, selection.pos = selection.position)
    seg1 = unname(convert.recomb.to.seg(g))
    
    sumstat1 = sum((seg.ratios - seg1)^2)
    
    final = c(sumstat1, selection.position, selection.strength)
    names(final) = c("sumstat1", "selection.position", "selection.strength")
    final
}, mc.cores = 50)



#### test code ####

# count1 = make.counter()
# lapply(1:8, function(x){
# 	i = count1()
# 	seg.ratios = convert.recomb.to.seg(init.simulations[[i]])

# 	abc.sim1 = lapply(1:20000, function(x){
# 		selection.strength = round(rexp(1, 0.2)) + 4
# 		selection.position = round(runif(1, 0, 225))

# 		genotype.data = gen.data(225, init.sim.parameters$pop.size[[i]], recombination.profiles = recombination.profiles1a, selection = selection.strength, selection.pos = selection.position)
# 		seg1 = unname(convert.recomb.to.seg(g))

# 		sumstat1 = sum((seg.ratios - seg1)^2)

# 		final = c(genotype.data, sumstat1, selection.position, selection.strength)
# 		names(final) = c("genotype.data", "sumstat1", "selection.position", "selection.strength")
# 		final
# 	})

# 	abc.sim1
# })
	


