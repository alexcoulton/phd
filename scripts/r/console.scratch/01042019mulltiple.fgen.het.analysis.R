evaluate.hets = function(df1){
	#gets the % number of heterozygotes out of the entire genotyping dataset
	#	args:
	#	df1 - genotyping dataframe
	num.hets = sum(unlist(lapply(df1, function(x){
		length(which(x == "H"))
	})))

	(num.hets / (nrow(df1) * ncol(df1))) * 100

}

load("rotation1scripts_v4/saved.objects/recomb.sim/fgen/simulations.f5.keep.intermediate.gen")

simulations = simulations.f5.keep.intermediate.gen

library(parallel)
hets = mclapply(simulations, function(x){
	unlist(lapply(x, evaluate.hets))
}, mc.cores = 60)


hets2 = combine.list.of.data.frames(hets)
lapply(hets2, mean)
lapply(hets2, sd)

