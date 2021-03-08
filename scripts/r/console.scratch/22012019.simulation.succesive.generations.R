source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")


produce.replicates = function(){
	selec.strengths = c(10, 8, 6, 4, 2)

	selec.data1 = lapply(selec.strengths, function(x){
		gen.data(225, 300, recombination.profiles = recombination.profiles1a, selection = x, selection.pos = 200)	
	})

	# no.selec = list(gen.data(225, 300, recombination.profiles = recombination.profiles1a))

	#generate 10 successive generations for all selection strengths
	generations1 = Map(function(data1, sel.strength){
		nextgen = gen.next.generation(data1, 225, sel.strength, 200, recombination.profiles = recombination.profiles1a)
		gen2 = list(nextgen)
		all.gens = lapply(1:10, function(x){
			nextgen <<- gen.next.generation(nextgen, 225, sel.strength, 200, recombination.profiles = recombination.profiles1a)
			nextgen
			})
		c(gen2, all.gens)
	}, selec.data1, selec.strengths)	

	generations1
}


produce.replicate.no.selec = function(){
	data1 = gen.data(225, 300, recombination.profiles = recombination.profiles1a)

	nextgen = gen.next.generation(data1, 225, recombination.profiles = recombination.profiles1a)
	gen2 = list(nextgen)
	all.gens = lapply(1:10, function(x){
		nextgen <<- gen.next.generation(nextgen, 225, recombination.profiles = recombination.profiles1a)
		nextgen
		})
	c(list(data1), gen2, all.gens)
	
}


replicates.no.selec = lapply(1:10, function(x){
	produce.replicate.no.selec()
	})



replicates = lapply(1:10, function(x){
	produce.replicates()
	})


library(ggplot2)
library(gridExtra)
g = make.plots(replicates[[1]][[1]])
pdf("rotation1scripts_v4/plots/simulation/successive.generations/plot1.pdf", 10, 10)
do.call(grid.arrange, g)
dev.off()


q = make.plots(replicates.no.selec[[1]])
pdf("rotation1scripts_v4/plots/simulation/successive.generations/plot1.no.selec.pdf", 10, 10)
do.call(grid.arrange, q)
dev.off()

#examine deviation from 1:1 ratio in successive generations
sum12 = lapply(replicates.no.selec[[1]], function(x){
	g = convert.recomb.to.seg(x)
	sum((g - 0.5)^2)

	#browser()
	})




seg.data = lapply(generations1, function(x){
	seg.data1 = lapply(x, function(y){
		seg1 = convert.recomb.to.seg(y)
		seg1
		})

	sumofs1 = lapply(seg.data1, function(x){
		sum((x - seg.data1[[11]])^2)
		})
	unlist(sumofs1)
	})


seg.data2 = as.data.frame(seg.data)
colnames(seg.data2) = c("10", "8", "6", "4", "2")


calculate.heterozygotes = function(geno.data1){
	#calculates the number of heterozygous markers in a genotyping dataset
	#args:
	#	geno.data1 - a genotyping data.frame in rQTL format
	sum(unlist(lapply(geno.data1, function(col1) length(which(col1 == "H")))))
}


#count number of heterozygotes in each generation
het.data.no.selec = lapply(replicates.no.selec, function(x){
	unlist(lapply(x, function(y){
		calculate.heterozygotes(y)
		}))
	})

het.data = lapply(replicates, function(x){
	lapply(x, function(y){
		unlist(lapply(y, function(z){
			calculate.heterozygotes(z)
			}))
		})
	})


het.data2 = lapply(1:5, function(x){
	lapply(het.data, function(y){
		y[[x]]
		})
	})

lapply(het.data2, function(x){
	g = as.data.frame(do.call(rbind, x))
	# browser()
	unlist(lapply(g, function(x) sum(x) / 10))
	})



lapply(het.data, function(x){

	})





loc200data = lapply(generations1, function(x){
	locdata1 = lapply(x, function(y){
		length(which(y[, 200] == "H"))
		})
	unlist(locdata1)
	})