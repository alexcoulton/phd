




#what number of simulations have distorted markers?
length(which(unlist(lapply(list.of.first.gen.recomb.5b.225.300, function(x){
	g = check.for.distorted.markers(x)
	if(length(g) > 0) return(T)
	return(F)
	}))))



#in those simulations that did have distorted markers, how many markers were distorted?
q2 = unlist(lapply(list.of.first.gen.recomb.5b.225.300, function(x){
	g = check.for.distorted.markers(x)
	length(g)
	}))

mean(q2[which(q2 != 0)])
sd(q2[which(q2 != 0)])


list.of.selec = list(list.of.first.gen.w.prob.w.selec, list.of.first.gen.w.prob.w.selec2, list.of.first.gen.w.prob.w.selec3, list.of.first.gen.w.prob.w.selec4, list.of.first.gen.w.prob.w.selec5, list.of.first.gen.w.prob.w.selec6, list.of.first.gen.w.prob.w.selec7, list.of.first.gen.w.prob.w.selec8)

#grab number of simulations in which the peak of segregation distortion is within 10 markers of the point of selection
num.within.10 = lapply(list.of.selec.1000pop, function(x){
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		mean(which(g1 == max(g1)))
	}))

	length(which(g > 190 & g < 210))
})

#check number of simulations with distorted markers for list.of.selec
num.with.dist.markers = lapply(list.of.selec.1000pop, function(y){
	length(which(unlist(lapply(y, function(x){
	g = check.for.distorted.markers(x)
	if(length(g) > 0) return(T)
	return(F)
	}))))

	})

#grab the mean and standard deviations of number of distorted markers for each group of 1000 simulations
mean.sd.num.dist = lapply(list.of.selec.1000pop, function(y){
	q2 = unlist(lapply(y, function(x){
		g = check.for.distorted.markers(x)
		length(g)
	}))
	
	g1 = list(mean(q2[which(q2 != 0)]), sd(q2[which(q2 != 0)]))
	names(g1) = c("mean", "sd")
	g1
	})



plots3 = make.plots(list.of.first.gen.w.prob.w.selec)
for(i in 1:1000){
	png(p("rotation1scripts_v4/plots/simulation/selec/plot", i, ".png"), 1000, 1000)
	plots3[i]
	dev.off()
}



#what is the mean peak of distortion and what is the sd of the peak of distortion between simulations
peak.dist = lapply(list.of.selec.1000pop, function(x){
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		mean(which(g1 == max(g1)))
	}))

	g
})


list.of.non.selec = list(list.of.first.gen.recomb.anywherepop300, list.of.first.gen.recomb.anywherepop1000, list.of.first.gen.recomb.anywherepop10000, all.recomb.prob2, all.recomb.prob.pop10000, list.of.first.gen.recomb.5b.225.300)


peak.dist2 = lapply(list.of.non.selec, function(x){
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		mean(which(g1 == max(g1)))
	}))

	g
})


#what is the magnitude of distortion in the simulations at the peak
magnitude.dist1000pop = lapply(list.of.selec.1000pop, function(x){
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		max(g1)
	}))

	g
})

s(magnitude.dist1000pop, "rotation1scripts_v4/saved.objects/magnitude.dist1000pop")


#what is the magnitude of distortion in the simulations at the peak
magnitude.dist10000pop = lapply(list.of.selec.10000pop, function(x){
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		max(g1)
	}))

	g
})

s(magnitude.dist10000pop, "rotation1scripts_v4/saved.objects/magnitude.dist10000pop")



#what is the magnitude of distortion in the simulations at the peak
magnitude.dist2 = lapply(list.of.non.selec, function(x){
	g = unlist(lapply(x, function(y){
		g = convert.recomb.to.seg(y)
		g1 = abs(g - 0.5)
		max(g1)
	}))

	g
})


list.of.selec.1000pop = list(list.of.first.gen.w.prob.w.selec.1000pop, list.of.first.gen.w.prob.w.selec.1000pop2, list.of.first.gen.w.prob.w.selec.1000pop3, list.of.first.gen.w.prob.w.selec.1000pop4, list.of.first.gen.w.prob.w.selec.1000pop5, list.of.first.gen.w.prob.w.selec.1000pop6, list.of.first.gen.w.prob.w.selec.1000pop7, list.of.first.gen.w.prob.w.selec.1000pop8)


list.of.selec.10000pop = list(list.of.first.gen.w.prob.w.selec.10000pop, list.of.first.gen.w.prob.w.selec.10000pop2, list.of.first.gen.w.prob.w.selec.10000pop3, list.of.first.gen.w.prob.w.selec.10000pop4, list.of.first.gen.w.prob.w.selec.10000pop5, list.of.first.gen.w.prob.w.selec.10000pop6, list.of.first.gen.w.prob.w.selec.10000pop7, list.of.first.gen.w.prob.w.selec.10000pop8)


file1 = list.files("rotation1scripts_v4/saved.objects/recomb.sim/wselection/")
file2 = paste0("rotation1scripts_v4/saved.objects/recomb.sim/wselection/", file1)

for(i in file2){
	load(i)
}

