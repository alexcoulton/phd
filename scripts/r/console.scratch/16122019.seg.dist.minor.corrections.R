setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")

source("rotation1scripts_v4/scripts/r/modelling/pedigree.sim.functions.R")



#### 16/12/2019 MINOR CORRECTIONS ####

#performing additional simulations

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop96.pos100"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 96)    
})


lapply(seq(12, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop96.pos100"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 96)    
})


lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop300.pos100"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 300)    
})

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop1000.pos100"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 1000)    
})

lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop10000.pos100v2"), 2, num.sims = 1000, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 10000)    
})


lapply(seq(2, 20, 2), function(x){
    run.1000.sims(p("axcf2selec", x, "pop10000.pos100.500sims"), 2, num.sims = 500, w.selec = T, selec.str = x, selec.pos = 100, pop.size = 10000)    
})



# THIS ONE WORKED - dont need to replace any
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec20pop10000.pos100v2")
axcf2selec20pop10000.pos100v2 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec20pop10000.pos100.500sims")
axcf2selec20pop10000.pos100.500sims = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec18pop10000.pos100v2")
axcf2selec18pop10000.pos100v2 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec18pop10000.pos100.500sims")
axcf2selec18pop10000.pos100.500sims = all.sim.files

null.coord1 = which(sapply(axcf2selec18pop10000.pos100v2, class) == "NULL")
axcf2selec18pop10000.pos100v2 = axcf2selec18pop10000.pos100v2[-null.coord1]

which(sapply(axcf2selec18pop10000.pos100.500sims, class) == "NULL")

num.replace1 = 1000 - length(axcf2selec18pop10000.pos100v2)

axcf2selec18pop10000.pos100v3 = c(axcf2selec18pop10000.pos100v2, axcf2selec18pop10000.pos100.500sims[1:num.replace1])
length(axcf2selec18pop10000.pos100v3)

analysis1 = perform.analysis(axcf2selec18pop10000.pos100v3, selection = T, selection.position = 100, selection.strength = 18, recombination.position = "AxC1A", F2 = T, fgen = 2)
write.csv(analysis1, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/axcf2selec18pop10000.pos100v3.csv", row.names = F)


fix.10kpops = function(xv2, x500, selec.str, project.name1){
	null.coord1 = which(sapply(xv2, class) == "NULL")
	if(length(null.coord1) > 0) xv2 = xv2[-null.coord1]

	if(length(xv2) != 1000){
		num.replace1 = 1000 - length(xv2)

		xv3 = c(xv2, x500[1:num.replace1])	

		analysis1 = perform.analysis(xv3, selection = T, selection.position = 100, selection.strength = selec.str, recombination.position = "AxC1A", F2 = T, fgen = 2)
		write.csv(analysis1, paste0("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/", project.name1, ".csv"), row.names = F)	
	}	
}




load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec16pop10000.pos100v2")
axcf2selec16pop10000.pos100v2 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec16pop10000.pos100.500sims")
axcf2selec16pop10000.pos100.500sims = all.sim.files

fix.10kpops(axcf2selec16pop10000.pos100v2, axcf2selec16pop10000.pos100.500sims, 16, "axcf2selec16pop10000.pos100v3")

# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec14pop10000.pos100v2")
# axcf2selec14pop10000.pos100v2 = all.sim.files
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec14pop10000.pos100.500sims")
# axcf2selec14pop10000.pos100.500sims = all.sim.files

fix.10kpops(axcf2selec14pop10000.pos100v2, axcf2selec14pop10000.pos100.500sims, 14, "axcf2selec14pop10000.pos100v3")

# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec12pop10000.pos100v2")
# axcf2selec12pop10000.pos100v2 = all.sim.files
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec12pop10000.pos100.500sims")
# axcf2selec12pop10000.pos100.500sims = all.sim.files

fix.10kpops(axcf2selec12pop10000.pos100v2, axcf2selec12pop10000.pos100.500sims, 12, "axcf2selec12pop10000.pos100v3")


# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec10pop10000.pos100v2")
# axcf2selec10pop10000.pos100v2 = all.sim.files
# load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec10pop10000.pos100.500sims")
# axcf2selec10pop10000.pos100.500sims = all.sim.files

fix.10kpops(axcf2selec10pop10000.pos100v2, axcf2selec10pop10000.pos100.500sims, 10, "axcf2selec10pop10000.pos100v3")




pos100selec.files = list.files("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses", pattern = "pos100", full.names = T)

pop96 = pos100selec.files[grep("pop96", pos100selec.files)][-11]
pop300 = pos100selec.files[grep("pop300", pos100selec.files)]
pop1000 = pos100selec.files[grep("pop1000\\.", pos100selec.files)]
pop10000 = pos100selec.files[grep("pop10000.*v2|pop10000.*v3", pos100selec.files)]

pop96.2 = lapply(pop96, read.csv)
pop96.3 = bind_rows(pop96.2)
pop96.4 = arrange(pop96.3, selection.strength)

pop300.2 = lapply(pop300, read.csv)
pop300.3 = bind_rows(pop300.2)
pop300.4 = arrange(pop300.3, selection.strength)

pop1000.2 = lapply(pop1000, read.csv)
pop1000.3 = bind_rows(pop1000.2)
pop1000.4 = arrange(pop1000.3, selection.strength)

pop10000.2 = lapply(pop10000, read.csv)
pop10000.3 = bind_rows(pop10000.2)
pop10000.4 = arrange(pop10000.3, selection.strength)

all.pop1 = bind_rows(pop96.4, pop300.4, pop1000.4, pop10000.4)

write.csv(all.pop1, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/all.pos.100.anal.csv", row.names = F)





##### PERFORM SIMULATIONS W NO SELECTION (SKELETON MARKERS ONLY) #### 




run.1000.sims("axcf2pop96.binnedmarkers", 2, num.sims = 1000, pop.size = 96, cm.pos = binned.markers1a)    

run.1000.sims("axcf2pop300.binnedmarkers", 2, num.sims = 1000, pop.size = 300, cm.pos = binned.markers1a)    

run.1000.sims("axcf2pop1000.binnedmarkers", 2, num.sims = 1000, pop.size = 1000, cm.pos = binned.markers1a)    

run.1000.sims("axcf2pop10000.binnedmarkers", 2, num.sims = 1000, pop.size = 10000, cm.pos = binned.markers1a)    



#### SUBSET EXISTING SIMULATIONS ####




load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop96")
axcf2.noselec.pop96 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop300")
axcf2.noselec.pop300 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop1000")
axcf2.noselec.pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop10000")
axcf2.noselec.pop10000 = all.sim.files

bin.coords1 = match(binned.markers1a, axc.chr1a.cm)

axcf2.noselec.pop96.bins = mclapply(axcf2.noselec.pop96, function(x){
	x[, bin.coords1]
	}, mc.cores = 20)

axcf2.noselec.pop300.bins = mclapply(axcf2.noselec.pop300, function(x){
	x[, bin.coords1]
	}, mc.cores = 20)

axcf2.noselec.pop1000.bins = mclapply(axcf2.noselec.pop1000, function(x){
	x[, bin.coords1]
	}, mc.cores = 20)

axcf2.noselec.pop10000.bins = mclapply(axcf2.noselec.pop10000, function(x){
	x[, bin.coords1]
	}, mc.cores = 20)

analysis1 = perform.analysis(axcf2.noselec.pop96.bins, recombination.position = "AxC1A.Bins", F2 = T, fgen = 2)
analysis2 = perform.analysis(axcf2.noselec.pop300.bins, recombination.position = "AxC1A.Bins", F2 = T, fgen = 2)
analysis3 = perform.analysis(axcf2.noselec.pop1000.bins, recombination.position = "AxC1A.Bins", F2 = T, fgen = 2)
analysis4 = perform.analysis(axcf2.noselec.pop10000.bins, recombination.position = "AxC1A.Bins", F2 = T, fgen = 2)

write.csv(analysis1, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/axcf2.noselec.pop96.bins.csv", row.names = F)
write.csv(analysis2, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/axcf2.noselec.pop300.bins.csv", row.names = F)
write.csv(analysis3, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/axcf2.noselec.pop1000.bins.csv", row.names = F)
write.csv(analysis4, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/axcf2.noselec.pop10000.bins.csv", row.names = F)


#### DETECTION OF SDR (REVIEWER 2) ####

#load pos 200 files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec2pop1000")
axcf2selec2pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec4pop1000")
axcf2selec4pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec6pop1000")
axcf2selec6pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec8pop1000")
axcf2selec8pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec10pop1000")
axcf2selec10pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec12pop1000")
axcf2selec12pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec14pop1000")
axcf2selec14pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec16pop1000")
axcf2selec16pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec18pop1000")
axcf2selec18pop1000 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec20pop1000")
axcf2selec20pop1000 = all.sim.files



load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec2pop1000.pos100")
axcf2selec2pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec4pop1000.pos100")
axcf2selec4pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec6pop1000.pos100")
axcf2selec6pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec8pop1000.pos100")
axcf2selec8pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec10pop1000.pos100")
axcf2selec10pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec12pop1000.pos100")
axcf2selec12pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec14pop1000.pos100")
axcf2selec14pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec16pop1000.pos100")
axcf2selec16pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec18pop1000.pos100")
axcf2selec18pop1000.pos100 = all.sim.files
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2selec20pop1000.pos100")
axcf2selec20pop1000.pos100 = all.sim.files



allpos200 = list(axcf2selec2pop1000, axcf2selec4pop1000, axcf2selec6pop1000, axcf2selec8pop1000, axcf2selec10pop1000, axcf2selec12pop1000, axcf2selec14pop1000, axcf2selec16pop1000, axcf2selec18pop1000, axcf2selec20pop1000) 
allpos100 = list(axcf2selec2pop1000.pos100, axcf2selec4pop1000.pos100, axcf2selec6pop1000.pos100, axcf2selec8pop1000.pos100, axcf2selec10pop1000.pos100, axcf2selec12pop1000.pos100, axcf2selec14pop1000.pos100, axcf2selec16pop1000.pos100, axcf2selec18pop1000.pos100, axcf2selec20pop1000.pos100) 

find.sdr.new = function(x, use.fdr){
	if(use.fdr == T){
		g = check.for.distorted.markers.fdr(x, F2 = T)	
	} else {
		g = check.for.distorted.markers(x, thresholdlevel = 0.001, F2 = T)
	}
	
	diffs <- c(1, diff(g))
	start_indexes <- c(1, which(diffs > 1))
	end_indexes <- c(start_indexes - 1, length(g))
	coloned <- paste(g[start_indexes], g[end_indexes], sep=":")
	ranges1 = data.frame(coloned, stringsAsFactors = F)
	r1 = multi.str.split(ranges1[, 1], ":", 1)
	r2 = multi.str.split(ranges1[, 1], ":", 2)
	ranges1 = data.frame(r1, r2, stringsAsFactors = F)
	coords1 = which(ranges1[, 1] == ranges1[, 2])
	if(length(coords1) > 0) ranges1 = ranges1[-coords1, ]
	ranges1$r2 = as.numeric(ranges1$r2)
	ranges1$r1 = as.numeric(ranges1$r1)
	ranges1$sdr.length = (ranges1$r2 - ranges1$r1) + 1
	ranges1
}

library(dplyr)


all.sdr200 = lapply(allpos200, function(x){
		g = mclapply(x, function(y){
			find.sdr.new(y, use.fdr = T)
		}, mc.cores = 60)
	g2 = bind_rows(g)
	coord1 = which(is.na(g2$r1))
	if(length(coord1) > 0) g2 = g2[-coord1, ]
	list(g, g2)
	})

all.sdr100 = lapply(allpos100, function(x){
		g = mclapply(x, function(y){
			find.sdr.new(y, use.fdr = T)
		}, mc.cores = 60)
	g2 = bind_rows(g)
	coord1 = which(is.na(g2$r1))
	if(length(coord1) > 0) g2 = g2[-coord1, ]
	list(g, g2)
	})

all.sdr200.nofdr = lapply(allpos200, function(x){
	g = mclapply(x, function(y){
			find.sdr.new(y, use.fdr = F)
		}, mc.cores = 60)
	g2 = bind_rows(g)
	coord1 = which(is.na(g2$r1))
	if(length(coord1) > 0) g2 = g2[-coord1, ]
	list(g, g2)
	})

all.sdr100.nofdr = lapply(allpos100, function(x){
	g = mclapply(x, function(y){
			find.sdr.new(y, use.fdr = F)
		}, mc.cores = 60)
	g2 = bind_rows(g)
	coord1 = which(is.na(g2$r1))
	if(length(coord1) > 0) g2 = g2[-coord1, ]
	list(g, g2)
	})

s.s1 = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20)

all.sdr200.2 = Map(function(x, y){
	x[[2]]$selec.strength = y
	x[[2]]
	}, all.sdr200, s.s1)


all.sdr200.3 = bind_rows(all.sdr200.2)

all.sdr100.2 = Map(function(x, y){
	x[[2]]$selec.strength = y
	x[[2]]
	}, all.sdr100, s.s1)


all.sdr100.3 = bind_rows(all.sdr100.2)


lapply(all.sdr200, function(x){
	nrow(x[[2]])
	})

numsdr1 = as.data.frame(t(data.frame(lapply(all.sdr200, function(x){
	nrow(x[[2]])
	}))))



numsdr2 = as.data.frame(t(data.frame(lapply(all.sdr100, function(x){
	nrow(x[[2]])
	}))))


numsdr3 = as.data.frame(t(data.frame(lapply(all.sdr200.nofdr, function(x){
	nrow(x[[2]])
	}))))


numsdr4 = as.data.frame(t(data.frame(lapply(all.sdr100.nofdr, function(x){
	nrow(x[[2]])
	}))))


numsdr1$s.strength = s.s1
numsdr1$pos = 200
numsdr1$test = "FDR 0.05"
numsdr2$s.strength = s.s1
numsdr2$pos = 100
numsdr2$test = "FDR 0.05"
numsdr3$s.strength = s.s1
numsdr3$pos = 200
numsdr3$test = "0.001"
numsdr4$s.strength = s.s1
numsdr4$pos = 100
numsdr4$test = "0.001"



total.sdr = bind_rows(numsdr1, numsdr2, numsdr3, numsdr4)
colnames(total.sdr)[1] = "Number significant SDR"

num.sim.at.least1sdr.1 = sapply(all.sdr200, function(x){
	length(which(sapply(x[[1]], function(y){
		if(nrow(y) == 1 & is.na(y[1, 1])){
				return(F)
			} else {
				return(T)
			}
		})))
	})



num.sim.at.least1sdr.2 = sapply(all.sdr100, function(x){
	length(which(sapply(x[[1]], function(y){
		if(nrow(y) == 1 & is.na(y[1, 1])){
				return(F)
			} else {
				return(T)
			}
		})))
	})

num.sim.at.least1sdr.3 = sapply(all.sdr200.nofdr, function(x){
	length(which(sapply(x[[1]], function(y){
		if(nrow(y) == 1 & is.na(y[1, 1])){
				return(F)
			} else {
				return(T)
			}
		})))
	})


num.sim.at.least1sdr.4 = sapply(all.sdr100.nofdr, function(x){
	length(which(sapply(x[[1]], function(y){
		if(nrow(y) == 1 & is.na(y[1, 1])){
				return(F)
			} else {
				return(T)
			}
		})))
	})


num.sim.df1 = data.frame(num.sim.at.least1sdr.1)
num.sim.df1$pos = 200
num.sim.df1$test = "FDR 0.05"
num.sim.df1$selec.strength = s.s1
colnames(num.sim.df1)[1] = "num.sig"
num.sim.df2 = data.frame(num.sim.at.least1sdr.2)
num.sim.df2$pos = 100
num.sim.df2$test = "FDR 0.05"
num.sim.df2$selec.strength = s.s1
colnames(num.sim.df2)[1] = "num.sig"
num.sim.df3 = data.frame(num.sim.at.least1sdr.3)
num.sim.df3$pos = 200
num.sim.df3$test = "0.001"
num.sim.df3$selec.strength = s.s1
colnames(num.sim.df3)[1] = "num.sig"
num.sim.df4 = data.frame(num.sim.at.least1sdr.4)
num.sim.df4$pos = 100
num.sim.df4$test = "0.001"
num.sim.df4$selec.strength = s.s1
colnames(num.sim.df4)[1] = "num.sig"

all.at.least1.df = bind_rows(num.sim.df1, num.sim.df2, num.sim.df3, num.sim.df4)

write.csv(total.sdr, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/sdr_selecpos_1.csv", row.names = F)
write.csv(all.at.least1.df, "rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/sdr_selecpos_2.csv", row.names = F)