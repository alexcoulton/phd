#Coulton et al. (2020) - Segregation distortion: utilizing simulated genotyping data to evaluate statistical methods
#This script implements selection for PedigreeSim (Voorrips & Maliepaard, 2012) for a single seed descent population structure of a biparental cross.

#The script was written in R version 3.6.1 and has been tested with the compiled version of PedigreeSim that was supplied in the supplementary data of Voorrips & Maliepaard (2012).

#Instructions:
#You must have PedigreeSim installed to utilize this script. Analyses are run with the run.1000.sims() function. See the function notes for details on the 
#arguments that can be supplied, these include the filial generation, the number of simulations, whether to apply selection or not, the selection strength,
#selection position, population size and distribution of markers in centimorgans. Results can be analysed using the perform.analysis() function. The script
#will store the results of simulations and analysis of simulations in the PedigreeSim installation folder.

#Ensure that this script, as well as the accompanying functions script (pedigree.sim.functions.R) are located within the PedigreeSim installation folder

#The genotyping data resulting from simulations are stored as a list of data.frame objects in the "simresults" folder. This can be loaded into R using the
#load() function, and will be opened in the envrionment as an object called all.sim.files

#Example usage:
#run.1000.sims(project.name = "wheat.selection.1", f.gen = 5, num.sims = 1000, w.selec = T, selec.str = 10, selec.pos = 100, pop.size = 300, cm.pos = seq(0, 100, 5), perform.analysis.flag = T)



####################################################

#CHANGE THIS TO YOUR PEDIGREESIM INSTALLATION PATH (remember to include a trailing forward slash)
pedigreesim.path = "/path/to/pedigreesim/installation/"

#SET NUMBER OF CPU CORES FOR ANALYSIS
number.of.cores = 50

####################################################
####################################################
library(parallel)
full.path = pedigreesim.path
setwd(full.path)

#load functions 
source("./pedigree.sim.functions.R")
####################################################



#### SAMPLE CODE ####

#run a simulation
run.1000.sims(project.name = "newselect", f.gen = 5, num.sims = 30, w.selec = T, selec.str = 2, selec.pos = 5, pop.size = 300, cm.pos = seq(0, 100, 5), perform.analysis.flag = T)

#load the results of the simulation into R for further analysis
load(paste0(pedigreesim.path, "simresults/newselect"))
newselect = all.sim.files

#check the genotype ratios at every locus in the first simulation
lapply(newselect[[1]], function(x){
	a = length(which(x == "A"))
	h = length(which(x == "H"))
	b = length(which(x == "B"))
	g = c(a, h, b)
	names(g) = c("A", "H", "B")
	g
	})