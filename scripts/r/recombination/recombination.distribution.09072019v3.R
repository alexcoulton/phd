#### INITIALIZE LIBRARIES / FUNCTIONS #### 

library(RMySQL)
library(readxl)
library(ape)
library(grid)
library(viridis)
library(rmarkdown)
library(Biostrings)
library(captioner)
library(gridExtra)
library(readr)
library(ASMap)
library(bookdown)
library(cowplot)
library(data.table)
library(BatchMap)
#library(SNPolisher)
library(zoo) #for na.approx()
library(dplyr)
library(sp)
library(parallel)
library(ggplot2)
library(pander)
library(kableExtra)
library(ggrepel)
library(data.table)

#setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R")

#### centromeres ####

#taken from S11 of supplementary data from 2018 IWGSC paper in Science. Where more than one region has been defined,
#the region that overlaps the highest centromeric TE family (either Cereba or Quinta; figure S7) has been taken.
centromeres = read.table("rotation1scripts_v4/original_data/IWGSC/centromere.positions.txt", sep = " ", stringsAsFactors = F)
centromere.positions = sapply(centromeres$V2, function(x){
    g = strsplit(x, "-")
    g = as.numeric(g[[1]])
    mean(g)
})

centromere.positions.bp = centromere.positions * 1000000

centromeres$pos = centromere.positions
centromeres$pos.bp = centromere.positions.bp
colnames(centromeres)[1] = "chr"

#### chromosome counts ####

chromosomecounts <<- read.table("rotation1scripts_v4/original_data/iwgsc.chromosome.counts.txt", stringsAsFactors = F)
chromosomecounts2 = chromosomecounts
chromosomecounts2$V2 = multi.str.split(chromosomecounts2$V2, "chr", 2)



chromosomecounts.v2 = chromosomecounts[-22, ]
chromosomecounts.v2$centro = centromeres$pos.bp
chromosomecounts.v2$centro.mb = chromosomecounts.v2$centro / 1000000
chromosomecounts.v2$centro.per = (chromosomecounts.v2$centro / chromosomecounts.v2$V1) * 100
# phys.midpoint = chromosomecounts.v2[which(chromosomecounts.v2$V2 == paste0("chr", y)), ]$centro.per


#### GENE DISTRIBUTION / RECOMBINATION ANALYSIS #### 

gff1 = read.table("~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", sep = "\t", stringsAsFactors = F, header = T)
gff1 = gff1[, -1]
gff1 = reset.colnames(gff1)
gff1$V4 = gff1$V4 / 1000000

gff1$V1 = multi.str.split(gff1$V1, "chr", 2)
colnames(gff1)[1] = "chr"
gff1 = gff1[-which(gff1$chr == "Un"), ]
gff2 = split(gff1, gff1$chr)



gene.dist.hists = Map(function(x, y){
    
    g = make.ggplot.histogram(x$V4, num.bins = 1000, breaks = seq(1, max(x$V4), (max(x$V4) / 1000)), xlabel = "Physical Position", ylabel = "Frequency", plot.title = "Histogram")
        
    g + geom_vline(xintercept = y, color = "red", size = 1) + ggtitle(paste0("Gene distribution ", x$chr[[1]])) + theme_bw(base_size = 22)
    
}, gff2[1:21], centromeres$pos)

do.call(grid.arrange, gene.dist.hists)




#### PERFORM AXC ANALYSIS ####

#analyses MRD between replicate populations (distribution of recombination) as well as recombination frequency
source("rotation1scripts_v4/scripts/r/recombination/axc.analysis.R")

#### PERFORM AXP ANALYSIS #### 

#analyses MRD between temperature treatments (distribution of recombination) as well as recombination frequency
source("rotation1scripts_v4/scripts/r/recombination/axp.analysis.R")

#saving workspace
#save.image("~/project.phd.main/rotation1scripts_v4/scripts/r/recombination/workspace20052020.RData")





#### FURTHER RECOMBINATION ANALYSIS ####

source("rotation1scripts_v4/scripts/r/recombination/further.recombination.analysis.R")


#### PERFORM SIMULATION ANALYSIS ####

#source("rotation1scripts_v4/scripts/r/recombination/simulation.analysis.R")


##### PERFORM QTL ANALYSIS ####

source("rotation1scripts_v4/scripts/r/recombination/recomb.qtl.analysis.R")


#complete workspace for rmarkdown file rendering
save.image("~/project.phd.main/rotation1scripts_v4/scripts/r/recombination/workspace22052020.RData")


