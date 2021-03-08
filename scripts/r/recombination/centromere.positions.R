setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

gff1 = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.genes.onlyv2.gff3", sep = "\t", stringsAsFactors = F)
gff2 = split(gff1, gff1$V1)

centromeres = read.table("rotation1scripts_v4/original_data/IWGSC/centromere.positions.txt", sep = " ", stringsAsFactors = F)
#taken from S11 of supplementary data from 2018 IWGSC paper in Science. Where more than one region has been defined,
#the region that overlaps the highest centromeric TE family (either Cereba or Quinta; figure S7) has been taken.
centromere.positions = sapply(centromeres$V2, function(x){
    g = strsplit(x, "-")
    g = as.numeric(g[[1]])
    mean(g)
})

centromere.positions.bp = centromere.positions * 1000000

centromeres$pos = centromere.positions
centromeres$pos.bp = centromere.positions.bp

hist(gff2[[1]]$V4, breaks = 1000)
abline(v = centromeres$pos.bp[[1]])


Map(function(x, y){
    hist(x$V4, breaks = 1000)
    abline(v = y)
}, gff2, centromeres$pos.bp)



