# setwd("~/Google Drive/R/scripts/rotation1scripts/") #mac
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/") #windows

library(dunn.test)
library(foreign)
library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)

source("scripts/r/recombination_analysis_functions.R")
source("scripts/r/functions.R")

list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)

#load in mapping data
load.data = function(x) read.csv(paste("processed_data/genetic_maps/q", x, "_multipoint_ordercons_asmap_est_w_chr.csv", sep = ""),
                                                                 header = T, stringsAsFactors = F)

maps = lapply(1:4, load.data)

find.cm.total = function(LG, map){
    #args:
    #LG - integer linkage group number
    #map - genetic map file
    temp.map = filter(map, lg == paste("LG", LG, sep = ""))
    total.cm = max(as.numeric(temp.map$pos))
}

find.biggest.cm.gap = function(LG, map){
    #find biggest cM gaps between adjacent markers in map
    #args:
    #LG - integer linkage group number
    #map - genetic map file
    temp.map = filter(map, lg == paste("LG", LG, sep = ""))
    biggest.diff = max(diff(as.numeric(temp.map$pos), 1))
    return(biggest.diff)
}

total.cm = lapply(maps, function(x){
    lapply(list.of.full.lgs.to.test, find.cm.total, map = x)
})

biggest.gaps = lapply(maps, function(x){
    lapply(list.of.full.lgs.to.test, find.biggest.cm.gap, map = x)
})

total.cm.df = newdf(c("lg","quad","totalcm"))
total.cm.df = newdf(c("lg","quad","totalcm"))

calculate.mean.cm = function(quadno) mean(as.numeric(filter(total.cm.df, quad == quadno)$data))
lapply(1:4, calculate.mean.cm)

calculate.sd.cm = function(quadno) sd(as.numeric(filter(total.cm.df, quad == quadno)$data))
lapply(1:4, calculate.sd.cm)

    convert.l.o.l.to.df = function(list.of.lists){
    #convert a list of lists to dataframe
    wholedf = newdf(c("data", "quad"))
    for(i in 1:length(list.of.lists)){
        df = as.data.frame(t(as.data.frame(list.of.lists[[i]])))
        rownames(df) = NULL
        df$listno = i
        colnames(df) = c("data", "quad")
        wholedf = rbind(wholedf, df)
    }
    wholedf = wholedf[-1,]
    return(wholedf)
}

total.cm.df = convert.l.o.l.to.df(total.cm)
total.cm.df$lg = rep(list.of.full.lgs.to.test, 4)

biggest.gaps.df = convert.l.o.l.to.df(biggest.gaps)
biggest.gaps.df$lg = rep(list.of.full.lgs.to.test, 4)

comparison.files = list.files("processed_data/pairwiseLGcomparisons/q1_q2/")
comparison.files = comparison.files[-length(comparison.files)]
comparison.files = paste("processed_data/pairwiseLGcomparisons/q1_q2/", comparison.files, sep = "")
comparisons = lapply(comparison.files, read.csv, header = T, stringsAsFactors = F)

#grab physical distance percentages with updated function

#fix some erroneous markers
comparisons[[8]]$blast_percentage = genetictophysical("3A", filter(maps[[1]], lg == "LG3")$marker, markerformat = ".", 40)
comparisons[[3]]$blast_percentage = genetictophysical("3B", filter(maps[[1]], lg == "LG15")$marker, markerformat = ".", 40)
comparisons[[10]]$blast_percentage = genetictophysical("7B", filter(maps[[1]], lg == "LG7")$marker, markerformat = ".", 40)


#grab min and max blast percentages for all linkage groups
min.percentages = as.numeric(as.vector(lapply(comparisons, function(x) min(as.numeric(na.omit(x$blast_percentage))))))
max.percentages = as.numeric(as.vector(lapply(comparisons, function(x) max(as.numeric(na.omit(x$blast_percentage))))))
min.percentages = min.percentages[-6]
max.percentages = max.percentages[-6]

biggest.phys.gaps = lapply(comparisons, function(x) max(diff(na.omit(x$blast_percentage), 1)))




