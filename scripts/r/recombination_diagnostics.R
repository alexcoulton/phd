# setwd("~/Google Drive/R/scripts/rotation1scripts/") #mac
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/") #windows

library(ggsignif)
library(dunn.test)
library(foreign)
library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)

source("scripts/r/recombination_analysis_functions.R")
source("scripts/r/functions.R")

#load in genotype data for analysis
load.data = function(x) read.csv(paste("processed_data/genotypes/alexq", x, "_man_cur_sk_map_genogeno.csv", sep = ""),
                                                                 header = T, stringsAsFactors = F)

initial.data = lapply(1:4, load.data)

genotypelists = lapply(initial.data, makegenotypelist)

########testing total events
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

test.of.total.num.events = makehistlist(7, genotypelists[[1]], T)

events = as.numeric(lapply(unique(test.of.total.num.events$individual), function(x){
    fil.events = filter(test.of.total.num.events, individual == x)
    num.events = nrow(fil.events)
    return(num.events)
}))

Mode(events)
hist(events)
################



get.central.lg.position = function(lg){
    filtered.genotype.df=initial.data[[1]][1:3, which(initial.data[[1]][1,]==lg)]
    central.position=ncol(filtered.genotype.df)/2
    return(central.position)
}

list.of.central=1:35
list.of.central=as.numeric(as.character(lapply(list.of.central, get.central.lg.position)))

add.dist.col = function(histlist, lg){
    histlist$dist.from.centre = abs(as.numeric(list.of.central[lg])-as.numeric(histlist$position))
    return(histlist)
}


#     ____________________________________________________________________________
#     CHECK RECOMBINATION FREQUENCIES FOR NORMAL DIST.                                             ####

#rewriting this section with Map() and expand.grid()
list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)
list.of.chromosome.assignments.for.lgs = c("2A", "3A", "5A", "7B", "6B", "1A", "3B", "7A", "2D", "4A")

grab.chromosome.assignment.for.lg = function(lg){
    index = which(list.of.full.lgs.to.test == lg)
    list.of.chromosome.assignments.for.lgs[index]
}

grab.list.of.markers = function(lg){
    g = initial.data[[1]][, which(initial.data[[1]][1,] == lg)]
    colnames(g)
}


all.combos = expand.grid(list.of.full.lgs.to.test, 1:4)

temp.histlists = Map(function(x, y) makehistlist(linkagegroup = x, genolist2 = genotypelists[[y]], return_recombinationcounts = T),
        all.combos$Var1, all.combos$Var2)

events.list = lapply(1:40, function(x){
    unique.individuals = unique(temp.histlists[[x]]$individual)
    as.numeric(lapply(unique.individuals, function(y){
        fil.events = filter(temp.histlists[[x]], individual == y)
        if(nrow(fil.events) == 1 & fil.events$is.recomb.event[1] == "N"){
            num.events = 0
        } else num.events = nrow(fil.events)
        return(num.events)
    }))
})

events.list

indexer = list(c(1:10), c(11:20), c(21:30), c(31:40))

events.list.dfs = lapply(indexer, function(x){
    as.data.frame(events.list[x])
})

events.list.dfs = lapply(events.list.dfs, function(x){
    colnames(x) = paste("lg_", list.of.full.lgs.to.test, sep = "")
    return(x)
})

#checking whether all individuals are included in temp.histlists
Map(function(x, y) all(unique(as.numeric(temp.histlists[[x]]$individual)) == 2:nrow(initial.data[[y]])), c(1, 11, 21, 31), 1:4)

convert.event.list.df.into.logical = function(df, threshold){
    df[] = lapply(df, function(x){
        x = x > threshold
    })
    return(df)
}

events.list.dfs.logical[] = lapply(events.list.dfs, convert.event.list.df.into.logical, threshold = 8)

grab.individuals.over.threshold = function(df){
    allcoords = expand.grid(1:nrow(df), 1:ncol(df))
    listoftruecoords = list()
    listoftruecoords = Map(function(x, y){
        if(df[x, y] == T){
            return(c(x, y))
        } else return("")
    }, allcoords$Var1, allcoords$Var2)
    listoftruecoords = listoftruecoords[-c(which(listoftruecoords == ""))]
    
    dftrue = as.data.frame(t(as.data.frame(listoftruecoords)))
    dftrue = add_column(dftrue, lg = "")
    colnames(dftrue) = c("individual", "lg_index", "lg")
    dftrue$individual = (dftrue$individual + 1) #add 1 to make tally up with rqtl format (has an additional row at the top of the df specifying chromosome / lg)
    #populate lg column with actual integers of LGs from indexes
    dftrue$lg = lapply(dftrue$lg_index, function(x) list.of.full.lgs.to.test[x])
    
    return(dftrue)
}

ind.over.threshold = lapply(events.list.dfs.logical, grab.individuals.over.threshold)

#conversion to character need for writing
ind.over.threshold = lapply(ind.over.threshold, function(x) data.frame(lapply(x, as.character)))

#write individuals to remove to files
for(i in 1:4){
    write.csv(ind.over.threshold[[i]], paste("processed_data/listofindividualsover.recomb.thres.q", i, ".csv", sep = ""),
                            row.names = F)
}

filter.init.data = function(lg, quad){
g = initial.data[[1]][c(1, (ind.over.threshold[[quad]]$individual[ind.over.threshold[[quad]]$lg == lg] + 1)),
                                            c(1,2, which(initial.data[[quad]][1,] == lg))]
return(g)
}

identify.possible.erroneous.recombination.events = function(lg, quad, threshold){
    g = filter.init.data(lg, quad)
    g1 = makegenotypelist(g)
    g2 = makehistlist(lg, g1, T)
    
    
    phys.pos = replacena(genetictophysical(chromosome = grab.chromosome.assignment.for.lg(lg), listofmarkers = colnames(g)[3:ncol(g)], 
                                                                                 markerformat = ".", threshold = 70))
    g2$phys = lapply(as.numeric(g2$position), function(x) phys.pos[x])
    
    g2$possible.error = ""
    g2$phys = as.numeric(g2$phys)
    
    for(i in 2:(nrow(g2) - 1)){
        if(g2$phys[i] != g2$phys[i-1] & g2$phys[i] != g2$phys[i+1] & abs(g2$phys[i] - g2$phys[i-1]) < threshold & abs(g2$phys[i] - g2$phys[i+1]) < threshold){
            g2$possible.error[i] = "Y"
        }
    } 
    return(g2)
}

lapply(list.of.full.lgs.to.test[1:7], function(x){
    print(x)
    length(unique(identify.possible.erroneous.recombination.events(x, 1, 8)$individual))
})


extract.lg.from.initial.data = function(lg, quad){
    g = initial.data[[quad]][, c(1,2, which(initial.data[[quad]][1,] == lg))]
    return(g)
}

extracted.lgs = Map(extract.lg.from.initial.data, list.of.full.lgs.to.test, c(rep(1, 10), rep(2, 10), rep(3, 10), rep(4, 10)))
extracted.lgs = list(extracted.lgs[1:10], extracted.lgs[11:20], extracted.lgs[21:30], extracted.lgs[31:40])

#remove individuals over threshold.
for(quad in 1:4){
    for(i in unique(ind.over.threshold[[quad]]$lg_index)){
        to.remove = which(rownames(extracted.lgs[[quad]][[i]]) %in% ind.over.threshold[[quad]][ind.over.threshold[[quad]]$lg_index == i,]$individual)
        extracted.lgs[[quad]][[i]] = extracted.lgs[[quad]][[i]][-to.remove,]
    }
}

#remove NA values
for(q in 1:4){
        for(i in 1:10){
            extracted.lgs[[q]][[i]][1,1] = ""
    }
}


#write files
for(q in 1:4){    
    for(i in 1:10){
        write.csv(extracted.lgs[[q]][[i]], file = paste("original_data/genotype.data/seperated.lgs.w.individuals.over.threshold.recombination.count.removed/q",
                                                                                         q, "lg", list.of.full.lgs.to.test[[i]], ".csv", sep = ""),
                            row.names = F)
    }
}
    




hist(as.numeric(lapply(events.list, Mode)))
as.numeric(lapply(events.list, mean))

lapply(events.list, function(x){x > 7})

lapply(lapply(events.list, function(x){x > 7}), which)



#output pdf histograms of mean recombination count in each individual for each lg for each quad
Map(function(x, lg, quad){
    pdf(file = paste("plots/recombination_in_individuals/lg", lg, "quad_", quad, "recombination_events_in_individuals.pdf", sep = ""))
    #windows(2,2)
    # par(mfrow=c(2,2))
    hist(x, breaks = 15, xlim = c(0, 25), ylim = c(0, 30))
    dev.off()
}, events.list, all.combos$Var1, all.combos$Var2)


#     ____________________________________________________________________________
#     BIMODAL CHECK                                                                                                                     ####

# solution to the following section is to use Map() instead of nested lapply()
#breaking the code up into smaller loops makes debugging easier

combinations.of.quad.and.lg = expand.grid(list.of.full.lgs.to.test, 1:4)

temp.histlist2 = Map(function(x, y) makehistlist(linkagegroup = x, genolist2 = genotypelists[[y]], return_recombinationcounts = T),
                                         combinations.of.quad.and.lg$Var1, combinations.of.quad.and.lg$Var2)

temp.histlist2 = Map(function(x, y) add.dist.col(x, y), temp.histlist2, combinations.of.quad.and.lg$Var1)

temp.histlist.means = lapply(temp.histlist2, function(x){
    list.of.means.from.centre = vector()
    for(i in unique(x$individual)){
        fil.data = filter(x, individual == i)
        mean.dist.from.centre = mean(fil.data$dist.from.centre)
        list.of.means.from.centre = c(list.of.means.from.centre, mean.dist.from.centre)
    }
    return(list.of.means.from.centre)
})


hist_seq = function(x) as.numeric(rownames(combinations.of.quad.and.lg[combinations.of.quad.and.lg$Var1 == x,]))

hist.indexes = lapply(list.of.full.lgs.to.test, hist_seq)

for(i in 1:length(hist.indexes)){
    pdf(file = paste("plots/bimodal_dist_within_individuals_test/bimodal_test_lg", list.of.full.lgs.to.test[i], ".pdf", sep = ""))
    par(mfrow=c(2,2))
    lapply(temp.histlist.means[hist.indexes[[i]]], hist)
    dev.off()
}





