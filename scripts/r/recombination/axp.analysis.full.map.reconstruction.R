

library(rmarkdown)
library(Biostrings)
#library(captioner)
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

#setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R")

load.data = function(x) read.csv(paste("rotation1scripts_v4/processed_data/genotypes/alexq", x, "_man_cur_sk_map_genogeno.csv", sep = ""), 
                                                                 header = T, stringsAsFactors = F)


#put in order of temperature, q1 = 10(deg)C, q3 = 14(deg)C, q2 = 26(deg)C, q4 = 28(deg)C
initial.data = lapply(c(1, 3, 2, 4), load.data)

list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)
list.of.chromosome.assignments.for.lgs = c("2A", "3A", "5A", "7B", "6B", "1A", "3B", "7A", "2D", "4A")


#### make full map from skeleton ####

library(dplyr)

comb.data1 = do.call(bind_rows, initial.data)
comb.data1 = comb.data1[-c(87, 88, 89, 168, 169, 170, 244, 245, 246), ]

chromo.assigns1 = read_csv("rotation1scripts_v4/processed_data/genetic_maps/q1_multipoint_ordercons_asmap_est_w_chr.csv")
# comb.data1[1, 3:ncol(comb.data1)] = chromo.assigns1$chr_truncated



all.markers.axpf2 = read.csv("rotation1scripts_v4/original_data/genotype.data/manuallycuratedsnps_flipped.csv", stringsAsFactors = F)

al3 = all.markers.axpf2[, -(match(colnames(comb.data1)[3:ncol(comb.data1)], colnames(all.markers.axpf2)[3:ncol(all.markers.axpf2)]) + 2)]
al3 = al3[match(comb.data1$probeset_id, al3$probeset_id), ]

count1 = 1
marker.matches1 = lapply(al3, function(x){
    g = sapply(comb.data1, function(y){
        x2 = x
        x2 = x2[-1]
        y = y[-1]
        
        hyphen.coord = which(y == "-")
        if(length(hyphen.coord) > 0) x2[hyphen.coord] = "-"
        
        hyphen.coord2 = which(x2 == "-")
        if(length(hyphen.coord2) > 0) y[hyphen.coord2] = "-"
        
        all(x2 == y)
        
        
    })
    
    # which(g)
    search.name = colnames(al3)[count1]
    count1 <<- count1 + 1
    g2 = list(search.name, which(g))
    names(g2) = c("al3.marker", "comb.data1.match")
    g2
    
})

marker.matches1 = marker.matches1[-c(1:2)]

comb.data2 = as.data.frame(t(comb.data1))
comb.data2 = convert.to.character.data.frame(comb.data2)
comb.data3 = add.column.at.position(comb.data2, 0)
comb.data3$empty.col = rownames(comb.data2)
#splitting dataframe so coordinates (in marker.matches1) are not affected by insertion of new markers
comb.data4 = split.df.into.list.of.dfs.by.column(comb.data3, "empty.col", F)

#insert the markers
count1 = 1
lapply(marker.matches1, function(x){
    
    # if(length(x$comb.data1.match) > 1) browser()
    
    if(length(x$comb.data1.match) > 0){
        
        to.insert = al3[, which(colnames(al3) == x$al3.marker)]
        matching_0 = comb.data4[[max(x$comb.data1.match)]]
        matching = unname(unlist(comb.data4[[max(x$comb.data1.match)]][1, ]))
        to.insert = c(x$al3.marker, to.insert)
        to.insert[2] = matching[2]
        to.insert = as.data.frame(t(data.frame(to.insert)))
        colnames(to.insert) = colnames(matching_0)
        
        marker.insertion = bind_rows(matching_0, to.insert)
        comb.data4[[max(x$comb.data1.match)]] <<- marker.insertion
        
    }
    
    1
})

#bind everything back together and transpose into appropriate format
comb.data5 = bind_rows(comb.data4)
comb.data6 = convert.to.character.data.frame(as.data.frame(t(comb.data5)))

colnames(comb.data6) = comb.data6[1, ]
comb.data6 = comb.data6[-1, ]

#order chromosomes
# comb.data7 = put.geno.in.order.of.chromo(comb.data6)
comb.data7.map = extract.centimorgan.from.genotypes(comb.data6)

#order markers within bins
comb.data7.map$chromo = as.numeric(comb.data7.map$chromo)
comb.data7.map2 = split(comb.data7.map, comb.data7.map$chromo)
comb.data7.map3 = lapply(comb.data7.map2, function(x){
    g = split(x, x$cm)
    g = lapply(g, function(y){
        y[sort(y$marker, index.return = T)$ix, ]
    })
    bind_rows(g)
})

comb.data7.map4 = bind_rows(comb.data7.map3)


#assign chromosomes 

comb.data7.map4$chr1 = chromo.assigns1[match(comb.data7.map4$marker, chromo.assigns1$marker), ]$chr_truncated

c1 = split(comb.data7.map4, comb.data7.map4$chromo)
c2 = lapply(c1, function(x) unique(na.omit(x$chr1)))
c3 = Map(function(x, y){
    x$chr1 = y
    x
}, c1, c2)

c4 = bind_rows(c3)


c5 = c4[which(c4$chromo %in% list.of.full.lgs.to.test), ]
c5$marker = switch.affy.format(c5$marker)
#correct an incorrectly labelled LG
c5[which(c5$chromo == 15), ]$chr1 = "3B"


comb.data6.reordered = comb.data6[, c(1, 2, which(colnames(comb.data6) %in% switch.affy.format(c5$marker)))]
comb.data6.reordered = comb.data6.reordered[, c(1, 2, match(switch.affy.format(c5$marker), colnames(comb.data6.reordered)))]

comb.data6.reordered[1, 3:ncol(comb.data6.reordered)] = c5$chr1

comb.phys1 = get.phys.for.chromos(comb.data6.reordered, F, T, ".")

c5$phys = comb.phys1$phys.pos

c5 = split(c5, c5$chr1)

inverted.lgs = unlist(lapply(c5, function(x){
    # print(na.omit(x$phys.dist.bp))
    check.if.inverted(na.omit(x$phys))
    
}))

for(i in which(inverted.lgs)){
    c5[[i]] = reverse.genetic.map(c5[[i]], "cm")
}

#order by physical positon within bins
c5 = lapply(c5, function(x){
    g = split(x, x$cm)
    g = lapply(g, function(y){
        y[sort(y$phys, index.return = T, na.last = T)$ix, ]
    })
    bind_rows(g)
})




#### 10022020 NEW ORDER BY PHYS WITHIN BINS ####
c5 = bind_rows(c5)

comb.data6.reordered = comb.data6.reordered[, c(1, 2, match(switch.affy.format(c5$marker), colnames(comb.data6.reordered)))]
cc.new = comb.data6.reordered
q1 = cc.new[c(1, grep("Q1", cc.new$probeset_id)), ]
q2 = cc.new[c(1, grep("Q2", cc.new$probeset_id)), ]
q3 = cc.new[c(1, grep("Q3", cc.new$probeset_id)), ]
q4 = cc.new[c(1, grep("Q4", cc.new$probeset_id)), ]

initial.data = list(q1, q3, q2, q4)
initial.data2 = list(q1, q3, q2, q4)
