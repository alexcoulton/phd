#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                #####     

find.hotspot.and.coldspot.ranges = function(dataframe1, cm.density.threshold){
    #args:
    #cm.density.threshold: only return hotspots that are above this value of cm.density
    if(missing(cm.density.threshold)) cm.density.threshold = 0
    
    hotspotlengths = as.numeric()
    coldspotlengths    = as.numeric()
    hotspot.df = newdf(c("cm.start", "cm.end", "phys.pos.start", "phys.pos.end", "cm.density"), no.rows = T)
    coldspot.df = newdf(c("cm.start", "cm.end", "phys.pos.start", "phys.pos.end", "cm.density"), no.rows = T)
    for(i in 1:nunique(dataframe1$cM)){
        bincoord = which(dataframe1$cM == unique(dataframe1$cM)[i])
        first.row.in.block = dataframe1[bincoord, ][1, ]
        last.row.in.prev.block = dataframe1[(bincoord - 1), ][1, ]
        
        if(i != 1){
            
            hotspot.df = add_row(hotspot.df)
            
            hotspot.df[nrow(hotspot.df), ]$cm.start = last.row.in.prev.block$cM
            hotspot.df[nrow(hotspot.df), ]$cm.end = first.row.in.block$cM
            hotspot.df[nrow(hotspot.df), ]$phys.pos.start = last.row.in.prev.block$phys.dist.bp
            hotspot.df[nrow(hotspot.df), ]$phys.pos.end = first.row.in.block$phys.dist.bp
        }

        
        first.row.phys = first.row.in.block$phys.dist.bp
        last.row.phys = last.row.in.prev.block$phys.dist.bp
        
        # browser()
        
        hotspot.len = first.row.phys - last.row.phys
        hotspotlengths = c(hotspotlengths, hotspot.len)
        
        names(hotspotlengths) = c(names(hotspotlengths)[-length(hotspotlengths)], unique(dataframe1$cM)[i])
        
        if(length(bincoord) > 1){
            coldspot.df = add_row(coldspot.df)
            
            fir = dataframe1[bincoord, ][1, ]
            las = dataframe1[bincoord, ][length(bincoord), ]

            coldspot.df[nrow(coldspot.df), ]$cm.start = fir$cM
            coldspot.df[nrow(coldspot.df), ]$cm.end = las$cM
            coldspot.df[nrow(coldspot.df), ]$phys.pos.start = fir$phys.dist.bp
            coldspot.df[nrow(coldspot.df), ]$phys.pos.end = las$phys.dist.bp
            
            
            
            coldspot.len = las$phys.dist.bp - fir$phys.dist.bp
            coldspotlengths = c(coldspotlengths, coldspot.len)
            
            names(coldspotlengths) = c(names(coldspotlengths)[-length(coldspotlengths)], unique(dataframe1$cM)[i])
        }
        
        # browser()
    }
    
    coldspot.df[] = lapply(coldspot.df, as.numeric)
    hotspot.df[] = lapply(hotspot.df, as.numeric)
    
    hotspot.df$phys.length = as.numeric(hotspot.df$phys.pos.end) - as.numeric(hotspot.df$phys.pos.start)
    hotspot.df$cm.length = as.numeric(hotspot.df$cm.end) - as.numeric(hotspot.df$cm.start)
    coldspot.df$phys.length = as.numeric(coldspot.df$phys.pos.end) - as.numeric(coldspot.df$phys.pos.start)
    coldspot.df$cm.length = as.numeric(coldspot.df$cm.end) - as.numeric(coldspot.df$cm.start)
    hotspot.df$cm.density = (hotspot.df$cm.end - hotspot.df$cm.start) / (hotspot.df$phys.length / 1000000)
    coldspot.df$cm.density = (coldspot.df$cm.end - coldspot.df$cm.start) / (coldspot.df$phys.length / 1000000)
    
    if(cm.density.threshold != 0) hotspot.df = hotspot.df[which(hotspot.df$cm.density > cm.density.threshold), ]
    
    if(any(hotspot.df$phys.length < 0)){
        print(p("Removed ", length(which(hotspot.df$phys.length < 0)), " neg. hotspots"))
        hotspot.df = hotspot.df[-which(hotspot.df$phys.length < 0), ]
    }
    
    if(any(coldspot.df$phys.length < 0)){
        print(p("Removed ", length(which(coldspot.df$phys.length < 0)), " neg. coldspots"))
        coldspot.df = coldspot.df[-which(coldspot.df$phys.length < 0), ]
    }
    
    g1 = list(hotspot.df, coldspot.df)
    names(g1) = c("hotspot.df", "coldspot.df")
    return(g1)
}

library(IRanges)
find.introgression.recombination.overlaps = function(intro.df1, recombination.df1){
    
    intro.ranges = IRanges(intro.df1$intro.start, intro.df1$intro.end)
    
    if(any(recombination.df1$phys.length < 0)){
        recombination.df1 = recombination.df1[-which(recombination.df1$phys.length < 0), ]    
    }
    
    recombination.ranges = IRanges(recombination.df1$phys.pos.start, recombination.df1$phys.pos.end)
    
    overlaps = findOverlaps(intro.ranges, recombination.ranges)
    overlap.lengths = ranges(overlaps, intro.ranges, recombination.ranges)

    all.overlaps = cbind(as.data.frame(overlaps), as.data.frame(overlap.lengths))
    
    
    #add additional fields from introgression data based on subjectHits
    if(nrow(all.overlaps) == 0){
        all.overlaps$cm.density = as.numeric()
        all.overlaps$recomb.start = as.numeric()
        all.overlaps$recomb.end = as.numeric()
        all.overlaps$wild.rel = as.character()
        all.overlaps$samp.type = as.character()
        all.overlaps$phy.dist = as.numeric()
        all.overlaps$introgression.id = as.numeric()
        all.overlaps$snp.density = as.numeric()
    } else {
        all.overlaps$cm.density = ""
        all.overlaps$cm.density = as.numeric(all.overlaps$cm.density)    
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$cm.density[i] = recombination.df1$cm.density[all.overlaps$subjectHits[i]]
        }
        
        all.overlaps$recomb.start = ""
        all.overlaps$recomb.start = as.numeric(all.overlaps$recomb.start)
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$recomb.start[i] = recombination.df1$phys.pos.start[all.overlaps$subjectHits[i]]
        }
        
        all.overlaps$recomb.end = ""
        all.overlaps$recomb.end = as.numeric(all.overlaps$recomb.end)
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$recomb.end[i] = recombination.df1$phys.pos.end[all.overlaps$subjectHits[i]]
        }
        
        all.overlaps$wild.rel = ""
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$wild.rel[i] = intro.df1$wild.relative[all.overlaps$queryHits[i]]
        }
        
        all.overlaps$samp.type = ""
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$samp.type[i] = intro.df1$samp.type[all.overlaps$queryHits[i]]
        }
        
        all.overlaps$phy.dist = ""
        all.overlaps$phy.dist = as.numeric(all.overlaps$phy.dist)
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$phy.dist[i] = intro.df1$phy.dist.from.avalon[all.overlaps$queryHits[i]]
        }
        
        all.overlaps$introgression.id = ""
        all.overlaps$introgression.id = as.numeric(all.overlaps$introgression.id)
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$introgression.id[i] = intro.df1$X1[all.overlaps$queryHits[i]]
        }
        
        all.overlaps$snp.density = ""
        all.overlaps$snp.density = as.numeric(all.overlaps$snp.density)
        for(i in 1:nrow(all.overlaps)){
            all.overlaps$snp.density[i] = intro.df1$snp.density[all.overlaps$queryHits[i]]
        }
        
    }
    
    
    
    
    if(nrow(all.overlaps) > 0){
        all.overlaps$per.overlap = ""
        
        for(i in 1:nrow(all.overlaps)){
            int1 = as.data.frame(intro.ranges[all.overlaps$queryHits[i]])
            rec1 = as.data.frame(recombination.ranges[all.overlaps$subjectHits[i]])
            
            #find how much of the recombination interval (hot or cold) the introgression covers as a
            #percentage
            recombination.width = as.numeric(rec1$width)
            overlap.percentage = (all.overlaps[i, ]$width / recombination.width) * 100
            
            
            
            all.overlaps$per.overlap[i] = overlap.percentage
            # print(all.overlaps)
        }
        
    }
        
    
    
    return(all.overlaps)
}

total.num.hotspots = 0
total.num.unique.subjectHits = 0

find.all.overlaps = function(gen.map, introgression.df, cm.den.threshold){
    count1 = make.counter()
    if(missing(cm.den.threshold)) cm.den.threshold = 0
    g = Map(function(genetic.map, chromo1){
        print(count1())
        recomb.ranges = find.hotspot.and.coldspot.ranges(genetic.map, cm.den.threshold)
        #init
        overlaps.hot = ""
        overlaps.cold = ""
        #check if there are any introgressions this particular chromosome
        if(is.null(introgression.df[[chromo1]]) == F){
            
            overlaps.hot = find.introgression.recombination.overlaps(introgression.df[[chromo1]], recomb.ranges$hotspot.df)    
            overlaps.cold = find.introgression.recombination.overlaps(introgression.df[[chromo1]], recomb.ranges$coldspot.df) 
            overlaps.hot$per.overlap = as.numeric(overlaps.hot$per.overlap)
            overlaps.cold$per.overlap = as.numeric(overlaps.cold$per.overlap)
            
            overlaps.hot = filter(overlaps.hot, per.overlap > 60)
            overlaps.cold = filter(overlaps.cold, per.overlap > 60)
            
        }
        both.overlaps = list(overlaps.hot, overlaps.cold)
        names(both.overlaps) = c("hot", "cold")
        
        # total.num.hotspots <<- total.num.hotspots + nrow(recomb.ranges$hotspot.df)
        # total.num.unique.subjectHits <<- total.num.unique.subjectHits + length(unique(overlaps.hot$subjectHits))
        # 
        # print(p("recomb.ranges ", nrow(recomb.ranges$hotspot.df), ", ", length(unique(overlaps.hot$subjectHits))))
        # browser()
        return(both.overlaps)
        
    }, gen.map, names(gen.map))
    
    rm.coord = which(unlist(lapply(g, function(x) all(x == ""))))
    if(length(rm.coord) != 0) g = g[-rm.coord]
    g
    
} 

#PLOT - distribution of recombination for a particular chromosome in s.x.r
# ggplot(all.m.8.iwgsc.4b.rev[[2]][[7]], aes(x = phys.dist, y = cm.diff)) + geom_point() + geom_line()

find.shared.overlaps.closure = function(hot.or.cold2){
    find.shared.overlaps = function(all.over.df1, all.over.df2){
        qq1 = function(df1, hot.or.cold) lapply(df1, function(x) x[[hot.or.cold]])
        
        hot1 = qq1(all.over.df1, hot.or.cold2)
        hot2 = qq1(all.over.df2, hot.or.cold2)
        
        hot1 = hot1[-which(unlist(lapply(hot1, nrow)) == 0)]
        hot2 = hot2[-which(unlist(lapply(hot2, nrow)) == 0)]
        
        
        shared.chromo = which(names(hot1) %in% names(hot2))
        if(length(shared.chromo) > 0) {
            hot1 = hot1[shared.chromo]    
            # browser()
            lapply(names(hot1), function(x){
                hot1[[x]]$start[which(hot1[[x]]$start %in% hot2[[x]]$start)]
            })
            
            
        } else {
            return("No shared overlaps")
        }
        
    }
    
}

find.cold.overlaps = find.shared.overlaps.closure("cold")
find.hot.overlaps = find.shared.overlaps.closure("hot")

count.all.overlaps = function(all.overlaps.ob){
    g = lapply(all.overlaps.ob, function(x){
        hotspot = as.data.frame(x[[1]])
        coldspot = as.data.frame(x[[2]])
        h.h = nunique(hotspot$subjectHits)
        c.h = nunique(coldspot$subjectHits)
        list(h.h, c.h)
        # browser()
    })
    g1 = unlist(g)
    
    no.hotspot.overlaps = sum(g1[grep("1$", names(unlist(g)))])
    no.coldspot.overlaps = sum(g1[grep("2$", names(unlist(g)))])
    q = list(no.hotspot.overlaps, no.coldspot.overlaps)
    names(q) = c("hot", "cold")
    
    return(q)
}

grab.mean.phylogenetic.dist = function(all.overlaps.ob){
    all.hot = unlist(lapply(all.overlaps.ob, function(x){
        x[[1]]$phy.dist
    }))
    
    all.cold = unlist(lapply(all.overlaps.ob, function(x){
        x[[2]]$phy.dist
    }))
    
    g = list(all.hot, all.cold)
    g
}

remove.redundant.overlaps = function(all.overlaps.obj){
    g = list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))
    
    for(i in g){
        Map(function(chromosome.name){
            hot.redundant = which(all.overlaps.obj[[i[1]]][[chromosome.name]][["hot"]]$subjectHits %in% all.overlaps.obj[[i[2]]][[chromosome.name]][["hot"]]$subjectHits)
            all.overlaps.obj[[i[1]]][[chromosome.name]][["hot"]] <<- all.overlaps.obj[[i[1]]][[chromosome.name]][["hot"]][-hot.redundant, ]
            
            cold.redundant = which(all.overlaps.obj[[i[1]]][[chromosome.name]][["cold"]]$subjectHits %in% all.overlaps.obj[[i[2]]][[chromosome.name]][["cold"]]$subjectHits)
            all.overlaps.obj[[i[1]]][[chromosome.name]][["cold"]] <<- all.overlaps.obj[[i[1]]][[chromosome.name]][["cold"]][-cold.redundant, ]
            # browser()
        }, names(all.overlaps.obj[[i[1]]]))        
    }
    
    all.overlaps.obj
}

#     ____________________________________________________________________________
#     ANALYSIS OF OVERLAP (INTROGRESSIONS & HOT/COLDSPOTS)                                        ####

#GENERATE LIST OF OVERLAPS

cm.den = 0

qqq3 = c(1, 1, 2, 2, 3, 3, 4, 4)

#finds recombination and introgression overlaps given a genetic map and introgression dataframe
fo2 = function(map.no, int.df) find.all.overlaps(all.m.8.iwgsc.4b.rev[[map.no]], int.df, cm.den)

all.overlaps.list = Map(function(int1df, mapno){
    fo2(mapno, int1df)
}, all.intro4[[1]], qqq3)

all.overlaps.list2 = Map(function(int1df, mapno){
    fo2(mapno, int1df)
}, all.intro4[[2]], qqq3)

all.overlaps.list2.rm.redun = remove.redundant.overlaps(all.overlaps.list2)

#checking whether the introgressions that overlap with recombination areas are low in snp.density
pdf("rotation1scripts_v4/plots/snp.density.comparison/all.intro.vs.overlapping.intro.pdf", 7, 7)
par(mfrow = c(2, 1))
hist((all.intro4.comb[[2]]$snp.density*1000000), breaks = 1000, xlim = c(0, 50), main = "All Introgressions", xlab = "SNP Density (# SNPs / Mb)")
hist((best.hot4$snp.density*1000000), breaks = 1000, xlim = c(0, 50), main = "Overlapping Introgressions", xlab = "SNP Density (# SNPs / Mb)")
dev.off()

all.overlaps.list3 = Map(function(int1df, mapno){
    fo2(mapno, int1df)
}, all.intro4[[3]], qqq3)

overlap.phy.dist = lapply(all.overlaps.list2, function(x){
    lapply(x, function(q){
        g = lapply(q, function(q1){
            mean(q1$phy.dist)
        })
        g[[1]]/g[[2]]
        
    })
})

length(which(unlist(overlap.phy.dist) < 1))
length(which(unlist(overlap.phy.dist) > 1))

overlap.counts1 = lapply(1:8, function(x){
    count.all.overlaps(all.overlaps.list2.rm.redun[[x]])
})

#add up overlap counts from both parents in each respective cross
t.list = list(c(1, 2), c(3, 4), c(5, 6), c(7, 8))
num.overlaps2 = lapply(t.list, function(x){
    hot = overlap.counts1[[x[1]]][[1]] + overlap.counts1[[x[2]]][[1]]
    cold = overlap.counts1[[x[1]]][[2]] + overlap.counts1[[x[2]]][[2]]
    g = list(hot, cold)
    names(g) = c("hot", "cold")
    g
})


lapply(all.overlaps.list2, function(q){
    f1 = do.call(rbind, lapply(q, function(x){
        x[[1]]
    }))    
})




#     ____________________________________________________________________________
#     GERMPLASM INVESTIGATION                                                                                                 ####

germplasm = read_excel("rotation1scripts_v4/original_data/sacha.introgression.data/germplasm.list.xlsx")

germplasm$`Sample type`[which(germplasm$`Sample type` == "Progenitor")] = "progenitor"

#shows the number of overlaps for samples classed as progenitors and relatives respectively
overlap.by.sample.type = lapply(all.overlaps.list2, function(q){
    g.hot = table(unlist(lapply(q, function(x) x[[1]]$wild.rel)))
    g.hot.sp = table(germplasm[match(names(g.hot), germplasm$Filename), ]$`Sample type`)    
    g.hot.sp
})

g.cold = table(unlist(lapply(all.overlaps.av, function(x) x[[2]]$wild.rel)))
g.cold.sp = table(germplasm[match(names(g.cold), germplasm$Filename), ]$`Sample type`)

#     ____________________________________________________________________________
#     ANALYSE ALL COLD/HOTSPOTS                                                                                             ####

#gets all of the coldspots that overlap with introgressions
#for ease of running statistics on the whole dataset
best.cold = lapply(all.overlaps.list2.rm.redun, function(q){
    q1 = lapply(q, function(x){
        # browser()
        g = x$cold
        g3 = newdf(colnames(g), no.rows = T)
        if(nrow(x$cold) > 0){
            
            un = unique(g$subjectHits)
            
            for(i in un){
                g1 = g[which(g$subjectHits == i), ]
                g2 = g1[match(max(g1$per.overlap), g1$per.overlap), ]
                g3 = rbind(g3, g2)
            }
        }
        return(g3)
    })
    return(q1)
}) 

best.cold2 = lapply(best.cold, function(x) if(length(x) > 0) combine.list.of.data.frames(x))
best.cold3 = combine.list.of.data.frames(best.cold2)

#gets all of the hotspots that overlap with introgressions
#for ease of running statistics on the whole dataset
best.hot = lapply(all.overlaps.list2.rm.redun, function(q){
    q1 = lapply(q, function(x){
        # browser()
        g = x$hot
        g3 = newdf(colnames(g), no.rows = T)
        if(nrow(x$hot) > 0){
            
            un = unique(g$subjectHits)
            
            for(i in un){
                g1 = g[which(g$subjectHits == i), ]
                g2 = g1[match(max(g1$per.overlap), g1$per.overlap), ]
                g3 = rbind(g3, g2)
            }
        }
        return(g3)
    })
    return(q1)
}) 

best.hot2 = lapply(best.hot, function(x) if(length(x) > 0) do.call(rbind, x))
if(length(which(unlist(lapply(best.hot2, is.null)))) > 0) best.hot2 = best.hot2[-which(unlist(lapply(best.hot2, is.null)))]
best.hot3 = do.call(rbind, best.hot2)
best.hot3$origin = row.names(best.hot3)
best.hot4 = arrange(best.hot3, cm.density)
best.hot4 = best.hot4[nrow(best.hot4):1, ]

best.hot4$recomb.width = as.numeric(best.hot4$recomb.end) - as.numeric(best.hot4$recomb.start)

#get ranges for all hotspots and coldspots, regardless of whether they overlap with 
#introgressions or not
hotspot.and.coldspot.ranges = lapply(all.m.8.iwgsc.4b.rev, function(x) lapply(x, find.hotspot.and.coldspot.ranges))
names(hotspot.and.coldspot.ranges) = names(all.m.8.iwgsc.4b.rev)

#calculate mean, s.d. of physical lengths of hotspots and coldspots; also find the 
#median centiMorgan density for hotspots. each row represents a chromosome in
#the genetic map
stat = lapply(lapply(hotspot.and.coldspot.ranges, function(x){
    chromostat = lapply(x, function(q){
        hm = mean(q[[1]]$phys.length) / 1000000
        cm = mean(q[[2]]$phys.length) / 1000000
        hsd = sd(q[[1]]$phys.length / 1000000)
        csd = sd(q[[2]]$phys.length / 1000000)
        m.cm.den = median(q[[1]]$cm.density)
        
        g = data.frame(hm, hsd, cm, csd, m.cm.den)
        return(g)
    })
    return(chromostat)
}), combine.list.of.data.frames)

#see whether the phylogenetic distances differ between hotspots and coldspot introgression overlaps (they don't)
t.test(best.hot4$phy.dist, best.cold3$phy.dist)

#num. coldspots
n.cold = combine.list.of.data.frames(lapply(hotspot.and.coldspot.ranges, function(x){
    combine.list.of.data.frames(lapply(x, function(q) q$coldspot.df))
}))

n.hot = combine.list.of.data.frames(lapply(hotspot.and.coldspot.ranges, function(x){
    combine.list.of.data.frames(lapply(x, function(q) q$hotspot.df))
}))

find.total.hotspots.for.individual.map = function(map.number){
    g = lapply(all.m.8.iwgsc.4b.rev[[map.number]], find.hotspot.and.coldspot.ranges)
    number.hot = sum(unlist(lapply(g, function(x) nrow(x$hotspot.df))))
    number.cold = sum(unlist(lapply(g, function(x) nrow(x$coldspot.df))))
    q = list(number.hot, number.cold)
    names(q) =    c("hot", "cold")
    q
}

#grab the number of hot/coldspots in each individual map
tot.hotspots.ind.maps = lapply(1:4, find.total.hotspots.for.individual.map)

#grab the percentage of 
Map(function(overlaps, hotspots){
    hot = overlaps$hot / hotspots$hot
    cold = overlaps$cold / hotspots$cold
    g = list(hot, cold)
    names(g) = c("hot", "cold")
    g
}, num.overlaps2, tot.hotspots.ind.maps)


#prepare hotspot data for plot of all hotspots in terms of their cM density
#i.e. how hot are these hotspots, are the boiling or luke-warm?
#grab all of the cm.density values for hotspots
stat2 = lapply(lapply(hotspot.and.coldspot.ranges, function(x){
    chromostat = lapply(x, function(q){
        g = data.frame(q$hotspot.df$cm.density, q$hotspot.df$phys.pos.start, q$hotspot.df$phys.pos.end)
        return(g)
    })
    return(chromostat)
}), function(x) x)

#add chromosome column to dataframes in the list
stat2.1 = lapply(stat2, function(x){
    count2 = make.counter()
    lapply(x, function(q){
        q$chromo = names(x)[count2()]
        return(q)
    })
})

#combine chromosome dataframes within each genetic map
stat2.2 = lapply(stat2.1, combine.list.of.data.frames)

#add map column to dataframes in the list
count3 = make.counter()
stat2.3 = lapply(stat2.2, function(x){
    x$map = names(stat2.2)[count3()]
    return(x)
})

stat2.4 = combine.list.of.data.frames(stat2.3)

stat2.4$q.hotspot.df.cm.density = as.numeric(stat2.4$q.hotspot.df.cm.density)
stat2.5 = arrange(stat2.4, q.hotspot.df.cm.density)
stat2.6 = stat2.5[nrow(stat2.5):1, ]

#mark hotspots that overlap with introgressions
stat2.6$overlap = F
stat2.6$overlap[match(best.hot4$recomb.start, stat2.6$q.hotspot.df.phys.pos.start)] = T

stat2.6$recomb.width = as.numeric(stat2.6$q.hotspot.df.phys.pos.end) - as.numeric(stat2.6$q.hotspot.df.phys.pos.start)

#run some statistics

unlist(lapply(stat, function(x) mean(x$hm))) #get the mean of hotspot mean physical lengths
unlist(lapply(stat, function(x) mean(x$cm))) #get the mean of coldspot mean physical lengths
unlist(lapply(stat, function(x) mean(x$hsd))) #hotspot mean of s.ds
unlist(lapply(stat, function(x) mean(x$csd))) #coldspot mean of s.ds

#run t.test on mean coldspot length vs mean hotspot length
lapply(stat, function(x) t.test(x$cm, x$hm))

iwgsc.chromo.lengths = read.table("rotation1scripts_v4/original_data/iwgsc.chromosome.counts.txt")
iwgsc.chromo.lengths = iwgsc.chromo.lengths[, 2:1]
iwgsc.chromo.lengths[, 1] = substr(as.character(iwgsc.chromo.lengths[, 1]), 4, 5)
iwgsc.chromo.lengths[, 2] = iwgsc.chromo.lengths[, 2] / 1000000
iwgsc.chromo.lengths = reset.colnames(iwgsc.chromo.lengths)

#calculate marker density for each map
map.marker.densities = lapply(all.m.8.iwgsc.4b.rev, function(x){
    count1 = make.counter()
    mean(unlist(lapply(x, function(q){
        chromo = names(x)[count1()]
        length.chromo = iwgsc.chromo.lengths[grep(chromo, iwgsc.chromo.lengths$V1), 2]
        nrow(q) / length.chromo
    })))
})

plot((stat2.6$recomb.width / 1000000), log(stat2.6$q.hotspot.df.cm.density))
plot((best.hot4$recomb.width / 1000000), log(best.hot4$cm.density))

#     ____________________________________________________________________________
#     CM DENSITY PLOTTING                                                                                                         ####

library(gridExtra)
library(ggplot2)

#plot of cm.density as a function of width of hotspot
plot3 = ggplot(stat2.6, aes(x = (recomb.width / 1000000), y = log(q.hotspot.df.cm.density))) + 
    geom_point() + coord_cartesian(xlim = c(0, 100), ylim = c(-5, 10))
plot4 = ggplot(best.hot4, aes(x = (recomb.width / 1000000), y = log(cm.density))) + 
    geom_point() + coord_cartesian(ylim = c(-5, 10))

grid.arrange(plot3, plot4, ncol = 2)

#all hotspots and their cm.densities
plot1 = ggplot(stat2.6, aes(x = 1:nrow(stat2.6), y = log(stat2.6$q.hotspot.df.cm.density))) + geom_point() +
    theme_bw() + ylab("Ln(cM Density) (cM/Mb)") +
    xlab("Index") + scale_color_manual(values = c("#000000", "#fff200")) + scale_shape(solid = F) +
    coord_cartesian(xlim = c(0, 3303), ylim = c(-7, 10)) + geom_hline(yintercept = log(1), color = "red", size = 1) +
    ggtitle("(a)") + geom_hline(yintercept = 1.299048, color = "blue", size = 1)

#only hotspots which overlap with introgressions
plot2 = ggplot(best.hot4, aes(x = 1:nrow(best.hot4), y = log(cm.density))) + geom_point() + theme_bw() + xlab("Index") +
    ylab("Ln(cM Density) (cM/Mb)") +
    coord_cartesian(xlim = c(0, 3303), ylim = c(-7, 10)) + 
    geom_hline(yintercept = log(1), color = "red", size = 1) +
    ggtitle("(b)") + geom_hline(yintercept = 1.299048, color = "blue", size = 1)


pdf("rotation1scripts_v4/plots/hotspot.cm.density.plots/cmdensity.w.introgressionsv3.pdf", width = 8, height = 5)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

#     ____________________________________________________________________________
#     FUNCTIONS                                                                                                                             ####

grab.valid.markers = function(introgression.data, parent.num){
    #grab introgression markers that are in 1-25% or 75-100% of chromosome (not near centromere)
    #args:
    #introgression.data: object generated by gen.intr()
    #parent.num: integer, either 1 or 2, will grab markers for that parent
    q = lapply(introgression.data, function(x){
        g = x[[parent.num]][[2]]
        
        if(is.null(nrow(g)) == F){
            g2 = g[abs(g$phys.dist - 50) > 25, ]    
            return(g2)
        } else {
            return()
        }
        
    })
    q[-which(unlist(lapply(q, is.null)))]
}

bin.valid.markers = function(valid.markers){
    g = lapply(valid.markers, function(x) x$phys.dist)
    unlist(lapply(g, function(x) nunique(as.character(cut(x, seq(0, 100, 1)))))) #bins of 2% correspond to roughly 12Mb each
}

get.markers.w.no.recomb.and.no.intr = function(genetic.map, introgression.data){
    #returns the number of markers that don't have any recombination between them, 
    #but also don't map to any introgressions
    #args:
    #genetic.map: a dataframe containing the genetic map
    #introgression.data: introgression information from gen.intr()
    Map(function(gen.map, chromo.name){
        q = gen.map
        coord = which(q$cm.diff == 0)
        coord = coord[2:length(coord)]
        q2 = q[coord, ][which(q$phys.dist.bp[coord] - q$phys.dist.bp[(coord-1)] > 10000), ]
        
        q3 = q2[abs(q2$phys.dist - 50) > 25, ]
        
        coord2 = introgression.data[[chromo.name]][[1]][[4]]
        coord3 = introgression.data[[chromo.name]][[2]][[4]]
        q.no.intro = q[q$cm.diff == 0, ][coord2, ]
        q.no.intro2 = q[q$cm.diff == 0, ][coord3, ]
        
        q4 = q3[which(q3$marker %in% q.no.intro$marker), ]
        q5 = q4[which(q4$marker %in% q.no.intro2$marker), ]
        nrow(q5)
    }, genetic.map, names(genetic.map))
}

grab.geno.vector.from.rqtl.format = function(markername, rqtl.df){
    #NB. Assumes the marker names are the column names of the data frame
    g = unname(unlist(rqtl.df[, which(colnames(rqtl.df) == markername)]))
    g[2:length(g)]
}