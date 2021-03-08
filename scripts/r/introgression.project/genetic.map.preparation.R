# setwd("C:/Users/ac14037/project.phd.main/")
setwd("E:/phd.project.main/")
# setwd("~/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

#     ____________________________________________________________________________
#     FUNCTIONS                                                                                                                             ####

fix.flipped.chromo = function(df.to.check){
    flipped.chromosomes = lapply(unique(df.to.check$chr), function(x){
        print(x)
        g = filter(df.to.check, chr == x)
        start.mean = mean(g$phys.dist.bp[1:6])
        end.mean = mean(g$phys.dist.bp[(nrow(g)-5):nrow(g)])
        # browser()
        if(start.mean > end.mean){
            return(x)
        } else{
            return()
        }
    })
    
    flipped.chromosomes = unlist(flipped.chromosomes)
    print(p(flipped.chromosomes, " are flipped"))
    
    #reverse selected chromosomes, reset centiMorgan values
    flipped.chromosome.dfs = lapply(flipped.chromosomes, function(x){
        g = filter(df.to.check, chr == x)
        g = g[nrow(g):1, ]
        
        #reverse cm distance
        total.cm = max(as.numeric(g$cM))
        
        g$cM = abs(g$cM - total.cm)
        
        return(g)
    })
    
    names(flipped.chromosome.dfs) = flipped.chromosomes
    
    #replace flipped chromosomes with newly inverted ones
    for(i in names(flipped.chromosome.dfs)){
        df.to.check[which(df.to.check$chr == i), ] = flipped.chromosome.dfs[[i]]
    }
    
    return(df.to.check)
    
}

extract.monotonic.chromo = function(chromo, fil.df){
    if(missing(fil.df)) fil.df = cs.x.p.map.w.genocounts
    g = filter(fil.df, chr == chromo)
    
    g[match(na.omit(g$phys.dist.bp)[longest_subseq.R(na.omit(g$phys.dist.bp))], g$phys.dist.bp), ]
}

#     ____________________________________________________________________________
#     FURTHER SACHA GENETIC MAP PROCESSING (Other biparental crosses)                 ####

library(readxl)

cs.x.p.map = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", sheet = 6)
a.x.c.map = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", sheet = 2)
s.x.r.map = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", sheet = 3)
o.x.s.map = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", sheet = 4)
a.x.p.map = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", sheet = 5)
a.x.p.map = a.x.p.map[-which(a.x.p.map$chr == "5A2" | a.x.p.map$chr == "7A2" | a.x.p.map$chr == "7A3"), ] #remove additional LGs from a.x.p.map as these won't work with genetictophysical()

all.maps = list(a.x.c.map, s.x.r.map, o.x.s.map, a.x.p.map, cs.x.p.map)

all.maps.trunc.2 = lapply(all.maps, function(q){
    g = q[, 1:3]
    # browser()
    phys.distances = lapply(unique(q$chr), function(x){
        print(x)
        g1 = g[g$chr == x, 1]
        # browser()
        genetictophysical(x, unname(unlist(g1)), "-", 30, F, nrgene = F)
        
    })
    
    g$phys.dist.per = unlist(phys.distances)
    g$phys.dist.bp = unlist(lapply(phys.distances, function(x) attr(x, "bp")))
    
    
    
    return(g)
    
})

#inverse any chromosomes in sacha's data that have the cM distance in the wrong order
order.cm = function(data.frame1){
    if(data.frame1$cM[1] > data.frame1$cM[nrow(data.frame1)]){
        data.frame1 = data.frame1[nrow(data.frame1):1, ]
    }
    return(data.frame1)
}

all.maps.trunc2 = lapply(all.maps.trunc.2, function(x){
    g2 = subapply(x, "chr", order.cm)    
})

order.within.bins = function(df){
    df2 = df
    #within bins, order cs.x.p.map.w.genocounts by physical sequence
    unique.chromo.cm.combos = unique(df2[, 2:3])
    for(i in 1:nrow(unique.chromo.cm.combos)){
        g = df2[which(df2$chr == as.character(unique.chromo.cm.combos[i, 1]) & df2$cM == as.numeric(unique.chromo.cm.combos[i, 2])), ]
        sorted = sort(g$phys.dist.bp, na.last = F)
        g = g[match(sorted, g$phys.dist.bp), ]
        # browser()
        df2[which(df2$chr == as.character(unique.chromo.cm.combos[i, 1]) & df2$cM == as.numeric(unique.chromo.cm.combos[i, 2])), ] = g
    }
    return(df2)
}

all.maps.trunc.ord = lapply(all.maps.trunc2, order.within.bins)
all.maps.trunc.ord.wo.na = lapply(all.maps.trunc.ord, function(x){
    x = x[-which(is.na(x$phys.dist.bp)), ]
    return(x)
})

#NB. Requires fix.flipped.chromo() from historical.recombination.mapping.R

all.maps.trunc.ord.wo.na.flip = lapply(all.maps.trunc.ord.wo.na, fix.flipped.chromo)

check.inverse.cM = function(df.to.check, chr.to.fil){
    g = filter(df.to.check, chr == chr.to.fil)
    # browser()
    if(mean(g$cM[1:6]) > mean(g$cM[(nrow(g) - 5):nrow(g)])){
        print(p(chr.to.fil, "inverse cM "))
        return(chr.to.fil)
    } 
}

c4 = make.counter()

all.maps.fix.rev.cM = lapply(all.maps.trunc.ord.wo.na.flip, function(x){
    print(c4())
    list.chro = unique(x$chr)
    any.rev.cm = unlist(lapply(list.chro, function(y) check.inverse.cM(x, y)))
    if(!is.null(any.rev.cm)){
        lapply(any.rev.cm, function(chr.to.rev){
            temp.x = filter(x, chr == chr.to.rev)
            temp.x = temp.x[nrow(temp.x):1, ]
            x[which(x$chr == chr.to.rev), ] <<- temp.x
        })
    }
    return(x)
})

all.maps.fix.rev.cM2 = lapply(all.maps.fix.rev.cM, order.within.bins)

#remove small chromosomes and chromosomes with bad concordance to physical map

all.maps.good.conc = lapply(all.maps.fix.rev.cM2, function(x){
    check.condordance.to.phys = function(df1){
        g = extract.monotonic.chromo(df1$chr[1], df1)
        (nrow(g) / nrow(df1) * 100)
    }
    conc.stats = subapply(x, "chr", check.condordance.to.phys)
    
    to.keep = which(conc.stats[, 1] > 80)
    q = split.df.into.list.of.dfs.by.column(x, "chr")
    q2 = combine.list.of.data.frames(q[to.keep])
    return(q2)
})

find.max.gap = function(phys.positions){
    #examines how well distributed the SNPs are across the chromosome by physical position
    #the smaller the number the better (0 = perfect, 100 = worst)
    max(unlist(lapply(1:100, function(x){
        min(abs(phys.positions - x))
    })))
}

#remove any chromosomes with massive gaps w/ respect to physical assembly (or e.g. small chromosomes)
all.maps.good.conc2 = lapply(all.maps.good.conc, function(x){
    g = split.df.into.list.of.dfs.by.column(x, "chr")
    concordance = unlist(lapply(g, function(q) find.max.gap(q$phys.dist.per)))
    to.keep = which(concordance < 30)
    combine.list.of.data.frames(g[to.keep])
})

#remove markers that have large physical distances from both neighbouring markers
all.m.3 = lapply(all.maps.good.conc2, function(x){
    remove.large.discrepancies = function(df1){
        g = calculate.mean.diff.neighbours(df1$phys.dist.per)
        g.coord = which(g[2:(length(g)-1)] > 40) + 1
        print(g.coord)
        if(length(g.coord) > 0){
            df1 = df1[-g.coord, ]    
        }
        return(df1)
    }
    subapply(x, "chr", remove.large.discrepancies)
})

all.m.4 = lapply(all.m.3, function(x){
    che = function(df1){
        df1$cm.diff = calculate.cm.diff(df1$cM)
        return(df1)
    }
    subapply(x, "chr", che)
})

# s(all.m.4, "rotation1scripts_v4/saved.objects/all.m.4", "further.sacha.genetic.map.processing.R")

all.m.4.1 = lapply(all.m.4, function(x){
    colnames(x)[4] = "phys.dist"
    return(x)
})

all.m.5 = lapply(all.m.4.1, function(x) split.df.into.list.of.dfs.by.column(x, "chr"))
names(all.m.5) = c("a.x.c", "s.x.r", "o.x.s", "a.x.p", "cs.x.p")

#     ____________________________________________________________________________
#     CALCULATE LONGEST INCREASING SUBSEQUENCES                                                             ####

longest.subs = lapply(all.m.5, function(q) unlist(lapply(q, function(x) x$marker[longest_subseq.R(x$phys.dist.bp)])))

all.maps.longest.subseq = Map(function(map, longest.subseq){
    coords = match(longest.subseq, map$marker)
    map = map[coords, ]
    return(map)
}, all.maps, longest.subs)

all.maps.geno.data = lapply(all.maps.longest.subseq, function(x){
    g = as.data.frame(t(x))
    colnames(g) = unname(unlist(g[1, ]))
    g = g[-3, ]
    
    g = cbind(g[, 1], g)
    g[, 1] = rownames(g)
    g = reset.rownames(g)
    g[1, 1] = ""
    g[2, 1] = ""
    
    g = g[-1, ]
    colnames(g) = c("", colnames(g)[2:ncol(g)])
    return(g)
})

count1 = make.counter()
lapply(all.maps.geno.data, function(x){
    write.csv(x, p("rotation1scripts_v4/processed_data/genotypes/sacha.longest.incr.subsequence/map", count1(), ".csv"), row.names = F)
})

mapfiles = c("map1.csv", "map2.csv", "map3.csv", "map4.csv", "map5.csv")

library(ASMap)
library(qtl)
longest.incr.cross = lapply(mapfiles, function(x){
    sacha.cross = read.cross(format = "csv", dir = "./", file = p("rotation1scripts_v4/processed_data/genotypes/sacha.longest.incr.subsequence/", x), genotypes = c("A", "B"), estimate.map = F)
    sacha.cross = convert2riself(sacha.cross)
    sacha.map.reesti = quickEst(sacha.cross, map.function = "kosambi")
    sacha.map.reesti2 = pull.map(sacha.map.reesti, as.table = T)
    return(sacha.map.reesti2)
})

#keep chromosomes that are above 40 and below 300 centimorgans in length
longest.incr.cross2 = lapply(longest.incr.cross, function(x){
    mm = function(q){
        g = max(q$pos)
        if(g > 40 & g < 300) return(q)
    } 
    subapply(data.frame.to.edit = x, colname.to.split.df.by = "chr", mm)
})

lapply(longest.incr.cross2, function(x) unique(x$chr))

all.m.6 = lapply(all.m.5, combine.list.of.data.frames)

#match up newly calculated maps with previous maps - now phys order is perfectly monotonic
#NB. Sacha AxP map did not get calculated and so is excluded here
all.m.7 = Map(function(long.incr.markers, old.gen.map){
    coords = match(rownames(long.incr.markers), old.gen.map$marker)
    old.gen.map2 = old.gen.map[coords, ]
    old.gen.map2$cM = long.incr.markers$pos
    if(nrow(old.gen.map2) > 0){
        addcmdiff.col = function(df1){
            df1$cm.diff = calculate.cm.diff(df1$cM)
            return(df1)
        }
        subapply(data.frame.to.edit = old.gen.map2, colname.to.split.df.by = "chr", addcmdiff.col)
    }
    return(old.gen.map2)
}, longest.incr.cross2, all.m.6)

#NOW DO A COMPARISON BETWEEN OLD V. W/ REFSEQ V1.0 CF. NRGENE

#remove apogee x paragon cross
all.m.7.1 = all.m.7[-4]

all.m.8.iwgsc.4b.rev = lapply(all.m.7.1, split.df.into.list.of.dfs.by.column, colname.to.split = "chr")

names(all.m.8.iwgsc.4b.rev) = c("a.x.c", "s.x.r", "o.x.s", "cs.x.p")

# s(all.m.8.iwgsc.4b.rev, "rotation1scripts_v4/saved.objects/all.m.8.iwgsc.4b.rev", "genetic.map.preparation.R")
# s(all.m.8.iwgsc, "rotation1scripts_v4/saved.objects/all.m.8.iwgsc", "further.sacha.genetic.map.processing.R")
# s(all.m.7, "rotation1scripts_v4/saved.objects/all.m.7", "further.sacha.genetic.map.processing.R")
# s(all.maps.geno.data, "rotation1scripts_v4/saved.objects/all.maps.geno.data", "genetic.map.preparation.R")
