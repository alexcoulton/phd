#### INIT ####
library(qtl)
library(gridExtra)
library(ggplot2)
library(IRanges)
library(tibble)
library(dplyr)
setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/recombination_analysis_functions.R")

#### FUNCTIONS #####

convert.recomb.to.seg = function(recomb.data){
    seg.data = unlist(lapply(recomb.data, function(x){
        # browser()
        a = length(which(x == "A"))
        b = length(which(x == "B"))
        a / b
    }))    
}

convert.recomb.to.seg2 = function(recomb.data){
    seg.data = lapply(recomb.data, function(x){
        # browser()
        a = length(which(x == "A"))
        b = length(which(x == "B"))
        chisq.test(c(a, b))
    })
}

check.for.distorted.markers = function(recomb.data){
    g = convert.recomb.to.seg2(recomb.data)
    g2 = which(unname(unlist(lapply(g, function(x) x[3]))) < 0.05)
    return(g2)
}

prepare.gg.data = function(recomb.data){
    new.data = data.frame(convert.recomb.to.seg(recomb.data), "", stringsAsFactors = F)
    new.data[, 2][check.for.distorted.markers(recomb.data)] = T
    colnames(new.data) = c("seg.ratio", "sig")
    new.data
}

make.plots = function(list.of.recomb.data){
    lapply(list.of.recomb.data, function(x){
        g = prepare.gg.data(x)
        # browser()
        ggplot(g, aes(x = 1:nrow(g), y = seg.ratio, group = sig, color = sig)) + geom_point() + 
            coord_cartesian(ylim = c(0.5, 1.6)) + geom_hline(yintercept = 1)
    })
}

#make function to add columns at selected position in dataframe (add_column not working correctly)
add.column.at.position = function(dataframe1, after_pos){
    #adds a new column to a dataframe, populated 1/2 with "A" and 1/2 "B"
    #args:
    # dataframe1 - arbitrary dataframe to add column to
    # after_pos - column coordinate to insert new column after
    if((nrow(dataframe1) %% 2) == 0){
        a1 = rep("A", (nrow(dataframe1) / 2))
        b1 = rep("B", (nrow(dataframe1) / 2))
        empty.col = as.character(c(a1, b1))
    } else {
        a1 = rep("A", ceiling(nrow(dataframe1) / 2))
        b1 = rep("B", floor(nrow(dataframe1) / 2))
        empty.col = as.character(c(a1, b1))
    }
    
    empty.col = as.character(rep("A", nrow(dataframe1)))
    dataframe2 = cbind(dataframe1[, 1:after_pos], empty.col, dataframe1[, (after_pos+1):ncol(dataframe1)])
    dataframe2 = convert.to.character.data.frame(dataframe2)
    dataframe2
}

split.geno.df.into.list.by.chromosome = function(y){
    #split the main dataframe into a list of dataframes, each containing a single chromosomes worth of genotyping data
    #args:
    # y - genotyping dataframe
    lapply(chromos, function(x){
        #retain first two rows in every dataframe
        g = c(1, 2, which(y[1, ] == x))
        # browser()
        ynew = y[, c(2, which(y[1, ] == x))]
        ynew
    })
}

#this function is a duplicate of convert.recomb.to.seg() in recombination.modelling.R
calc.ratio.a.to.b.alleles = function(z){
    #calculate ratio of A alleles to B alleles for each dataframe in list of chromosome genotype dataframes
    #args:
    # z - a list of genotyping dataframes
    lapply(z, function(y){
        
        unlist(lapply(y[, 2:ncol(y)], function(x){
            a = length(which(x == "A"))
            b = length(which(x == "B"))
            
            a / (a + b)
        }))
    })
    
}

calc.chi.sq.a.b = function(z){
    #calculate chi square p value for ratio of A alleles to B alleles in each chromosome genotyping dataframe
    #args:
    # z - a list of genotyping dataframes
    lapply(z, function(y){
        
        unlist(lapply(y[, 2:ncol(y)], function(x){
            a = length(which(x == "A"))
            b = length(which(x == "B"))
            h = length(which(x == "H"))
            # browser()
            chisq.test(c(a, h, b), p = c(1/4, 1/2, 1/4))[[3]]
        }))
    })
}

make.ratio.df = function(seg.ratio1, seg.chi1){
    ratiodf = as.data.frame(list(seg.ratio1, seg.chi1))
    colnames(ratiodf) = c("seg.ratio1", "seg.chi1")
    ratiodf$axis = 1:nrow(ratiodf)
    ratiodf$chi1t = ""
    ratiodf$chi1t[which(ratiodf$seg.chi1 < 0.05)] = "T"
    ratiodf$chi1t[which(ratiodf$seg.chi1 < 0.01)] = "T2"
    ratiodf$chi1t[which(ratiodf$seg.chi1 < 0.001)] = "T3"
    ratiodf
}

make.cake.plot = function(sim.data, title, facet1){
    #facet1 indicates whether a list of genotyping dataframes have been supplied.
    if(missing(facet1)) facet1 = F
    if(missing(title)) title = ""
    
    if(facet1 == T){
        alldat1 = lapply(sim.data, function(y){
            # browser()
            sim.data.f2.bar = as.data.frame(t(sapply(y, geno.count)))
            sim.data.f2.bar = as.data.frame(cbind(rownames(sim.data.f2.bar), sim.data.f2.bar))
            # browser()
            s2 = melt(sim.data.f2.bar)
            colnames(s2) = c("marker", "Genotype", "count")
            s2$percent = (s2$count / nrow(y)) * 100
            s2$marker = rep(1:ncol(y), 3)
            s2$seg = ""
            s2$sig.text = ""
            
            seg.ratios1 = sapply(convert.recomb.to.seg.f2(y), function(x) x[[3]])
            seg.marker1 = which(seg.ratios1 < 0.05)
            s2$seg[which(s2$marker %in% seg.marker1)] = T
            s2$sig.text[which(s2$seg == T & s2$Genotype == "A")] = "*"
            s2
        })
        
        alldat2 = Map(function(x, y){
            x$num1 = y
            x
        }, alldat1, c("(a)", "(b)", "(c)", "(d)"))
        
        alldat3 = combine.list.of.data.frames(alldat2)
        
        # browser()
        # browser()
        
        ggplot(alldat3, aes(x = marker, y = percent, fill = Genotype, alpha = seg )) +
            geom_bar(stat = "identity", width = 1) +
            theme_bw() + scale_fill_manual(values = c("#EE6C4D", "#231942", "#98C1D9"), labels = c("AA", "AB", "BB")) + 
            scale_alpha_discrete(range = c(1, 1)) + guides(alpha = F) + coord_cartesian(ylim = c(-4, 100)) + 
            scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
            geom_hline(yintercept = 25) + geom_hline(yintercept = 75) + ggtitle(title) +
            geom_text(aes(y = -4.5, label = sig.text)) + ylab("Percentage of Genotypes (%)") + xlab("Marker") + facet_grid(. ~ num1) +
            theme(axis.text.x = element_text(angle = 45), strip.text = element_text(size = 12), strip.background = element_rect(colour = "white", fill = "#ffffff"),
                        panel.grid = element_blank(), strip.placement = "outside")
        
    } else {
        sim.data.f2.bar = as.data.frame(t(sapply(sim.data, geno.count)))
        sim.data.f2.bar = as.data.frame(cbind(rownames(sim.data.f2.bar), sim.data.f2.bar))
        # browser()
        s2 = melt(sim.data.f2.bar)
        colnames(s2) = c("marker", "Genotype", "count")
        s2$percent = (s2$count / nrow(sim.data)) * 100
        s2$marker = rep(1:ncol(sim.data), 3)
        s2$seg = ""
        s2$sig.text = ""
        
        seg.ratios1 = sapply(convert.recomb.to.seg.f2(sim.data), function(x) x[[3]])
        seg.marker1 = which(seg.ratios1 < 0.05)
        s2$seg[which(s2$marker %in% seg.marker1)] = T
        s2$sig.text[which(s2$seg == T & s2$Genotype == "A")] = "*"
        
        ggplot(s2, aes(x = marker, y = percent, fill = Genotype, alpha = seg )) +
            geom_bar(stat = "identity", width = 1) +
            theme_classic() + scale_fill_manual(values = c("#EE6C4D", "#231942", "#98C1D9"), labels = c("AA", "AB", "BB")) + 
            scale_alpha_discrete(range = c(1, 1)) + guides(alpha = F) + coord_cartesian(ylim = c(-4, 110)) + 
            scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
            geom_hline(yintercept = 25) + geom_hline(yintercept = 75) + ggtitle(title) +
            geom_text(aes(y = -3, label = sig.text)) + ylab("Percentage of Genotypes (%)") + xlab("Marker")
    }
    
    

    # browser()
    # browser()
    
    
        
    #todo: need to add centimorgan positions to x-axis / keep only skeleton markers
}

#### FUNCTIONS - AXC SDR ANALYSIS ####

#how many of the markers exhibiting segregation distortion have neighbours that also exhibit seg. dist.? - i.e. segregation distortion regions
#first replicate

make.morestat1 = function(chicolname){
    nummore1 = lapply(ratiodfs2.3, function(x){
        g1 = grab.all1(x, chicolname)
        
        # mean(rle(diff(g1))$lengths)
        
        num1 = sum(which(rle(diff(g1))$lengths > 1)) #number of seg. dist. markers that have a neighbour also exhibiting seg. dist.
        morethan1marker = length(which(rle(diff(g1))$lengths > 1)) #number of seg. dist. regions (SDR)
        regions = rle(diff(g1))$lengths
        regionlength = regions[which(regions > 1)]
        #the coordinates made by the following line only have relevance within a single chicolname, they are not transferable between,
        #I've made a different function, grab.sdrs(), to investigate segregation distortion regions instead
        # marker.coords.in.sdr = which(diff(g1) == 1)
        
        list(num1, morethan1marker, regionlength) #marker.coords.in.sdr)
        
    })
    
    q1 = sum(unlist(lapply(nummore1, function(x) x[[1]])))
    q2 = sum(unlist(lapply(nummore1, function(x) x[[2]])))
    q3 = mean(na.omit(unlist(lapply(nummore1, function(x) x[[3]])))) # mean length of SDR
    q4 = sd(na.omit(unlist(lapply(nummore1, function(x) x[[3]])))) # mean length of SDR
    # browser()
    # q5 = lapply(nummore1, function(x) x[[4]])
    results1 = list(q1, q2, q3, q4) #, q5)
    names(results1) = c("num.markers.w.neighbour", "num.sdr", "mean.len.sdr", "sd.len.sdr") #, "marker.coords.in.sdr")
    results1
}


grab.sdrs = function(chicolname){
    i = 1
    lapply(ratiodfs2.3, function(x){
        lapply(x, function(y){
            # if(count1() == 11) browser()
            g = grab.all1(y, chicolname)
            # print(i)
            # if(i == 3) browser()
            
            
            #split numeric vector into list of numeric vectors, each element containing only consecutive integers
            diffs <- c(1, diff(g))
            start_indexes <- c(1, which(diffs > 1))
            end_indexes <- c(start_indexes - 1, length(g))
            coloned <- paste(g[start_indexes], g[end_indexes], sep=":")
            ranges1 = data.frame(coloned, stringsAsFactors = F)
            r1 = multi.str.split(ranges1[, 1], ":", 1)
            r2 = multi.str.split(ranges1[, 1], ":", 2)
            ranges1 = data.frame(r1, r2, stringsAsFactors = F)
            coords1 = which(ranges1[, 1] == ranges1[, 2])
            if(length(coords1) != 0) ranges1 = ranges1[-coords1, ]
            
            # if(count1() == 3) browser()
            
            if(nrow(ranges1) == 1 & ranges1[1, 1] == "NA"){
                return(NA)
            }
            #in some cases, only a single seg. dist. marker is found, which means ranges1 will be an empty df with no rows and no NAs
            #but we still want to remove this case as it is not an SDR
            if(nrow(ranges1) == 0) return(NA)
            
            
            ranges1[] = lapply(ranges1, as.numeric)
            
            
            
            #put into IRanges object
            IRanges(ranges1$r1, ranges1$r2)
        })
        # i <<- i + 1
    })
}

sdr.stat = function(sdr.df){
    #grab the mean and sd of the width of segregation distortion regions (defined as 2 consecutive markers exhibiting seg dist.)
    sdr.df.comb = sdr.df[-which(is.na(sdr.df))]
    sdr.df.comb = lapply(sdr.df.comb, as.data.frame)
    sdr.df.comb = combine.list.of.data.frames(sdr.df.comb)
    g = c(nrow(sdr.df.comb), sum(sdr.df.comb$width), mean(sdr.df.comb$width), sd(sdr.df.comb$width))
    names(g) = c("num.sdr.regions", "total.num.markers", "mean.width", "sd.width")
    g
}

numsegdistmarkers = function(chicolname){
    sum(unlist(lapply(ratiodfs2.3, function(x){
        length(grab.all1(x, chicolname))
    })))    
}

#### FUNCTIONS - RECOMBIATION PROFILES #### 

make.recombination.profiles = function(genotypelist, chromosome, halve.events, threshold){
    #args:
    # genotypelist - genotype dataframe made with makegenotypelist()
    # chromosome - string indicating the chromosome to make recombination profiles for
    # halve.events - boolean flag indicating whether to reduce the number of events in each recombination profile
    # threshold - integer, profiles with a higher number of events than this are removed
    
    if(missing(halve.events)) halve.events = T
    
    histlist1 = makehistlist(chromosome, genotypelist, T, homozygous = T)
    
    
    
    recomb.positions = lapply(unique(histlist1$individual), function(x){
        g = filter(histlist1, individual == x)
        g1 = as.numeric(g$position)
        g1
    })
    
    #remove profiles with over 5 events
    if(length(which(unlist(lapply(recomb.positions, function(x) length(x) > 5 )))) > 0){
        recomb.positions2 = recomb.positions[-which(unlist(lapply(recomb.positions, function(x){
            length(x) > 5
        })))]    
    } else {
        recomb.positions2 = recomb.positions
    }
    
    
    #remove ~half of the events from each recombination profile
    if(halve.events == T){
        recomb.positions3 = lapply(recomb.positions2, function(x){
            diffs = diff(x)
            diffs.sort = sort(diffs, index.return = T)$ix
            # browser()
            g = x
            if(length(x) > 1){
                g = x[-diffs.sort[1:round((length(x) / 2))]]
                browser()
            }
            g
        })    
    } else {
        recomb.positions3 = recomb.positions2
    }
    
    
    # mean(unlist(lapply(recomb.diffs, function(x) length(x))))
    recomb.positions3
    
    
}

#### PROCESSING ####

cxaf2asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa_noaxc_f2.csv", stringsAsFactors = F, header = T)
cxaf2.c651.asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa651-redone-map.csv")
cxaf2.c651.asmap$X = switch.affy.format(cxaf2.c651.asmap$X)
cxaf2.c651.asmap$marker = switch.affy.format(cxaf2.c651.asmap$marker)

cxaf2.c651.asmap.list = split.df.into.list.of.dfs.by.column(cxaf2.c651.asmap, "chromo", do.sort = F)

c.x.af2 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.csv", stringsAsFactors = F)

c.x.af2.no.axc = c.x.af2[c(1, 2, 3, grep("65-1|65-3", c.x.af2$probeset_id)), ]
c.x.af2.no.axc[1, 1] = ""
# write.csv(c.x.af2.no.axc, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.noaxc.csv", row.names = F)

#genotyping data with asmap genetic map order
c.x.af2.asmap.ord = c.x.af2[, c(1, 2, match(cxaf2.c651.asmap$marker, colnames(c.x.af2)))]
c.x.af2.asmap.ord[1, 3:ncol(c.x.af2.asmap.ord)] = cxaf2.c651.asmap$chromo

c.x.af2v2 = c.x.af2.asmap.ord

c.x.af2v3 = c.x.af2v2[c(1, grep("65-1|65-3", c.x.af2v2$probeset_id)), ]

# write.csv(c.x.af2v3, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.af2_onlycxa.geneticmaporder.csv", row.names = F)

#process a.x.c DH cross
a.x.cdh = read.csv("rotation1scripts_v4/processed_data/genotypes/avalon.x.cadenza/a.x.c.dh.flipped.probeset.matching.axcf2.csv", stringsAsFactors = F, header = T)

c.x.af2v2[1, ][1:10]

a.x.cdh2 = a.x.cdh[, na.omit(match(colnames(c.x.af2v2), colnames(a.x.cdh)))]
a.x.cdh2[1, ] = c.x.af2v2[1, na.omit(match(colnames(a.x.cdh2), colnames(c.x.af2v2)))]

nacolcoords = which(is.na(match(colnames(c.x.af2v2), colnames(a.x.cdh2)))) - 1

#add NA columns where missing probes from the DH population are
a.x.cdh3 = a.x.cdh2
for(i in nacolcoords){
    a.x.cdh3 = add.column.at.position(a.x.cdh3, i)
    # browser()
}

#update chromosome information
a.x.cdh3[1, ] = c.x.af2v2[1, ]

#process a.x.c f2 crosses

a.x.c512 = read.csv("rotation1scripts_v4/processed_data/genotypes/a.x.c.f2/axc51_2.csv", stringsAsFactors = F, header = T)
a.x.c611 = read.csv("rotation1scripts_v4/processed_data/genotypes/a.x.c.f2/axc61_1.csv", stringsAsFactors = F, header = T)

a.x.c512 = a.x.c512[, match(colnames(c.x.af2v2), colnames(a.x.c512))]
a.x.c512[1, ] = c.x.af2v2[1, ]

a.x.c611 = a.x.c611[, match(colnames(c.x.af2v2), colnames(a.x.c611))]
a.x.c611[1, ] = c.x.af2v2[1, ]

c.x.a651 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa651_redone.csv", stringsAsFactors = F, header = T)
c.x.a653 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa653_redone.csv", stringsAsFactors = F, header = T)

#add chromosome assignments to genotype data
c.x.a651 = c.x.a651[, match(colnames(c.x.af2v2), colnames(c.x.a651))]
c.x.a651[1, ] = c.x.af2v2[1, ]

#add chromosome assignments again
c.x.a653 = c.x.a653[, match(colnames(c.x.af2v2), colnames(c.x.a653))]
c.x.a653[1, ] = c.x.af2v2[1, ]

write.csv(c.x.a651, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa651_redone_w_chromos.csv", row.names = F)
write.csv(c.x.a653, "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa653_redone_w_chromos.csv", row.names = F)


c.x.acomb = rbind(c.x.a651, c.x.a653[3:nrow(c.x.a653), ])

#main analysis
#seperate out the two cadenza x avalon populations from the main dataframe
c651 = c.x.af2v2[c(1, grep("65-1", c.x.af2v2$probeset_id)), ]
c653 = c.x.af2v2[c(1, grep("65-3", c.x.af2v2$probeset_id)), ]

#get list of chromosomes in c.x.a cross
chromos = unique(unlist(as.list(c.x.a653[1, ])))[-(1:2)]

list.of.geno.dfs.each.pop = list(c.x.a651, c.x.a653, a.x.c512, a.x.c611, a.x.cdh3)

combined.cad.x.ava = list(rbind(list.of.geno.dfs.each.pop[[1]], list.of.geno.dfs.each.pop[[2]][2:nrow(list.of.geno.dfs.each.pop[[2]]), ]))
combined.ava.x.cad = list(rbind(list.of.geno.dfs.each.pop[[3]], list.of.geno.dfs.each.pop[[4]][2:nrow(list.of.geno.dfs.each.pop[[4]]), ]))

all.combined.geno.data = list(rbind(list.of.geno.dfs.each.pop[[1]], list.of.geno.dfs.each.pop[[2]][2:nrow(list.of.geno.dfs.each.pop[[2]]), ], list.of.geno.dfs.each.pop[[3]][2:nrow(list.of.geno.dfs.each.pop[[3]]), ], list.of.geno.dfs.each.pop[[4]][2:nrow(list.of.geno.dfs.each.pop[[4]]), ]))

list.of.geno.dfs.each.pop = c(list.of.geno.dfs.each.pop, all.combined.geno.data, combined.cad.x.ava, combined.ava.x.cad)

list.of.pop.names = c("c.x.a651", "c.x.a653", "a.x.c512", "a.x.c611", "a.x.cdh3", "all.comb", "cxa.comb", "axc.comb")

geno.dfs.2 = lapply(list.of.geno.dfs.each.pop, split.geno.df.into.list.by.chromosome)

geno.dfs.ratio = lapply(geno.dfs.2, calc.ratio.a.to.b.alleles)

length(geno.dfs.ratio[[1]][[1]])
length(geno.dfs.ratio[[5]][[1]])

geno.dfs.chi = lapply(geno.dfs.2, calc.chi.sq.a.b)

ratiodfs2 = Map(function(x, y){
    Map(function(a, b){
        make.ratio.df(a, b)
    }, x, y)
}, geno.dfs.ratio, geno.dfs.chi)

ratiodfs2.1 = lapply(1:15, function(y){
    lapply(ratiodfs2, function(x){
        x[[y]]
    })    
})

ratiodfs2.2 = lapply(ratiodfs2.1, function(x){
    lapply(x, function(y){
        y$fdr = p.adjust(y$seg.chi1, "BH")
        y    
        y$chifdr = ""
        y$chifdr[which(y$fdr < 0.05)] = "T"
        y$chifdr[which(y$fdr < 0.01)] = "T1"
        y$chifdr[which(y$fdr < 0.001)] = "T2"
        y
    })
    
})

ratiodfs2.2.2 = lapply(ratiodfs2.2, function(x){
    lapply(x, function(y){
        y$bonf = p.adjust(y$seg.chi1, "bonferroni")
        y    
        y$chibonf = ""
        y$chibonf[which(y$bonf < 0.05)] = "T"
        y$chibonf[which(y$bonf < 0.01)] = "T1"
        y$chibonf[which(y$bonf < 0.001)] = "T2"
        y$chimixed = ""
        y$chimixed[which(y$chi1t != "")] = "T"
        y$chimixed[which(y$chifdr != "")] = "T1"
        y$chimixed[which(y$chibonf != "")] = "T2"
        y
    })
    
})



#add centimorgan positions to dataframe
ratiodfs2.3 = Map(function(x, y){
    lapply(x, function(q){
        q$cm = y$pos
        q
    })
}, ratiodfs2.2.2, cxaf2.c651.asmap.list)




#### EMPIRICAL CAKE PLOTS CXA #### 

trunc.geno.dfs = geno.dfs.2[1:4]

trunc.geno.dfs2 = transpose.list.of.lists(trunc.geno.dfs)

q = trunc.geno.dfs2[[3]][[3]]
make.cake.plot(q[2:nrow(q), 2:ncol(q)])

i = 1
cake.plots1 = lapply(trunc.geno.dfs2, function(x){
    g = make.cake.plot(x[[1]][2:nrow(x[[1]]), 2:ncol(x[[1]])], p(chromos[i], " C x A Replicate 1"))
    g1 = make.cake.plot(x[[2]][2:nrow(x[[2]]), 2:ncol(x[[2]])], "C x A Replicate 2")
    g2 = make.cake.plot(x[[3]][2:nrow(x[[3]]), 2:ncol(x[[3]])], "A x C Replicate 1")
    g3 = make.cake.plot(x[[4]][2:nrow(x[[4]]), 2:ncol(x[[4]])], "A x C Replicate 2")
    
    
    i <<- i + 1
    list(g, g1, g2, g3)
})

#### CAKE PLOT PUBLICATION FIGURE #### 
#just 6b
q = trunc.geno.dfs2[[3]]
g = make.cake.plot(q[[1]][2:nrow(q[[1]]), 2:ncol(q[[1]])], "(a)") + theme(legend.position = "none")
g1 = make.cake.plot(q[[2]][2:nrow(q[[2]]), 2:ncol(q[[2]])], "(b)")
g2 = make.cake.plot(q[[3]][2:nrow(q[[3]]), 2:ncol(q[[3]])], "(c)") + theme(legend.position = "none")
g3 = make.cake.plot(q[[4]][2:nrow(q[[4]]), 2:ncol(q[[4]])], "(d)") + theme(legend.position = "none")

q2 = Map(function(x, y){
    x$plotnum = y
    x
}, q, 1:4)

allqdat = list(q[[1]][2:nrow(q[[1]]), 2:ncol(q[[1]])],
                             q[[2]][2:nrow(q[[2]]), 2:ncol(q[[2]])], 
                             q[[3]][2:nrow(q[[3]]), 2:ncol(q[[3]])], 
                             q[[4]][2:nrow(q[[4]]), 2:ncol(q[[4]])])

final.cake = make.cake.plot(allqdat, "", T)
ggsave("rotation1scripts_v4/plots/simulation/pedsim/Fig8v2.eps", final.cake, width = 17.4, height = 8, units = "cm", device = cairo_pdf)



q[[1]]


grid.arrange(g, g1, g2, g3)

cak2 = plot_grid(g, g1, g2, g3, ncol = 2, align = "v", axis = "r")

tiff("rotation1scripts_v4/plots/simulation/pedsim/cake.plot.tiff", units = "cm", width = 17.4, height = 7.3, res = 1000)
plot_grid(g, g1, g2, g3, ncol = 2, align = "v", axis = "r")
dev.off()


plot_grid(g, g1, g2, g3, ncol = 2)

save_plot("rotation1scripts_v4/plots/simulation/pedsim/cakeplot5.eps", cak2, ncol = 2, device = cairo_pdf, base_width = 6.850394, base_height = 7.2)

q1 = arrangeGrob(g, g1, g2, g3)



do.call(grid.arrange, cake.plots1[[3]])

lapply(1:15, function(x){
    pdf(p("rotation1scripts_v4/plots/cake.plots.axc.f2/", chromos[x], ".cakeplot.pdf"), 8, 8)
    do.call(grid.arrange, cake.plots1[[x]]) 
    dev.off()
})



q1 = trunc.geno.dfs2[[7]][[1]]
q2 = make.cake.plot(q1[2:nrow(q1), 2:ncol(q1)])
pdf(p("rotation1scripts_v4/plots/cake.plots.axc.f2/caketest.pdf"), 8, 8)
q2
dev.off()



#plot dataframes
yscale1 = c(0, 1)
library(gridExtra)
c1 = make.counter()    
lapply(ratiodfs2.3, function(xx){
    
    count = c1()
    plots1 = lapply(xx[1:4], function(y){
        ggplot(y, aes(x = cm, y = seg.ratio1, group = chimixed, color = chimixed, shape = chimixed)) + 
            geom_point() + geom_hline(yintercept = 0.5) + coord_cartesian(ylim = yscale1) + ggtitle(chromos[count]) + 
            theme_classic() + scale_color_manual(values=c("#492de5", "#999999", "#E69F00", "#56B4E9")) #+ geom_path()
    })

    pdf(p("rotation1scripts_v4/plots/c.x.af2/asmaporder/allpops/paperfig/just1to4_", chromos[count], ".pdf"), 15, 3)
    do.call(grid.arrange, c(plots1, list(ncol = 4)))
    dev.off()
})

p.adjust(ratiodfs2.2[[1]][[3]]$seg.chi1[1:60], "BH")

segratios1a = ratiodfs[[15]]
# s(segratios1a, "rotation1scripts_v4/saved.objects/segratios1a", "c.x.aprocessing2.R")

segratios5b = ratiodfs[which(chromos == "5B")]
# s(segratios5b, "rotation1scripts_v4/saved.objects/segratios5b", "c.x.aprocessing2.R")


#### PLOT SEGREGATION DISTORTION ####

#remove 4B and 5A as these are messed up in one of the a.x.c replicates

names(ratiodfs2.3) = chromos


bad.chromos1 = c(which(chromos == "4B"), which(chromos == "5A"))
ratiodfs2.3 = ratiodfs2.3[-c(bad.chromos1)]

#### FIRST REPLICATE STATS (CADENZA X AVALON) ####

#make a nice function to save some code (grab all 3 significance categories from chi-square test)
grab.all1 = function(df1, c.name){
    which(df1[[c.name]] == "T" | df1[[c.name]] == "T2" | df1[[c.name]] == "T3")
} 

sdr.regions1 = grab.sdrs("chi1t")

sdr.regions2 = lapply(sdr.regions1, function(x) x[1:4])

#transpose / invert a list of lists
# lapply(1:4, function(i) lapply(sdr.regions2, "[[", i))

sdr.regions3 = lapply(1:4, function(y){
    lapply(sdr.regions2, function(x){
        x[[y]]
    })
})

sdr.stat(sdr.regions3[[1]])
sdr.stat(sdr.regions3[[2]])
sdr.stat(sdr.regions3[[3]])
sdr.stat(sdr.regions3[[4]])


sdr.regions3 = lapply(sdr.regions3, function(x){
    names(x) = chromos[-bad.chromos1]
    x
})

#find overlapping SDRs
overlaps1 = Map(function(x, y){
    if(is.na(x[1]) | is.na(y[1])){
        NA
    } else {
        overlapsRanges(x, y)    
    }
}, sdr.regions3[[1]], sdr.regions3[[2]])

overlaps1 = overlaps1[!is.na(overlaps1)]
overlaps1 = overlaps1[sapply(overlaps1, function(x) length(x) > 0)]

overlaps2 = Map(function(x, y){
    if(is.na(x[1]) | is.na(y[1])){
        NA
    } else {
        overlapsRanges(x, y)    
    }
}, sdr.regions3[[3]], sdr.regions3[[4]])

overlaps2 = overlaps2[!is.na(overlaps2)]
overlaps2 = overlaps2[sapply(overlaps2, function(x) length(x) > 0)]






overlaps1







make.morestat1("chi1t")


#between the two replicates, how many markers show significant segregation distortion?
sum(unlist(lapply(ratiodfs2.3, function(x){
    # browser()
    g = grab.all1(x[[1]], "chi1t")
    g1 = grab.all1(x[[2]], "chi1t")
    length(unique(c(g, g1)))
})))

#how many of the markers that show significant seg. dist. in first pop also show it in second? 
sum(unlist(lapply(ratiodfs2.3, function(x){
    g = grab.all1(x[[1]], "chi1t")
    g1 = grab.all1(x[[2]], "chi1t")
    
    length(which(g %in% g1))
})))

#get names of markers that are distorted in both replicates
lapply(ratiodfs2.3, function(x){
    g = grab.all1(x[[1]], "chi1t")
    g1 = grab.all1(x[[2]], "chi1t")
    
    # browser()
    
    names1 = rownames(x[[1]][g[which(g %in% g1)], ])
    names1
    
})

#when both datasets are combined, how many markers show significant segregation distortion?
length(unlist(lapply(ratiodfs2.3, function(x){
    g = grab.all1(x[[7]], "chi1t")
    g
})))

#of those markers that are distorted in the combined dataframe, how many are distorted in the original dataframe and how many are new distortions?
length(unlist(lapply(ratiodfs2.3, function(x){
    g1 = grab.all1(x[[1]], "chi1t")
    g1.1 = grab.all1(x[[2]], "chi1t")
    g2 = grab.all1(x[[7]], "chi1t")
    
    which(g2 %in% g1 | g2 %in% g1.1)
})))

#how many seg dist markers in the first replicate
length(unlist(lapply(ratiodfs2.3, function(x){
    grab.all1(x[[1]], "chi1t")
})))

#how many seg dist markers in the second replicate
length(unlist(lapply(ratiodfs2.3, function(x){
    grab.all1(x[[2]], "chi1t")
})))


#### SECOND REPLICATE STATS (AVALON X CADENZA) ####

sum(unlist(lapply(ratiodfs2.3, function(x){
    g = grab.all1(x[[3]], "chi1t")
    g1 = grab.all1(x[[4]], "chi1t")
    
    length(unique(c(g, g1)))
})))

#how many of the markers that show significant seg. dist. in first pop also show it in second? 
sum(unlist(lapply(ratiodfs2.3, function(x){
    g = grab.all1(x[[3]], "chi1t")
    g1 = grab.all1(x[[4]], "chi1t")
    
    length(which(g %in% g1))
})))

#when both datasets are combined, how many markers show significant segregation distortion?
length(unlist(lapply(ratiodfs2.3, function(x){
    g = grab.all1(x[[8]], "chi1t")
    g
})))

#of those markers that are distorted in the combined dataframe, how many are distorted in the original dataframe and how many are new distortions?
length(unlist(lapply(ratiodfs2.3, function(x){
    g1 = grab.all1(x[[3]], "chi1t")
    g1.1 = grab.all1(x[[4]], "chi1t")
    g2 = grab.all1(x[[8]], "chi1t")
    
    which(g2 %in% g1 | g2 %in% g1.1)
})))

#how many seg dist markers in the first replicate
length(unlist(lapply(ratiodfs2.3, function(x){
    grab.all1(x[[3]], "chi1t")
})))

#how many seg dist markers in the second replicate
length(unlist(lapply(ratiodfs2.3, function(x){
    grab.all1(x[[4]], "chi1t")
})))






#### RECOMBINATION PROCESSING ####

all.geno.lists1 = lapply(list.of.geno.dfs.each.pop, makegenotypelist)

load("rotation1scripts_v4/saved.objects/all.m.8.iwgsc.4b.rev")

a.x.c.sacha = combine.list.of.data.frames(all.m.8.iwgsc.4b.rev[[1]])


c.x.af2.gen.ord = c.x.af2[, c(1, 2, na.omit(match(switch.affy.format(a.x.c.sacha$marker), colnames(c.x.af2))))]
c.x.af2.gen.ord[1, 1:ncol(c.x.af2.gen.ord)] = a.x.c.sacha$chr[match(colnames(c.x.af2.gen.ord), switch.affy.format(a.x.c.sacha$marker))]


cxa.gen.651 = c.x.af2.gen.ord[c(1, grep("CxA 65-1", c.x.af2.gen.ord$probeset_id)), ]
cxa.gen.653 = c.x.af2.gen.ord[c(1, grep("CxA 65-3", c.x.af2.gen.ord$probeset_id)), ]

cxa.genord.geno = makegenotypelist(cxa.gen.651)
cxa.genord.geno2 = makegenotypelist(cxa.gen.653)

chr.to.check = "1A"
cxa.genord.hist = makehistlist(chr.to.check, cxa.genord.geno, T)
cxa.genord.hist2 = makehistlist(chr.to.check, cxa.genord.geno, F)
cxa.genord.hist3 = makehistlist(chr.to.check, cxa.genord.geno2, T)
cxa.genord.hist3 = makehistlist(chr.to.check, cxa.genord.geno2, F)


# par(mfrow = c(2, 1))
hist(as.numeric(cxa.genord.hist2), breaks = 300)
hist(as.numeric(cxa.genord.hist3), breaks = 300)


v(cxa.genord.hist)
library(dplyr)
pos.diffs = lapply(unique(cxa.genord.hist$individual), function(x){
    g = filter(cxa.genord.hist, individual == x)
    g1 = as.numeric(g$position)
    g1
    
    
})

#remove individuals with over 5 events
pos.diffs2 = pos.diffs[-which(unlist(lapply(pos.diffs, function(x){
    length(x) > 5
})))]

length(which(c.x.af2v2[1, ] == "1A"))

#save recombination profiles for use in simulation
recombination.profiles = pos.diffs2
#s(recombination.profiles, "rotation1scripts_v4/saved.objects/recombination.profiles", "c.x.aprocessing2.R")

recombination.profiles1a = make.recombination.profiles(geno12, "1A")
recombination.profiles1a.full = make.recombination.profiles(geno12, "1A", halve.events = F)

halved.profile1a = table(sapply(recombination.profiles1a, length))
full.profile1a = table(sapply(recombination.profiles1a.full, length))




cxa.new.geno1 = makegenotypelist(c.x.acomb)
recomb1a.v2 = make.recombination.profiles(cxa.new.geno1, "1A")

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

recomb1a.v2 %in% recombination.profiles1a

recombination.profiles1a.full = make.recombination.profiles(geno12, "6A", F)
r1 = unlist(recombination.profiles1a.full)
length(r1)
which(duplicated(r1))
lapply(recombination.profiles1a.full, length)


geno12 = makegenotypelist(c.x.af2v2)

recombination.profiles5b = make.recombination.profiles(geno12, "5B")
# s(recombination.profiles1a, "rotation1scripts_v4/saved.objects/recombination.profiles1a", "c.x.aprocessing2.R")
# s(recombination.profiles5b, "rotation1scripts_v4/saved.objects/recombination.profiles5b", "c.x.aprocessing2.R")

#making a reduced recombination profile with evenly spaced markers along 1A
c.x.af2v2.1a = c.x.af2v2[, c(1, 2, which(c.x.af2v2[1, ] == "1A"))]
c.x.af2v2.1a.reduced = c.x.af2v2.1a[, c(1, 2, seq(3, 226, 27))]
geno12.1a.red = makegenotypelist(c.x.af2v2.1a.reduced)
recombination.profiles1a.reduced = make.recombination.profiles(geno12.1a.red, "1A")
s(recombination.profiles1a.reduced, "rotation1scripts_v4/saved.objects/recombination.profiles1a.reduced", "c.x.aprocessing2.R")


c.x.af2v2.1a = c.x.af2v2[, c(1, 2, which(c.x.af2v2[1, ] == "1A"))]
c.x.af2v2.1a.reduced = c.x.af2v2.1a[, c(1, 2, seq(3, 226, 70))]
geno12.1a.red = makegenotypelist(c.x.af2v2.1a.reduced)
recombination.profiles1a.reduced.further = make.recombination.profiles(geno12.1a.red, "1A")

s(recombination.profiles1a.reduced.further, "rotation1scripts_v4/saved.objects/recombination.profiles1a.reduced.further", "c.x.aprocessing2.R")

ncol(c.x.af2v2.1a.reduced)


cxa.genord.hist$individual = as.numeric(cxa.genord.hist$individual)

cxa.split = split.df.into.list.of.dfs.by.column(cxa.genord.hist, "individual")


sortbymean = sort(unlist(lapply(pos.diffs, function(x) mean(x))), index.return = T)$ix
posw10 = which(unlist(lapply(pos.diffs, function(x){
    g = which(x == 10)
    if(length(g) == 0) return(F)
    return(T)
})))

cxa.split2 = cxa.split[posw10]
cxa.unsplit = combine.list.of.data.frames(cxa.split2)
cxa.unsplit$position = as.numeric(cxa.unsplit$position)

ggplot(cxa.unsplit, aes(x = individual, y = position, group = individual, color = individual)) + geom_point()



#do recombination analysis / plots

# c651 = c651[, -1]
# c653 = c653[, -1]

c651 = convert.to.character.data.frame(c651)
c653 = convert.to.character.data.frame(c653)


c651.geno = makegenotypelist(c651)
c653.geno = makegenotypelist(c653)


#analyse frequency of recombination events
g = lapply(asmaplgs, function(x){
    g = make.recombination.profiles(c651.geno, x, F)
    mean(unlist(lapply(g, length)))
})

mean(na.omit(unlist(g)))
sd(na.omit(unlist(g)))

#analyse frequency for second CxA replicate
g2 = lapply(asmaplgs, function(x){
    g = make.recombination.profiles(c653.geno, x, F)
    mean(unlist(lapply(g, length)))
})

mean(na.omit(unlist(g2)))
sd(na.omit(unlist(g2)))

make.recombination.profiles(c651.geno, "1A", F)
make.recombination.profiles(c651.geno, "1A", F)


c651.hist = makehistlist("5B", c651.geno, T)
c653.hist = makehistlist("1B", c653.geno, T)

c651.hist = makehistlist("5B", c651.geno, F)
c653.hist = makehistlist("1B", c653.geno, F)

lgs = paste("chr", listofwheatchromosomes, sep = "")
asmaplgs = unique(cxaf2.c651.asmap$chromo)


c65.all = lapply(asmaplgs, function(x){
    lapply(all.geno.lists1, function(y){
        as.numeric(makehistlist(x, y, F))
    })
    # c651hist = as.numeric(makehistlist(x, c651.geno, F))
    # c653hist = as.numeric(makehistlist(x, c653.geno, F))
    # list(c651hist, c653hist)
})




#### PLOTTING RECOMBINATION ####
count1 = make.counter()
hists1 = lapply(c65.all, function(x){
    chromo = asmaplgs[count1()]    
    count2 = make.counter()
    lapply(x, function(y){
        make.ggplot.histogram(y, 100, "Marker No.", p(list.of.pop.names[count2()], " ", chromo))
    })
})

c1 = make.counter()
lapply(hists1, function(x){
    g = c1()
    pdf(p("rotation1scripts_v4/plots/c.x.af2/recombination/asmapord/allpops/gg", asmaplgs[g], "recomb.pdf"), 20, 13)
    do.call(grid.arrange, x)
    dev.off()
})

do.call(grid.arrange, hists1)

q = as.data.frame(c65.all[[1]][[1]])
ggplot(q, aes(q[, 1])) + geom_histogram(bins = 200) + labs(title = "chocho")

#recombination plots
c1 = make.counter()
lapply(c65.all, function(x){
    g = c1()
    pdf(p("rotation1scripts_v4/plots/c.x.af2/recombination/asmapord/onlycxa/gg", asmaplgs[g], "recomb.pdf"), 20, 13)
    q1 = as.data.frame(x[[1]])
    q2 = as.data.frame(x[[2]])
    plot1 = ggplot(q1, aes(q1[, 1])) + geom_histogram(bins = 200) + labs(title = asmaplgs[g])
    plot2 = ggplot(q2, aes(q2[, 1])) + geom_histogram(bins = 200) + labs(title = asmaplgs[g])    
    grid.arrange(plot1, plot2)
    
    # par(mfrow = c(2, 1))
    # hist(as.numeric(x[[1]]), breaks = 1000, main = p(asmaplgs[g], " c.x.a 65-1"))
    # hist(as.numeric(x[[2]]), breaks = 1000, main = p(asmaplgs[g], " c.x.a 65-3"))
    dev.off()
})

par(mfrow = c(2, 1))
hist(as.numeric(c651.hist), breaks = 1000, main = "slkdjf")
hist(as.numeric(c653.hist), breaks = 1000)

library(gridExtra)
grid.arrange(hist(as.numeric(c651.hist), breaks = 1000), hist(as.numeric(c653.hist), breaks = 1000))


#### DOUBLE HAPLOID ####


load("rotation1scripts_v4/saved.objects/all.m.8.iwgsc.4b.rev")
a.x.c.dh = all.m.8.iwgsc.4b.rev[[1]]

a.x.c.dh = read.csv("rotation1scripts_v4/processed_data/genotypes/avalon.x.cadenza/a.x.c.flipped.csv", stringsAsFactors = F)
a.x.c.dh2 = a.x.c.dh[, na.omit(match(colnames(c.x.af2v2), colnames(a.x.c.dh)))]

a.x.c.dh2[1, ] = as.character(c.x.af2v2[1, match(colnames(a.x.c.dh2), colnames(c.x.af2v2))])

v(a.x.c.dh2)

a.x.c.dh2.geno = makegenotypelist(a.x.c.dh2)
dh2.1 = makehistlist("chr5A", a.x.c.dh2.geno, "T")

f2.1 = makehistlist("chr5A", c651.geno, "T")

par(mfrow = c(2, 1))
hist(as.numeric(dh2.1), breaks = 1000)
hist(as.numeric(f2.1), breaks = 1000)

v(c.x.af2v2)

c.chr1 = c.x.af2v2[-1, which(c.x.af2v2[1, ] == "chr1A")]
apply(c.chr1, 1, function(x){
    g = x[1]
    pos = as.numeric()
    for(i in 1:length(x)){
        if(x[i] != g){
            pos = c(pos, i)
            g = x[i]
        }
    }
    pos
    browser()
})


#### RECOMBINATION STATS ####

names(c65.all) = asmaplgs

par(mfrow = c(2, 2))
lgnum = 6
hist(sort(c65.all[[lgnum]][[1]]), breaks = 100)
hist(sort(c65.all[[lgnum]][[2]]), breaks = 100)
hist(sort(c65.all[[lgnum]][[3]]), breaks = 100)
hist(sort(c65.all[[lgnum]][[4]]), breaks = 100)

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/recombination.distribution.test.6b")
recombination.distribution.test = all.sim.files
recombination.distribution.test6b = all.sim.files

process.sim.recomb = function(geno){
    rd1 = geno
    rd1 = as.data.frame(rbind(rd1[1, ], rd1))
    rd1 = convert.to.character.data.frame(rd1)
    rd1[1, ] = 1
    
    rd1 = cbind(rd1[, 1:2], rd1)
    rd1 = as.data.frame(rd1)
    rd1 = convert.to.character.data.frame(rd1)
    rd1[, 2] = rownames(rd1)
    rd1[, 1] = c("", 1:(nrow(rd1) - 1))
    
    rd1[1, 2] = ""
    rd1[2, 2] = "F2_1_1"
    
    rd1.geno = makegenotypelist(rd1)
    rd1.hist = sort(makehistlist(1, rd1.geno, F, F, F))
    rd1.hist
}

rd1 = process.sim.recomb(recombination.distribution.test6b[[1]])
rd2 = process.sim.recomb(recombination.distribution.test6b[[2]])
rd3 = process.sim.recomb(recombination.distribution.test6b[[3]])
rd4 = process.sim.recomb(recombination.distribution.test6b[[4]])

rd.all = list(rd1, rd2, rd3, rd4)

rd.all = lapply(recombination.distribution.test6b, process.sim.recomb)

hist(rd1, breaks = 100)
hist(rd2, breaks = 100)
hist(rd3, breaks = 100)
hist(rd4, breaks = 100)

#calculate distance of recombination positions from central marker
all.mann.test = lapply(list(rd.all), function(y){
    #get center from first population
    center = round(mean(y[[1]]))
    center = 100
    manntest1 = lapply(y, function(x){
        g = sort(x)
        
        abs(g - center)
        
    })
    manntest1
})


cxamann = lapply(all.mann.test, function(x){
    # browser()
    wilcox.test(x[[1]], x[[2]])
})

lapply(c65.all[[8]], function(x) list(mean(x), sd(x)))
tester123 = c65.all[[1]]
tester123 = lapply(tester123, function(x) sort(x))

plot(tester123[[1]], tester123[[2]])

#do some comparison of recombination profiles

#how many unique hotspot differences are there between the two populations?
lapply(c65.all, function(x) length(which(!unique(x[[1]]) %in% unique(x[[2]]))))


c65allstat = lapply(c65.all, function(x){
    
    data.frame(mean(x[[1]]), sd(x[[1]]), mean(x[[2]]), sd(x[[2]]))
    
    # list(mean(x[[1]]), mean(x[[2]]), sd(x[[1]]), sd(x[[2]]))
})

c65allstat = combine.list.of.data.frames(c65allstat)
rownames(c65allstat) = asmaplgs
c65allstat[] = lapply(c65allstat, round, digits = 2)


lapply(c65.all[[1]], function(x){
    g = sort(x)
    mean(g)
    g - 110
})

par(mfrow = c(1, 1))
boxplot(c65.all[[13]])

#calculate distance of recombination positions from central marker
all.mann.test = lapply(c65.all, function(y){
    #get center from first population
    center = round(mean(y[[1]]))
    manntest1 = lapply(y, function(x){
        g = sort(x)
        
        abs(g - center)
        
    })
    manntest1
})


cxamann = lapply(all.mann.test, function(x){
    # browser()
    wilcox.test(x[[1]], x[[2]])
})

cxamann.p = unlist(lapply(cxamann, function(x){
    x[3]
}))

names(cxamann.p) = multi.str.split(names(cxamann.p), "\\.", 1)
round(cxamann.p, 2)

3 / 15

#c x a unique events analysis
cxauniqueevents = unlist(lapply(c65.all, function(x){
    length(which(!unique(x[[1]]) %in% unique(x[[2]])))
}))


#apogee x paragon temp shift experiment - from recombination_analysisv3.R

apogeno1 = lapply(genotypelists, function(x){
    lapply(unique(x$lg), function(y){
        makehistlist(y, x, F)
    })
})

#analyse frequency of recombination events in AxP temp. shit populations
lapply(genotypelists[c(1, 3, 2, 4)], function(x){
    #analyse frequency of recombination events
    g = lapply(list.of.full.lgs.to.test, function(y){
        g = make.recombination.profiles(x, y, F)
        mean(unlist(lapply(g, length)))
    })
    
    g1 = c(mean(na.omit(unlist(g))), sd(na.omit(unlist(g))))
    names(g1) = c("mean", "sd")
    g1
})


apogeno1.2 = lapply(apogeno1, function(x) x[list.of.full.lgs.to.test])


apogenostat = lapply(apogeno1.2, function(x){
    g = lapply(x, function(y) data.frame(mean(y), sd(y)))
    combine.list.of.data.frames(g)
})

#put in order of temperature, q1 = 10(deg)C, q3 = 14(deg)C, q2 = 26(deg)C, q4 = 28(deg)C
apogeno1.2 = apogeno1.2[c(1, 3, 2, 4)]

apogeno1.3 = apogeno1.2[c(2, 3)]
apogeno1.4 = Map(function(x, y){
    list(x, y)
}, apogeno1.3[[1]], apogeno1.3[[2]])


#compare unique recombination events
apo.unique.events = unlist(lapply(apogeno1.4, function(x){
    length(which(!unique(x[[1]]) %in% unique(x[[2]])))
}))

names(apo.unique.events) = list.of.corresponding.chromosomes.to.lgs

all.mann.test2 = lapply(apogeno1.4, function(y){
    #get center from first population
    center = round(mean(y[[1]]))
    manntest1 = lapply(y, function(x){
        g = sort(x)
        
        abs(g - center)
        
    })
    manntest1
})


apomann = lapply(all.mann.test2, function(x){
    wilcox.test(x[[1]], x[[2]])
})

apomann.p = unlist(lapply(apomann, function(x){
    x[3]
}))

names(apomann.p) = list.of.corresponding.chromosomes.to.lgs
round(apomann.p, 2)

which(apomann.p < 0.05)

apogenostat2 = do.call(cbind, apogenostat)

apogenostat3 = apogenostat2[, c(3, 4, 5, 6)]

apogenostat3[] = lapply(apogenostat3, round, digits = 2)


list.of.corresponding.chromosomes.to.lgs = c("2A", "3A", "5A", "7B", "6B", "1A", "3B", "7A", "2D", "4A")
list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)

rownames(apogenostat3) = list.of.corresponding.chromosomes.to.lgs

g = apogenostat3$mean.y. - apogenostat3$mean.y..1
names(g) = list.of.corresponding.chromosomes.to.lgs
g2 = c65allstat$mean.x..1... - c65allstat$mean.x..2...
names(g2) = asmaplgs


par(mfrow = c(2, 1))
lapply(apogeno1.4, function(x){
    hist(x[[1]], breaks = 100)
    hist(x[[2]], breaks = 100)
})

temphist = apogeno1.4[[1]][[1]]

chromo1 = 1
pdf("rotation1scripts_v4/plots/c.x.af2/recombination/2abirminghamAxP.pdf", 20, 13)
grid.arrange(make.ggplot.histogram(apogeno1.4[[chromo1]][[1]], xlabel = "Locus", plot.title = "AxP 14(deg)C 2A"), make.ggplot.histogram(apogeno1.4[[chromo1]][[2]], xlabel = "Locus", plot.title = "AxP 26(deg)C"))
dev.off()





#### MISC ####

#checking proportion of cadenza alleles on chromosome 5B
g5 = c.x.af2v2[, which(c.x.af2v2[1, ] == "chr5B")]

g6 = unlist(lapply(c.x.af2v2, function(x){
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    a / b
}))

length(g6)
length(which(g6 < 1))
length(which(g6 > 1))



#testing magnitude of distortion vs. chi-square test 


g = c(rep("A", 34), rep("B", 20))
g1 = length(which(g == "A"))
g2 = length(which(g == "B"))
chisq.test(c(g1, g2))

chisq.test(c(10000, 10034))

((10034 - 10017)^2) / 10017 + ((10000 - 10017)^2) / 10017
289 / 10017

chisq.test(c(34, 20))

((34 - 27)^2) / 27

lapply(list.of.geno.dfs.each.pop, function(x){
    x[, which(colnames(x) == "AX.94789875")]
})




v(ratiodfs2.3[[2]][[3]])

plot(ratiodfs2.3[[2]][[3]]$cm[100:130], ratiodfs2.3[[2]][[3]]$seg.ratio1[100:130])


max.distortion = function(chromosome1){
    #find the four markers exhibiting the highest segregation distortion in the reciprocal c.x.a crosses for a particular chromosome. shows only the markers that are in the top 4 for more than one of the reciprocal cross populations
    #args:
    # chromosome1 - a string specifying the chromosome to search, e.g. "1A"
    chromosome1 = which(chromos == chromosome1)
    top.markers = lapply(1:4, function(x){
        coords1 = sort(abs(ratiodfs2.3[[chromosome1]][[x]]$seg.ratio1 - 0.5), decreasing = T, index.return = T)$ix
        ratiodfs2.3[[chromosome1]][[x]][coords1[1:4], ]
    })
    
    # browser()
    
    top.markers2 = lapply(top.markers, function(x){
        x$marker = rownames(x)
        x$dist = abs(x$seg.ratio1 - 0.5)
        x[which(x$seg.chi1 < 0.05), ]
    })
    
    top.markers3 = combine.list.of.data.frames(top.markers2)
    sorted.markers = sort(table(top.markers3$marker), decreasing = T)
    sorted.markers2 = sorted.markers[which(sorted.markers > 1)]    
    names(sorted.markers2)
}

max.distortion("1A")

distorters = sapply(chromos, max.distortion)

#make clusterplot directories for all of the segregation distortion markers
lapply(unlist(distorters), function(x){
    dir.create(p("rotation1scripts_v4/plots/clusterplots/", x))
})

#move the files generated by the autohotkey script aas_win_activate.ahk to their respective directories
lapply(unlist(distorters), function(x){
    basepath = "rotation1scripts_v4/plots/clusterplots/"
    files1 = list.files("rotation1scripts_v4/plots/clusterplots/")
    # browser()
    g = strsplit(x, "\\.")
    files2 = files1[grep(g[[1]][2], files1)]
    dir1 = files2[grep("AX", files2)]
    files2 = files2[-grep("AX", files2)]
    
    lapply(files2, function(y){
        file.rename(p(basepath, y), p(basepath, dir1, "/", y))
    })    
})


# basepath = "rotation1scripts_v4/plots/clusterplots/"
# lapply(list.files(basepath), function(x){
#     files1 = list.files(p(basepath, x))
#     file.remove(p(basepath, x, "/", files1[grep("combined", files1)]))
# })

#combine all four png files for each marker into a single png
library(png)
library(grid)
library(gridExtra)
basepath = "rotation1scripts_v4/plots/clusterplots/"
lapply(list.files(basepath), function(x){
    plots1 = lapply(list.files(p(basepath, x)), function(y){
        pic1 = readPNG(p(basepath, x, "/", y))
        rasterGrob(pic1)
    })
    # browser()
    png(p(basepath, x, "combined.png"), 1000, 1000)
    do.call(grid.arrange, plots1)
    dev.off()
})



unlist(distorters)

ratiodfs2.3[[15]][[2]]$seg.ratio1[37]

ratiodfs2.3[[15]][[2]][37, ]

ratiodfs2.3[[2]][[1]][119, ]
ratiodfs2.3[[2]][[4]][119, ]

test1 = lapply(ratiodfs2.3[[1]], function(x){
    rownames(x)
})

lapply(1:4, function(x){
    all(test1[[1]] == test1[[x]])    
})

#### F2 CAKE PLOT ####

load("rotation1scripts_v4/saved.objects/recomb.sim/recomb.profile1a.no.selec.sims")
library(reshape2)

cake1 = make.cake.plot(recomb.profile1a.no.selec.sims[[12]])


cake.max.dist.sim = make.cake.plot(list.of.first.gen.recomb1a.pop96[[950]])
cake.min.dist.sim = make.cake.plot(list.of.first.gen.recomb1a.pop96[[857]])
cake.real = make.cake.plot(realdat1.geno)

grid.arrange(cake.max.dist.sim, cake.min.dist.sim, cake.real, plot4, plot5, ncol = 1)

generate.cm.for.sim.data(recomb.profile1a)


#test without alpha desgination of significance
ggplot(s2, aes(x = marker, y = count, fill = Genotype)) + geom_bar(stat = "identity", width = 1) +
    theme_classic() + scale_fill_manual(values = c("#EE6C4D", "#231942", "#98C1D9")) + guides(alpha = F) + 
    geom_text(aes(label = sig.text), nudge_y = 210)



c("#EE6C4D", "#E8EDDF", "#98C1D9")

#### COMPARE REAL DATA TO SIMULATED DATA ####

# load("rotation1scripts_v4/saved.objects/recomb.sim/recomb.profile1a.no.selec.sims")
# 
# load("rotation1scripts_v4/saved.objects/recomb.sim/fgen/f6.no.selec.300pop")
# 
# load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb1a.pop96")
# 
# load("rotation1scripts_v4/saved.objects/recomb.sim/fgen/f6.no.selec.96pop")

source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

generate.cm.for.sim.data = function(sim1.rqtl, fgen){
    # browser()
    #prepare data.frame for use with ASMap
    
    if(missing(fgen)) fgen = 2
    
    sim1.rqtl2 = sim1.rqtl
    sim1.rqtl2 = convert.to.character.data.frame(sim1.rqtl2)
    sim1.rqtl2[sim1.rqtl2 == "H"] = "X"
    sim1.rqtl2 = as.data.frame(t(sim1.rqtl2))
    sim1.rqtl2 = convert.to.character.data.frame(sim1.rqtl2)
    sim1.rqtl3 = mstmap.data.frame(sim1.rqtl2, p("RIL", fgen), "kosambi", p.value = 4)
    
    
    sim1.cm = pull.map(sim1.rqtl3, as.table = T)
    # browser()
    rownames(sim1.cm) = as.numeric(multi.str.split(rownames(sim1.cm), "marker", 2))
    #fix internal order of bins for each skeleton marker
    
    sim1.cm = sim1.cm[sort(as.numeric(rownames(sim1.cm)), decreasing = T, index.return = T)$ix, ]
    
    sim1.cm = sim1.cm[nrow(sim1.cm):1, ]
    
    #reverse the cM order 
    
    sim1.cm$pos = reverse.cm.distances(sim1.cm$pos)
    
    sim1.cm
    
    
}

prepare.sim.data.for.plot = function(sim1, fgen){
    # browser()
    sim1.seg = convert.recomb.to.seg(sim1)
    
    if(missing(fgen)) fgen = 2
    
    if(fgen == 2){
        sim1.seg.chi = convert.recomb.to.seg.f2(sim1)    
    } else {
        sim1.seg.chi = convert.recomb.to.seg2(sim1)
    }
    
    sim1.seg.chi = unlist(lapply(sim1.seg.chi, function(x) x[[3]]))
    
    #remove final data point to match length of real 1A data
    sim1.seg = sim1.seg[-225]
    sim1.seg.chi = sim1.seg.chi[-225]
    
    sim1.dat = make.ratio.df(sim1.seg, sim1.seg.chi)
    sim.cm = generate.cm.for.sim.data(sim1, fgen)
    sim1.dat$cm = sim.cm$pos
    sim1.dat
}

generate.cm.for.sim.data(recomb.profile1a.no.selec.sims[[4]])

plots5 = make.plots(recomb.profile1a.no.selec.sims)


# sim1.dat = prepare.sim.data.for.plot(recomb.profile1a.no.selec.sims[[1]])
# sim2.dat = prepare.sim.data.for.plot(recomb.profile1a.no.selec.sims[[7]])
# sim3.dat = prepare.sim.data.for.plot(recomb.profile1a.no.selec.sims[[12]])

#load pedigreeSim data
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop96")

axcf2.noselec.pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6.noselec.pop96")

axcf6.noselec.pop96 = all.sim.files

min.pvalues.f2 = lapply(axcf2.noselec.pop96, function(x){
    min(convert.recomb.to.seg.f2(x, T))
})

f2.sim.max.distortion = axcf2.noselec.pop96[[which(unlist(min.pvalues.f2) == min(unlist(min.pvalues.f2)))]]
f2.sim.min.distortion = axcf2.noselec.pop96[[which(unlist(min.pvalues.f2) == max(unlist(min.pvalues.f2)))]]


sim1.min = prepare.sim.data.for.plot(f2.sim.min.distortion, 2)
sim1.max = prepare.sim.data.for.plot(f2.sim.max.distortion, 2)

realdat1 = ratiodfs2.3[[15]][[1]]

realdat1.geno = geno.dfs.2[[1]][[15]]
realdat1.geno = realdat1.geno[-(1:3), -1]

min.pvalues.f6 = lapply(axcf6.noselec.pop96, function(x){
    min(convert.recomb.to.seg2(x, T))
})

f6.sim.max.distortion = axcf6.noselec.pop96[[which(unlist(min.pvalues.f6) == min(unlist(min.pvalues.f6)))]]
f6.sim.min.distortion = axcf6.noselec.pop96[[which(unlist(min.pvalues.f6) == max(unlist(min.pvalues.f6)))]]

sim1.f6.min = prepare.sim.data.for.plot(f6.sim.min.distortion, 6)
sim1.f6.max = prepare.sim.data.for.plot(f6.sim.max.distortion, 6)

yscale1 = c(0.35, 0.65)
yscale1 = c(0.25, 0.75)

plotter1 = function(plotdata, title){
    ggplot(plotdata, aes(x = cm, y = seg.ratio1)) + 
        geom_path() + geom_point(aes(color = chi1t, shape = chi1t), size = 1) + geom_hline(yintercept = 0.5) + coord_cartesian(ylim = yscale1) + ggtitle(title) + 
        theme_classic() + 
        scale_color_manual(values=c("#492de5", "#999999", "#E69F00", "#56B4E9"), labels = c("N.S.", "0.05", "0.01", "0.001")) +
        scale_shape_manual(values = c(16, 15, 17, 18), labels = c("N.S.", "0.05", "0.01", "0.001")) +
        xlab("Genetic distance (cM)") + ylab("Segregation ratio") +
        guides(color = guide_legend(title = "Sig."), shape = guide_legend(title = "Sig.")) #+
        #geom_hline(yintercept = 0.44) + geom_hline(yintercept = 0.56)
}

srdat1 = list(realdat1[, which(colnames(realdat1) %in% c("seg.ratio1", "seg.chi1", "axis", "chi1t", "cm", "num1"))], sim1.min, sim1.max, sim1.f6.min, sim1.f6.max)

srdat2 = Map(function(x, y){
    x$num1 = y
    x
}, srdat1, c("(a)", "(b)", "(c)", "(d)", "(e)"))

srdat3 = combine.list.of.data.frames(srdat2)

plotter1.5 = function(plotdata, title){
    ggplot(plotdata, aes(x = cm, y = seg.ratio1)) + 
        geom_path() + geom_point(aes(color = chi1t, shape = chi1t), size = 1.3) + geom_hline(yintercept = 0.5) + coord_cartesian(ylim = yscale1) + ggtitle(title) + 
        theme_bw() + 
        scale_color_manual(values=c("#492de5", "#999999", "#E69F00", "#56B4E9"), labels = c("N.S.", "0.05", "0.01", "0.001")) +
        scale_shape_manual(values = c(16, 15, 17, 18), labels = c("N.S.", "0.05", "0.01", "0.001")) +
        xlab("Genetic distance (cM)") + ylab("Segregation ratio") +
        guides(color = guide_legend(title = "Sig."), shape = guide_legend(title = "Sig.")) + facet_grid(num1 ~ .) +
        theme(strip.text.y = element_text(angle = 0, size = 12), strip.background = element_rect(color = "white", fill = "#ffffff"), strip.placement = "outside",
                    panel.grid = element_blank())
            #+
    #geom_hline(yintercept = 0.44) + geom_hline(yintercept = 0.56)
}

#### sim vs real seg PUBLICATION FIGURE #### 
seg.comp.plot = plotter1.5(srdat3, "")
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig2v2_plos.eps", seg.comp.plot, width = 11, height = 19, units = "cm", device = cairo_pdf)

plot1 = plotter1(realdat1, "Real Data 1A")
plot2 = plotter1(sim1.min, "Simulated 1A Min")
plot3 = plotter1(sim1.max, "Simulated 1A Max")

plot4 = plotter1(sim1.f6.min, "Simulated 1A F6 Min")
plot5 = plotter1(sim1.f6.max, "Simulated 1A F6 Max")






# plots.combined = grid.arrange(plot1, plot2, plot3, plot4, plot5, ncol = 1)

compgrob = arrangeGrob(plot1, plot2, plot3, plot4, plot5, ncol = 1)
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/real.sim.seg.comparison.eps", compgrob, width = 8.4, height = 23.4, units = "cm", device = cairo_pdf)


pdf("rotation1scripts_v4/plots/simulation/pedsim/paper.comp.pdf", 5, 10)
grid.arrange(plot1, plot2, plot3, plot4, plot5, ncol = 1)
dev.off()





pdf("rotation1scripts_v4/plots/simulation/real.sim.comparison.pdf", 15, 10)
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
dev.off()

load("rotation1scripts_v4/saved.objects/list.of.no.selec2")

no.selec.sim.dat = lapply(list.of.no.selec2, function(x){
    prepare.sim.data.for.plot(x[[6]])
})

yscale1 = c(0.33, 0.67)
p1 = Map(function(x, y){
    plotter1(x, title = y)
}, no.selec.sim.dat, c("Population size: 96",
                                             "Population size: 300",
                                             "Population size: 1000",
                                             "Population size: 10000"))

pdf("rotation1scripts_v4/plots/simulation/pop.size.no.selection.comparison.pdf", 10, 10)
do.call(grid.arrange, p1)
dev.off()




#### HEATMAP FOR SIM. DATA #### 
gen.asmap.for.sim = function(sim.dat1, p.value1){
    if(missing(p.value1)) p.value1 = 4
    sim1.rqtl2 = sim.dat1
    sim1.rqtl2[sim1.rqtl2 == "H"] = "X"
    sim1.rqtl2 = as.data.frame(t(sim1.rqtl2))
    sim1.rqtl2 = convert.to.character.data.frame(sim1.rqtl2)
    sim1.rqtl3 = mstmap.data.frame(sim1.rqtl2, "RIL2", "kosambi", p.value = p.value1)
    sim1.rqtl3
}

#make heatmap of simulated data for rf and lod scores 
asmap1 = gen.asmap.for.sim(list.of.first.gen.recomb1a.pop96[[100]])

g = convert.sim.data.to.rqtl.format(list.of.first.gen.recomb1a.pop96[[100]])
write.csv(g, "rotation1scripts_v4/temp/test1.csv", row.names = F)
g1 = rqtl.read("rotation1scripts_v4/temp/test1.csv")
# g2 = est.rf(g1)
est.map(g1, map.function = "haldane")
pull.map(g1, as.table = T)

heatMap(g1, what = "rf")


pull.map(asmap1, as.table = T)


asmap.rf = est.rf(asmap1)
asmap.rf = pull.rf(asmap.rf)
asmap.rf = as.data.frame(asmap.rf)
asmap.rf[] = lapply(asmap.rf, as.numeric)
sapply(as.data.frame(asmap.rf), function(x) max(na.omit(x)))

heatMap(asmap1)




#### EXAMINE RECOMBINATION DATA #### 

#genetic map of only cxa651
# cxa651redonewchromos = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa651_redone_w_chromos.csv", stringsAsFactors = F, header = T)
#genetic map of both cxa651 and cxa653
cxa651redonewchromos = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.af2_onlycxa.geneticmaporder.csv", stringsAsFactors = F, header = T)
cxa651.geno1 = makegenotypelist(cxa651redonewchromos)
cxa651.hist = makehistlist("1A", cxa651.geno1, F)
cxa651.hist2 = makehistlist("1A", cxa651.geno1, T)

real.tdat = unlist(lapply(unique(cxa651.hist2$individual), function(x){
    nrow(cxa651.hist2[which(cxa651.hist2$individual == x), ])
}))


convert.sim.data.to.rqtl.format = function(sim1){
    # browser()
    sim2 = convert.to.character.data.frame(sim1)
    sim2 = rbind(sim2[1, ], sim2)
    sim2[1, ] = "1"
    
    colnames(sim2) = paste0("marker", 1:ncol(sim2))
    sim2 = cbind(sim2[, 1], sim2)
    
    sim2 = convert.to.character.data.frame(sim2)
    sim2[2:nrow(sim2), 1] = paste0("ind", 1:(nrow(sim2) - 1))
    sim2[1, 1] = ""
    
    sim2 = cbind(sim2[, 1], sim2)
    sim2 = convert.to.character.data.frame(sim2)
    sim2[2:nrow(sim2), 1] = 1:(nrow(sim2) - 1)
    colnames(sim2)[1:2] = c("v1", "v2")
    
    sim2
}




sim.dat96pop1 = list.of.first.gen.recomb1a.pop96[[1]]
sim.dat.96pop2 = convert.sim.data.to.rqtl.format(sim.dat96pop1)

sim.geno1 = makegenotypelist(sim.dat.96pop2)

sim.geno2 = makehistlist(1, sim.geno1, F)
sim.geno3 = makehistlist(1, sim.geno1, T)

sim.tdat = unlist(lapply(unique(sim.geno3$individual), function(x){
    nrow(sim.geno3[which(sim.geno3$individual == x), ])
}))

t.test(real.tdat, sim.tdat)

wilcox.test(real.tdat, sim.tdat)


par(mfrow = c(2, 1))
hist(cxa651.hist, breaks = seq(0, 250, 5))
hist(sim.geno2, breaks = seq(0, 250, 5))



hist1 = make.ggplot.histogram(cxa651.hist, breaks = seq(0, 220, 2), xlabel = "Locus", ylabel = "Number of Recombination Events", plot.title = "Real Data Recombination Distribution") + theme_classic()
hist2 = make.ggplot.histogram(sim.geno2, breaks = seq(0, 220, 2), xlabel = "Locus", ylabel = "Number of Recombination Events", plot.title    = "Simulated Data Recombination Distribution") + theme_classic()

grid.arrange(hist1, hist2)


pdf("rotation1scripts_v4/plots/simulation/recomb.distribution.real.vs.simulated.data.pdf", 10, 5)
grid.arrange(hist1, hist2)
dev.off()


load("rotation1scripts_v4/simulation/pedigreesim/simresults/axcf2noselecpop300.het.test")


# validation - test het levels over subsequent filial generations

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2noselecpop300.het.test")
axcf2noselecpop300.het.test = all.sim.files
axcf2noselecpop300.het.test.t = sapply(axcf2noselecpop300.het.test.t, evaluate.hets)
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf3noselecpop300.het.test")
axcf3noselecpop300.het.test = all.sim.files
axcf3noselecpop300.het.test.t = sapply(axcf3noselecpop300.het.test.t, evaluate.hets)
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf4noselecpop300.het.test")
axcf4noselecpop300.het.test = all.sim.files
axcf4noselecpop300.het.test.t = sapply(axcf4noselecpop300.het.test.t, evaluate.hets)
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf5noselecpop300.het.test")
axcf5noselecpop300.het.test = all.sim.files
axcf5noselecpop300.het.test.t = sapply(axcf5noselecpop300.het.test.t, evaluate.hets)
load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf6noselecpop300.het.test")
axcf6noselecpop300.het.test = all.sim.files
axcf6noselecpop300.het.test.t = sapply(axcf6noselecpop300.het.test.t, evaluate.hets)


#### countXO analysis####
#compare recombination events in real data to simulated data
cross1 = read.cross("csv", "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/", "cxa651_redone_w_chromos.csv", estimate.map = F)

crossovers = countXO(cross1, chr = "1A")

sim1 = list.of.first.gen.recomb1a.pop96[[1]]
sim1 = convert.sim.data.to.rqtl.format(sim1)

write.csv(sim1, "rotation1scripts_v4/processed_data/genotypes/simulated/sim1.rqtl.csv", row.names = F)

cross2 = read.cross("csv", "rotation1scripts_v4/processed_data/genotypes/simulated/", "sim1.rqtl.csv", estimate.map = F)
crossovers2 = countXO(cross2)


par(mfrow = c(2, 1))
hist(crossovers, breaks = seq(-1, 30, 1), ylim = c(0, 60))
hist(crossovers2, breaks = seq(-1, 30, 1), ylim = c(0, 60))

make.ggplot.histogram(crossovers, breaks = seq(-1, 30, 1), ylab = "count")
make.ggplot.histogram(crossovers2, breaks = seq(-1, 30, 1), ylab = "count")

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop96")
axcf2.noselec.pop96 = all.sim.files


g = lapply(axcf2.noselec.pop96, function(x){
    g1 = convert.to.character.data.frame(x)
    gen.asmap.for.sim(g1)
    })


g2 = sapply(g, function(x){
    q = pull.map(x, as.table = T)
    max(q$pos)
})


axc.psim1 = convert.sim.data.to.rqtl.format(axcf2.noselec.pop96[[1]])

write.csv(axc.psim1, "rotation1scripts_v4/processed_data/genotypes/simulated/pedsim1.rqtl.csv", row.names = F)

cross3 = rqtl.read("rotation1scripts_v4/processed_data/genotypes/simulated/pedsim1.rqtl.csv")

crossovers3 = countXO(cross3)

mean(crossovers3)

sd(crossovers3)

crossovers
wilcox.test(crossovers, crossovers3)

library(ASMap)

#prepare data.frame for use with ASMap
sim1.rqtl2 = sim1.rqtl
sim1.rqtl2[sim1.rqtl2 == "H"] = "X"
sim1.rqtl2 = sim1.rqtl2[-1, ]
sim1.rqtl2 = as.data.frame(t(sim1.rqtl2))
sim1.rqtl2 = convert.to.character.data.frame(sim1.rqtl2)
sim1.rqtl2 = sim1.rqtl2[-1, ]
sim1.rqtl3 = mstmap.data.frame(sim1.rqtl2, "RIL2", "kosambi", p.value = 4)

sim1.cm = pull.map(sim1.rqtl3, as.table = T)
sim1.cm = as.numeric(sim1.cm)

sim1.cm$pos
reverse.cm.distances(sim1.cm$pos)

rownames(sim1.cm) = as.numeric(multi.str.split(rownames(sim1.cm), "V", 2))
#fix internal order of bins for each skeleton marker
sim1.cm = sim1.cm[sort(as.numeric(rownames(sim1.cm)), decreasing = T, index.return = T)$ix, ]
sim1.cm = sim1.cm[nrow(sim1.cm):1, ]
#reverse the cM order 
sim1.cm$pos = reverse.cm.distances(sim1.cm$pos)



#recombination profile for chinese spring x paragon
#NB. This won't work, see note "Comparing simulated to real data in OneNote (22/03/2019)"
library(readxl)
setwd("E:/phd.project.main/")
cs.x.p.map = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic maps.xlsx", sheet = 6)
cs.x.p.map.list = split.df.into.list.of.dfs.by.column(cs.x.p.map, "chr")
cs.x.p1a.phys = genetictophysical("1A", cs.x.p.map.list[[1]]$marker, "-")

cs.x.p.map.1a = cs.x.p.map.list[[1]]

#get chinese spring x paragon map in an appropriate format for makegenotypelist()
cs = cs.x.p.map.1a
cs = as.data.frame(t(cs))
cs = cbind(cs[, 1], cs)
cs[, 1] = rownames(cs)
cs = cs[-3, ]
cs = cbind(cs[, 1], cs)
cs[, 1] = 1:nrow(cs)
cs[1:2, 1] = ""
cs[2, 2] = ""
cs = convert.to.character.data.frame(cs)
colnames(cs) = cs[1, ]
cs = cs[-1, ]

cs.geno1 = makegenotypelist(cs)
library(tibble)
library(dplyr)
cs.recomb1 = make.recombination.profiles(cs.geno1, "1A")
cs.recomb1.full = make.recombination.profiles(cs.geno1, "1A", F)



#### NO SELEC. COMPARISON SIMULATED DATA ####

# load("rotation1scripts_v4/saved.objects/pop.stats")

#get pedigreesim pop stats made with 15042019.make.seg.dist.pedsim.plot.R
pop.stats2 = read.csv("rotation1scripts_v4/temp/all.stats2.csv", stringsAsFactors = F)

pop.stats2 = as.data.frame(pop.stats2)

pop.stats2[] = lapply(pop.stats2, as.numeric)

# v(pop.stats2)
# 
plotter2 = function(pop.data, means, sds, title, ylimits, hline, bounds, facets1){
    if(missing(ylimits)) ylimits = c(0.4, 0.6)
    if(missing(hline)) hline = F
    if(missing(bounds)) bounds = 0
    if(missing(facets1)) facets1 = F
    # browser()
    if(facets1 == T){
        lower.bounds1 = data.frame(c(96, 300, 1000, 10000), c((36 / (36 + 56)), (128 / (128 + 163)), (453 / (453 + 514)), 4751 / (4944 + 4751)))
        colnames(lower.bounds1) = c("all1", "lower_bound")
        upper.bounds1 = data.frame(c(96, 300, 1000, 10000), c((56 / (36 + 56)), (163 / (128 + 163)), (514 / (453 + 514)), 4944 / (4944 + 4751)))
        colnames(upper.bounds1) = c("all1", "upper_bound")    
        
        #facets plot1
        plot1 = ggplot(data = pop.data, aes(x = markernum, y = means)) + geom_line(size = 1) +
            geom_ribbon(aes(ymin = (means - sds), ymax = (means + sds)), alpha = 0.25, linetype = 3) +
            coord_cartesian(ylim = ylimits) + xlab("Marker Number") + ylab("Mean segregation ratio") +
            ggtitle(title) + theme_bw() + theme(axis.text.x = element_text(angle = 45)) + facet_grid(. ~ all1) +
            geom_hline(data = lower.bounds1, aes(yintercept = lower_bound), linetype = 2) +
            geom_hline(data = upper.bounds1, aes(yintercept = upper_bound), linetype = 2) +
            geom_hline(yintercept = 0.5, linetype = 3) +
            theme(strip.text = element_text(size = 12), strip.background = element_rect(colour = "white", fill = "#ffffff"),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")
        
        plot1
    } else {
        
        #old plot1
        plot1 = ggplot(data = pop.data, aes(x = 1:nrow(pop.data), y = means)) + geom_line(size = 1) +
            geom_ribbon(aes(ymin = (means - sds), ymax = (means + sds)), alpha = 0.25, linetype = 3) +
            coord_cartesian(ylim = ylimits) + xlab("Marker Number") + ylab("Mean segregation ratio") +
            ggtitle(title) + theme_classic() + theme(axis.text.x = element_text(angle = 45))
        
        
        
        
        
        if(bounds != 0){
            plot1 = plot1 + 
                geom_hline(yintercept = bounds[1] / (bounds[1] + bounds[2]), linetype = 2) + 
                geom_hline(yintercept = bounds[2] / (bounds[1] + bounds[2]), linetype = 2)
        } 
        if(hline == T) plot1 = plot1 + geom_hline(yintercept = 0.5, linetype = 3)
        
        plot1
    }
}

#### F6 SELECTION PLOT ####

plot96 = plotter2(pop.stats2, pop.stats2$axcf6selec20.pop96.seg2.mean, pop.stats2$axcf6selec20.pop96.seg2.sd, "Pop.: 96", bounds = c(36, 56), hline = T, ylimits = c(0.39, 0.61)) 
plot97 = plotter2(pop.stats2, pop.stats2$axcf6selec20.pop300.seg2.mean, pop.stats2$axcf6selec20.pop300.seg2.sd, "Pop.: 300", bounds = c(128, 163), hline = T, ylimits = c(0.39, 0.61))
plot98 = plotter2(pop.stats2, pop.stats2$axcf6selec20.pop1000.seg2.mean, pop.stats2$axcf6selec20.pop1000.seg2.sd, "Pop.: 1000", bounds = c(453, 514), hline = T, ylimits = c(0.39, 0.61))
plot99 = plotter2(pop.stats2, pop.stats2$axcf6selec20.pop10000.seg2.mean, pop.stats2$axcf6selec20.pop10000.seg2.sd, "Pop.: 10000", bounds = c(4944, 4751), hline = T, ylimits = c(0.39, 0.61))

v(pop.stats2)
all.pop20 = list(pop.stats2[, c(1, 5)], pop.stats2[, c(2, 6)], pop.stats2[, c(3, 7)], pop.stats2[, c(4, 8)])

all.pop20 = lapply(all.pop20, function(x){
    colnames(x) = c("mean", "sd")
    x
})

all.pop201 = Map(function(x, y){
    x$all1 = y
    x$markernum = 1:nrow(x)
    x
}, all.pop20, c(96, 300, 1000, 10000))

all.pop202 = combine.list.of.data.frames(all.pop201)

all.pop.selec.20.plot = plotter2(all.pop202, all.pop202$mean, all.pop202$sd, "", facets1 = T) + ggtitle("(a)") 
all.pop.selec.20.plot

grid.arrange(plot96, plot97, plot98, plot99)

pdf("rotation1scripts_v4/plots/simulation/pop.size.noise.comparison.pdf", 7, 7)
grid.arrange(plot96, plot97, plot98, plot99)
dev.off()


load("rotation1scripts_v4/saved.objects/all.pop.peak.pedsim")



#### F6 NO SELEC PLOT ####

noselec.data = read.csv("rotation1scripts_v4/temp/noselec.data.csv", stringsAsFactors = F)

p1 = plotter2(noselec.data, noselec.data$axcf6.noselec.pop96.seg2.mean, pop.stats2$axcf6selec20.pop96.seg2.sd, "Population Size: 96", bounds = c(36, 56), hline = T, ylimits = c(0.39, 0.61))
p2 = plotter2(noselec.data, noselec.data$axcf6.noselec.pop300.seg2.mean, pop.stats2$axcf6selec20.pop300.seg2.sd, "Population Size: 300", bounds = c(128, 163), hline = T, ylimits = c(0.39, 0.61))
p3 = plotter2(noselec.data, noselec.data$axcf6.noselec.pop1000.seg2.mean, pop.stats2$axcf6selec20.pop1000.seg2.sd, "Population Size: 1000", bounds = c(453, 514), hline = T, ylimits = c(0.39, 0.61))
p4 = plotter2(noselec.data, noselec.data$axcf6.noselec.pop10000.seg2.mean, pop.stats2$axcf6selec20.pop10000.seg2.sd, "Population Size: 10000", bounds = c(4944, 4751), hline = T, ylimits = c(0.39, 0.61))


nsdat1 = noselec.data[, c(1, 5)]
nsdat2 = noselec.data[, c(2, 6)]
nsdat3 = noselec.data[, c(3, 7)]
nsdat4 = noselec.data[, c(4, 8)]

nsdatall = list(nsdat1, nsdat2, nsdat3, nsdat4)

nsdatall = lapply(nsdatall, function(x){
    colnames(x) = c("mean", "sd")
    x
})

nsdatall1 = Map(function(x, y){
    x$all1 = y
    x$markernum = 1:nrow(x)
    x
}, nsdatall, c(96, 300, 1000, 10000))

nsdatall2 = combine.list.of.data.frames(nsdatall1)

#### PUBLICATION FIG 3 ####
fig3.facets = plotter2(nsdatall2, nsdatall2$mean, nsdatall2$sd, "", bounds = c(36, 56), hline = T, ylimits = c(0.39, 0.61), T)
fig3.facets

ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig3v2.eps", fig3.facets, width = 17.4, height = 8, units = "cm", device = cairo_pdf)


pdf("rotation1scripts_v4/plots/simulation/pedsim/f6.noselec.seg.pdf", 8, 8)
grid.arrange(p1, p2, p3, p4)
dev.off()


figure3 = arrangeGrob(p1, p2, p3, p4)
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/fig3.eps", figure3, width = 17.4, height = 17.4, units = "cm", device = cairo_pdf)


#### ALL SIMS IN ONE DISTORTION PLOT ####

load("rotation1scripts_v4/saved.objects/pop96.anal2")
pop96.anal2 = reset.colnames(pop96.anal2)
pop96.anal2.melt = melt(pop96.anal2)
pop96.anal2.melt = convert.to.character.data.frame(pop96.anal2.melt)
pop96.anal2.melt$xaxis = as.numeric(multi.str.split(pop96.anal2.melt$variable, "V", 2))
pop96.anal2.melt$group1 = as.character(rep(1:1000, 225))
pop96.anal2.melt$value = as.numeric(pop96.anal2.melt$value)

plot1 = ggplot(pop96.anal2.melt, aes(x = xaxis, y = value, group = group1)) + geom_line() + theme_classic()
plot1 + geom_hline(yintercept = (36 / (36 + 56)), color = "red", size = 2) + geom_hline(yintercept = (56 / (36 + 56)), color = "red", size = 2)







#### PEAK OF DIST. HISTOGRAMS ####
pop96.anal3 = as.data.frame(t(pop96.anal2))
nrow(pop96.anal3)

hist(sapply(pop96.anal3, function(x){
    x2 = abs(x - 0.5)
    g = which(x2 == max(x2))
    g[ceiling(length(g) / 2)]
    # if(length(g[(length(g) / 2)]) == 0)
}))


load("rotation1scripts_v4/saved.objects/all.pop.peak.pedsim")
library(scales)
peak1 = make.ggplot.histogram(all.pop.peak.pedsim[[1]], breaks = seq(1, 225, 1), ylabel = "Count") + 
    theme_classic() + coord_cartesian(ylim = c(0, 800)) + scale_y_sqrt(breaks = seq(50, 850, 100), labels = function(x) paste0(x/10, "%")) +
    xlab("Marker number") + ylab("% Sim. w/ peak of dist. at X") +    theme(axis.text.x = element_text(angle = 45))
peak2 = make.ggplot.histogram(all.pop.peak.pedsim[[2]], breaks = seq(1, 225, 1), ylabel = "Count") + 
    theme_classic() + coord_cartesian(ylim = c(0, 800)) + scale_y_sqrt(breaks = seq(50, 850, 100), labels = function(x) paste0(x/10, "%")) +
    xlab("Marker number") + ylab("% Sim. w/ peak of dist. at X") +    theme(axis.text.x = element_text(angle = 45))
peak3 = make.ggplot.histogram(all.pop.peak.pedsim[[3]], breaks = seq(1, 225, 1), ylabel = "Count") + 
    theme_classic() + coord_cartesian(ylim = c(0, 800)) + scale_y_sqrt(breaks = seq(50, 850, 100), labels = function(x) paste0(x/10, "%")) +
    xlab("Marker number") + ylab("% Sim. w/ peak of dist. at X") +    theme(axis.text.x = element_text(angle = 45))
peak4 = make.ggplot.histogram(all.pop.peak.pedsim[[4]], breaks = seq(1, 225, 1), ylabel = "Count") + 
    theme_classic() + coord_cartesian(ylim = c(0, 800)) + scale_y_sqrt(breaks = seq(50, 850, 100), labels = function(x) paste0(x/10, "%")) +
    xlab("Marker number") + ylab("% Sim. w/ peak of dist at X") +    theme(axis.text.x = element_text(angle = 45))



all.pop.peak.pedsim2 = Map(function(x, y){
    data.frame(x, y)
}, all.pop.peak.pedsim, c(96, 300, 1000, 10000))
all.pop.peak.pedsim3 = combine.list.of.data.frames(all.pop.peak.pedsim2)

histograms1 = ggplot(all.pop.peak.pedsim3, aes(all.pop.peak.pedsim3$x)) + geom_histogram() + facet_grid(. ~ y) + theme_bw() +
    theme(strip.text = element_text(size = 12), strip.background = element_rect(colour = "white", fill = "#ffffff"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Marker Number") + ylab("Num. sim. w/ peak dist. at marker") + 
    ggtitle("(b)") + theme(axis.text.x = element_text(angle = 45), strip.placement = "outside")

grid.arrange(all.pop.selec.20.plot, histograms1)

#### PUBLICATION FIGURE 6 ####
fig6 = arrangeGrob(all.pop.selec.20.plot, histograms1)
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig6v2.eps", fig6, width = 17.4, height = 15, units = "cm", device = cairo_pdf)



grid.arrange(peak1, peak2, peak3, peak4)

pdf("rotation1scripts_v4/plots/simulation/f6.distortion.w.selec.peak.dist.hist.pdf", 8, 8)
grid.arrange(peak1, peak2, peak3, peak4)
dev.off()

pdf("rotation1scripts_v4/plots/simulation/pedsim/f6.distortion.w.selec.seg.and.peak.dist.hist.pdf", 18, 9)
grid.arrange(plot96, plot97, plot98, plot99, peak1, peak2, peak3, peak4, ncol = 4)
dev.off()

grid.arrange(all.pop.selec.20.plot, peak1, peak2, peak3, peak4, ncol = 4)




twopanelplot = arrangeGrob(plot96, plot97, plot98, plot99, peak1, peak2, peak3, peak4, ncol = 4)
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig6.eps", twopanelplot, width = 17.4, height = 15, units = "cm", device = cairo_pdf)


#### probability vector ####

prob.vec2 = makehistlist("1A", geno12, F)
prob.vec2 = table(factor(prob.vec2, levels = 1:max(prob.vec2)))
prob.vec3 = as.numeric(prob.vec2) / sum(as.numeric(prob.vec2))
s(prob.vec3, "rotation1scripts_v4/saved.objects/prob.vec3", "c.x.aprocessing2.R")


prob.vec2 = table(prob.vec2) / sum(table(prob.vec2))
all.prob.vec = data.frame(1:225, 1:225)

attributes(na.omit(match(1:225, as.data.frame(table(prob.vec2))$prob.vec2)))$na.action


gen.hap(225, prob.vec = prob.vec3)

gen1.pv = gen.data(225, 1000, prob.vec = prob.vec3)

gen1.pv.rqtl = convert.sim.data.to.rqtl.format(gen1.pv)
write.csv(gen1.pv.rqtl, "rotation1scripts_v4/temp/gen1.pv.rqtl.csv", row.names = F)



gen1.rqtl = read.cross("csv", "rotation1scripts_v4/temp/", "gen1.pv.rqtl.csv", genotypes=c("A", "H", "B"),estimate.map=F)

table(countXO(gen1.rqtl))

par(mfrow = c(1, 1))
heatMap(gen1.rqtl, what = "rf")

gen1.rf = as.data.frame(pull.rf(gen1.rqtl))
sapply(gen1.rf, function(x) max(na.omit(x)))





acf(ratiodfs2.3[["1A"]][[1]]$seg.ratio1, lag.max = 1, plot = F)

sim.autocorrelations = lapply(axcf2.noselec.pop96, function(x){
    acf(convert.recomb.to.seg(x), plot = F, lag.max = 1)    
})



#### PREPARE SUPPLEMENTARY FILES ####

axc512 = read_delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/axc51-2.all.snps.txt", delim = "\t", comment = "##")
axc611 = read_delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/axc61-1.all.snps.txt", delim = "\t", comment = "##")
cxa651 = read_delim("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/c.x.a65-1.all.snps.txt", delim = "\t", comment = "##")
cxa653 = read_delim("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/c.x.a65-3.all.snps.txt", delim = "\t", comment = "##")

colnames(axc512)
colnames(axc611)
colnames(cxa651)
colnames(cxa653)


library(RMySQL)
    
mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')
    
rs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Cadenza%'") #grab Avalon genotyping data for 35k array
rs2 = fetch(rs, n = -1)


qs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Avalon%'") #grab Avalon genotyping data for 35k array
qs2 = fetch(qs, n = -1)


rs2 = rs2[which(rs2$Var_col == "Cadenza"), ]
qs2 = qs2[which(qs2$Var_col == "Avalon"), ]

all(match(rs2$Probe_row, axc512$probeset_id) == 1:nrow(rs2))
all(match(qs2$Probe_row, axc512$probeset_id) == 1:nrow(qs2))


g = cbind(axc512, rs2$Matrix_value, qs2$Matrix_value)
colnames(g)[(ncol(g) - 1):ncol(g)] = c("Cadenza", "Avalon")

all1 = cbind(g, axc611[, 2:ncol(axc611)], cxa651[, 2:ncol(cxa651)], cxa653[, 2:ncol(cxa653)])
which(colnames(all1) %in% c("Cadenza", "Avalon"))
all2 = all1[, c(1, 99, 100, 2:98, 101:ncol(all1))]

all3 = as.data.frame(t(all2))
all3 = cbind(all3[, 1], all3)
all3[, 1] = rownames(all3)
all3 = cbind(all3[, 1], all3)
all3[, 1] = ""
all3 = reset.rownames(all3)
all3 = convert.to.character.data.frame(all3)
colnames(all3) = all3[1, ]
all3[1, ] = ""

all3[all3 == "AA"] = "A"
all3[all3 == "AB"] = "H"
all3[all3 == "BB"] = "B"
all3[all3 == "NoCall"] = "-"





all4 = cleanup(all3, 2:3)
all4 = convert.to.character.data.frame(all4)
all4 = assign.parental.genotypes(all4, 2:3)

# write.csv(all4, "rotation1scripts_v4/temp/cxageno.csv")


rqtl1 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/axc51-2.all.snps.txt")

all5 = as.data.frame(t(all4))
all5 = convert.to.character.data.frame(all5)
parent.coords = which(sapply(all5, function(x){
    
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    
    if(a > 4000 | b > 4000){
        T
    } else {
        F
    }
}))

other.coords = 1:nrow(all4)
other.coords = other.coords[-c(1, parent.coords)]


all6 = all4[c(1, parent.coords, other.coords), ]

nrow(cxaf2.c651.asmap)

all7 = as.data.frame(t(all6))
all8 = all7[-1, -1]
all8 = convert.to.character.data.frame(all8)
all8[1, 3:6] = c("Cadenza2", "Avalon2", "Cadenza3", "Avalon3")
all8 = reset.colnames(all8)

all9 = all8[c(1, na.omit(match(switch.affy.format(cxaf2.c651.asmap$marker), rownames(all8)))), ]
all9 = cbind(cxaf2.c651.asmap, all9)


missingmarkers = switch.affy.format(cxaf2.c651.asmap$marker)[which(!switch.affy.format(cxaf2.c651.asmap$marker) %in% rownames(all8))]

missingmarkers2 = all3[c("probeset_id", missingmarkers)]
g = missingmarkers2[c(1, parent.coords, other.coords), ]


consensus.cadenza = g[c(2, 4, 6, 3, 5, 7), ]
consensus.cadenza2 = consensus.cadenza


consensus.cadenza2[] = lapply(consensus.cadenza2, function(x){
    # browser()
    cad = table(x[1:3])
    ava = table(x[4:6])
    
    cad = cad[which(names(cad) == "A" | names(cad) == "B")]
    ava = ava[which(names(ava) == "A" | names(ava) == "B")]
    
    cad = names(cad[which(cad == max(cad))])
    ava = names(ava[which(ava == max(ava))])
    
    if(length(unique(cad)) == 1){
        x[1:3] = cad
    }
    
    if(length(unique(ava)) == 1){
        x[4:6] = ava
    }
    x
    
})

which(!sapply(consensus.cadenza2, function(x){
    if(length(unique(x[1:3])) == 1 & length(unique(x[4:6])) == 1) return(T)
    F
}))

consensus.cadenza2[, 288] = c("B", "B", "B", "A", "A", "A")

new1 = all3[c(1, parent.coords, other.coords), ]
new1[1:7, ] = new1[c(1, 2, 4, 6, 3, 5, 7), ]


v(new1)
new1[2:7, match(colnames(consensus.cadenza2), colnames(new1))] = consensus.cadenza2

new2 = cleanup(new1, c(2, 5))
new2 = assign.parental.genotypes(new2, c(2, 5))

new2 = as.data.frame(t(new2), stringsAsFactors = F)
new2 = new2[-1, -1]
new2 = reset.colnames(new2)
colnames(new2) = new2[1, ]




new3 = new2[match(switch.affy.format(cxaf2.c651.asmap$marker), rownames(new2)), ]

new4 = as.data.frame(cbind(cxaf2.c651.asmap, new3), stringsAsFactors = F)
new4 = new4[, -c(4, 5, 9)]

colnames(new4)[1:12] = c("Marker", "Linkage Group", "Genetic Position (cM)", "Chromosome", "Physical Distance (%)", "Physical Distance (BP)",
                                                 "Cadenza1", "Cadenza2", "Cadenza3", "Avalon1", "Avalon2", "Avalon3")
write.csv(new4, "E:/ac14037.backup/Google Drive/University/PhD/Seg dist simulation/supplementary.files/genetic.map.csv", row.names = F)

colnames(new4)
#cadenza population sizes
length(grep("65\\-3", colnames(new4)))
length(grep("65\\-1", colnames(new4)))
#avalon population sizes
length(grep("61\\-1", colnames(new4)))
length(grep("51\\-2", colnames(new4)))

dir.exists("rotation1scripts_v4/")

#### All figures ####

all.rf.plot2
seg.comp.plot #pass to grid.arrange
fig3.facets
plotheat
plot1
plot2
grid.arrange(genetic.rf.plot1)
final.cake

main_plots = list(all.rf.plot2, seg.comp.plot, fig3.facets, plotheat, plot1, plot2, genetic.rf.plot1, final.cake)

#### Supplementary figures ####

fig6
recombinationplot1
selectionpos.plot
sdrplot.all
recombinationplot3

supp_plots = list(fig6, recombinationplot1, selectionpos.plot, sdrplot.all, recombinationplot3)

s(main_plots, '~/project.phd.main/rotation1scripts_v4/saved.objects/main_plots', 'c.x.aprocessing2.R')
s(supp_plots, '~/project.phd.main/rotation1scripts_v4/saved.objects/supp_plots', 'c.x.aprocessing2.R')


#### PLOTTING UPDATES ####

main_plots[[2]] = main_plots[[2]] +
    geom_path() + geom_point(aes(color = chi1t, shape = chi1t), size = 2) + geom_hline(yintercept = 0.5) + coord_cartesian(ylim = yscale1) +
    theme_bw(base_size = 22) +
    scale_color_manual(values=c("#492de5", "#999999", "#E69F00", "#56B4E9"), labels = c("N.S.", "0.05", "0.01", "0.001")) +
    scale_shape_manual(values = c(16, 15, 17, 18), labels = c("N.S.", "0.05", "0.01", "0.001")) +
    xlab("Genetic distance (cM)") + ylab("Segregation ratio") +
    guides(color = guide_legend(title = "Sig."), shape = guide_legend(title = "Sig.")) + facet_grid(num1 ~ .) +
    theme(strip.text.y = element_text(angle = 0, size = 12), strip.background = element_rect(color = "white", fill = "#ffffff"), strip.placement = "outside",
                panel.grid = element_blank())






main_plots[[3]] + geom_line(size = 1) +
    theme_bw(base_size = 22) + theme(axis.text.x = element_text(angle = 45)) + 
    theme(strip.text = element_text(size = 20), strip.background = element_rect(colour = "white", fill = "#ffffff"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.placement = "outside")
