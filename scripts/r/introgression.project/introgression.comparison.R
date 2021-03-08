#comparative introgression analysis
setwd("C:/Users/ac14037/project.phd.main/")
library(readr)
library(dplyr)

intro.liberal = read_csv("rotation1scripts_v4/processed_data/introgressions/full.comp.4b.rev/all.intro.4b.rev.no.overlaps.csv")
intro.bar.cons = read_csv("rotation1scripts_v4/processed_data/introgressions/full.comp.4b.rev.order.cons.bar/all.intro.bar.cons.no.overlaps")

load("rotation1scripts_v4/saved.objects/rows.to.rm.same.intro")

intro.bar.cons = intro.bar.cons[-rows.to.rm.same.intro, ]

#checking if intro.bar.cons has any introgressions that span the same width but are from different wild relatives
# intro.bar.cons2 = intro.bar.cons[, c(2, 4, 5, 6, 7, 8, 9)]
# nrow(intro.bar.cons) / nrow(unique(intro.bar.cons2))

intro.ultra.cons = read_csv("rotation1scripts_v4/processed_data/introgressions/ultra.cons.bar.cons/all.intro.ult.bar.cons.no.overlaps.csv")

all.intro = list(intro.liberal, intro.bar.cons, intro.ultra.cons)
names(all.intro) = c("intro.liberal", "intro.bar.cons", "intro.ultra.cons")

load("rotation1scripts_v4/saved.objects/germplasm")

#perform some filtering (minimum number.snps and intro.length); assign species and sample types from germplasm df
all.intro2 = lapply(all.intro, function(x){
    grab.germplasm.field = function(name.of.field){
        m.coord = match(x$wild.relative, germplasm$Name)
        germplasm[[name.of.field]][m.coord]
    }
    x$species = grab.germplasm.field("Species")
    x$samp.type = grab.germplasm.field("Sample type")
    x$intro.length.mb = x$intro.length / 1000000
    x$phy.dist.from.avalon = unlist(as.list(tree.dist3[1, match(substr(x$wild.relative, 1, 10), colnames(tree.dist3))])) #REQUIRES tree.dist3 FROM phylogenetic.distances.R
    x$snp.density = x$snp.density * 1000000 #convert snp density from snps / bp to snps / Mb
    x = filter(x, snp.density > 5)
    x = filter(x, intro.length.mb > 0.05)
    x = filter(x, number.snps > 3)
    x$chromosome = unlist(lapply(x$chromosome, function(q1){ #ultra cons dataset doesn't need substr, chromosomes already in proper format
        if(nchar(q1) > 2){
            q2 = substr(q1, 4, 5)
        } else {
            q2 = q1
        }
        q2
    }))
    # if(nchar(x$chromosome) > 2) x$chromosome = substr(x$chromosome, 4, 5)
    x = x[-which(x$elite == "Synthetic_39"), ]
    x
})

#split into list of dfs by elite and chromosomes
all.intro3 = lapply(all.intro2, function(x){
    g = split.df.into.list.of.dfs.by.column(x, "elite")
    g = g[c(1, 2, 6, 7, 4, 8, 3, 5)]
    g = lapply(g, function(q) arrange(q, chromosome, intro.start))
    g = lapply(g, function(q) split.df.into.list.of.dfs.by.column(q, "chromosome"))
    g
})

#rm redundant introgressions w/ 90 percent overlap
all.intro4 = lapply(all.intro3, function(x){
    lapply(x, function(q){
        lapply(q, function(p){
            tryCatch(remove.smaller.in.90.per.overlap(p), error = function(e) p) #in case of an error, just return the chromosome without the function applied to it
        })
    })
})

all.intro4.comb = lapply(all.intro4, function(x){
    combine.list.of.data.frames(lapply(x, combine.list.of.data.frames))
})

#     ____________________________________________________________________________
#     STATISTICS AND PLOTS                                                                                                        ####

#check whether elite varieties in the same cross have the same introgressions
check.for.same.introgressions = function(parent1, parent2){
    g = filter(all.intro4.comb[[2]], elite == parent1)
    g1 = filter(all.intro4.comb[[2]], elite == parent2)
    
    g.1 = g[, -(1:2)]
    g1.1 = g1[, -(1:2)]
    
    g.rows = do.call(paste0, g.1)
    g1.rows = do.call(paste0, g1.1)
    
    g[which(g.rows %in% g1.rows), ]$X1
}

a.x.c.same = check.for.same.introgressions("Avalon", "Cadenza")
r.x.s.same = check.for.same.introgressions("Rialto", "Savannah")
o.x.s.same = check.for.same.introgressions("Opata", "Synthetic")
c.x.p.same = check.for.same.introgressions("ChineseSpring", "Paragon")

# rows.to.rm.same.intro = c(a.x.c.same, r.x.s.same, o.x.s.same, c.x.p.same)

# s(rows.to.rm.same.intro, "rotation1scripts_v4/saved.objects/rows.to.rm.same.intro", "introgression.comparison.R")

#examine the frequency of introgressions in terms of the patristic distance from Avalon
png("rotation1scripts_v4/plots/introgression.histogram.patristic.distance/int.hist.1.png", 1000, 1000)
hist(all.intro4.comb[[2]]$phy.dist.from.avalon, breaks = 100, xlab = "Patristic Distance from Avalon", main = "Introgression Histogram")
dev.off()



lapply(all.intro4.comb, function(x){
    g = table(x$samp.type)
    g[[1]]/g[[2]]
})

sp.freqs1 = sort(table(all.intro4.comb[[2]]$species))
p.dist1 = all.intro4.comb[[2]]$phy.dist.from.avalon[match(names(sort(table(all.intro4.comb[[2]]$species))), all.intro4.comb[[2]]$species)]
p.dist.df = data.frame(1:44, (p.dist1*15000))
colnames(p.dist.df) = c("V1", "V2")

#barplots showing number of introgressions detected for all species (one for each dataset)
pdf("rotation1scripts_v4/plots/introgression.hist/sp.hist7.pdf", 15, 7)
# png("rotation1scripts_v4/plots/introgression.hist/sp.hist3.png", 1400, 800)
ggplot(as.data.frame(sort(table(all.intro4.comb[[2]]$species))), aes(x = Var1, y = Freq)) +
        geom_bar(stat = "identity") + theme_bw() +
        theme(axis.text.x = element_text(angle = 80, hjust = 1, size = 12), legend.position = "none") + xlab("Species") +
        ylab("Frequency") + 
        geom_point(data = p.dist.df, aes(x = V1, y = V2, color = "red"), stat = "identity") +
        geom_line(data = p.dist.df, aes(x = V1, y = V2, color = "red"), stat = "identity")
dev.off()



#examine the distribution of introgressions using histograms
par(mfrow = c(2, 2))

breaks2 = c(breaks1, histinfo3$breaks[53:69])

histinfo = hist(all.intro4[[1]][[1]][[1]]$intro.start, xlim = c(0, 7e8), breaks = 100)
breaks1 = histinfo$breaks
histinfo2 = hist(all.intro4[[2]][[1]][[1]]$intro.start, xlim = c(0, 7e8), breaks = breaks1)
histinfo3 = hist(all.intro4[[3]][[1]][[1]]$intro.start, xlim = c(0, 7e8), breaks = breaks1)


#circos plot experiment

data(UCSC.HG19.Human.CytoBandIdeogram)
RCircos.Set.Core.Components(UCSC.HG19.Human.CytoBandIdeogram)

RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Histogram.Data)
data.col <- 4
track.num <- 8
side <- "in"
pdf("rotation1scripts_v4/plots/circos.plot.test/plot1.pdf", 5, 5)

#change color of histogram on circos plot
rc.param = RCircos.Get.Plot.Parameters()
rc.param$track.background = "white"
rc.param$hist.color = "black"
RCircos.Reset.Plot.Parameters(rc.param)

RCircos.Set.Plot.Area()
RCircos.Histogram.Plot(RCircos.Histogram.Data, data.col, track.num, side)
dev.off()

#     ____________________________________________________________________________
#     OLD INTROGRESSION ANALYSIS CODE                                                                                 ####

#examine the contribution of introgressions under 10 SNPs
check.dist.of.small.ints = function(all.introgressions, elite.var, chromo, sample.type){
    #generates list of dataframes, one df per wheat relative. 
    #columns of df:
    #diff - the difference, in number of SNPs, between this run of consecutive SNPs
    #and the previous run of consecutive SNPs for this pariticular wheat relative
    #length - number of consecutive identical SNPs
    #snp.start - start position of this particular run of consecutive SNPs
    
    #args:
    #elite.var - character string, name of elite variety to analyse
    #chromo - character string, name of chromosome to analyse
    #sample.type (optional) - character string, either "relative" or "progenitor"
    if(missing(sample.type)) sample.type = F
    if(sample.type != F){
        w.rel.list = sort(table(filter(all.introgressions, chromosome == chromo, samp.type == sample.type, elite == elite.var)$wild.relative))
        
        w.rel.max = names(w.rel.list[length(w.rel.list)])
        avg.snp.dist = lapply(names(w.rel.list), function(x){
            q = c(0, diff(arrange(filter(all.introgressions, chromosome == chromo, samp.type == sample.type, elite == elite.var, wild.relative == x), snp.start)$snp.start))
            q1 = arrange(filter(all.introgressions, chromosome == chromo, samp.type == sample.type, elite == elite.var, wild.relative == x), snp.start)$number.snps
            q2 = arrange(filter(all.introgressions, chromosome == chromo, samp.type == sample.type, elite == elite.var, wild.relative == x), snp.start)$snp.start
            
            # q1 = list(mean(q), sd(q))
            # names(q1) = c("mean", "sd")
            
            q3 = list(q, q1, q2)
            names(q3) = c("diff", "length", "snp.start")
            q3 = as.data.frame(q3)
            q3
            
        })
    } else {
        w.rel.list = sort(table(filter(all.introgressions, chromosome == chromo, elite == elite.var)$wild.relative))
        
        w.rel.max = names(w.rel.list[length(w.rel.list)])
        avg.snp.dist = lapply(names(w.rel.list), function(x){
            q = c(0, diff(arrange(filter(all.introgressions, chromosome == chromo, elite == elite.var, wild.relative == x), snp.start)$snp.start))
            q1 = arrange(filter(all.introgressions, chromosome == chromo, elite == elite.var, wild.relative == x), snp.start)$number.snps
            q2 = arrange(filter(all.introgressions, chromosome == chromo, elite == elite.var, wild.relative == x), snp.start)$snp.start
            
            q3 = list(q, q1, q2)
            names(q3) = c("diff", "length", "snp.start")
            q3 = as.data.frame(q3)
            q3
            
        })
    }
    names(avg.snp.dist) = names(w.rel.list)
    avg.snp.dist
    
}

q = check.dist.of.small.ints(all.intro4.comb[[2]], "Paragon", "2A", "relative")

q1 = q[[length(q)]]
q1 = q[[8]]

ggplot(q1, aes(x = snp.start, y = length)) + geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = 10, color = "red", linetype = "dashed")

q = check.dist.of.small.ints("Paragon", "2A")

q1 = all.introgressions[which(all.introgressions$samp.type == "progenitor"), ]
ggplot(q1, aes(x = 1:nrow(q), y = Freq)) + geom_bar(stat = "identity") 
# coord_cartesian(xlim = c(0, 30), ylim = c(0, 13000))

#introgression distribution plot
png("rotation1scripts_v4/plots/introgression.hist/distribution.hist.png", 1000, 500)
par(mfrow = c(1, 2))
hist(filter(final.blast3, sseqid == "chr1A")$sstart, breaks = 100, xlab = "Physical position (bp)", main = "(a) 2A SNPs")
hist(filter(a6, chromosome == "2A")$intro.start, breaks = 100, xlab = "Physical position (bp)", main = "(b) 2A Introgressions")
dev.off()

#varieties in the genotype file but not in the germplasm list
unique(a6$wild.relative)[which(!unique(a6$wild.relative) %in% unique(germplasm$Name))]

#     ____________________________________________________________________________
#     FUNCTIONS                                                                                                                             ####

remove.smaller.in.90.per.overlap = function(intro.df1){
    #if two introgressions partially overlap, and the overlap contains 90% of the range of the smaller of the two introgressions,
    #the smaller of the two introgressions is discarded.
    g1 = intro.df1
    g2 = IRanges(g1$intro.start, g1$intro.end)
    g3 = findOverlaps(g2, drop.self = T, drop.redundant = T)
    
    
    g4 = cbind(as.data.frame(g3), as.data.frame(ranges(g3, g2, g2)))
    g4$query.start = ""
    g4$query.end = ""
    g4$query.length = ""
    g4$subject.start = ""
    g4$subject.end = ""
    g4$subject.length = ""
    g4[] = lapply(g4, as.numeric)
    
    g2.1 = as.data.frame(g2)
    for(i in 1:nrow(g4)){
        g4$query.start[i] = g2.1$start[g4$queryHits[i]]
        g4$query.end[i] = g2.1$end[g4$queryHits[i]]
        g4$query.length[i] = g4$query.end[i] - g4$query.start[i]        
        g4$subject.start[i] = g2.1$start[g4$subjectHits[i]]
        g4$subject.end[i] = g2.1$end[g4$subjectHits[i]]
        g4$subject.length[i] = g4$subject.end[i] - g4$subject.start[i]
    }
    
    g4$per.cov.min = ""
    g4$per.cov.min = as.numeric(g4$per.cov.min)
    g4$q.or.s.smaller = ""
    for(i in 1:nrow(g4)){
        q = min(g4$query.length[i], g4$subject.length[i])
        q.or.s = which(c(g4$query.length[i], g4$subject.length[i]) == q)
        
        if(q.or.s == 1){
            g4$q.or.s.smaller[i] = "q"
        } else {
            g4$q.or.s.smaller[i] = "s"
        }
        g4$per.cov.min[i] = (g4$width[i] / q) * 100
    }
    
    g5 = g4[which(g4$per.cov.min > 90), ]    
    rows.to.rm = as.numeric()
    for(i in 1:nrow(g5)){
        # print(g5$q.or.s.smaller[i])
        if(length(g5$q.or.s.smaller[i]) > 0){
            if(!is.na(g5$q.or.s.smaller[i])){
                if(g5$q.or.s.smaller[i] == "s"){
                    rows.to.rm = c(rows.to.rm, g5$subjectHits[i])        
                } else {
                    rows.to.rm = c(rows.to.rm, g5$queryHits[i])
                }
            }
        }
    }
    if(length(rows.to.rm) > 0) intro.df1 = intro.df1[-rows.to.rm, ]
    return(intro.df1)
    
}

#this function is for the Wilkins instance of R
rm.redundant.overlaps = function(intro.df1){
    ranges1 = IRanges(intro.df1$intro.start, intro.df1$intro.end)
    ranges.overlap = findOverlaps(ranges1, type = "within", ignoreSelf = T, ignoreRedundant = T) #NB. the names of these arguments has changed in the latest version if IRanges
    ro.df = as.data.frame(ranges.overlap)
    
    q1 = intro.df1[-ro.df$queryHits, ]
    q1
}

split.and.arrange = function(intro.df1){
    intro.apo2 = split.df.into.list.of.dfs.by.column(intro.df1, "chr")
    intro.apo3 = lapply(intro.apo2, function(x){
        arrange(x, bp.start)
    })
    return(intro.apo3)
}