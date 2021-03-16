#modified verison of axc.analysis.R
#testing MRD in F5 Chinese Spring X Paragon population

#### SETUP REPLICATE AxC ####

#original processing of genotyping data from geno calls done separately in AAS
# source("rotation1scripts_v4/scripts/r/recombination/axc.analysis/axc.initial.data.processing.R")

#new processing of geno data with geno calls done together in AAS 
source("rotation1scripts_v4/scripts/r/recombination/axc.analysis/axc.initial.data.processing.all.geno.comb.R")

axc.pops = c("cxa651", "cxa653", "axc512", "axc611")

l.geno.dfs3.w.change = list.of.geno.dfs.each.pop3


cs.map = allen.maps.formatted[[1]][[1]]

samplepool = 2:nrow(cs.map)
cs.sample1 = sample(samplepool, 70)

samplepool = samplepool[-which(samplepool %in% cs.sample1)]


cs.sample2 = sample(samplepool, 70)

cs.map.split1 = cs.map[c(1, cs.sample1), ]
cs.map.split2 = cs.map[c(1, cs.sample2), ]

cs.map.split1$probeset_id[2:nrow(cs.map.split1)] = paste0(cs.map.split1$probeset_id[2:nrow(cs.map.split1)], "_split1")
cs.map.split2$probeset_id[2:nrow(cs.map.split2)] = paste0(cs.map.split2$probeset_id[2:nrow(cs.map.split2)], "_split2")


asmaplgs = unique(axc.avg2$chromo)


l.geno.dfs3.w.change = list(cs.map.split1, cs.map.split2)

#need to take longest increasing subsequence of combined list of genotypes as there are potentially more cM bins
lgeno3comb = combinelist.of.genos(l.geno.dfs3.w.change)
lgeno3comb.centi = extract.centimorgan.from.genotypes(lgeno3comb)

convert.map.format = function(genetic.map){
    colnames(genetic.map) = c("cM", "chr", "marker")
    genetic.map[, c(3, 2, 1)]
}

lgeno3comb.centi = convert.map.format(lgeno3comb.centi)

axc.liss = qtlanalysis.all(lgeno3comb, lgeno3comb.centi, "RIL2", T, "-")

#update genotyping df to reflect liss
l.geno.dfs3.w.change = lapply(l.geno.dfs3.w.change, function(x){
    # browser()
    x[, c(1, 2, match(axc.liss[[2]]$marker, colnames(x)))]
})


recomb4 = lapply(l.geno.dfs3.w.change, detect.recombination)

booboo1 = list()

geno.axc = lapply(l.geno.dfs3.w.change, detect.recombination, conservative = T, liss.df = axc.liss)

lapply(geno.axc, nrow)

#### REMOVE DOUBLE EVENTS ####

booboo2 = do.call(bind_rows, booboo1)

g = unlist(lapply(1:nrow(booboo2), function(x){
    paste0(booboo2[x, 1], booboo2[x, 2], booboo2[x, 3])
}))



#2 events in each gamete #class 1
g1 = c(
    which(g == "ABA"),
    which(g == "BAB")
)

#2 events in one gamete #class 2
g2 = c(
    which(g == "AHA"),
    which(g == "BHB")
)

#3 events #class 3
g3 = c(
    which(g == "ABH"),
    which(g == "BAH"),
    which(g == "HBA"),
    which(g == "HAB")
)

#2 events, one in each gamete #class 4
g4 = c(
    which(g == "HAH"),
    which(g == "HBH"),
    which(g == "AHB"),
    which(g == "BHA")
)


booboo2$recomb.class = ""
booboo2$recomb.class[g1] = 1
booboo2$recomb.class[g2] = 2
booboo2$recomb.class[g3] = 3
booboo2$recomb.class[g4] = 4

#set MB distance between markers - setting it to 30
booboo.rm = booboo2[which(booboo2$dist.double < 40), ]
booboo.rm$quad = ""
booboo.rm$quad[grep("split1", booboo.rm$individual)] = 1
booboo.rm$quad[grep("split2", booboo.rm$individual)] = 2

#remove double events from the genotyping dataframes
lapply(1:nrow(booboo.rm), function(x){
    quad = as.numeric(booboo.rm[x, ]$quad)
    
    marker.name = booboo.rm[x, ]$probe.after.transition
    individual.name = booboo.rm[x, ]$individual
    
    l.geno.dfs3.w.change[[quad]][which(l.geno.dfs3.w.change[[quad]]$probeset_id == individual.name), which(colnames(l.geno.dfs3.w.change[[quad]]) == marker.name)] <<- "-"
    
})




geno.axc2 = lapply(l.geno.dfs3.w.change, detect.recombination, conservative = F, liss.df = axc.liss)

lapply(geno.axc, nrow)

#### PERFORM ANALYSIS REPLICATE AxC ####

#geno.centi1 = lapply(l.geno.dfs3.w.change, extract.centimorgan.from.genotypes)

#lower.qc.phys = get.phys.for.chromos(l.geno.dfs3.w.change[[1]], F, T, "-")

#geno.centi.len = lapply(geno.centi1, function(x){
    #x$chromo = as.factor(x$chromo)
    #levels(x$chromo) = unique(x$chromo)
    #aggregate(x$cm, list(x$chromo), max)
#})



##get total map length for all genetic chromosomes for each population
#total.map.length = lapply(geno.centi.len, function(x){
    #sum(x$x)
#})

#total.map.length2 = lapply(geno.centi.len, function(x){
    #sum(x[-12, ]$x)
#})



#gc.len2 = comb.treatments(geno.centi.len)
#colnames(gc.len2) = c("Chromosome", "cm", "Treatment")
#gc.len2$Treatment = as.factor(gc.len2$Treatment)

#gc.len2$Treatment = factor(gc.len2$Treatment, levels = unique(gc.len2$Treatment))
#aggregate(gc.len2$cm, list(gc.len2$Treatment), max)

#library(viridis)

#chromo.maplengths.plot = ggplot(gc.len2, aes(x = Chromosome, y = cm, group = Treatment, fill = Treatment)) + geom_bar(stat = 'identity', position = "dodge") +
    #scale_fill_viridis(discrete = T, labels = axc.pops, name = "Population") + theme_bw() + ylab("Length of chromosome (cM)")

#chromo.maplengths.plot

#plots1 = prepare.genetic.map.comparison.plots(l.geno.dfs3.w.change, geno.centi1)








#### distribution analysis 2 ####

#add physical positions to recombination dataframe
geno.axc2 = lapply(geno.axc2, function(x){
    x$phys.pos = lower.qc.phys[match(x$probe.before.transition, lower.qc.phys$marker), ]$phys.pos
    x
})

axc.ind.averages = lapply(geno.axc2, ind.avg, use.centromere = T, scale.by.arm = T, manual.center = 50)

axc.avg2 = comb.treatments(axc.ind.averages, 2)
levels(axc.avg2$chromo) = unique(axc.avg2$chromo)

#axc.avg3 = axc.avg2[-which(axc.avg2$treatment == 3), ]

#get mean MRD for each treatment
axc.avg2$treatment = factor(axc.avg2$treatment, levels = unique(axc.avg2$treatment))
aggregate(axc.avg2$phys.dist, list(axc.avg2$treatment), function(x) mean(na.exclude(x)))


mrd.long1 = aggregate(axc.avg2$phys.dist.long, list(axc.avg2$treatment), function(x) mean(na.exclude(x)))
aggregate(axc.avg2$phys.dist.long, list(axc.avg2$treatment), function(x) sd(na.exclude(x)))
mrd.short1 = aggregate(axc.avg2$phys.dist.short, list(axc.avg2$treatment), function(x) mean(na.exclude(x)))
aggregate(axc.avg2$phys.dist.short, list(axc.avg2$treatment), function(x) sd(na.exclude(x)))

#mrd.long1$Group.1 = axc.pops
#mrd.short1$Group.1 = axc.pops


axc.kruskal.all.long = kruskal.test(axc.avg2$phys.dist.long ~ axc.avg2$treatment)
axc.kruskal.all.long
axc.kruskal.all.short = kruskal.test(axc.avg2$phys.dist.short ~ axc.avg2$treatment)
axc.kruskal.all.short

#kruskal.test(axc.avg3$phys.dist ~ axc.avg3$treatment)



#remove.treatment = function(nvector1, arm1, mrd.avg.df){
    #if(missing(arm1)) arm1 = "none"
    #if(missing(mrd.avg.df)) mrd.avg.df = avg.r1
    #mrd.avg.df = mrd.avg.df[-which(mrd.avg.df$treatment %in% nvector1), ]
    #if(arm1 == "none") return(kruskal.test(mrd.avg.df$phys.dist ~ mrd.avg.df$treatment))
    #if(arm1 == "short") return(kruskal.test(mrd.avg.df$phys.dist.short ~ mrd.avg.df$treatment))
    #if(arm1 == "long") return(kruskal.test(mrd.avg.df$phys.dist.long ~ mrd.avg.df$treatment))
#}

##(3)

#lapply(1:4, function(x) remove.treatment(x, "short", axc.avg2))
#lapply(1:4, function(x) remove.treatment(x, "long", axc.avg2))

#axc.1.rm.p.short = round(p.adjust(sapply(1:4, function(x) unclass(remove.treatment(x, "short", axc.avg2))$p.value), method = "bonf"), 4)
#axc.1.rm.p.long = round(p.adjust(sapply(1:4, function(x) unclass(remove.treatment(x, "long", axc.avg2))$p.value), method = "bonf"), 4)

##get all pairwise combinations of numbers 1 to 4
#combinations.to.rm = as.data.frame(combn(1:4, 2))
#combinations.to.rm = lapply(combinations.to.rm, function(x) x)
#combinations.to.rm = c(list(1, 2, 3, 4), combinations.to.rm)


#long.kruskal = lapply(combinations.to.rm, function(x){
    #remove.treatment(x, "long", axc.avg2)
#})

#axcg1 = round(sapply(long.kruskal, function(x) unclass(x)$p.value), digits = 5)
#axcg1.1 = round(sapply(long.kruskal, function(x) unclass(x)$statistic), digits = 5)
#axcg2 = round(p.adjust(sapply(long.kruskal, function(x) unclass(x)$p.value), "bonf"), digits = 5)

#short.kruskal = lapply(combinations.to.rm, function(x){
    #remove.treatment(x, "short", axc.avg2)
#})

#axcg3 = round(sapply(short.kruskal, function(x) unclass(x)$p.value), digits = 5)
#axcg3.1 = round(sapply(short.kruskal, function(x) unclass(x)$statistic), digits = 5)
#axcg4 = round(p.adjust(sapply(short.kruskal, function(x) unclass(x)$p.value), "bonf"), digits = 5)

#axcg4.2 = unlist(lapply(combinations.to.rm, function(x) paste(axc.pops[x], collapse = ", ")))

##### TABLE 1 ####
#axcg5 = data.frame(axcg4.2, axcg1, axcg1.1, axcg2, axcg3, axcg3.1, axcg4)



#colnames(axcg5) = c("Treatments removed", "p.long", "chi.long", "p.long.bonf", "p.short", "chi.short", "p.short.bonf")











#histogram.average.distances.of.recomb = ggplot(axc.avg2, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment)) + theme_bw() + xlab("Physical Distance From Center of Chromosome (%)") +
    #ylab("Frequency")

#means.axc.phys = as.data.frame(aggregate(axc.avg2$phys.dist, list(axc.avg2$treatment), function(x) mean(na.exclude(x))))
#colnames(means.axc.phys) = c("Population", "Mean")
#sd.axc.phys = as.data.frame(aggregate(axc.avg2$phys.dist, list(axc.avg2$treatment), function(x) sd(na.exclude(x))))
#colnames(sd.axc.phys) = c("Population", "SD")

#hist.per.chromo = ggplot(axc.avg2, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("Physical Distance From Center of Chromosome (%)") +
    #ylab("Frequency")


#axc.split1 = asmaplgs[1:8]
#axc.split2 = asmaplgs[9:15]

#axc.avg2.2 = axc.avg2[which(axc.avg2$chromo %in% axc.split1), ]
#axc.avg2.3 = axc.avg2[which(axc.avg2$chromo %in% axc.split2), ]




#hist.per.chromo2 = ggplot(axc.avg2.2, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("") +
    #ylab("Frequency")

#hist.per.chromo3 = ggplot(axc.avg2.3, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("Physical Distance From Center of Chromosome (%)") +
    #ylab("Frequency")

#axc.avg2$chromo = factor(axc.avg2$chromo, levels = asmaplgs)
#means.per.chromo = as.data.frame(aggregate(axc.avg2$phys.dist, list(axc.avg2$treatment, axc.avg2$chromo), function(x) mean(na.exclude(x))))
#means.per.chromo2 = dcast(means.per.chromo, Group.2 ~ Group.1)

#apply(means.per.chromo2, 1, function(x){
    #which(x[2:5] == max(x[2:5]))
#})

#length(which(apply(means.per.chromo2, 1, function(x){
    #which(x[2:5] == max(x[2:5]))
#}) == 3))



####p.values2####
asmaplgs = unique(lower.qc.phys$chromo)

perform.chromo.kruskal = function(ind.avg.input, var.to.test, correction1){
    if(missing(correction1)) correction1 = F
    p.val1 = round(unlist(sapply(asmaplgs, function(x){
        a3 = ind.avg.input[which(ind.avg.input$chromo == x), ]
        
        kruskal.test(a3[[var.to.test]] ~ a3$treatment)[3]
        
    })), digits = 5)
    
    chi.val1 = round(unlist(sapply(asmaplgs, function(x){
        a3 = ind.avg.input[which(ind.avg.input$chromo == x), ]
        
        kruskal.test(a3[[var.to.test]] ~ a3$treatment)[1]
        
    })), digits = 3)
    
    df.val1 = round(unlist(sapply(asmaplgs, function(x){
        a3 = ind.avg.input[which(ind.avg.input$chromo == x), ]
        
        kruskal.test(a3[[var.to.test]] ~ a3$treatment)[2]
        
    })), digits = 3)
    
    
    
    names(p.val1) = asmaplgs
    names(chi.val1) = asmaplgs
    names(df.val1) = asmaplgs
    
    if(correction1 == T){
        p.val1 = p.adjust(p.val1, method = "bonf")
    }
    
    
    data.frame(p.val1, chi.val1, df.val1)
    
}

#all treatments

make.mrd.chromo.df = function(mrd.df1, pop.rm, pop.names){
    if(missing(pop.rm)) pop.rm = as.numeric()
    if(missing(pop.names)) pop.names = "None"
    
    if(length(pop.rm) > 0){
        mrd.df1 = mrd.df1[-which(mrd.df1$treatment %in% pop.rm), ]
        pop.names = paste0(pop.names[pop.rm], collapse = ", ")
    } 
    
    chromo.long.mrd = perform.chromo.kruskal(mrd.df1, "phys.dist.long", T)
    chromo.short.mrd = perform.chromo.kruskal(mrd.df1, "phys.dist.short", T)
    chromo.mrd1 = bind_cols(chromo.long.mrd, chromo.short.mrd)
    chromo.mrd1$Chromosome = asmaplgs
    chromo.mrd1$pop.rm = pop.names
    colnames(chromo.mrd1) = c("p SA", "chi SA", "df SA", "p LA", "chi LA", "df LA", "Chromosome", "Population removed")
    
    chromo.mrd1 = put.col.first(chromo.mrd1, "Population removed")
    chromo.mrd1 = put.col.first(chromo.mrd1, "Chromosome")    
    chromo.mrd1
}

#### TABLE 2 ####

make.mrd.chromo.df(axc.avg2)
make.mrd.chromo.df(axc.avg2, 3, axc.pops)
make.mrd.chromo.df(axc.avg2, c(3:4), axc.pops)
make.mrd.chromo.df(axc.avg2, 1:2, axc.pops)

#### EXAMINE RECOMBINATION FREQUENCY #### 

recomb.freq1 = lapply(geno.axc2, function(x){
    # browser()
    x = arrange(x, individual)
    aggregate(x$num.events, list(factor(x$individual, levels = unique(x$individual))), sum)
})

recomb.freq2 = comb.treatments(recomb.freq1)

# recomb.freq2$treatment = axc.pops[recomb.freq2$treatment]

hist(recomb.freq1[[1]]$x)

anova.freq = aov(recomb.freq2$x ~ as.factor(recomb.freq2$treatment))
anova.freq.tuk = TukeyHSD(anova.freq)






library(broom)
aft2 = tidy(anova.freq.tuk)
aft2.p = round(aft2$adj.p.value[c(1, 3, 4, 6)], digits = 4)

aft.means = aggregate(recomb.freq2$x, list(factor(recomb.freq2$treatment, levels = unique(recomb.freq2$treatment))), mean)
colnames(aft.means) = c("treatment", "mean")
aft.sds = aggregate(recomb.freq2$x, list(factor(recomb.freq2$treatment, levels = unique(recomb.freq2$treatment))), sd)
colnames(aft.sds) = c("treatment", "sd")

aft3 = bind_cols(aft.means, aft.sds)

aft3 = relabel.treat.w.pop.names(aft3, "treatment", axc.pops, F)
aft3 = relabel.treat.w.pop.names(aft3, "treatment1", axc.pops, F)
aft3$mean = round(aft3$mean, digits = 3)
aft3$sd = round(aft3$sd, digits = 3)
# install.packages("ggsignif")
library(ggsignif)

recomb.freq.plot1 = ggplot(aft3, aes(x = treatment, y = mean)) + geom_bar(stat = 'identity', fill = "#e6e6e6", color = "#000000") + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.5) + 
    xlab("Population") + ylab("Mean recombination events")

recomb.freq.plot1 = add.comparison.bars(recomb.freq.plot1, anova.freq.tuk, 30, axc.pops)

recomb.freq.plot1

#calculate recombination frequency in individuals
recomb.events = lapply(geno.axc2, function(x){
    x$chromo = factor(x$chromo, levels = unique(x$chromo))
    aggregate(x$num.events, list(x$chromo), function(y) sum(y))
})

par(mfrow = c(2, 2))
barplot(recomb.events[[1]]$x)
barplot(recomb.events[[2]]$x)
barplot(recomb.events[[3]]$x)
barplot(recomb.events[[4]]$x)

recomb.events2 = comb.treatments(recomb.events)
colnames(recomb.events2) = c("Chromosome", "num.events", "Treatment")

frequency.plot1 = ggplot(recomb.events2, aes(x = Chromosome, y = num.events)) + geom_bar(stat = "identity") + facet_grid(rows = vars(Treatment))


#see if number of recombination events correlates to chromosome length
iwgsc.chromo.lengths = read.delim("rotation1scripts_v4/original_data/iwgsc.chromosome.counts.txt", sep = " ", header = F, stringsAsFactors = F)
iwgsc.chromo.lengths$V2 = multi.str.split(iwgsc.chromo.lengths$V2, "chr", 2)

iwgsc.chr3 = iwgsc.chromo.lengths[match(unique(recomb.events2$Chromosome), iwgsc.chromo.lengths$V2), ]
frequency.plot2 = ggplot(iwgsc.chr3, aes(x = V2, y = V1)) + geom_bar(stat = 'identity')

plot_grid(frequency.plot1, frequency.plot2, rel_heights = c(0.7, 0.3), ncol = 1, align = "v", axis = "l")




low.2 = split.df.into.list.of.dfs.by.column(lower.qc.phys, "chromo")
coverage = sapply(low.2, function(x){
    100 - mean(diff(na.omit(x$phys.pos)[longest_subseq.R(na.omit(x$phys.pos))]))
})

#grab closest marker to center of chromosome by physical position
lapply(low.2, function(x){
    g = abs(na.omit(x$phys.pos) - 50)
    g2 = na.omit(x$phys.pos)[which(g == min(g))]
    x[which(x$phys.pos == g2), ]
})




recomb.events2$Chromosome = factor(recomb.events2$Chromosome, levels = unique(recomb.events2$Chromosome))
recomb.summed = aggregate(recomb.events2$num.events, list(recomb.events2$Chromosome), sum)

cor.df1 = bind_cols(recomb.summed, iwgsc.chr3)
colnames(cor.df1) = c("Chromosome", "num.events", "Physical.Length", "Chromosome2")

corplot1 = ggplot(cor.df1, aes(Physical.Length, num.events)) + geom_point() + xlab("Length of chromosome (bp)") + ylab("Number of recombination events")
corplot2 = ggplot(cor.df1, aes(coverage, num.events)) + geom_point() + xlab("Chromosome marker coverage") + ylab("Number of recombination events")

grid.arrange(corplot1, corplot2, ncol = 2)

cor.df1$coverage = coverage




cor.result1 = cor.test(cor.df1$num.events, cor.df1$Physical.Length)
unclass(cor.result1)[1:4]

cor.result2 = cor.test(cor.df1$num.events, cor.df1$coverage, method = "spearman")

cor.test(recomb.summed$x, coverage, method = "spearman")

# nlcor(recomb.summed$x, coverage)

report_cor = function(x){
    x = unclass(x)
    x = lapply(x[1:5], round, digits = 4)
    paste0("t = ", x[[1]], ", df = ", x[[2]], ", p = ", x[[3]], ", cor = ", x[[4]])
}

report_cor_spear = function(x){
    x = unclass(x)
    x = lapply(x[c(1, 3, 4)], round, digits = 4)
    paste0("S = ", x[[1]], ", p = ", x[[2]], ", rho = ", x[[3]])
}

report_cor(cor.result1)
report_cor_spear(cor.result2)







#### GENE DISTRIBUTION ANALYSIS ####

individual.genetic.maps = lapply(l.geno.dfs3.w.change, function(x){
    g = extract.centimorgan.from.genotypes(x)
    g$phys2 = lower.qc.phys$phys.bp
    g$phys3 = lower.qc.phys$phys.pos
    g
}) 

indi.maps.test = split(individual.genetic.maps[[1]], individual.genetic.maps[[1]]$chromo)

lapply(indi.maps.test, function(x){
    sequence1 = seq(0, 100, 25)
    
    distances1 = sapply(sequence1, function(z){
        min(abs(z - x$phys3))
    })
    
    length(which(distances1 > 25))
})



count2 = 1
templabels = axc.pops

axc.genedistplotsv2 = lapply(individual.genetic.maps, function(y){
    add.label1 = T
    cslm2 = split(y, y$chromo)
    cslm2 = lapply(cslm2, function(x){
        x$phys.mb = x$phys2 / 1000000
        x
    })
    
    
    cslm3 = lapply(cslm2, function(x){
        x$cmdiff = c(0, diff(x$cm))
        x$physdiff = c(0, diff(x$phys.mb))
        x$cmmb = x$cmdiff / x$physdiff
        x
    })
    
    count1 = 1
    plots1 = lapply(cslm3, function(x){
        # browser()
        
        x2.backup = x
        g5 = max(na.omit(x2.backup$cmmb)) + 2
        # x2.backup$cmmb = g5
        x2.backup$cmmb = -2
        x2.backup$backup = T
        x$backup = F
        
        x3 = bind_rows(x, x2.backup)
        
        # plot1 = ggplot(x3, aes(x = phys.mb, y = cmmb, group = backup)) + geom_line() + ggtitle(paste0(unique(individual.genetic.maps[[1]]$chromo)[count1], " ", templabels[count2])) + 
        #     theme_classic() + ylab("cM / Mb") + xlab("Physical position (Mb)") + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(-50, 900), ylim = c(0, (g5 + 5))) + geom_point() +
        #     geom_text_repel(aes(label = ifelse(backup == F, marker, '')), ylim = c((g5), (g5 + 10)), segment.color = "grey80", angle = 320, direction = "x")
        # 
        
        if(add.label1 == T){
            if(count2 == 4){
                plot1 = ggplot(x, aes(x = phys.mb, y = cmmb)) + ggtitle(paste0(unique(individual.genetic.maps[[1]]$chromo)[count1], " ", templabels[count2])) + 
                    theme_classic() + ylab("cM / Mb") + xlab("Physical position (Mb)") + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(0, 773), ylim = c(0, 40)) +
                    geom_text_repel(aes(label = marker), ylim = c(20, 40), xlim = c(0, 750), segment.color = "grey80", angle = 320, size = 2.5, direction = "x", nudge_x = 5) + geom_line() + geom_point() +
                    scale_y_continuous(breaks = c(0, 5, 10, 15), labels = c(0, 5, 10, 15))
            } else {
                plot1 = ggplot(x, aes(x = phys.mb, y = cmmb)) + geom_line() + ggtitle(paste0(unique(individual.genetic.maps[[1]]$chromo)[count1], " ", templabels[count2])) + 
                    theme_classic() + ylab("cM / Mb") + xlab("Physical position (Mb)") + scale_x_continuous(expand = c(0.01, 0)) + geom_point() + coord_cartesian(xlim = c(0, 773), ylim = c(0, 17))
            }
        } else {
            plot1 = ggplot(x, aes(x = phys.mb, y = cmmb)) + geom_line() + ggtitle(paste0(unique(individual.genetic.maps[[1]]$chromo)[count1], " ", templabels[count2])) + 
                theme_classic() + ylab("cM / Mb") + xlab("Physical position (Mb)") + scale_x_continuous(expand = c(0.01, 0)) + geom_point() + coord_cartesian(xlim = c(0, 773), ylim = c(0, 17))
        }
        
        
        
        # 
        # geom_text_repel(aes(label = ifelse(variable == "10(deg)C", marker, '')), color = "black", 
        #                                 hjust = 1, direction = "y", nudge_x = -0.1, size = 2, xlim = c(0, 1.5), ylim = c(-buffer.val1, (max(na.omit(x2$value)) + buffer.val1)), 
        #                                 segment.size = 0.2, segment.color = "grey50") + ylim(-buffer.val1, (max(na.omit(x2$value)) + buffer.val1))
        
        
        
        count1 <<- count1 + 1
        plot1
    })
    count2 <<- count2 + 1
    # plots1
    cslm3
    
})

y = axc.genedistplotsv2[[4]][[2]]


axc.genedistplotsv2[[1]][[1]]

#transpose list of lists
axc.genedistplotsv3 = lapply(1:length(axc.genedistplotsv2[[1]]), function(y){
    lapply(axc.genedistplotsv2, function(x){
        x[[y]]
    })
})

axcavg.r1 = axc.avg2

calculate.on.chromosome = function(avg.mrd.df, calc.function){
    avg.mrd.df$chromo = factor(avg.mrd.df$chromo, levels = unique(avg.mrd.df$chromo))
    axcmeans.per.chromo.long = as.data.frame(aggregate(avg.mrd.df$phys.dist.long, list(avg.mrd.df$treatment, avg.mrd.df$chromo), calc.function))
    axcmeans.per.chromo.long2 = dcast(axcmeans.per.chromo.long, Group.2 ~ Group.1)
    
    axcmeans.per.chromo.short = as.data.frame(aggregate(avg.mrd.df$phys.dist.short, list(avg.mrd.df$treatment, avg.mrd.df$chromo), calc.function))
    axcmeans.per.chromo.short2 = dcast(axcmeans.per.chromo.short, Group.2 ~ Group.1)
    
    axcmeans.per.chromo.long3 = melt(axcmeans.per.chromo.long2)
    axcmeans.per.chromo.short3 = melt(axcmeans.per.chromo.short2)
    
    axcmeans.per.chromo.long3$arm = "Long"
    axcmeans.per.chromo.short3$arm = "Short"
    bind_rows(axcmeans.per.chromo.long3, axcmeans.per.chromo.short3)
}



axc.mrd.means = calculate.on.chromosome(axcavg.r1, function(x) mean(na.exclude(x)))
axc.mrd.sds = calculate.on.chromosome(axcavg.r1, function(x) sd(na.exclude(x)))

axc.mrd.means$sd = axc.mrd.sds$value

axcavg.r1$chromo = factor(axcavg.r1$chromo, levels = unique(axcavg.r1$chromo))
axcmeans.per.chromo.long = as.data.frame(aggregate(axcavg.r1$phys.dist.long, list(axcavg.r1$treatment, axcavg.r1$chromo), function(x) mean(na.exclude(x))))
axcmeans.per.chromo.long2 = dcast(axcmeans.per.chromo.long, Group.2 ~ Group.1)

axcmeans.per.chromo.short = as.data.frame(aggregate(axcavg.r1$phys.dist.short, list(axcavg.r1$treatment, axcavg.r1$chromo), function(x) mean(na.exclude(x))))
axcmeans.per.chromo.short2 = dcast(axcmeans.per.chromo.short, Group.2 ~ Group.1)

axcmeans.per.chromo.long3 = melt(axcmeans.per.chromo.long2)
axcmeans.per.chromo.short3 = melt(axcmeans.per.chromo.short2)

axcmeans.per.chromo.long3$arm = "Long"
axcmeans.per.chromo.short3$arm = "Short"

axc.mrd.both.arms = bind_rows(axcmeans.per.chromo.long3, axcmeans.per.chromo.short3)

axc.plot.mrd = ggplot(axc.mrd.means, aes(x = Group.2, y = value, fill = variable)) + 
    geom_bar(stat = 'identity', position = 'dodge') +
    scale_fill_viridis(discrete = T, labels = axc.pops) + xlab("Chromosome") + 
    ylab("Mean MRD") + labs(fill = "Population") +
    geom_errorbar(aes(ymin = value, ymax = (value + sd), group = variable), position = 'dodge', size = 0.2) +
    facet_grid(rows = vars(arm))

axc.plot.mrd


#### RECOMBINATION DISTRIBUTION PLOTS ####
all.plotsaxc1 = lapply(axc.genedistplotsv3, function(x){
    marginval = 0.1
    chromo = x[[1]]$chromo[1]
    
    plots1 = Map(function(y, z){
        ggplot(y, aes(x = phys.mb, y = cmdiff)) + geom_point() + geom_line() + 
            ylim(c(0, max(sapply(x, function(y) max(y$cmdiff))))) + xlim(c(0, (max(x[[1]]$phys.mb) + 20))) + 
            ylab("Recombination") + xlab("") + ylab("") + theme_bw() + ggtitle(z) + theme(plot.margin = unit(c(marginval, marginval, marginval, marginval), "cm")) +
            theme(plot.title = element_text(margin = margin(t = 0, b = -20), hjust = 0.01))
    }, x, axc.pops)
    
    #make plot showing marker names
    gdf1 = x[[1]]
    labeldf1 = data.frame(x[[1]]$marker, seq(0, max(x[[1]]$phys.mb) + 20, (max(x[[1]]$phys.mb) + 20) / nrow(gdf1))[1:nrow(gdf1)])
    labeldf1$y = 22
    colnames(labeldf1) = c("marker", "x", "y")
    labeldf2 = gdf1[, c(3, 6, 7)]
    colnames(labeldf2) = c("marker", "x", "y")
    labeldf1$group = 1:nrow(labeldf1)
    labeldf2$group = 1:nrow(labeldf2)
    labeldf2$y = 17
    labeldf2$y = gdf1[, c(3, 6, 7)]$cmdiff
    labeldf2$y = cummax(gdf1[, c(3, 6, 7)]$cmdiff) + 2
    
    # labeldf2$y[(which(c(0, diff(labeldf2$y)) > 0) + 1)] = labeldf2$y[(which(c(0, diff(labeldf2$y)) > 0) + 1)] + 2
    labeldf2$y = cummax(labeldf2$y)
    
    labeldf1.1 = labeldf1
    labeldf1.1$x = labeldf1.1$x + 13 #line label end nudge
    labeldf1.1$y = 15 #max line coordinate
    
    labeldf3 = bind_rows(labeldf1.1, labeldf2)
    
    label.subplot = x[[1]]
    label.subplot$y = 0
    label.subplot.line = labeldf3
    
    label.subplot.line$y[1:max(labeldf3$group)] = 10
    label.subplot.line$x[1:max(labeldf3$group)] = label.subplot.line$x[1:max(labeldf3$group)] - 2 # x buffer value
    label.subplot.line$y[(max(labeldf3$group) + 1):nrow(label.subplot.line)] = 0.2
    
    label.subplot.text = labeldf1
    label.subplot.text$y = 14
    
    
    label.subplot2 = ggplot(label.subplot, aes(x = phys.mb, y = y)) + geom_line() + ylim(c(0, 17)) + 
        xlim(c(0, max(x[[1]]$phys.mb) + 20)) + ylab("Recombination") + xlab("") + ylab("") + theme_classic() + ggtitle(paste0(chromo, " Markers")) +
        geom_line(data = label.subplot.line, aes(x = x, y = y, group = group), colour = "#a3a3a3", size = 0.1) + #angled lines
        #geom_line(data = labeldf4, aes(x = x, y = y, group = group), colour = "#a3a3a3", size = 0.1) + # vertical lines
        geom_text(data = label.subplot.text, aes(x = x, y = y, label = marker), angle = 330, size = 2) + geom_point() + 
        theme(plot.margin = unit(c(marginval, marginval, marginval, 0.55), "cm")) +
        theme(plot.title = element_text(margin = margin(t = 0, b = 0), hjust = 0.01)) + 
        theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    genedisthist = gene.dist.hists[[which(listofwheatchromosomes == chromo)]] + xlim(c(0, max(x[[1]]$phys.mb) + 20)) + 
        theme(axis.title.x = element_blank()) +
        theme(plot.margin = unit(c(marginval, marginval, marginval, marginval), "cm"))
    
    c(list(genedisthist), plots1, list(label.subplot2))
    
})

num1 = 2
grid.arrange(arrangeGrob(all.plotsaxc1[[num1]][[6]], 
                                                 all.plotsaxc1[[num1]][[5]], 
                                                 all.plotsaxc1[[num1]][[4]], 
                                                 all.plotsaxc1[[num1]][[3]], 
                                                 all.plotsaxc1[[num1]][[2]],
                                                 all.plotsaxc1[[num1]][[1]],
                                                 left = textGrob("Recombination (cM Diff)", rot = 90, vjust = 2),
                                                 bottom = textGrob("Physical Position (Mb)", vjust = -1),
                                                 ncol = 1))



#### MRD DISTRIBUTION FIGURE ####

axc.avg2.rename = axc.avg2
axc.avg2.rename$treatment = as.numeric(as.character(axc.avg2.rename$treatment))
axc.avg2.rename[which(axc.avg2.rename$treatment == 1), ]$treatment = "cxa651"
axc.avg2.rename[which(axc.avg2.rename$treatment == 2), ]$treatment = "cxa653"
axc.avg2.rename[which(axc.avg2.rename$treatment == 3), ]$treatment = "axc512"
axc.avg2.rename[which(axc.avg2.rename$treatment == 4), ]$treatment = "axc611"

axc.avg2.rename.rev = axc.avg2.rename[nrow(axc.avg2.rename):1, ]
axc.avg2.rename.rev$treatment = factor(axc.avg2.rename.rev$treatment, levels = c(axc.pops))
colnames(axc.avg2.rename.rev)[4:5] = c("Short arm", "Long arm")
axc.avg2.rename.rev2 = melt(axc.avg2.rename.rev, id = c("marker.dist", "phys.dist", "chromo", "individual", "treatment"))


histogram.average.distances.of.recomb.axc = ggplot(axc.avg2.rename.rev2, aes(x = value)) + geom_histogram(bins = 100) + 
    facet_grid(rows = vars(treatment), cols = vars(variable)) + theme_bw() + xlab("Physical distance along chromosome arm (%)") +
    ylab("Frequency")

histogram.average.distances.of.recomb.axc


#### COMPARE GENETIC MAP WITHOUT LISS ####


geno.centi.all = lapply(list.of.geno.dfs.each.pop3, extract.centimorgan.from.genotypes)

geno.centi.all = lapply(geno.centi.all, function(x){
    x$liss.keep = F
    x$liss.keep[which(x$marker %in% geno.centi1[[1]]$marker)] = T
    x
})

count2 = 1
templabels = axc.pops

geno.centi.all2 = lapply(geno.centi.all, function(y){
    add.label1 = T
    cslm2 = split(y, y$chromo)
    
    cslm3 = lapply(cslm2, function(x){
        x$cmdiff = c(0, diff(x$cm))
        x$xaxis = 1:nrow(x)
        x
    })
    
    cslm3
})


geno.centi.all3 = transpose.list.of.lists(geno.centi.all2)

geno.centi.all2.bound = lapply(geno.centi.all2, bind_rows)


allcmdiff = as.data.frame(lapply(geno.centi.all2.bound, function(x){
    x$cmdiff
}))

allcmdiff = as.data.frame(t(reset.colnames(allcmdiff)))

allcmdiff.dist = dist(allcmdiff)
allcmdiff.dist2 = as.data.frame(as.matrix(allcmdiff.dist))

fit <- cmdscale(allcmdiff.dist, eig=TRUE, k=2) # k is the number of dim
mds.df = as.data.frame(fit$points)
rownames(mds.df) = axc.pops
mds.df$pop = axc.pops

mds.all.plot.axc = ggplot(mds.df, aes(x = V1, y = V2)) + geom_point() + geom_text_repel(aes(label = pop)) + theme_bw()
mds.all.plot.axc

all.chromo.mds.plots = lapply(geno.centi.all3, function(z){
    allcmdiff = as.data.frame(lapply(z, function(x){
        x$cmdiff
    }))
    
    allcmdiff = as.data.frame(t(reset.colnames(allcmdiff)))
    
    allcmdiff.dist = dist(allcmdiff)
    
    fit <- cmdscale(allcmdiff.dist, eig=TRUE, k=2) # k is the number of dim
    mds.df = as.data.frame(fit$points)
    rownames(mds.df) = axc.pops
    mds.df$pop = axc.pops
    
    ggplot(mds.df, aes(x = V1, y = V2)) + geom_point() + geom_text_repel(aes(label = pop)) + theme_bw()
    
})



geno.centi.all.plots = lapply(geno.centi.all3, function(x){
    plots1 = lapply(x, function(y){
        ggplot(y, aes(x = xaxis, y = cmdiff, color = liss.keep)) + geom_point() + geom_line() + theme_bw() + scale_color_manual(values = c("#AAAAAA", "#000000"))
    })
    
})

num1 = 2
grid.arrange(geno.centi.all.plots[[num1]][[1]],
                         geno.centi.all.plots[[num1]][[2]],
                         geno.centi.all.plots[[num1]][[3]],
                         geno.centi.all.plots[[num1]][[4]],
                         ncol = 1)


do.call(grid.arrange, geno.centi.all.plots[[1]])


g = lapply(list.of.geno.dfs.each.pop3, detect.recombination, conservative = F)

v1 = sapply(g, function(x) sum(x$num.events))
axc.pops
v2 = sapply(geno.axc2, nrow)


Map(function(x, y){
    (y / x) * 100
}, v1, v2)



