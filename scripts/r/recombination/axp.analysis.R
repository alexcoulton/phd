#to be used in conjunction with recombination.distribution.09072019v2.R

#### RECOMBINATION TEMP SHIFT AxP ####

axptemps1 = c("10(deg)C", "14(deg)C", "26(deg)C", "28(deg)C")

source("rotation1scripts_v4/scripts/r/recombination/axp.analysis.full.map.reconstruction.R")

##axc##
# initial.data2 = list.of.geno.dfs.each.pop3
# initial.data2 = lapply(initial.data2, function(x){
#     colnames(x)[3:ncol(x)] = switch.affy.format(colnames(x)[3:ncol(x)])
#     x
# })
##/axc##







initial.data.comb = combinelist.of.genos(initial.data2)

initial.data.comb.map = extract.centimorgan.from.genotypes(initial.data.comb)
initial.data.comb.map = convert.map.format(initial.data.comb.map)

phys1 = get.phys.for.chromos(initial.data2[[1]], F, T, ".")

initial.data.comb.map$phys = phys1$phys.pos
initial.data.comb.map$phys.bp = phys1$phys.bp

initial.data.comb.map.split = split(initial.data.comb.map, initial.data.comb.map$chr)

lapply(initial.data.comb.map.split, function(x) which(duplicated(x$cM)))


#do marker filtering by longest increasing subsequence / chromosome coverage
axp.liss = qtlanalysis.all(initial.data.comb, initial.data.comb.map, "RIL2", T, ".")


initial.data2 = lapply(initial.data2, function(x){
    # browser()
    x[, c(1, 2, match(axp.liss[[2]]$marker, colnames(x)))]
})


geno.centi.temp = lapply(initial.data2, extract.centimorgan.from.genotypes)


geno.centi.temp7b = geno.centi.temp[[2]][which(geno.centi.temp[[2]]$chromo == "7B"), ]
                                                                                 
gc7b.bp = attributes(genetictophysical("7B", geno.centi.temp7b$marker, "."))$bp

geno.centi.temp7b$bp = gc7b.bp

geno.centi.temp7b[which(geno.centi.temp7b$bp == 253910467), ]


geno.centi.temp3b = geno.centi.temp[[2]][which(geno.centi.temp[[2]]$chromo == "3B"), ]

gc3b.bp = attributes(genetictophysical("3B", geno.centi.temp3b$marker, "."))$bp

geno.centi.temp3b$bp = gc3b.bp

geno.centi.temp3b[which(geno.centi.temp3b$bp == 442638399), ]

geno.centi.len.axp = lapply(geno.centi.temp, function(x){
    x$chromo = as.factor(x$chromo)
    levels(x$chromo) = unique(x$chromo)
    aggregate(x$cm, list(x$chromo), max)
})

chromosome.names.axp = unlist(lapply(geno.centi.temp, function(x) unique(x$chromo))[[1]])

total.map.length.axp = sapply(geno.centi.len.axp, function(x){
    sum(x$x)
})

#### GENETIC MAP FIGURE ####
plots2 = prepare.genetic.map.comparison.plots(initial.data2, geno.centi.temp, list.of.populations = axptemps1, labels.on = T)

plots2[[10]]
plots2[[5]]

gen.plot1 = grid.arrange(plots2[[10]], plots2[[5]], ncol = 2)
ggsave("rotation1scripts_v4/plots/recombination.paper/figure3v3.eps", gen.plot1, width = 18, height = 11.5, units = "cm", device = cairo_pdf)



axp.gc.len2 = comb.treatments(geno.centi.len.axp)
colnames(axp.gc.len2) = c("Chromosome", "cm", "Treatment")
axp.gc.len2$Treatment = as.factor(axp.gc.len2$Treatment)

axp.gc.len2$Treatment = factor(axp.gc.len2$Treatment, levels = unique(axp.gc.len2$Treatment))
aggregate(axp.gc.len2$cm, list(axp.gc.len2$Treatment), max)

axp.gc.len2 = relabel.treat.w.pop.names(axp.gc.len2, "Treatment", axptemps1, F)

#### PAPER FIGURE 4 MAP LENGTH ####
axp.chromo.maplengths.plot = ggplot(axp.gc.len2, aes(x = Chromosome, y = cm, group = Treatment, fill = Treatment)) + geom_bar(stat = 'identity', position = "dodge") +
    scale_fill_viridis(discrete = T) + theme_bw(base_size = 22) + ylab("Length of chromosome (cM)")



#get physical data for chromosomes

iniphys2 = get.phys.for.chromos(initial.data2[[2]])

#### DETECT RECOMBINATION ####

# recomb.before.liss.axp = lapply(initial.data, detect.recombination, conservative = F)
# lapply(recomb.before.liss.axp, nrow)


booboo1 = list()

geno.temp = lapply(initial.data2, detect.recombination, conservative = T, liss.df = axp.liss)

lapply(geno.temp, nrow)


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
booboo.rm = booboo2[which(booboo2$dist.double < 30), ]


if(length(grep("65-1", booboo.rm$individual)) > 10){
    #axc analysis
    booboo.rm$quad = 0
    booboo.rm$quad[grep("65-1", booboo.rm$individual)] = 1
    booboo.rm$quad[grep("65-3", booboo.rm$individual)] = 2
    booboo.rm$quad[grep("51-2", booboo.rm$individual)] = 3
    booboo.rm$quad[grep("61-1", booboo.rm$individual)] = 4
} else {
    #axp analysis
    booboo.rm$quad = multi.str.split(booboo.rm$individual, "_", 2)
    booboo.rm$quad = as.numeric(substring(booboo.rm$quad, 2, 2))
    
    booboo.rm$quad[which(booboo.rm$quad == 3)] = 7
    booboo.rm$quad[which(booboo.rm$quad == 2)] = 3
    booboo.rm$quad[which(booboo.rm$quad == 7)] = 2
}




#remove double events from the genotyping dataframes
lapply(1:nrow(booboo.rm), function(x){
    quad = booboo.rm[x, ]$quad
    marker.name = booboo.rm[x, ]$probe.after.transition
    individual.name = booboo.rm[x, ]$individual
    
    initial.data2[[quad]][which(initial.data2[[quad]]$probeset_id == individual.name), which(colnames(initial.data2[[quad]]) == marker.name)] <<- "-"
    
})




geno.temp = lapply(initial.data2, detect.recombination, conservative = T, liss.df = axp.liss)

lapply(geno.temp, nrow)


#### processing ####
lg.lengths = table(as.character(initial.data2[[1]][1, ]))



#add physical positions to recombination dataframe
geno.temp2 = lapply(geno.temp, function(x){
    phys.pos.marker1 = iniphys2[match(x$probe.before.transition, iniphys2$marker), ]$phys.pos
    phys.pos.marker1.bp = iniphys2[match(x$probe.before.transition, iniphys2$marker), ]$phys.bp
    
    phys.pos.marker2 = iniphys2[match(x$probe.after.transition, iniphys2$marker), ]$phys.pos
    phys.pos.marker2.bp = iniphys2[match(x$probe.after.transition, iniphys2$marker), ]$phys.bp
    
    phys.pos.midpoint = (phys.pos.marker1 + phys.pos.marker2) / 2
    phys.bp.midpoint = (phys.pos.marker1.bp + phys.pos.marker2.bp) / 2
    
    x$phys.pos = phys.pos.midpoint
    x$phys.bp = phys.bp.midpoint
    x
})

gg3 = lapply(geno.temp2, function(x){
    x[which(x$chromo == "2A"), ]$phys.bp
})

par(mfrow = c(4, 1))
hist(gg3[[4]], breaks = 1000, ylim = c(0, 25))
hist(gg3[[3]], breaks = 1000, ylim = c(0, 25))
hist(gg3[[2]], breaks = 1000, ylim = c(0, 25))
hist(gg3[[1]], breaks = 1000, ylim = c(0, 25))

# average.recomb.positions = lapply(geno.temp2, ind.avg)
# average.recomb.positions = lapply(geno.temp2, ind.avg, use.centromere = T)
average.recomb.positions = lapply(geno.temp2, ind.avg, use.centromere = T, scale.by.arm = T, manual.center = 50)

#add treatment number to each individual dataframe
for(i in 1:4){
    average.recomb.positions[[i]]$treatment = i
}

#bind dataframes for statistical analysis
avg.r1 = bind_rows(average.recomb.positions)

g = aov(avg.r1$phys.dist.short ~ as.factor(avg.r1$treatment))
summary(g)
TukeyHSD(g)


summary(aov(avg.r1$phys.dist.long ~ avg.r1$treatment))

avg.lm = lm(log(avg.r1$phys.dist.long) ~ avg.r1$treatment)
hist(resid(avg.lm))


axpavg.r1 = avg.r1
axpavg.r1$chromo = factor(axpavg.r1$chromo, levels = unique(axpavg.r1$chromo))
axpmeans.per.chromo.long = as.data.frame(aggregate(axpavg.r1$phys.dist.long, list(axpavg.r1$treatment, axpavg.r1$chromo), function(x) mean(na.exclude(x))))
axpmeans.per.chromo.long2 = dcast(axpmeans.per.chromo.long, Group.2 ~ Group.1)

axpmeans.per.chromo.short = as.data.frame(aggregate(axpavg.r1$phys.dist.short, list(axpavg.r1$treatment, axpavg.r1$chromo), function(x) mean(na.exclude(x))))
axpmeans.per.chromo.short2 = dcast(axpmeans.per.chromo.short, Group.2 ~ Group.1)

#### mean mrd per chromo ####
axpmeans.per.chromo.short2
axpmeans.per.chromo.long2

axp.split1 = unique(avg.r1$chromo)[1:5]
axp.split2 = unique(avg.r1$chromo)[6:10]

avg.r1.2 = avg.r1[which(avg.r1$chromo %in% axp.split1), ]
avg.r1.3 = avg.r1[which(avg.r1$chromo %in% axp.split2), ]


hist.per.chromoaxp2 = ggplot(avg.r1.2, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("") +
    ylab("Frequency")

hist.per.chromoaxp3 = ggplot(avg.r1.3, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("") +
    ylab("Frequency") + xlab("Physical Distance From Center of Chromosome (%)")

grid.arrange(hist.per.chromoaxp2, hist.per.chromoaxp3, ncol = 1)




# avg.r1 = relabel.treat.w.pop.names(avg.r1, "treatment", axptemps1, reverse.treatment.order = T)
#### PAPER FIGURE 1 ####
avg.r1[which(avg.r1$treatment == 1), ]$treatment = "10(deg)C"
avg.r1[which(avg.r1$treatment == 2), ]$treatment = "14(deg)C"
avg.r1[which(avg.r1$treatment == 3), ]$treatment = "26(deg)C"
avg.r1[which(avg.r1$treatment == 4), ]$treatment = "28(deg)C"

avg.r1.rev = avg.r1[nrow(avg.r1):1, ]
avg.r1.rev$treatment = factor(avg.r1.rev$treatment, levels = c("28(deg)C", "26(deg)C", "14(deg)C", "10(deg)C"))
colnames(avg.r1.rev)[4:5] = c("Short arm", "Long arm")
avg.r1.rev2 = melt(avg.r1.rev, id = c("marker.dist", "phys.dist", "chromo", "individual", "treatment"))


histogram.average.distances.of.recomb.axp = ggplot(avg.r1.rev2, aes(x = value)) + geom_histogram(bins = 100) + 
    facet_grid(rows = vars(treatment), cols = vars(variable)) + theme_bw() + xlab("Physical distance along chromosome arm (%)") +
    ylab("Frequency")

histogram.average.distances.of.recomb.axp


#### kruskal test ####
# (1)

avg.r1$treatment[which(avg.r1$treatment == "10(deg)C")] = 1
avg.r1$treatment[which(avg.r1$treatment == "14(deg)C")] = 2
avg.r1$treatment[which(avg.r1$treatment == "26(deg)C")] = 3
avg.r1$treatment[which(avg.r1$treatment == "28(deg)C")] = 4

kruskal.test(avg.r1$phys.dist.short ~ avg.r1$treatment)
kruskal.test(avg.r1$phys.dist.long ~ avg.r1$treatment)

remove.treatment = function(nvector1, arm1, mrd.avg.df){
    if(missing(arm1)) arm1 = "none"
    if(missing(mrd.avg.df)) mrd.avg.df = avg.r1
    mrd.avg.df = mrd.avg.df[-which(mrd.avg.df$treatment %in% nvector1), ]
    if(arm1 == "none") return(kruskal.test(mrd.avg.df$phys.dist ~ mrd.avg.df$treatment))
    if(arm1 == "short") return(kruskal.test(mrd.avg.df$phys.dist.short ~ mrd.avg.df$treatment))
    if(arm1 == "long") return(kruskal.test(mrd.avg.df$phys.dist.long ~ mrd.avg.df$treatment))
}

#(3)

lapply(1:4, function(x) remove.treatment(x, "short"))
lapply(1:4, function(x) remove.treatment(x, "long"))
round(p.adjust(sapply(1:4, function(x) unclass(remove.treatment(x, "short"))$p.value), method = "bonf"), 4)
round(p.adjust(sapply(1:4, function(x) unclass(remove.treatment(x, "long"))$p.value), method = "bonf"), 4)

#get all pairwise combinations of numbers 1 to 4
combinations.to.rm = as.data.frame(combn(1:4, 2))

long.kruskal = lapply(combinations.to.rm, function(x){
    remove.treatment(x, "long")
})

g1 = round(sapply(long.kruskal, function(x) unclass(x)$p.value), digits = 5)
g1.1 = round(sapply(long.kruskal, function(x) unclass(x)$statistic), digits = 5)
g2 = round(p.adjust(sapply(long.kruskal, function(x) unclass(x)$p.value), "bonf"), digits = 5)

short.kruskal = lapply(combinations.to.rm, function(x){
    remove.treatment(x, "short")
})

g3 = round(sapply(short.kruskal, function(x) unclass(x)$p.value), digits = 5)
g3.1 = round(sapply(short.kruskal, function(x) unclass(x)$statistic), digits = 5)
g4 = round(p.adjust(sapply(short.kruskal, function(x) unclass(x)$p.value), "bonf"), digits = 5)

g4.2 = unlist(lapply(combinations.to.rm, function(x) paste(axptemps1[x], collapse = ", ")))

#### TABLE 1 ####
g5 = data.frame(g4.2, g1, g1.1, g2, g3, g3.1, g4)



colnames(g5) = c("Treatments removed", "p.long", "chi.long", "p.long.bonf", "p.short", "chi.short", "p.short.bonf")






all.kruskal = lapply(combinations.to.rm, function(x){
    remove.treatment(x, "none")
})

sapply(all.kruskal, function(x) unclass(x)$p.value)


avg.lm = lm(avg.r1$phys.dist ~ avg.r1$treatment)
hist(resid(avg.lm))

qqnorm(resid(avg.lm))
qqline(resid(avg.lm))

#get p values for all chromosomes
all.chromop1.long = unlist(sapply(unique(avg.r1$chromo), function(x){
    a3 = avg.r1[which(avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist.long ~ a3$treatment)[3]
    
}))

e1 = round(p.adjust(all.chromop1.long, "bonf"), digits = 5)

all.chromop1.short = unlist(sapply(unique(avg.r1$chromo), function(x){
    a3 = avg.r1[which(avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist.short ~ a3$treatment)[3]
    
}))

all.chromop1.botharms = unlist(sapply(unique(avg.r1$chromo), function(x){
    a3 = avg.r1[which(avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist ~ a3$treatment)[3]
    
}))


axp.long.and.short = c(all.chromop1.long, all.chromop1.short)

round(p.adjust(axp.long.and.short, method = "bonf"), digits = 5)
round(p.adjust(all.chromop1.botharms, method = "bonf"), digits = 5)

e2 = round(p.adjust(all.chromop1.short, "bonf"), digits = 5)

#get chi-square values for all chromosomes
e3 = unlist(sapply(unique(avg.r1$chromo), function(x){
    a3 = avg.r1[which(avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist.long ~ a3$treatment)[1]
    
}))

e4 = unlist(sapply(unique(avg.r1$chromo), function(x){
    a3 = avg.r1[which(avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist.short ~ a3$treatment)[1]
    
}))

e5 = unlist(sapply(unique(avg.r1$chromo), function(x){
    a3 = avg.r1[which(avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist.short ~ a3$treatment)[2]
    
}))

#### TABLE 2 ####
e6 = data.frame(e1, e3, e2, e4, e5)
colnames(e6) = c("p.long", "chi.long", "p.short", "chi.short", "df")
e6$chr = multi.str.split(rownames(e6), "\\.", 1)

#### mean mrd all chromo ####

mean.mrd.axp.long = aggregate(avg.r1$phys.dist.long, list(factor(avg.r1$treatment, levels = unique(avg.r1$treatment))), function(x) mean(na.omit(x)))
sd.mrd.axp.long = aggregate(avg.r1$phys.dist.long, list(factor(avg.r1$treatment, levels = unique(avg.r1$treatment))), function(x) sd(na.omit(x)))

mean.mrd.axp.short = aggregate(avg.r1$phys.dist.short, list(factor(avg.r1$treatment, levels = unique(avg.r1$treatment))), function(x) mean(na.omit(x)))
sd.mrd.axp.short = aggregate(avg.r1$phys.dist.short, list(factor(avg.r1$treatment, levels = unique(avg.r1$treatment))), function(x) sd(na.omit(x)))

stats.mrd.axp = bind_cols(mean.mrd.axp.long, sd.mrd.axp.long)
stats.mrd.axp = stats.mrd.axp[, -3]
colnames(stats.mrd.axp) = c("treatment", "mean", "sd")

stats.mrd.rep = apply(stats.mrd.axp, MARGIN = 1, function(x){
    paste(round(as.numeric(x[2]), digits = 2), "?", round(as.numeric(x[3]), digits = 2))
})
stats.mrd.rep2 = paste0(paste(stats.mrd.rep[1:3], collapse = ", "), " and ", stats.mrd.rep[4])


#(2)

# as.data.frame(round(p.adjust(c(all.chromop1, all.chromop1.cons), "bonf"), digits = 5))
# 
# all.chromop2 = as.data.frame(round(p.adjust(all.chromop1, "bonf"), digits = 5))
# all.chromop.fdr = as.data.frame(round(p.adjust(all.chromop1, "fdr"), digits = 5))
# colnames(all.chromop2) = "p-value"
# rownames(all.chromop2) = multi.str.split(rownames(all.chromop2), ".p.value", 1)

#do some stats
kruskal.test(list(average.recomb.positions[[1]]$phys.dist, average.recomb.positions[[2]]$phys.dist, average.recomb.positions[[3]]$phys.dist, average.recomb.positions[[4]]$phys.dist))

wilcox.test(average.recomb.positions[[2]]$phys.dist, average.recomb.positions[[4]]$phys.dist)



lapply(average.recomb.positions, function(x){
    mean(na.omit(x$phys.dist))
})

lapply(average.recomb.positions, function(x){
    sd(na.omit(x$phys.dist))
})


hist(average.recomb.positions[[1]]$phys.dist)

ks.test(average.recomb.positions[[1]]$phys.dist, "pnorm")

ks.test(average.recomb.positions[[1]]$marker.dist, "pnorm")


hist(average.recomb.positions[[1]]$phys.dist)

# qqplot(average.recomb.positions[[1]]$phys.dist)


hist((average.recomb.positions[[1]]$phys.dist))

hist(average.recomb.positions[[1]]$marker.dist)





histlists.temp = lapply(geno.temp2, makehistlist.new, physical = T)


perform.test(histlists.temp, c(1, 2), T)

round(p.adjust(perform.test(histlists.temp, c(1, 2), T), "bonf"), digits = 5)
round(p.adjust(perform.test(histlists.temp, c(2, 3), T), "bonf"), digits = 5)
round(p.adjust(perform.test(histlists.temp, c(2, 4), T), "bonf"), digits = 5)


par(mfrow = c(4, 1))
num1 = 2
hist(histlists.temp[[4]][[num1]], breaks = 1000, ylim = c(0, 30))
hist(histlists.temp[[3]][[num1]], breaks = 1000, ylim = c(0, 30))
hist(histlists.temp[[2]][[num1]], breaks = 1000, ylim = c(0, 30))
hist(histlists.temp[[1]][[num1]], breaks = 1000, ylim = c(0, 30))


#### recombination frequency ####

axp.recomb.freq1 = lapply(geno.temp2, function(x){
    # browser()
    x = arrange(x, individual)
    aggregate(x$num.events, list(factor(x$individual, levels = unique(x$individual))), sum)
})

axp.recomb.freq1 = lapply(geno.temp2, function(x){
    # browser()
    x = arrange(x, individual)
    aggregate(x$num.events, list(factor(x$individual, levels = unique(x$individual))), sum)
})



axp.recomb.freq2 = comb.treatments(axp.recomb.freq1)

# hist(axp.recomb.freq1[[1]]$x)

axp.anova.freq = aov(axp.recomb.freq2$x ~ as.factor(axp.recomb.freq2$treatment))
axp.anova.freq.tuk = TukeyHSD(axp.anova.freq)
library(broom)
axp.aft2 = tidy(axp.anova.freq.tuk)
axp.aft2.p = round(axp.aft2$adj.p.value[c(1, 3, 4, 6)], digits = 5)

#### RECOMB FREQ ALL INDIVIDUALS ####
mean(axp.recomb.freq2$x)
sd(axp.recomb.freq2$x)

axp.aft.means = aggregate(axp.recomb.freq2$x, list(factor(axp.recomb.freq2$treatment, levels = unique(axp.recomb.freq2$treatment))), mean)
colnames(axp.aft.means) = c("treatment", "mean")
axp.aft.sds = aggregate(axp.recomb.freq2$x, list(factor(axp.recomb.freq2$treatment, levels = unique(axp.recomb.freq2$treatment))), sd)
colnames(axp.aft.sds) = c("treatment", "sd")

axp.aft3 = bind_cols(axp.aft.means, axp.aft.sds)

axp.aft3 = relabel.treat.w.pop.names(axp.aft3, "treatment", axptemps1, F)
axp.aft3 = relabel.treat.w.pop.names(axp.aft3, "treatment1", axptemps1, F)

#### PAPER FIGURE 5 ####
library(ggsignif)
axp.recomb.freq.plot1 = ggplot(axp.aft3, aes(x = treatment, y = mean)) + geom_bar(stat = 'identity', fill = "#e6e6e6", color = "#000000") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.5) + 
    xlab("Population") + ylab("Recombination frequency") + theme_bw(base_size = 22)

axp.recomb.freq.plot1 = add.comparison.bars(axp.recomb.freq.plot1, axp.anova.freq.tuk, 30, axptemps1)
axp.recomb.freq.plot1

#### GENE DISTRIBUTION ANALYSIS ####

individual.genetic.maps = lapply(initial.data2, function(x){
    g = extract.centimorgan.from.genotypes(x)
    g$phys2 = iniphys2$phys.bp
    g$phys3 = iniphys2$phys.pos
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
templabels = c("10(deg)C", "14(deg)C", "26(deg)C", "28(deg)C")

axp.genedistplotsv2 = lapply(individual.genetic.maps, function(y){
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

y = axp.genedistplotsv2[[4]][[2]]


axp.genedistplotsv2[[1]][[1]]

#transpose list of lists
axp.genedistplotsv3 = lapply(1:length(axp.genedistplotsv2[[1]]), function(y){
    lapply(axp.genedistplotsv2, function(x){
        x[[y]]
    })
})

axpmeans.per.chromo.long3 = melt(axpmeans.per.chromo.long2)
axpmeans.per.chromo.short3 = melt(axpmeans.per.chromo.short2)

axpmeans.per.chromo.long3$arm = "Long"
axpmeans.per.chromo.short3$arm = "Short"

axp.mrd.both.arms = bind_rows(axpmeans.per.chromo.long3, axpmeans.per.chromo.short3)

plot.longarmsmrd = ggplot(axpmeans.per.chromo.long3, aes(x = Group.2, y = value, fill = variable)) + geom_bar(stat = 'identity', position = 'dodge') +
    scale_fill_viridis(discrete = T, labels = axptemps1) + xlab("Chromosome") + ylab("Mean MRD") + labs(fill = "Temperature") + ggtitle("Long arms") + theme_bw(base_size = 22)


plot.shortarmsmrd = ggplot(axpmeans.per.chromo.short3, aes(x = Group.2, y = value, fill = variable)) + geom_bar(stat = 'identity', position = 'dodge') +
    scale_fill_viridis(discrete = T, labels = axptemps1) + xlab("Chromosome") + ylab("Mean MRD") + labs(fill = "Temperature") + ggtitle("Short arms") + theme_bw(base_size = 22)

#### MRD PLOTS ####
grid.arrange(plot.longarmsmrd, plot.shortarmsmrd)


#### RECOMBINATION DISTRIBUTION PLOTS ####
all.plots1 = lapply(axp.genedistplotsv3, function(x){
    marginval = 0.1
    chromo = x[[1]]$chromo[1]
    xlimbuffer = 30
    
    plots1 = Map(function(y, z){
        ggplot(y, aes(x = phys.mb, y = cmdiff)) + geom_point() + geom_line() + 
            ylim(c(0, max(sapply(x, function(y) max(y$cmdiff))))) + xlim(c(0, (max(x[[1]]$phys.mb) + xlimbuffer))) + 
            ylab("Recombination") + xlab("") + ylab("") + theme_bw(base_size = 22) + ggtitle(z) + theme(plot.margin = unit(c(marginval, marginval, marginval, marginval), "cm")) +
            theme(plot.title = element_text(margin = margin(t = 0, b = -20), hjust = 0.01, vjust = -1, size = 17))
    }, x, axptemps1)
    
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
    labeldf1.1$x = labeldf1.1$x + 15 #line label end nudge
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
        xlim(c(0, max(x[[1]]$phys.mb) + xlimbuffer)) + ylab("Recombination") + xlab("") + ylab("") + theme_classic(base_size = 22) + ggtitle(paste0(chromo, " Markers")) +
        geom_line(data = label.subplot.line, aes(x = x, y = y, group = group), colour = "#a3a3a3", size = 0.1) + #angled lines
        #geom_line(data = labeldf4, aes(x = x, y = y, group = group), colour = "#a3a3a3", size = 0.1) + # vertical lines
        geom_text(data = label.subplot.text, aes(x = x, y = y, label = marker), angle = 330, size = 2) + geom_point() + 
        theme(plot.margin = unit(c(marginval, marginval, marginval, 0.55), "cm")) +
        theme(plot.title = element_text(margin = margin(t = 0, b = 0), hjust = 0.01)) + 
        theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    genedisthist = gene.dist.hists[[which(listofwheatchromosomes == chromo)]] + xlim(c(0, max(x[[1]]$phys.mb) + xlimbuffer)) + 
        theme(axis.title.x = element_blank()) +
        theme(plot.margin = unit(c(marginval, marginval, marginval, marginval), "cm")) 
    
    c(list(genedisthist), plots1, list(label.subplot2))
    
})

all.plots1[[1]][[1]]
all.plots1[[1]][[2]]
all.plots1[[1]][[6]]
num1 = 1
grid.arrange(arrangeGrob(all.plots1[[num1]][[6]], 
                                                 all.plots1[[num1]][[5]], 
                                                 all.plots1[[num1]][[4]], 
                                                 all.plots1[[num1]][[3]], 
                                                 all.plots1[[num1]][[2]],
                                                 all.plots1[[num1]][[1]],
                                                 left = textGrob("Recombination (cM Diff)", rot = 90, vjust = 2),
                                                 bottom = textGrob("Physical Position (Mb)", vjust = -1),
                                                 ncol = 1))






genedisthist6b = gene.dist.hists[[17]] + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(0, 730)) + 
    xlab("Physical Position (Mb)") + ggtitle("6B HC Gene Distribution")
genedisthist1a = gene.dist.hists[[1]] + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(0, 600)) + 
    xlab("Physical Position (Mb)") + ggtitle("1A HC Gene Distribution")
genedisthist1a
genedisthist3b = gene.dist.hists[[8]] + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(-50, 900)) + 
    xlab("Physical Position (Mb)") + ggtitle("3B HC Gene Distribution")
genedisthist3b

genedisthist2a = gene.dist.hists[[4]] + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(0, 773)) + 
    xlab("Physical Position (Mb)") + ggtitle("2A HC Gene Distribution")
genedisthist2a

length(list.of.chromosome.assignments.for.lgs)





#### prepare full genetic map comparison ####

comb.data.reduced = comb.data6[, c(1, 2, which(colnames(comb.data6) %in% switch.affy.format(c5$marker)))]
comb.data.reduced = comb.data.reduced[, c(1, 2, na.omit(match(switch.affy.format(c5$marker), colnames(comb.data.reduced))))]
comb.data.reduced[1, 3:ncol(comb.data.reduced)] = c5$chr1



full.map.geno1 = comb.data.reduced[c(1, 2, 3, grep("Q1", comb.data.reduced$probeset_id)), ]
full.map.geno2 = comb.data.reduced[c(1, 2, 3, grep("Q3", comb.data.reduced$probeset_id)), ]
full.map.geno3 = comb.data.reduced[c(1, 2, 3, grep("Q2", comb.data.reduced$probeset_id)), ]
full.map.geno4 = comb.data.reduced[c(1, 2, 3, grep("Q4", comb.data.reduced$probeset_id)), ]

full.map.geno.all = list(full.map.geno1, full.map.geno2, full.map.geno3, full.map.geno4)

full.geno.centi.temp = lapply(full.map.geno.all, extract.centimorgan.from.genotypes)

plots.full2 = prepare.genetic.map.comparison.plots(full.map.geno.all, full.geno.centi.temp, list.of.populations = axptemps1)
test123 = prepare.genetic.map.comparison.plots(full.map.geno.all, full.geno.centi.temp, list.of.populations = axptemps1)
test123$`7B`
plots.full2[[1]]
plots.full2$`3B`
plots.full2$`7B`

plots.full2$`6B`

plots.full2$`2A`
grid.arrange(plots2$`6B`, plots.full2$`6B`, ncol = 2)


grid.arrange(plots2$`6B`, plots.full2$`6B`, ncol = 2)


#### MAP COMPARISON ####

allen.axp = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic map AxP.csv", stringsAsFactors = F)
allen.axp2 = allen.axp[which(allen.axp$chr %in% unique(c5$chr1)), ]

#compare total shared markers (before genetic mapping)
all.markers.axpf2v1 = read.csv("rotation1scripts_v4/original_data/genotype.data/manuallycuratedsnps.csv", stringsAsFactors = F)
allmf2 = switch.affy.format(colnames(all.markers.axpf2v1))[-(1:2)]
length(which(allmf2 %in% allen.axp$marker))

#markers that were in my map but not Sacha's
markers.not.shared = allmf2[which(!allmf2 %in% allen.axp$marker)]

#markers that were in Sacha's map but not mine
markers.not.shared2 = allen.axp$marker[which(!allen.axp$marker %in% allmf2)]

allmf2

length(which(c5$marker %in% allen.axp$marker))


axpf2.genotyping = read.table("rotation1scripts_v4/original_data/recombination.distribution.chap.results/genotyping.comparison/axpf2.results.txt", sep = "\t", header = T)
axpf5.genotyping = read.table("rotation1scripts_v4/original_data/recombination.distribution.chap.results/genotyping.comparison/axpf5.results.txt", sep = "\t", header = T)

table(axpf5.genotyping[match(markers.not.shared, axpf5.genotyping$probeset_id), ]$ConversionType)
table(axpf2.genotyping[match(allmf2, axpf2.genotyping$probeset_id), ]$ConversionType)

table(axpf5.genotyping[match(allen.axp$marker, axpf5.genotyping$probeset_id), ]$ConversionType)



length(which(allen.axp == "AB"))
length(which(allen.axp == "AA"))
length(which(allen.axp == "BB"))

which(is.na(match(allen.axp$marker, axpf5.genotyping$probeset_id)))

allen.genotyping.info = axpf5.genotyping[na.omit(match(allen.axp$marker, axpf5.genotyping$probeset_id)), ]

table(axpf5.genotyping[na.omit(match(allen.axp$marker, axpf5.genotyping$probeset_id)), ]$ConversionType)

allen.genotyping.info[which(allen.genotyping.info$ConversionType == "NoMinorHom"), ][1:10, 1:4]
allen.genotyping.info[which(allen.genotyping.info$ConversionType == "CallRateBelowThreshold"), ][1:10, 1:4]
allen.other.markers1 = allen.genotyping.info[which(allen.genotyping.info$ConversionType == "Other"), ]$probeset_id


allen.axp[which(allen.axp$marker == "AX-94866192"), ]

switch.affy.format(axp.liss[[2]]$marker) %in% allen.axp$marker


ini.v2 = split(c5, c5$chr1)
allen.v2 = split(allen.axp2, allen.axp2$chr)


length(which(switch.affy.format(axp.liss[[2]]$marker) %in% allen.axp$marker))
length(axp.liss[[2]]$marker)

#compare length of maps 

sapply(ini.v2, function(x) max(x$cm))
sapply(allen.v2, function(x) max(x$cM))


map.comparison1 = Map(function(x, y){
    g = data.frame(na.omit(match(x$marker, y$marker)))
    colnames(g) = "num1"
    ggplot(g, aes(x = 1:length(num1), y = num1)) + geom_point()
}, ini.v2, allen.v2)






#### compare marker order ####


filtered.map1 = axp.liss[[2]]
filtered.map1 = filtered.map1[, c(3, 2, 1, 2)]
colnames(filtered.map1) = c("cm", "chromo", "marker", "chr1")
filtered.map1$marker = switch.affy.format(filtered.map1$marker)
filtered.map2 = split(filtered.map1, filtered.map1$chr1)

#compare clustering between my map and Allen et al SSD F5 A x P map
lapply(filtered.map2, function(x){
    as.character(na.omit(unique(allen.axp[match(x$marker, allen.axp$marker), ]$chr)))
})

map.comp2 = lapply(filtered.map2, function(x){
    g = allen.axp[match(x$marker, allen.axp$marker), ]$chr
    g[which(is.na(g))] = "NA"
    table(g)
})



#does the filtered map contain the same markers as the allen map
sapply(filtered.map2, function(x) x$chr1[1])
sapply(allen.v2, function(x) x$chr[1])

#examine colinearity of LISS vs Allen
count1 = 1
marker.order.plots2 = Map(function(x, y){
    y2 = y[, 1:3]
    x2 = x
    y3 = y2
    x2$map = "axpf2"
    y3$map = "axpf5"
    y3 = y3[, c(3, 2, 1, 4)]
    colnames(y3) = c("cm", "chromo", "marker", "map")
    x2 = x2[, c(1, 4, 3, 5)]
    colnames(x2)[2] = "chromo"
    
    if(length(which.max(na.omit(match(y3$marker, x2$marker)))) == 0) browser()
    
    if(which.max(na.omit(match(y3$marker, x2$marker))) < 30){
        y4 = y3
    } else {
        y4 = y3
    }
    
    
    #mark skeleton markers
    y5 = split(y4, y4$cm)
    y6 = lapply(y5, function(x){
        #sort markers within bins
        x = x[sort(x$marker, index.return = T)$ix, ]
        x$ske = ""
        x$ske[1] = T
        
        x
    })
    
    y7 = bind_rows(y6)
    
    x3 = split(x2, x2$cm)
    x4 = lapply(x3, function(y){
        #sort markers within bins
        y = y[sort(y$marker, index.return = T)$ix, ]
        y$ske = ""
        y$ske[1] = T
        
        y
    })
    
    x5 = bind_rows(x4)
    
    #only label markers at particular cM intervals
    y7$ske = ""
    g = seq(0, max(y7$cm), 5)
    y7$ske[sapply(g, function(x){
        g2 = abs(y7$cm - x)
        which.min(g2)
    })] = T
    
    y7[which(y7$ske == T), ]$ske[which(c(diff(y7[which(y7$ske == T), ]$cm), 0) < 2)] = F
    
    x5$ske = ""
    
    count1 <<- count1 + 1
    
    marker.matches1 = na.omit(match(x5$marker, y7$marker))
    
    if(sum(cumsum(diff(marker.matches1))) < -100){
        print("rev")
    
        x5 = reverse.genetic.map(x5, "cm")
    }
    
    y8 = y7
    y8 = split(y8, y8$cm)
    
    
    
    #need to perform matches with respect to the bins in the sacha map, not the strict order of markers themselves as markers within bins will not be ordered correctly.
    marker.matches1 = lapply(x5$marker, function(x){
        g = sapply(y8, function(y){
            length(which(x == y$marker))
        })
        
        g2 = which(g == 1)
        if(length(g2) == 0){
            matching.bin = NA
        } else {
            matching.bin = y8[[g2]]$cm[1]
        }
        
        data.frame(x, matching.bin)
    })
    
    chromo1 = x$chr1[1]
    
    marker.matches2 = bind_rows(marker.matches1)
    marker.matches2$axpf2cm = x5$cm
    marker.matches2$chromo = chromo1
    num.missing = length(which(is.na(marker.matches2$matching.bin)))
    
    
    marker.matches3 = marker.matches2[-which(is.na(marker.matches2$matching.bin)), ]
    marker.matches3$num.missing = num.missing
    marker.matches3$total.markers = nrow(marker.matches3)
    #normalize cM values?    
    marker.matches3$matching.bin = marker.matches3$matching.bin - min(marker.matches3$matching.bin)
    marker.matches3$axpf2cm = marker.matches3$axpf2cm - min(marker.matches3$axpf2cm)
    
    marker.matches3$nonmonotonic = F
    
    nonmono = which(!marker.matches3$matching.bin == cummax(marker.matches3$matching.bin))

    if(length(nonmono) > 0){
        marker.matches3$nonmonotonic[nonmono] = T
    }
    
    slope1 = (max(marker.matches3$matching.bin) - min(marker.matches3$matching.bin)) / (max(marker.matches3$axpf2cm) - min(marker.matches3$axpf2cm))
    
    yint1 = marker.matches3$matching.bin[1] - (slope1 * marker.matches3$axpf2cm[1])
    
    marker.matches3$slope = slope1
    marker.matches3$intercept = yint1
    
    marker.matches3
    
}, filtered.map2, allen.v2)
# }, ini.v2, allen.v2) #using map that contains all markers

marker.order.plots2

inversion.distances = lapply(marker.order.plots2, function(x){
    g = which(x$nonmonotonic == T)
    g2 = lapply(g, function(z){
        c(z-1, z)
    })
    
    sapply(g2, function(z){
        diff(x[z, ]$matching.bin)    
    })
    
    
})

abs(sd(unlist(inversion.distances)))

rsquared = sapply(marker.order.plots2, function(x){
    # cor.test(x$matching.bin, x$axpf2cm)
    g = lm(matching.bin ~ axpf2cm, x)
    unclass(summary(g))$r.squared
})

pvals1 = sapply(marker.order.plots2, function(x){
    # cor.test(x$matching.bin, x$axpf2cm)
    g = lm(matching.bin ~ axpf2cm, x)
    unclass(summary(g))$coefficients[2, 4]
})

res.hists = lapply(marker.order.plots2, function(x){
    # cor.test(x$matching.bin, x$axpf2cm)
    g = lm(matching.bin ~ axpf2cm, x)
    g2 = unclass(summary(g))$residuals
    g2 = data.frame(g2)
    ggplot(g2, aes(x = g2)) + geom_histogram()
})


p.adjust(pvals1, "bonf")

rsquared2 = data.frame(rsquared)
rsquared2$chromo = rownames(rsquared2)
rsquared2$nonmonotonic = F
rsquared2$rsquared = round(rsquared2$rsquared, digits = 3)

m1 = bind_rows(marker.order.plots2)

m2 = m1[match(unique(m1$slope), m1$slope), ]

#### NEW MAP COMP FIG ####
new.map.comp.fig1 = ggplot(m1, aes(x = axpf2cm, y = matching.bin, shape = nonmonotonic, color = nonmonotonic)) + geom_point() + scale_color_manual(values = c("#000000", "#a8a8a8")) +
    theme_bw() +
    ylab("A x P F5 Genetic distance (cM)") + xlab("A x P F2 Genetic distance (cM)") + facet_wrap(~ chromo) +
    geom_abline(data = m2, aes(slope = slope, intercept = intercept)) + 
    geom_text(data = rsquared2, x = 20, y = 270, color = "black", aes(label = rsquared)) + labs(color = "Non-monotonic", shape = "Non-monotonic")

new.map.comp.fig1

#### prepare cm/mb plot for Allen map ####

#convert allen map into rqtl format

#convert allen data into rQTL format
allen.axp.only.chr = allen.axp[-which(!allen.axp$chr %in% listofwheatchromosomes), ]
allen2 = as.data.frame(t(allen.axp.only.chr))
allen2 = convert.to.character.data.frame(allen2)

colnames(allen2) = allen2[1, ]
allen2 = allen2[-1, ]

allen2 = add.column.at.position(allen2, 0)
allen2$empty.col = rownames(allen2)

allen2 = add.column.at.position(allen2, 0)
allen2 = allen2[-2, ]
allen2[1, 2] = ""
allen2[2:nrow(allen2), 1] = 1:(nrow(allen2) - 1)
colnames(allen2)[1:2] = c("X", "probeset_id")
allen2[allen2 == "AA"] = "A"
allen2[allen2 == "AB"] = "H"
allen2[allen2 == "BB"] = "B"

source("rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R")

allen2.map = extract.centimorgan.from.genotypes(allen2)
allen2.map2 = convert.map.format(allen2.map)

allen2.map.split = split(allen2.map2, allen2.map2$chr)

allen.liss = qtlanalysis.all(allen2, allen2.map2, "RIL5", T)

allen3 = allen2[, c(1, 2, match(allen.liss[[2]]$marker, colnames(allen2)))]


#### attempt 2 prepare cm/mb plot for allen map ####

allen.iniphys2 = get.phys.for.chromos(allen3, marker.format = "-", interpolate.na = F)

allen.individual.genetic.maps = lapply(list(allen3), function(x){
    g = extract.centimorgan.from.genotypes(x)
    g$phys2 = allen.iniphys2$phys.bp
    g$phys3 = allen.iniphys2$phys.pos
    g
}) 

allen.indi.maps.test = split(allen.individual.genetic.maps[[1]], allen.individual.genetic.maps[[1]]$chromo)

lapply(allen.indi.maps.test, function(x){
    sequence1 = seq(0, 100, 25)
    
    distances1 = sapply(sequence1, function(z){
        min(abs(z - x$phys3))
    })
    
    length(which(distances1 > 25))
})



count2 = 1
templabels = c("10(deg)C", "14(deg)C", "26(deg)C", "28(deg)C")
allen.axp.genedistplots = lapply(allen.individual.genetic.maps, function(y){
    cslm2 = split(y, y$chromo)
    cslm2 = lapply(cslm2, function(x){
        x$phys.mb = x$phys2 / 1000000
        x
    })
    
    
    cslm3 = lapply(cslm2, function(x){
        x$cmdiff = c(diff(x$cm), 0)
        x$physdiff = c(diff(x$phys.mb), 0)
        x$cmmb = x$cmdiff / x$physdiff
        x
    })
    
    count1 = 1
    plots1 = lapply(cslm3, function(x){
        
        plot1 = ggplot(x, aes(x = phys.mb, y = cmmb)) + geom_line() + ggtitle(paste0(unique(allen.individual.genetic.maps[[1]]$chromo)[count1], " ", "A x P F5")) + 
            theme_bw() + ylab("cM / Mb") + xlab("Physical position (Mb)") + scale_x_continuous(expand = c(0.01, 0)) + coord_cartesian(xlim = c(0, 850)) + geom_point()
        count1 <<- count1 + 1
        plot1
    })
    count2 <<- count2 + 1
    plots1
    
})




#### ASSESS MRD FOR F5 POPULATION ####


#prepare allen marker set that only contains markers in the F2 marker set



allen4 = allen3[, (which(switch.affy.format(colnames(allen3)) %in% colnames(initial.data2[[1]])))]


allen4.map = extract.centimorgan.from.genotypes(allen4)
allen4.map2 = convert.map.format(allen4.map)

allen.liss = qtlanalysis.all(allen4, allen2.map2, "RIL5", T)

allen5 = allen4[, c(1, 2, match(allen.liss[[2]]$marker, colnames(allen4)))]




init.shared.order = initial.data2[[4]][, which(switch.affy.format(colnames(initial.data2[[2]])) %in% colnames(allen5))]

colnames(allen5)[3:ncol(allen5)] = switch.affy.format(colnames(allen5)[3:ncol(allen5)])

all(match(switch.affy.format(colnames(init.shared.order)), colnames(allen5)) == 1:135) # colnames are already in exactly the same order after selecting shared markers. 

initial.data.shared = lapply(initial.data2, function(x){
    x[, which(colnames(x) %in% colnames(allen5))]
})

#dataframe for comparison
initial.data.shared2 = c(initial.data.shared, list(allen5))


shared1.geno.centi.temp = lapply(initial.data.shared2, extract.centimorgan.from.genotypes)

shared1.geno.centi.len.axp = lapply(shared1.geno.centi.temp, function(x){
    x$chromo = as.factor(x$chromo)
    levels(x$chromo) = unique(x$chromo)
    aggregate(x$cm, list(x$chromo), max)
})

shared1.total.map.length.axp = sapply(shared1.geno.centi.len.axp, function(x){
    sum(x$x)
})


#get physical data for chromosomes

shared1.iniphys2 = get.phys.for.chromos(initial.data.shared2[[1]])

# SET LIBERAL / CONSERVATIVE RECOMBINATION DETECTION
shared1.geno.temp = lapply(initial.data.shared2, detect.recombination, conservative = F)

shared1.lg.lengths = table(as.character(initial.data.shared2[[1]][1, ]))

#add physical positions to recombination dataframe
shared1.geno.temp2 = lapply(shared1.geno.temp, function(x){
    x$phys.pos = shared1.iniphys2[match(x$probe.before.transition, shared1.iniphys2$marker), ]$phys.pos
    x
})



shared1.average.recomb.positions = lapply(shared1.geno.temp2, ind.avg, use.centromere = T)

#add treatment number to each individual dataframe
for(i in 1:5){
    shared1.average.recomb.positions[[i]]$treatment = i
}

#bind dataframes for statistical analysis
shared1.avg.r1 = bind_rows(shared1.average.recomb.positions)

shared1.axpavg.r1 = shared1.avg.r1
shared1.axpavg.r1$chromo = factor(shared1.axpavg.r1$chromo, levels = unique(shared1.axpavg.r1$chromo))
shared1.axpmeans.per.chromo = as.data.frame(aggregate(shared1.axpavg.r1$phys.dist, list(shared1.axpavg.r1$treatment, shared1.axpavg.r1$chromo), function(x) mean(na.exclude(x))))
shared1.axpmeans.per.chromo2 = dcast(shared1.axpmeans.per.chromo, Group.2 ~ Group.1)






shared1.axp.split1 = unique(shared1.avg.r1$chromo)[1:5]
shared1.axp.split2 = unique(shared1.avg.r1$chromo)[6:10]

shared1.avg.r1.2 = shared1.avg.r1[which(shared1.avg.r1$chromo %in% shared1.axp.split1), ]
shared1.avg.r1.3 = shared1.avg.r1[which(shared1.avg.r1$chromo %in% shared1.axp.split2), ]

# avg.r1.2 = relabel.treat.w.pop.names(avg.r1.2, "treatment", axptemps1, T)
# avg.r1.3 = relabel.treat.w.pop.names(avg.r1.3, "treatment", axptemps1, T)


shared1.hist.per.chromoaxp2 = ggplot(shared1.avg.r1.2, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("") +
    ylab("Frequency")

shared1.hist.per.chromoaxp3 = ggplot(shared1.avg.r1.3, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment), cols = vars(chromo)) + theme_bw() + xlab("") +
    ylab("Frequency") + xlab("Physical Distance From Center of Chromosome (%)")

grid.arrange(shared1.hist.per.chromoaxp2, shared1.hist.per.chromoaxp3, ncol = 1)




# avg.r1 = relabel.treat.w.pop.names(avg.r1, "treatment", axptemps1, reverse.treatment.order = T)

shared1.histogram.average.distances.of.recomb.axp = ggplot(avg.r1, aes(x = phys.dist)) + geom_histogram(bins = 100) + facet_grid(rows = vars(treatment)) + theme_bw() + xlab("Physical Distance From Center of Chromosome (%)") +
    ylab("Frequency")





# (1)
kruskal.test(shared1.avg.r1$phys.dist ~ shared1.avg.r1$treatment)

shared1.remove.treatment = function(nvector1){
    avg.r2 = shared1.avg.r1[-which(shared1.avg.r1$treatment %in% nvector1), ]
    kruskal.test(avg.r2$phys.dist ~ avg.r2$treatment)    
}

#### F5 significant for all chromo ####
#(3)
shared1.remove.treatment(c(3, 4))
shared1.remove.treatment(c(1, 2))[[3]] * 2 #bonferroni correction for two statistical tests



remove.treatment(c(3, 4))


shared1.avg.lm = lm(shared1.avg.r1$phys.dist ~ shared1.avg.r1$treatment)
hist(resid(shared1.avg.lm))

qqnorm(resid(avg.lm))
qqline(resid(avg.lm))


shared1.all.chromop1 = unlist(sapply(unique(shared1.avg.r1$chromo), function(x){
    a3 = shared1.avg.r1[which(shared1.avg.r1$chromo == x), ]
    
    kruskal.test(a3$phys.dist ~ a3$treatment)[3]
    
}))

shared1.mean.mrd.axp = aggregate(avg.r1$phys.dist, list(factor(avg.r1$treatment, levels = unique(avg.r1$treatment))), function(x) mean(na.omit(x)))
shared1.sd.mrd.axp = aggregate(avg.r1$phys.dist, list(factor(avg.r1$treatment, levels = unique(avg.r1$treatment))), function(x) sd(na.omit(x)))
shared1.stats.mrd.axp = bind_cols(mean.mrd.axp.long, sd.mrd.axp.long)
shared1.stats.mrd.axp = stats.mrd.axp[, -3]
colnames(stats.mrd.axp) = c("treatment", "mean", "sd")

shared1.stats.mrd.rep = apply(stats.mrd.axp, MARGIN = 1, function(x){
    paste(round(as.numeric(x[2]), digits = 2), "?", round(as.numeric(x[3]), digits = 2))
})
shared1.stats.mrd.rep2 = paste0(paste(stats.mrd.rep[1:3], collapse = ", "), " and ", stats.mrd.rep[4])


#(2)
# shared1.all.chromop2 = as.data.frame(round(p.adjust(all.chromop1, "bonf"), digits = 5))
# shared1.all.chromop.fdr = as.data.frame(round(p.adjust(all.chromop1, "fdr"), digits = 5))
# colnames(all.chromop2) = "p-value"
# rownames(all.chromop2) = multi.str.split(rownames(all.chromop2), ".p.value", 1)

#do some stats
kruskal.test(list(average.recomb.positions[[1]]$phys.dist, average.recomb.positions[[2]]$phys.dist, average.recomb.positions[[3]]$phys.dist, average.recomb.positions[[4]]$phys.dist))

wilcox.test(average.recomb.positions[[2]]$phys.dist, average.recomb.positions[[4]]$phys.dist)



lapply(average.recomb.positions, function(x){
    mean(na.omit(x$phys.dist))
})

lapply(average.recomb.positions, function(x){
    sd(na.omit(x$phys.dist))
})


hist(average.recomb.positions[[1]]$phys.dist)

ks.test(average.recomb.positions[[1]]$phys.dist, "pnorm")

ks.test(average.recomb.positions[[1]]$marker.dist, "pnorm")


hist(average.recomb.positions[[1]]$phys.dist)

# qqplot(average.recomb.positions[[1]]$phys.dist)


hist((average.recomb.positions[[1]]$phys.dist))

hist(average.recomb.positions[[1]]$marker.dist)





shared1.histlists.temp = lapply(geno.temp2, makehistlist.new, physical = T)


perform.test(histlists.temp, c(1, 2), T)

round(p.adjust(perform.test(histlists.temp, c(1, 2), T), "bonf"), digits = 5)
round(p.adjust(perform.test(histlists.temp, c(2, 3), T), "bonf"), digits = 5)
round(p.adjust(perform.test(histlists.temp, c(2, 4), T), "bonf"), digits = 5)


par(mfrow = c(2, 2))
num1 = 6
hist(histlists.temp[[1]][[num1]], breaks = 100)
hist(histlists.temp[[2]][[num1]], breaks = 100)
hist(histlists.temp[[3]][[num1]], breaks = 100)
hist(histlists.temp[[4]][[num1]], breaks = 100)






#### marker distribution fig setup #### 

marker.dist1 = split(axp.liss[[2]], axp.liss[[2]]$chr)
md2 = axp.liss[[2]]
md2$phys.mb = md2$phys.bp / 1000000

md2.chromo = chromosomecounts2[which(chromosomecounts2$V2 %in% md2$chr), ]
md2.chromo$start = 0
colnames(md2.chromo) = c("end", "chr", "start")
md2.chromo2 = melt(md2.chromo)
md2.chromo2$value = md2.chromo2$value / 1000000
colnames(md2.chromo2) = c("chr", "variable", "phys.mb")

md2.centromere = centromeres
colnames(md2.centromere) = c("chr", "v2", "phys.mb", "phys.bp")
md2.centromere = md2.centromere[which(md2.centromere$chr %in% md2.chromo2$chr), ]


#get the mean distance between markers

md3 = split(md2, md2$chr)

lapply(md3, function(x){
    g = mean(diff(x$phys.mb))
    g2 = sd(diff(x$phys.mb))
    c(g, g2)
})

marker.dist.plot1 = ggplot(md2, aes(x = chr, y = phys.mb)) + geom_point(data = md2) + geom_line(data = md2.chromo2) + geom_point(data = md2.centromere, shape = 95, size = 15, color = "red") +
    xlab("Chromosome") + ylab("Physical position (Mb)") + theme_classic() + ggtitle("b") + scale_y_reverse()

#jpeg("rotation1scripts_v4/plots/recombination.distribution11072019/marker.dist.plot.jpg", 1000, 1000, res = 150)
#marker.dist.plot1
#dev.off()



#### marker density for full map ####

tt1 = get.phys.for.chromos(comb.data.reduced, marker.format = ".", interpolate.na = F)
nrow(tt1) - length(which(is.na(tt1$phys.pos)))


tt2 = bind_cols(full.geno.centi.temp[[1]], tt1)
tt3 = tt2[, -4]


full.marker.dist1 = split(tt3, tt3$chromo)
full.md2 = tt3
full.md2$phys.mb = full.md2$phys.bp / 1000000
full.md2 = full.md2[, c(3, 2, 1, 4, 5, 6)]
colnames(full.md2) = c("marker", "chr", "cM", "phys", "phys.bp", "phys.mb")

full.md2.chromo = chromosomecounts2[which(chromosomecounts2$V2 %in% full.md2$chr), ]
full.md2.chromo$start = 0
colnames(full.md2.chromo) = c("end", "chr", "start")
full.md2.chromo2 = melt(full.md2.chromo)
full.md2.chromo2$value = full.md2.chromo2$value / 1000000
colnames(full.md2.chromo2) = c("chr", "variable", "phys.mb")


full.md2.centromere = centromeres
colnames(md2.centromere) = c("chr", "v2", "phys.mb", "phys.bp")
full.md2.centromere = md2.centromere[which(md2.centromere$chr %in% full.md2.chromo2$chr), ]

marker.dist.plot.full = ggplot(full.md2, aes(x = chr, y = phys.mb)) + geom_point(data = full.md2) + geom_line(data = full.md2.chromo2) + geom_point(data = full.md2.centromere, shape = 95, size = 15, color = "red") +
    xlab("Chromosome") + ylab("Physical position (Mb)") + theme_classic() + ggtitle("a") + scale_y_reverse()

length(na.omit(full.md2$phys))

full.md3 = split(full.md2, full.md2$chr)

lapply(full.md3, function(x){
    g = mean(diff(na.omit(x$phys.mb)))
    g2 = sd(diff(na.omit(x$phys.mb)))
    c(g, g2)
})


#violin plot
ggplot(full.md2, aes(x = chr, y = phys.mb)) + geom_violin(data = full.md2) + geom_line(data = full.md2.chromo2) + geom_point(data = full.md2.centromere, shape = 95, size = 15, color = "red") +
    xlab("Chromosome") + ylab("Physical position (Mb)") + theme_classic()

#### MARKER DIST. FIGURE ####
comb.marker.dist.plot = grid.arrange(marker.dist.plot.full, marker.dist.plot1, ncol    = 2)

# jpeg("rotation1scripts_v4/plots/recombination.distribution11072019/comb.marker.dist.plotv2.jpg", 1000, 700, res = 150)
# grid.arrange(marker.dist.plot.full, marker.dist.plot1, ncol    = 2)
# dev.off()

#publication EPS figure
# ggsave("rotation1scripts_v4/plots/recombination.distribution11072019/comb.marker.dist.plotv2.eps", comb.marker.dist.plot, width = 18, height = 12, units = "cm", device = cairo_pdf)


#### GLOBAL HOTSPOT ANALYSIS ####


axp.genedistplotsv2[[4]][[1]]

hotspot.coords1 = lapply(axp.genedistplotsv2, function(x){
    coldspot.coord = lapply(x, function(y){
        which(y$cmdiff > 0)
    })
    coldspot.coord
})

#invert list of lists
hotspot.coords2 = lapply(1:10, function(x){
    lapply(hotspot.coords1, function(y){
        y[[x]]
    })
})

coldspot.coords1 = lapply(axp.genedistplotsv2, function(x){
    coldspot.coord = lapply(x, function(y){
        which(y$cmdiff == 0)
    })
    coldspot.coord
})

#invert list of lists
coldspot.coords2 = lapply(1:10, function(x){
    lapply(coldspot.coords1, function(y){
        y[[x]]
    })
})

temp.dep.hotspots = Map(function(x, y){
    low.temp.cold = x[[1]][which(x[[1]] %in% x[[2]])]
    high.temp.hot = y[[3]][which(y[[3]] %in% y[[4]])]
    
    which(high.temp.hot %in% low.temp.cold)
    
}, coldspot.coords2, hotspot.coords2)





high.hotspot = lapply(hotspot.coords2, function(x){
    low.hotspots = x[[1]][which(x[[1]] %in% x[[2]])]
    high.hotspots = x[[3]][which(x[[3]] %in% x[[4]])]
    high.hotspots[which(!high.hotspots %in% low.hotspots)]
    # low.hotspots[which(!low.hotspots %in% high.hotspots)]
})


Map(function(x, y){
    x[y, ]
}, axp.genedistplotsv2[[2]], temp.dep.hotspots)



axp.genedistplotsv2[[2]][[1]][12, ]


axp.invert = lapply(1:10, function(x){
    lapply(axp.genedistplotsv2, function(y){
        y[[x]][, c(2, 3, 4, 5, 7)]
    })
})


axp.invert2 = lapply(axp.invert, function(x){
    bind_cols(x)
})


axpi3 = lapply(axp.invert2, function(x){
    marker1 = x[which(x$cmdiff == 0 & x$cmdiff1 == 0 & x$cmdiff2 > 0 & x$cmdiff3 > 0),
        c(1, 2, 3, 4, 5, 10, 15, 20) ]
    
    marker.before = x[(which(x$cmdiff == 0 & x$cmdiff1 == 0 & x$cmdiff2 > 0 & x$cmdiff3 > 0) - 1),
        c(2, 3, 4) ]
    colnames(marker.before) = c("marker.before", "phys.bp.marker.before", "phys.per.marker.before")
    
    bind_cols(marker1, marker.before)
    
})

axpi4 = bind_rows(axpi3)

axpi4$phys.diff = axpi4$phys2 - axpi4$phys.bp.marker.before
axpi4$phys.diff.mb = axpi4$phys.diff / 1000000

hcgff = read.delim("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", sep = "\t", stringsAsFactors = F)
hcgff = hcgff[which(hcgff$V2 == "chr1A"), ]
hcgff2 = hcgff[which(hcgff$V5 > 368198237 & hcgff$V6 < 370097158), ]

axpi4.genes = lapply(1:nrow(axpi4), function(x){
    hcgff2 = hcgff[which(hcgff$V2 == paste0("chr", axpi4$chromo[x])), ]
    hcgff2[which(hcgff2$V5 > axpi4$phys.bp.marker.before[x] & hcgff2$V6 < axpi4$phys2[x]), ]
})

axpi4$genes = sapply(axpi4.genes, nrow)

axpi5 = axpi4[sort(axpi4$cmdiff2, index.return = T, decreasing = T)$ix, ]

axpi5$phys2 = round(axpi5$phys2 / 1000000, digits = 2)
axpi5$phys.bp.marker.before = round(axpi5$phys.bp.marker.before / 1000000, digits = 2)

axpi5$centro = chromosomecounts.v2$centro.mb[match(paste0("chr", axpi5$chromo), chromosomecounts.v2$V2)]
axpi5$midpoint = (axpi5$phys2 + axpi5$phys.bp.marker.before) / 2
axpi5$dist.from.centro = axpi5$midpoint - axpi5$centro

#### TABLE 3 ####

axpi5

#### FUNCTIONAL ANNO ####

# functional.anno = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", sep = "\t", header = T, fill = T, stringsAsFactors = F)
# functional.anno$gene2 = multi.str.split(as.character(functional.anno$Gene.ID), "\\.", 1)
# 

# 
# genes3b1 = multi.str.split(axpi4.genes[[9]]$V10, ";", 1)
# genes3b2 = multi.str.split(genes3b1, "=", 2)
# gene.coord = unlist(sapply(genes3b2, function(x){
#     grep(x, functional.anno$Gene.ID)
# }))
# 
# v(functional.anno[gene.coord, ])
# 
# 
# 
# axp.genes2 = bind_rows(axpi4.genes)
# 
# library(Biostrings)
# 
# cds = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")
# 
# gene1 = multi.str.split(axp.genes2$V10, ";", 1)
# gene2 = multi.str.split(gene1, "=", 2)
# (370097158 - 368198237) / 1000000
# 
# cds2 = cds[unlist(lapply(gene2, function(x) grep(x, names(cds))))]
# 
# cds3 = cds2[-which(duplicated(multi.str.split(names(cds2), "gene=", 2)))]
# 
# writeLines(names(cds3), "rotation1scripts_v4/original_data/recombination.distribution.chap.results/allhotspotgenenames.txt")
# 
# writeXStringSet(cds3, "rotation1scripts_v4/original_data/recombination.distribution.chap.results/all.temp.dependent.hotspot.genes.fa")
# 
# 
# 
# 
# library(parallel)
# test.1 = mclapply(gene2, function(x){
#     nrow(functional.anno[grep(x, functional.anno$Gene.ID), ])
# }, mc.cores = 50)
# 
# 
# grep(gene2[[3]], functional.anno$Gene.ID)
# 
# 
# 
# go.enrich = read.delim("rotation1scripts_v4/original_data/recombination.distribution.chap.results/go_enrich/functional_annotation.txt", sep = "\t")
# go.enrich2 = read.delim("rotation1scripts_v4/original_data/recombination.distribution.chap.results/go_enrich/go_enrichment.txt", sep = "\t")
# 
# 
# go.enrich = read.delim("rotation1scripts_v4/original_data/recombination.distribution.chap.results/go_enrich/omicsbox_table.txt", sep = "\t")
# go.enrich2 = go.enrich[-which(is.na(go.enrich$X.GO)), ]
# 
# #grab only seq names that have go IDs
# writeLines(as.character(go.enrich2$SeqName), "rotation1scripts_v4/original_data/recombination.distribution.chap.results/seq.names.fisher.txt")
# 
# go.enricnrow(go.enrich2)


#### COMPARE GENETIC MAP WITHOUT LISS ####

make.mds.plot = function(geno.df1){
    geno.centi.all = lapply(geno.df1, extract.centimorgan.from.genotypes)
    
    
    geno.centi.all = lapply(geno.centi.all, function(x){
        x$liss.keep = F
        x$liss.keep[which(x$marker %in% geno.centi1[[1]]$marker)] = T
        x
    })
    
    count2 = 1
    templabels = axptemps1
    
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
    mds.df$pop = axptemps1
    
    ggplot(mds.df, aes(x = V1, y = V2)) + geom_point() + geom_text_repel(aes(label = pop)) + theme_bw()

    
}

axp.mds = make.mds.plot(initial.data)


recomb.noliss = lapply(initial.data, detect.recombination, conservative = F)

make.recomb.freq.plot = function(recomb.df, pop.names1, y.adjust){
    if(missing(y.adjust)) y.adjust = 45
    
    axp.recomb.freq1 = lapply(recomb.df, function(x){
        x = arrange(x, individual)
        aggregate(x$num.events, list(factor(x$individual, levels = unique(x$individual))), sum)
    })
    
    
    axp.recomb.freq2 = comb.treatments(axp.recomb.freq1)
    
    axp.anova.freq = aov(axp.recomb.freq2$x ~ as.factor(axp.recomb.freq2$treatment))
    axp.anova.freq.tuk = TukeyHSD(axp.anova.freq)
    library(broom)
    axp.aft2 = tidy(axp.anova.freq.tuk)
    axp.aft2.p = round(axp.aft2$adj.p.value[c(1, 3, 4, 6)], digits = 5)
    
    
    mean(axp.recomb.freq2$x)
    sd(axp.recomb.freq2$x)
    
    axp.aft.means = aggregate(axp.recomb.freq2$x, list(factor(axp.recomb.freq2$treatment, levels = unique(axp.recomb.freq2$treatment))), mean)
    colnames(axp.aft.means) = c("treatment", "mean")
    axp.aft.sds = aggregate(axp.recomb.freq2$x, list(factor(axp.recomb.freq2$treatment, levels = unique(axp.recomb.freq2$treatment))), sd)
    colnames(axp.aft.sds) = c("treatment", "sd")
    
    axp.aft3 = bind_cols(axp.aft.means, axp.aft.sds)
    
    axp.aft3 = relabel.treat.w.pop.names(axp.aft3, "treatment", pop.names1, F)
    axp.aft3 = relabel.treat.w.pop.names(axp.aft3, "treatment1", pop.names1, F)
    
    library(ggsignif)
    axp.recomb.freq.plot1 = ggplot(axp.aft3, aes(x = treatment, y = mean)) + geom_bar(stat = 'identity', fill = "#e6e6e6", color = "#000000") + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.5) + 
        xlab("Population") + ylab("Recombination frequency")
    
    axp.recomb.freq.plot1 = add.comparison.bars(axp.recomb.freq.plot1, axp.anova.freq.tuk, y.adjust, pop.names1)
    axp.recomb.freq.plot1
}

axp.recomb.freq.noliss = make.recomb.freq.plot(recomb.noliss, axptemps1)

grid.arrange(axp.recomb.freq.plot1, axp.recomb.freq.noliss, ncol = 2)



