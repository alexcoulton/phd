allpos100 = read.csv("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/all.pos.100.anal.csv")


sim.stats1 = read.csv("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/all.analysesv2.csv")
sim.noselec = sim.stats1[1:4, ]

sim.noselec$selection.strength = 0
sim.noselec$selection = F
sim.noselec$selection.position = 0
sim.noselec$nsp = 0
allpos100$selection.strength = 1 / allpos100$selection.strength


allpos1002 = bind_rows(allpos100, sim.noselec)

allpos1002$pop.size = factor(allpos1002$pop.size, levels = c(96, 300, 1000, 10000))



#### plot - comparison of population sizes #### 

plotpos100.1 = ggplot(allpos1002, aes(x = selection.strength, y = ns.fdr.0.05, group = pop.size, color = pop.size, shape = pop.size)) + geom_point(size = 2, alpha = 0.7) + 
    geom_line(size = 0.7) + geom_abline(intercept = 950, size = 1) + xlab("Selection strength") + ylab("Proportion of sim. w/ sig. distorted markers (%)") +
    guides(shape = guide_legend(title = "Pop. Size"), color = guide_legend(title = "Pop. Size")) + theme_classic() + ggtitle("(b)") +
    geom_hline(yintercept = 95)



library(reshape2)

allpos1002.1000pop = allpos1002[which(allpos1002$pop.size == 1000), c(6, 10:14)]
allpos100v2 = melt(allpos1002.1000pop, id.vars = "selection.strength")

#### plot - comparison of statsitical tests ####

plotpos100.2 = ggplot(allpos100v2, aes(x = selection.strength, y = value, group = variable, color = variable, shape = variable)) + scale_x_continuous(labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25"), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) + geom_line(size = 0.7, alpha = 0.8) + theme_classic() + xlab("Selection strength") +
    ylab("Proportion of sim. w/ sig. distorted markers (%)") + scale_color_manual(labels = c("0.05", "0.01", "0.001", "0.05 FDR", "0.05 Bonf."), values = c("#050517", "#FF8C42", "#FF3C38", "#372549", "#A23E48")) + scale_shape_manual(labels = c("0.05", "0.01", "0.001", "0.05 FDR", "0.05 Bonf."), values = c(15, 16, 17, 18, 3)) + coord_cartesian(xlim = c(0, 0.25)) + guides(color = guide_legend(title = "Threshold"), shape = guide_legend(title = "Threshold")) + geom_point(size = 2) + ggtitle("(a)")

grid.arrange(plotpos100.2, plotpos100.1, ncol = 2)

#### Comparison between selection position ####
POPSIZE = c(96, 300, 1000, 10000)
allpos1002.1000pop = allpos1002[which(allpos1002$pop.size %in% POPSIZE), c(2, 6, 10:14)]
allpos100v2 = melt(allpos1002.1000pop, id.vars = c("selection.strength", "pop.size"))

allpos200.data = sim.stats3[c("pop.size", "selec.str", "ns.0.05",
                                                            "ns.0.01", "ns.0.001", "fdr0.05", 
                                                            "ns.bon.0.05"
                                                            )]

allpos200.data2 = allpos200.data[which(allpos200.data$pop.size %in% POPSIZE), ]
p.data2 = melt(allpos200.data2, id.vars = c("selec.str", "pop.size"))


# allpos100v2$value = (allpos100v2$value / 1000) * 100
colnames(allpos100v2) = c("selec.str", "pop.size", "variable", "value")
allpos100v2$variable = p.data2$variable

allpos100v2$s.pos = 100
p.data2$s.pos = 200

allpos.comb1 = bind_rows(allpos100v2, p.data2)

allpos.comb2 = allpos.comb1#[which(allpos.comb1$var == "fdr0.05"), ]
# allpos.comb2 = allpos.comb1[which(allpos.comb1$var == "0.05 Bonferroni"), ]
# allpos.comb2 = allpos.comb1[which(allpos.comb1$var == "0.05"), ]
# allpos.comb2 = allpos.comb1[which(allpos.comb1$var == "0.01"), ]

allpos.comb2$s.pos = factor(allpos.comb2$s.pos, levels = unique(allpos.comb2$s.pos))

unique(allpos.comb2)

allpos.comb2$value = (allpos.comb2$value / 1000) * 100

#which(allpos.comb2$variable)
levels(allpos.comb2$variable) = c("0.05", "0.01", "0.001", "FDR 0.05", "Bonferroni 0.05")
selectionpos.plot = ggplot(allpos.comb2, aes(x = selec.str, y = value, group = s.pos, color = s.pos)) + geom_point() + geom_line(aes(group = s.pos)) + 
    facet_grid(rows = vars(pop.size), cols = vars(variable)) +
    xlab("Selection strength") + ylab("Simulations containing sig. dist. markers (%)") +
    labs(color = "S. Pos.")


ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig.s2.pdf", selectionpos.plot, width = 22, height = 17.4, units = "cm", device = cairo_pdf)


#### cm / mb plot ####
#i subsequently decided that calculating cm / mb makes no sense here

genodata1 = read.csv("~/project.phd.main/rotation1scripts_v4/original_data/S6 Coulton_et_al_Genotyping_Data.csv", stringsAsFactors = F)

genodata1 = genodata1[-(1:5), ]
colnames(genodata1) = genodata1[1, ]
genodata1 = genodata1[-1, ]


genodata2 = genodata1[which(genodata1$Chromosome == "1A"), ]
plot(genodata2$`Genetic Position (cM)`)
genodata2$`Genetic Position (cM)` = as.numeric(genodata2$`Genetic Position (cM)`)
recombinationplot1 = ggplot(genodata2, aes(x = 1:nrow(genodata2), y = `Genetic Position (cM)`)) + geom_line() +
    geom_point() + xlab("Marker position")

ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig.s1.pdf", recombinationplot1, width = 10, height = 10, units = "cm", device = cairo_pdf)

genodata3 = genodata1[which(genodata1$Chromosome == "6B"), ]
plot(genodata3$`Genetic Position (cM)`)
genodata3$`Genetic Position (cM)` = as.numeric(genodata3$`Genetic Position (cM)`)
recombinationplot3 = ggplot(genodata3, aes(x = 1:nrow(genodata3), y = `Genetic Position (cM)`)) + geom_line() +
    geom_point() + xlab("Marker position")

ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig.s3.pdf", recombinationplot3, width = 10, height = 10, units = "cm", device = cairo_pdf)



cm.diff1 = as.data.frame(diff(as.numeric(genodata2$`Genetic Position (cM)`)))
colnames(cm.diff1) = "val"
ggplot(cm.diff1, aes(x = 1:nrow(cm.diff1), y = val)) + geom_line()




plot(diff(as.numeric(genodata2$`Genetic Position (cM)`)) / diff(luzie.axiom.1a), type = "b")

genodata2[194:208, 1:5]

luzie.axiom = read.csv("rotation1scripts_v4/original_data/IWGSC/Luzie.probes.JIC_RefSeq1_axiom_order_Nov18.csv", stringsAsFactors = F)

luzie.axiom.1a = luzie.axiom[match(switch.affy.format(genodata2$Marker), luzie.axiom$SNP.id), ]
luzie1abp = luzie.axiom.1a$WGA.bp

luzie1asub = luzie.axiom.1a[match(na.omit(luzie1abp)[longest_subseq.R(na.omit(luzie1abp))], luzie.axiom.1a$WGA.bp), ]

c.x.af2.1asub = c.x.af2[, c(1, 2, match(luzie1asub$SNP.id, switch.affy.format(colnames(c.x.af2))))]

source("rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R")
library(data.table)
cm1asub = extract.centimorgan.from.genotypes(c.x.af2.1asub)

lu1aper = (luzie1asub$WGA.bp / 604003765) * 100

plot(luzie1asub$WGA.bp, cm1asub$cm)

g = (c(diff(cm1asub$cm), 0) / c(diff(lu1aper), 0))

g.coord = which(g > 400)
g2 = g[-g.coord]
lu1aper2 = lu1aper[-g.coord]

plot(lu1aper2, g2)

plot(g2)





#### SDR ANALYSIS ####

sdr.analysis1 = read.csv("~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/sdr_selecpos_1.csv")
sdr.analysis1$pos = as.character(sdr.analysis1$pos)
colnames(sdr.analysis1)[1:2] = c("num.sig", "selec.strength")
sdr.analysis1$selec.strength = 1 / sdr.analysis1$selec.strength
sdr.analysis2 = read.csv("~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/sdr_selecpos_2.csv")

sdr.analysis2$pos = as.character(sdr.analysis2$pos)
sdr.analysis2$selec.strength = 1 / sdr.analysis2$selec.strength

#total number of SDR among all 1000 simulations
sdrplot1 = ggplot(sdr.analysis1, aes(x = selec.strength, y = num.sig, group = pos, color = pos, shape = pos)) + geom_point() + facet_grid(cols = vars(test)) +
    ylab("Total number of SDR") + xlab("Selection strength") +
    labs(color = "Position", shape = "Position") + geom_line(aes(group = pos)) +
    ggtitle("(a)") + theme_gray(base_size = 22)

#number of simulations with at least one SDR
sdrplot2 = ggplot(sdr.analysis2, aes(x = selec.strength, y = num.sig, group = pos, color = pos, shape = pos)) + geom_point() + facet_grid(cols = vars(test)) +
    ylab("Num. sim. with ≥ 1 SDR") + xlab("Selection strength") +
    labs(color = "Position", shape = "Position") + geom_line(aes(group = pos)) +
    ggtitle("(b)") + theme_gray(base_size = 22)



sdrplot.all = grid.arrange(sdrplot1, sdrplot2)

save(sdrplot.all, file = '~/project.phd.main/rotation1scripts_v4/saved.objects/sdrplot.all')



ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig.s4.pdf", sdrplot.all, width = 15, height = 15, units = "cm", device = cairo_pdf)



