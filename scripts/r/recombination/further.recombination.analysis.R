
#### AXP DQC ANALYSIS ####

axp.dqc.info = read_delim("rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.alex.f2/dqc.info.txt", "\t", comment = "#")
lowqc.met3 = read_delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/all.replicates.comb.qc.metric/axc.all.replicates.txt",
                                                delim = "\t") 

axp.dqc1 = axp.dqc.info[grep("Q1", axp.dqc.info$`Sample Filename`), ]
axp.dqc2 = axp.dqc.info[grep("Q2", axp.dqc.info$`Sample Filename`), ]
axp.dqc3 = axp.dqc.info[grep("Q3", axp.dqc.info$`Sample Filename`), ]
axp.dqc4 = axp.dqc.info[grep("Q4", axp.dqc.info$`Sample Filename`), ]

all.axp.dqc.list = list(axp.dqc1, axp.dqc2, axp.dqc3, axp.dqc4)

lapply(all.axp.dqc.list, function(x) mean(x$DQC))

lapply(all.axp.dqc.list, function(x) mean(na.exclude(x$DQC)))
lapply(all.axp.dqc.list, function(x) mean(na.exclude(x$`QC call_rate`)))
lapply(all.axp.dqc.list, function(x) mean(na.exclude(x$call_rate)))
lapply(all.axp.dqc.list, function(x) mean(na.exclude(x$`QC het_rate`)))

#lapply(qc.metrics.files, function(x) mean(x$DQC))
#lapply(qc.metrics.files, function(x) mean(x$`QC call_rate`))
#lapply(qc.metrics.files, function(x) mean(x$call_rate))
#lapply(qc.metrics.files, function(x) mean(x$`QC het_rate`))




perform.dqc.analysis(axc.avg2, "phys.dist.short", lowqc.met3)
perform.dqc.analysis(axc.avg2, "phys.dist.long", lowqc.met3)






#### PREPARE HIGHER QC GENOTYPING DATA FOR AxC / CxA ####

#examine lower qc sample metrics





perform.dqc.analysis(axc.avg2, "phys.dist.long", lowqc.met3)
perform.dqc.analysis(axc.avg2, "phys.dist.short", lowqc.met3)

axc.temp.avg1 = axc.avg2[which(axc.avg2$treatment == 3), ]

# perform.dqc.analysis(axc.temp.avg1, "phys.dist.long", lowqc.met3)
# perform.dqc.analysis(axc.temp.avg1, "phys.dist.short", lowqc.met3)

perform.dqc.analysis(axpavg.r1, "phys.dist.long", axp.dqc.info)
perform.dqc.analysis(axpavg.r1, "phys.dist.short", axp.dqc.info)



# count1 = 1
# dqc.mrd.plot1 = lapply(axc.avg.chromo3, function(x){
#     g = ggplot(x, aes(x = phys.dist, y = DQC)) + geom_point() + theme_bw() + ggtitle(asmaplgs[count1]) + xlab("MRD")
#     count1 <<- count1 + 1
#     g
# })
# 
# count1 = 1
# qc.mrd.plot2 = lapply(axc.avg.chromo3, function(x){
#     g = ggplot(x, aes(x = phys.dist, y = `QC call_rate`)) + geom_point() + theme_bw() + ggtitle(asmaplgs[count1]) + xlab("MRD") + ylab("QC Call Rate")
#     count1 <<- count1 + 1
#     g
# })



# axcc2 = lapply(axc.avg.chromo2, compare.qc.phys.dist, qc.type = 2)
# 
# do.call(grid.arrange, axcc1)
# 
# do.call(grid.arrange, axcc2)



lowqc.met3$mrd = ""

# lowqc.met4 = lowqc.met3[-which(!lowqc.met3$`Sample Filename` %in% mrd.avg$Group.1), ]
# 
# lowqc.met5 = bind_cols(lowqc.met4, mrd.avg[match(lowqc.met4$`Sample Filename`, mrd.avg$Group.1), ])
# 
# 
# hist(lowqc.met5$x)
# 
# lm1 = lm(lowqc.met5$x ~ lowqc.met5$DQC)
# summary(lm1)
# 
# lm2 = lm(lowqc.met5$x ~ lowqc.met5$`QC call_rate`)
# summary(lm2)
# 
# 
# 
# 
# 
# lm2 = lm(lowqc.met5$x ~ lowqc.met5$`QC call_rate`)
# summary(lm2)
# 



#examine sample QC metrics
# 
# qc.metrics.files = read_delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/all.replicates.comb.qc.metric/axc.all.replicates.txt",
#                                                             delim = "\t") 
# 
# qccxa.651 = qc.metrics.files[grep("65-1", qc.metrics.files$`Sample Filename`), ]
# qccxa.653 = qc.metrics.files[grep("65-3", qc.metrics.files$`Sample Filename`), ]
# qcaxc.512 = qc.metrics.files[grep("51-2", qc.metrics.files$`Sample Filename`), ]
# qcaxc.611 = qc.metrics.files[grep("61-1", qc.metrics.files$`Sample Filename`), ]
# 
# qc.metrics.files = list(qccxa.651, qccxa.653, qcaxc.512, qcaxc.611)
# 
# # geno.order = c("c.x.a651", "c.x.a653", "a.x.c512", "a.x.c611")
# 
# lapply(qc.metrics.files, function(x) mean(x$DQC))
# lapply(qc.metrics.files, function(x) mean(x$`QC call_rate`))
# lapply(qc.metrics.files, function(x) mean(x$call_rate))
# lapply(qc.metrics.files, function(x) mean(x$`QC het_rate`))
# 
# 
# 
# 
# result.qc.met3 = sapply(qc.met2, function(x){
#     dqc = c(mean(x$DQC), sd(x$DQC))
#     qc.cr = c(mean(x$`QC call_rate`), sd(x$`QC call_rate`))
#     
#     qc.pass = mean(x[which(x$`Pass/Fail` == "Pass"), ]$`QC call_rate`)
#     
#     data.frame(dqc[1], dqc[2], qc.cr[1], qc.cr[2], qc.pass)
#     
# })
# 
# colnames(result.qc.met3) = geno.order
# rownames(result.qc.met3) = c("dqc.mean", "dqc.sd", "qc.cr.mean", "qc.cr.sd", "avg.qc.for.passing.samples")
# 
# 
# 
# 
# qc.met4 = comb.treatments(qc.met2)
# qc.met4$treatment = as.factor(qc.met4$treatment)
# 
# 
# 
# dqc.hist = ggplot(qc.met4, aes(x = DQC)) + geom_histogram(bins = 100) + facet_grid(vars(treatment)) + theme_bw()
# 
# printplot("dqc.hist2", dqc.hist)
# 
# exportresults(result.qc.met3, "result.qc.met3")
# 
# dqc.anova = aov(`QC call_rate` ~ treatment, data = qc.met4)
# summary(dqc.anova)
# 
# TukeyHSD(dqc.anova)
# 
# 
# par(mfrow = c(4, 1))
# hist(qc.met2[[1]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
# hist(qc.met2[[2]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
# hist(qc.met2[[3]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
# hist(qc.met2[[4]]$DQC, ylim = c(0, 25), xlim = c(0.86, 1), breaks = 50)
# 
# 
# 



#are dqc values significantly different between treatments?
# kruskal.test(sapply(qc.met2, function(x) x$DQC))
# 
# wilcox.test(sapply(qc.met2, function(x) x$DQC)[[1]], sapply(qc.met2, function(x) x$DQC)[[3]])


#examine SNP QC metrics
# 
# snp.qc.metrics.files = paste0("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.probe.metrics/", 
#                                                             list.files("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.probe.metrics//", pattern = ".txt"))
# 
# snp.qc.metrics.files = snp.qc.metrics.files[c(3, 4, 1, 2)]
# 
# snp.met2 = lapply(snp.qc.metrics.files, read.table, sep = "\t", header = T)
# 
# g = read.table("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/qc.probe.metrics/axc512.txt", sep = "\t")
# 
# g[which(g$V1 == "AX-94796264"), ]
# 
# snp.met3 = lapply(snp.met2, function(x){
#     x[which(x$probeset_id %in% highly.conservative.marker.list), ]
# })
# 
# lapply(snp.met3, function(x){
#     mean(na.exclude(x$HomFLD))
# })
# 
# lapply(snp.met3, function(x){
#     mean(na.exclude(x$FLD))
# })
# 
# lapply(snp.met3, function(x){
#     mean(na.exclude(x$HetSO))
# })
# 
# lapply(snp.met3, function(x){
#     mean(na.exclude(x$HomRO))
# })


#make random mixtures of populations
# 
# c.x.a.comb = as.data.frame(rbind(c.x.a651[2:nrow(c.x.a651), ], c.x.a653[2:nrow(c.x.a653), ]))
# 
# vec1 = 1:190
# 
# samp1 = sample(190, 95)
# samp2 = vec1[-match(samp1, vec1)]
# 
# c.x.a.random1 = as.data.frame(rbind(c.x.a651[1, ], c.x.a.comb[samp1, ]))
# c.x.a.random2 = as.data.frame(rbind(c.x.a651[1, ], c.x.a.comb[samp2, ]))
# 
# a.x.c.comb = as.data.frame(rbind(a.x.c512[2:nrow(a.x.c512), ], a.x.c611[2:nrow(a.x.c611), ]))
# 
# vec1 = 1:190
# 
# samp1 = sample(190, 95)
# samp2 = vec1[-match(samp1, vec1)]
# 
# a.x.c.random1 = as.data.frame(rbind(a.x.c512[1, ], a.x.c.comb[samp1, ]))
# a.x.c.random2 = as.data.frame(rbind(a.x.c512[1, ], a.x.c.comb[samp2, ]))
# 
# list.of.geno.dfs.each.pop.random = list(c.x.a.random1, c.x.a.random2, a.x.c.random1, a.x.c.random2)
# 
# #remove individuals with excessive (erroneous) number of recombination events
# list.of.geno.dfs.each.pop.random = lapply(list.of.geno.dfs.each.pop.random, function(x){
#     to.rm1 = which(x$probeset_id %in% indivi.lar1)
#     
#     if(length(to.rm1) > 0){
#         x = x[-to.rm1, ]
#     }
#     
#     x
# })
# 
# 
# #setting all linkage groups to the same to see if any individuals are skewing recombination distribution by a large amount
# list.of.geno.dfs.each.pop.random.all.one.lg = lapply(list.of.geno.dfs.each.pop.random, function(x){
#     x[1, 3:ncol(x)] = "1"
#     x
# })

#### PREPARE COMPARISON OF PARENTS ON AxC ARRAY TO OTHER AxC GENOTYPES ####

#first do comparison for the axc611 array

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'ac14037', password = 'mnb56ghj', db = 'alexcereals', host = '127.0.0.1')
avalon = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Avalon'")
avalon2 = fetch(avalon, n=-1)

cadenza = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Cadenza'")
cadenza2 = fetch(cadenza, n=-1)


avalon651 = read.delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/cxa651.txt", sep = "\t", comment.char = "#")
avalon651.2 = avalon651[, c(1, grep("O23", colnames(avalon651)))]

cadenza651.2 = avalon651[, c(1, grep("M23", colnames(avalon651)))]


avalon611 = read.delim("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/higher.qc/axc611.low.qc.txt", sep = "\t", comment.char = "#")
avalon611.2 = avalon611[, c(1, grep("O23", colnames(avalon611)))]

cadenza611.2 = avalon611[, c(1, grep("M23", colnames(avalon611)))]

# match(avalon2$Probe_row, avalon611.2$probeset_id)

# length(which(avalon2$Matrix_value == avalon611.2$O23_AxC.61.1.H12.CEL_call_code)) / 35143

# length(which(cadenza2$Matrix_value == cadenza611.2$M23_AxC.61.1.G12.CEL_call_code)) / 35143

# length(which(avalon2$Matrix_value == avalon651.2$O23_CxA.65.1.H12.CEL_call_code)) / 35143

# length(which(cadenza2$Matrix_value == cadenza651.2$M23_CxA.65.1.G12.CEL_call_code)) / 35143

t11 = as.data.frame(table(cadenza2$Matrix_value))
t11$var = "cadenza_cereals"
t12 = as.data.frame(table(avalon2$Matrix_value))
t12$var = "avalon_cereals"

t13 = as.data.frame(table(cadenza611.2$M23_AxC.61.1.G12.CEL_call_code))
t13$var = "cadenza_axc_array"
t14 = as.data.frame(table(avalon611.2$O23_AxC.61.1.H12.CEL_call_code))
t14$var = "avalon_axc_array"

t15 = as.data.frame(table(cadenza651.2$M23_CxA.65.1.G12.CEL_call_code))
t15$var = "cadenza_cxa_array"

t16 = as.data.frame(table(avalon651.2$O23_CxA.65.1.H12.CEL_call_code))
t16$var = "avalon_cxa_array"

geno.freq.comp = bind_rows(t11, t12, t13, t14, t15, t16)
colnames(geno.freq.comp) = c("Genotype", "Frequency", "Sample")
geno.freq.comp = dcast(geno.freq.comp, Genotype ~ Sample, value.var = "Frequency")
geno.freq.comp = geno.freq.comp[, c(1, 3, 6, 2, 4, 5, 7)]
colnames(geno.freq.comp) = c("Genotype", "Avalon Cereals", "Cadenza Cereals", "Avalon AxC", "Avalon CxA", "Cadenza AxC", "Cadenza CxA")


exportresults(bind_rows(t11, t12, t13, t14, t15, t16), "genotype.frequency.comparison")


#do comparison only for highly.conservative.marker.list
# 
# avalon651.3 = avalon651.2[which(avalon651.2$probeset_id %in% highly.conservative.marker.list), ]
# cadenza651.3 = cadenza651.2[which(cadenza651.2$probeset_id %in% highly.conservative.marker.list), ]
# 
# avalon3 = avalon2[which(avalon2$Probe_row %in% highly.conservative.marker.list), ]
# avalon611.3 = avalon611.2[which(avalon611.2$probeset_id %in% highly.conservative.marker.list), ]




# 
# qq1 = length(which(avalon3$Matrix_value == avalon611.3$O23_AxC.61.1.H12.CEL_call_code))
# qq2 = length(highly.conservative.marker.list)
# 
# qq1 / qq2
# 
# cadenza3 = cadenza2[which(cadenza2$Probe_row %in% highly.conservative.marker.list), ]
# cadenza611.3 = cadenza611.2[which(cadenza611.2$probeset_id %in% highly.conservative.marker.list), ]
# 
# qq3 = length(which(cadenza3$Matrix_value == cadenza611.3$M23_AxC.61.1.G12.CEL_call_code))
# qq4 = length(highly.conservative.marker.list)
# 
# qq3 / qq4
# 
# 
# mismatch.coords1 = which(cadenza3$Matrix_value != cadenza611.3$M23_AxC.61.1.G12.CEL_call_code)
# 
# cad3.1 = cadenza3[mismatch.coords1, ]
# cad3.2 = cadenza611.3[mismatch.coords1, ]
# cad3.3 = bind_cols(cad3.1, cad3.2)
# t22 = as.data.frame(table(cad3.3$Matrix_value))
# t22$var = "cad_cerealsdb"
# t23 = as.data.frame(table(cad3.3$M23_AxC.61.1.G12.CEL_call_code))
# t23$var = "cad_axc_array"
# 
# mismatch.coords2 = which(avalon3$Matrix_value != avalon611.3$O23_AxC.61.1.H12.CEL_call_code)
# 
# ava3.1 = avalon3[mismatch.coords2, ]
# ava3.2 = avalon611.3[mismatch.coords2, ]
# ava3.3 = bind_cols(ava3.1, ava3.2)
# t24 = as.data.frame(table(ava3.3$Matrix_value))
# t24$var = "ava_cerealsdb"
# t25 = as.data.frame(table(ava3.3$O23_AxC.61.1.H12.CEL_call_code))
# t25$var = "ava_axc_array"
# 
# 
# 
# 
# qq5 = length(which(avalon3$Matrix_value == avalon651.3$O23_CxA.65.1.H12.CEL_call_code))
# 
# qq5 / qq4
# 
# qq6 = length(which(cadenza3$Matrix_value == cadenza651.3$M23_CxA.65.1.G12.CEL_call_code))
# 
# qq6 / qq4
# 
# length(which(cadenza651.3$M23_CxA.65.1.G12.CEL_call_code == cadenza611.3$M23_AxC.61.1.G12.CEL_call_code)) / nrow(avalon651.3)
# 
# length(which(avalon651.3$O23_CxA.65.1.H12.CEL_call_code == avalon611.3$O23_AxC.61.1.H12.CEL_call_code)) / nrow(cadenza651.3)
# 
# 
# 
# mismatch.coords3 = which(avalon3$Matrix_value != avalon651.3$O23_CxA.65.1.H12.CEL_call_code)
# 
# ava3.11 = avalon3[mismatch.coords3, ]
# ava3.12 = avalon651.3[mismatch.coords3, ]
# ava3.13 = bind_cols(ava3.11, ava3.12)
# 
# t26 = as.data.frame(table(ava3.13$O23_CxA.65.1.H12.CEL_call_code))
# t26$var = "ava_cxa_array"
# 
# 
# mismatch.coords4 = which(cadenza3$Matrix_value != cadenza651.3$M23_CxA.65.1.G12.CEL_call_code)
# 
# cad3.11 = cadenza3[mismatch.coords4, ]
# cad3.12 = cadenza651.3[mismatch.coords4, ]
# cad3.13 = bind_cols(cad3.11, cad3.12)
# 
# t28 = as.data.frame(table(cad3.13$M23_CxA.65.1.G12.CEL_call_code))
# t28$var = "cad_cxa_array"
# 
# geno.freq.comp2 = bind_rows(t22, t24, t23, t25, t26, t28)
# colnames(geno.freq.comp2) = c("Genotype", "Frequency", "Sample")
# library(reshape2)
# geno.freq.comp2 = dcast(geno.freq.comp2, Genotype ~ Sample, value.var = "Frequency")
# geno.freq.comp2 = geno.freq.comp2[, c(1, 3, 6, 2, 4, 5, 7)]
# colnames(geno.freq.comp2) = c("Genotype", "Avalon Cereals", "Cadenza Cereals", "Avalon AxC", "Avalon CxA", "Cadenza AxC", "Cadenza CxA")
# 
# exportresults(geno.freq.comp2, "comparison.genotypes.mismatch.markers.highly.cons")
# 
# 


#compare chinese spring samples for reference

cs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Chinese%'")
cs2 = fetch(cs, n=-1)

unique(cs2$Var_col)

cs3 = cs2[which(cs2$Var_col == "Chinese_Spring"), ]

cs4 = cs2[which(cs2$Var_col == "192_Chinese_Spring"), ]


length(which(cs3$Matrix_value == cs4$Matrix_value)) / 35143



##### FUNCTIONAL ANNOTATION ANALYSIS ####


# functional.anno = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", sep = "\t", header = T, fill = T, stringsAsFactors = F)

#### PHYLOGENETIC ANALYSIS ####

library(readxl)
alab1 = read_excel("rotation1scripts_v4/original_data/recombination.gene.analysis/alabdullah.genes1.xlsx")
alab1 = alab1[-1, ]
colnames(alab1) = alab1[1, ]
alab1 = alab1[-1, ]
alab2 = alab1[, which(colnames(alab1) %in% c("Gene ID (v1.0)", "Gene name", "Expected function in wheat"))]
alab3 = alab2[match(unique(alab2$`Gene name`), alab2$`Gene name`), ]

alab.subset = alab2[match(unique(alab2$`Expected function in wheat`), alab2$`Expected function in wheat`), ]


alab.tree = readLines("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/alab.genes.tree.lengths.txt")
#alab.tree2 = multi.str.split(alab.tree, ": ", 2)
alab.tree2 = as.numeric(alab.tree)

random.tree = readLines("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/random.genes.tree.lengths.txt")
#random.tree2 = multi.str.split(random.tree, ": ", 2)
random.tree2 = as.numeric(random.tree)
report_t(t.test(alab.tree2, random.tree2))

par(mfrow = c(1, 1))
boxplot(alab.tree2, random.tree2, names = c("Meiotic genes", "Random gene set"))

# ------ WITHOUT BARLEY

alab.tree.nobar = readLines("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/alab.genes.tree.lengths.nobarley.txt")
#alab.tree.nobar2 = multi.str.split(alab.tree.nobar, ": ", 2)
alab.tree.nobar2 = as.numeric(alab.tree.nobar)

random.tree.nobar = readLines("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/random.genes.no.barley.tree.lengths.txt")
#random.tree.nobar2 = multi.str.split(random.tree.nobar, ": ", 2)
random.tree.nobar2 = as.numeric(random.tree.nobar)

report_t(t.test(alab.tree.nobar2, random.tree.nobar2))
report_mann(wilcox.test(alab.tree.nobar2, random.tree.nobar2))

boxplot(alab.tree.nobar2, random.tree.nobar2, names = c("Meiotic genes", "Random gene set"))

# ------ MORE PHYLOGENETIC ANALYSIS


library(ape)
#library(adephylo)

#read in trees calculated from superalignments
alab.tree = read.tree("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/supermatrices/alab.gene.supermatrix.fa.treefile")
random.tree = read.tree("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/supermatrices/random.gene.supermatrix.fa.treefile")


alab.tree.nobar = read.tree("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/supermatrices/alab.gene.supermatrix.no.barley.fa.treefile")
random.tree.nobar = read.tree("~/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/supermatrices/random.gene.supermatrix.nobarley.fa.treefile")


alab.tree$tip.label = c("Chinese Spring", "Cadenza", "Robigus", "Claire", "Paragon", "Barley")
alab.tree.nobar$tip.label = c("Chinese Spring", "Cadenza", "Robigus", "Claire", "Paragon")
random.tree$tip.label = c("Chinese Spring", "Cadenza", "Paragon", "Robigus", "Claire", "Barley")
random.tree.nobar$tip.label = c("Chinese Spring", "Cadenza", "Paragon", "Robigus", "Claire")

alabdist = as.matrix(cophenetic.phylo(alab.tree))
alabdist[upper.tri(alabdist)] <- 0
alabdist = as.data.frame(alabdist)
alabdist[alabdist == "0"] = ""

randomdist = as.matrix(cophenetic.phylo(random.tree))[c(1, 2, 4, 5, 3, 6), c(1, 2, 4, 5, 3, 6)]
randomdist[upper.tri(randomdist)] <- 0
randomdist = as.data.frame(randomdist)
randomdist[randomdist == "0"] = ""



par(mfrow = c(1, 2))
plot(alab.tree)
add.scale.bar(x = 0.2, y = 2.5, length = 0.05)
title("Meiosis Gene Set")
plot(random.tree)
add.scale.bar(x = 0.2, y = 2.5, length = 0.05)
title("Random Gene Set")

plot(unroot(alab.tree.nobar), type = "unrooted", cex = 0.8)
title("Meiosis Gene Set")
add.scale.bar(length = 0.0005)
plot(unroot(random.tree.nobar), type = "unrooted", cex = 0.8)
title("Random Gene Set")
add.scale.bar(length = 0.0005)

#### CHINESE SPRING RECOMBINATION ####

#chinese spring recombination plots

csxp.map.orig = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.cs.x.p.map.csv", stringsAsFactors = F)
csxp.geno = prepare.allen.map(csxp.map.orig)
csxp.map1 = get.allen.genetic.map(csxp.map.orig, csxp.geno)
csxp.liss.map = qtlanalysis.all(csxp.geno, csxp.map1, geno.prob.step = 0, skip.qtl.anal = T)

cslm1 = csxp.liss.map[[4]]
cslm2 = split(cslm1, cslm1$chr)
cslm2 = lapply(cslm2, function(x){
    x$phys.mb = x$phys.bp / 1000000
    x
})


cslm3 = lapply(cslm2, function(x){
    x$cmdiff = c(diff(x$cM), 0)
    x$physdiff = c(diff(x$phys.mb), 0)
    x$cmmb = x$cmdiff / x$physdiff
    x
})

cslm4 = unsplit(cslm3, cslm1$chr)
#cslm4 = cslm4[-which(cslm4$chr == "1B" | cslm4$chr == "5A" | cslm4$chr == "7B"), ]
#cslm4 = cslm4[which(cslm4$chr == "3B"), ]
centromeres.cs = centromeres[which(centromeres$chr %in% unique(cslm4$chr)), ]

chinese.spring.recombination = ggplot(cslm4, aes(x = phys.mb, y = cmmb)) + geom_line() + facet_wrap(. ~ chr, scales = "free_y", ncol = 1) + 
    geom_vline(aes(xintercept = pos), data = centromeres.cs) + geom_point(shape = 16, size = 1) +
    ylab("Recombination (cM / Mb)") + xlab("Physical position (Mb)")


cslm2b = cslm4[which(cslm4$chr == "2B"), ]

#png("E:/ac14037.backup/Google Drive/University/PhD/presentations/round table talk/figures/cs.recomb.png", 2500, 1500,
        #res = 400)
#ggplot(cslm2b, aes(x = phys.mb, y = cmmb)) + geom_line() +
    #geom_vline(aes(xintercept = 347.85)) + geom_point(shape = 16, size = 1) +
    #ylab("Recombination (cM / Mb)") + xlab("Physical position (Mb)") + ggtitle("Chromosome 2B, Chinese Spring X Paragon")
#dev.off()




plots111 = lapply(cslm3, function(x){
    ggplot(x, aes(x = phys.mb, y = cmmb)) + geom_line()
})
do.call(grid.arrange, plots111)

#### COMPARISON OF ARRAY TO GENOME ####

probes35k = readDNAStringSet("rotation1scripts_v4/original_data/fasta/35kprobes.fa")
probes820k = readDNAStringSet("rotation1scripts_v4/original_data/fasta/820kprobes.fa")

luzie.probes = read.csv("rotation1scripts_v4/original_data/IWGSC/Luzie.probes.JIC_RefSeq1_axiom_order_Nov18.csv", stringsAsFactors = F)

unique(luzie.probes$WGA.Chr)


luzie.probes = luzie.probes[-which(luzie.probes$WGA.Chr == "unknown"), ]
luzie.probes$mb = luzie.probes$WGA.bp / 1000000
colnames(luzie.probes)[2] = "chr"
luzie.probes2 = split(luzie.probes, luzie.probes$WGA.Chr)


probe.dist.hists = Map(function(x, y){
    g = make.ggplot.histogram(x$mb, num.bins = 1000, breaks = seq(1, max(na.omit(x$mb)), (max(na.omit(x$mb)) / 1000)), xlabel = "Physical Position", ylabel = "Frequency", plot.title = "Histogram") + 
        coord_cartesian(ylim = c(0, 50))
    g + geom_vline(xintercept = y) + ggtitle(x$WGA.Chr[[1]]) + theme_bw()
}, luzie.probes2, centromeres$pos)

gene.dist.hists[[1]]
probe.dist.hists[[1]]


all.genes.plot = ggplot(gff1, aes(x = V4)) + geom_histogram(binwidth = 1) + facet_grid(rows = vars(chr)) + geom_vline(aes(xintercept = pos), data = centromeres) + coord_cartesian(ylim = c(0, 60)) +
    xlab("Physical Position (Mb)") + ylab("Gene Frequency")
all.probes.plot = ggplot(luzie.probes, aes(x = mb)) + geom_histogram(binwidth = 1) + facet_grid(rows = vars(chr)) + geom_vline(aes(xintercept = pos), data = centromeres) + coord_cartesian(ylim = c(0, 60)) +
    xlab("Physical Position (Mb)") + ylab("Probe Frequency")

grid.arrange(all.genes.plot, all.probes.plot, ncol = 2)


#### CORRELATION CHROMOSOME LENGTH AND MRD ####


#axc analysis
# chromosomecounts <<- read.table("rotation1scripts_v4/original_data/iwgsc.chromosome.counts.txt", stringsAsFactors = F)
# chromosomecounts2 = chromosomecounts
# chromosomecounts2$V2 = multi.str.split(chromosomecounts2$V2, "chr", 2)
# 
# chromomrd2 = cbind(chromosomecounts2[match(means.per.chromo2$Chromosome, chromosomecounts2$V2), ], means.per.chromo2)
# 
# chromomrd2$V1 = chromomrd2$V1 / 1000000
# 
# cor.test(chromomrd2$V1, chromomrd2$`C x A 65-1`)
# cor.test(chromomrd2$V1, chromomrd2$`C x A 65-3`)
# cor.test(chromomrd2$V1, chromomrd2$`A x C 51-2`)
# cor.test(chromomrd2$V1, chromomrd2$`A x C 61-1`)
# 
# #axp analysis
# chromosomecounts2 = chromosomecounts
# chromosomecounts2$V2 = multi.str.split(chromosomecounts2$V2, "chr", 2)
# 
# chromomrd2 = cbind(chromosomecounts2[match(axpmeans.per.chromo2$Group.2, chromosomecounts2$V2), ], axpmeans.per.chromo2)
# 
# chromomrd2$V1 = chromomrd2$V1 / 1000000
# 
# #(5)
# cor.test(chromomrd2$V1, chromomrd2$`1`)
# cor.test(chromomrd2$V1, chromomrd2$`2`)
# cor.test(chromomrd2$V1, chromomrd2$`3`)
# cor.test(chromomrd2$V1, chromomrd2$`4`)
# 
# 
# 
# 
# 
# plot(chromomrd2$V1, chromomrd2$`1`)
# plot(chromomrd2$V1, chromomrd2$`2`)
# plot(chromomrd2$V1, chromomrd2$`3`)
# plot(chromomrd2$V1, chromomrd2$`4`)
