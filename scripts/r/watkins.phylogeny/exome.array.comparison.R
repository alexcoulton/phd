source("rotation1scripts_v4/scripts/r/functions.R")

#exome data
load("rotation1scripts_v4/saved.objects/integrated.snps.no.indels")
exome.geno = integrated.snps.no.indels
#preparing array data
library(RMySQL)

# mydb = dbConnect(MySQL(), user = 'ac14037', password = 'mnb56ghj', db = 'alexcereals', host = '127.0.0.1')

rs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Watkins%'") #grab Avalon genotyping data for 35k array
rs2 = fetch(rs, n=-1)

samples1 = list.files("rotation1scripts_v4/original_data/ensembl-vep/output.files/", pattern = "*.output.sift.txt")

samples2 = multi.str.split(samples1, "\\.", 1)
s.codes = multi.str.split(samples2, "_", 3)

s.codes = unique(s.codes)

s.codes.trunc = multi.str.split(s.codes, "1190", 2)

#run on server
# library(parallel)
# geno.data = mclapply(s.codes.trunc, function(x){
#     geno.data = rs2[grep(x, rs2$Var_col), ]
#     print("done")
#     geno.data
# }, mc.cores = 50)

load("rotation1scripts_v4/saved.objects/watkins.array.geno.data")

geno.data = watkins.array.geno.data
geno.data[[25]] = geno.data[[25]][which(geno.data[[25]]$Var_col == "Watkins_1190254"), ]

unlist(lapply(geno.data, function(x) unique(x$Var_col)))





library(dplyr)
geno.data2 = bind_rows(geno.data)

#fix up labels
geno.data2$Var_col = gsub("Watkins_", "Watkins_1190", geno.data2$Var_col)
geno.data2$Var_col = gsub("\\.1", "", geno.data2$Var_col)
geno.data2$Var_col = gsub("\\.2", "", geno.data2$Var_col)
geno.data2$Var_col = gsub("\\.3", "", geno.data2$Var_col)
geno.data2$Var_col[which(geno.data2$Var_col == "Watkins_11901190254")] = "Watkins_1190254"

library(tibble)
library(reshape2)

#### ADD BASE INFORMATION TO MARKERS ####
axiom.probes = read.csv("rotation1scripts_v4/original_data/axiom.probe.assignments/Axiom_WhtBrd-1_Annotation.r3.csv", comment.char = "#", stringsAsFactors = F, header = T)

#roughly 16% heterozygotes in the Watkins array data
# n.h = length(which(array.data == 1))
# n.a = length(which(array.data == 0))
# n.b = length(which(array.data == 2))
# (n.h / (n.h + n.a + n.b))* 100

axiom.probes$Allele.H = ""

r.code = c("A", "G")
m.code = c("A", "C")
w.code = c("A", "T")
s.code = c("C", "G")
y.code = c("C", "T")
k.code = c("G", "T")

all.codes = list(r.code, m.code, w.code, s.code, y.code, k.code)
code.letters = c("R", "M", "W", "S", "Y", "K")
for(i in 1:nrow(axiom.probes)){
    both.alleles = c(axiom.probes$Allele.A[i], axiom.probes$Allele.B[i])
    axiom.probes$Allele.H[i] = code.letters[which(sapply(all.codes, function(x){
        all(x %in% both.alleles)
    }))]
}

geno.data3 = geno.data2
geno.data3$base1 = ""
geno.data3$base2 = ""

#AA
geno.data3$base1[which(geno.data3$Matrix_value == "AA")] = axiom.probes$Allele.A[match(geno.data3$Probe_row[which(geno.data3$Matrix_value == "AA")], axiom.probes$Probe.Set.ID)]
geno.data3$base2[which(geno.data3$Matrix_value == "AA")] = axiom.probes$Allele.A[match(geno.data3$Probe_row[which(geno.data3$Matrix_value == "AA")], axiom.probes$Probe.Set.ID)]

#BB
geno.data3$base1[which(geno.data3$Matrix_value == "BB")] = axiom.probes$Allele.B[match(geno.data3$Probe_row[which(geno.data3$Matrix_value == "BB")], axiom.probes$Probe.Set.ID)]
geno.data3$base2[which(geno.data3$Matrix_value == "BB")] = axiom.probes$Allele.B[match(geno.data3$Probe_row[which(geno.data3$Matrix_value == "BB")], axiom.probes$Probe.Set.ID)]

#AB
geno.data3$base1[which(geno.data3$Matrix_value == "AB")] = axiom.probes$Allele.A[match(geno.data3$Probe_row[which(geno.data3$Matrix_value == "AB")], axiom.probes$Probe.Set.ID)]
geno.data3$base2[which(geno.data3$Matrix_value == "AB")] = axiom.probes$Allele.B[match(geno.data3$Probe_row[which(geno.data3$Matrix_value == "AB")], axiom.probes$Probe.Set.ID)]

geno.data3$base1[which(geno.data3$Matrix_value == "NoCall")] = 10000
geno.data3$base2[which(geno.data3$Matrix_value == "NoCall")] = 10000

#sort exome data 

g = multi.str.split(unique(geno.data3$Var_col), "Watkins_", 2)
exome.geno2 = exome.geno[, c(1, 2, 3, unlist(sapply(g, function(x) grep(x, colnames(exome.geno)))), 108)]

count1 = 1
snp.verify = apply(exome.geno2[, 4:97], 1, function(x){
    g = length(which(x != exome.geno2$ref[count1]))
    count1 <<- count1 + 1
    g
})

exome.geno3 = exome.geno2[which(snp.verify != 0), ]

ncol(exome.geno3)

exome.geno4 = exome.geno3[, -(1:3)]


#### FURTHER ARRAY DATA SORTING #### 

array.data = geno.data3

array.data2 = split(array.data, array.data$Var_col)

array.data3 = lapply(array.data2, function(x){
    g = as.data.frame(t(x))
    g = g[-(1:3), ]
    g
})


array.data1 = dcast(geno.data3, Var_col ~ Probe_row, value.var = "base1")
array.data2 = dcast(geno.data3, Var_col ~ Probe_row, value.var = "base2")

new.array.data = as.data.frame(matrix(nrow = (nrow(array.data1) * 2), ncol = (ncol(array.data1))))

new.array.data[seq(1, nrow(new.array.data), 2), ] = array.data1
new.array.data[seq(2, nrow(new.array.data), 2), ] = array.data2

new.array.data

new.array.data[new.array.data == "A"] = 1
new.array.data[new.array.data == "T"] = 2
new.array.data[new.array.data == "C"] = 3
new.array.data[new.array.data == "G"] = 4

# array.data[array.data == "AA"] = 0
# array.data[array.data == "AB"] = 1
# array.data[array.data == "BB"] = 2
# array.data[array.data == "NoCall"] = 10000


new.array.data2 = cbind(new.array.data[, 1], new.array.data)
new.array.data2 = cbind(new.array.data2[, 1], new.array.data2)

new.array.data2[, 1] = new.array.data2[, 3]
new.array.data2[, 2] = 0
new.array.data2[, 3] = 0

pop1 = as.numeric(as.factor(wilkins.data2$Region))

pop1 = unlist(lapply(pop1, function(x) rep(x, 2)))

new.array.data2[, 2] = pop1
new.array.data2[, 1] = multi.str.split(new.array.data2$`new.array.data2[, 1]`, "Watkins_1190", 2)

new.array.data2 = new.array.data2[, -3]

write.table(new.array.data2, "software/strauto/array.input_diploidv3.txt", sep = " ", row.names = F, quote = F, col.names = F)

#### exome structure input ####
exome.str = read.table("rotation1scripts_v4/original_data/STRUCTURE/exome.input.diploid.str", sep = " ", stringsAsFactors = F)
exome.str$V1 = paste0("Watkins_", multi.str.split(exome.str$V1, "_", 3))
exome.str$V1 = multi.str.split(exome.str$V1, "Watkins_1190", 2)

exome.str2 = exome.str[which(exome.str$V1 %in% new.array.data2[, 1]), ]
exome.str2 = exome.str2[-which(exome.str2$V1 == "224")[3:4], ]

g = match(new.array.data2$`new.array.data2[, 1]`, exome.str2$V1)

g[seq(2, 186, 2)] = g[seq(2, 186, 2)] + 1

exome.str3 = exome.str2[g, ]
exome.str3 = exome.str3[, -3]
exome.str3[, 2] = pop1

write.table(exome.str3, "software/strauto/exome.input.diploidv3.str", sep = " ", row.names = F, quote = F, col.names = F)




#### PREPARE ARRAY PCO PLOT ####

v(array.data)

array.data.calls = dcast(geno.data3, Var_col ~ Probe_row, value.var = "Matrix_value")
s(array.data.calls, "rotation1scripts_v4/saved.objects/array.data.calls", "exome.array.comparison.R")


array.data.calls.noind = array.data.calls[, -1]

nrow(array.data.calls)

combinations1 = as.data.frame(combn(93, 2))

count1 = 1
pairwise.dist = mclapply(combinations1, function(x){
    a = as.character(array.data.calls.noind[x[1], ])
    b = as.character(array.data.calls.noind[x[2], ])
    coord.nocall = c(which(a == "NoCall"), which(b == "NoCall"))

    a = a[-coord.nocall]
    b = b[-coord.nocall]
    
    
    g = length(which(a == b)) / ncol(array.data.calls.noind)
    print(count1)
    count1 <<- count1 + 1
    g
}, mc.cores = 50)

library(dplyr)

p1 = bind_rows(combinations1, pairwise.dist)
p1 = as.data.frame(t(p1))

#add self comparisons
p1$V3 = p1$V3 * 100
p1$V3 = 100 - p1$V3
p1.1 = bind_rows(p1, reset.colnames(data.frame(1:93, 1:93, 0)))

write.csv(p1.1, "rotation1scripts_v4/original_data/STRUCTURE/array.pairwise.comparisons.csv", row.names = F)

library(reshape2)
p2 = dcast(p1.1, V1 ~ V2, value.var = "V3")
array.sim.matrix = p2
s(array.sim.matrix, "rotation1scripts_v4/saved.objects/array.sim.matrix", "exome.array.comparison.R")

load("rotation1scripts_v4/saved.objects/array.sim.matrix")

array.sim.matrix = array.sim.matrix[, -1]

#make matrix symmetrical
array.sim.matrix[lower.tri(array.sim.matrix)] = t(array.sim.matrix)[lower.tri(t(array.sim.matrix))]
array.sim.matrix = as.matrix(array.sim.matrix)


as.df = as.data.frame(array.sim.matrix)

distances1 = sapply(as.df, function(x) sum(x))
distances2 = sort(distances1, index.return = T)$ix
wilkins.data2[distances2[1], ]

a1 = as.dist(array.sim.matrix)

plot(cmdscale(a1, k = 4)[, c(1, 2)])

cmdscale(a1)

array.pco = cmdscale(array.sim.matrix)
plot(array.pco)

array.pco = as.data.frame(array.pco)

#### add geo data ####

library(readxl)
wilkins.data = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 2)
wilkins.data$`Watkins Number` = gsub("Watkins_", "", wilkins.data$`Watkins Number`)
wilkins.data$`Watkins Number` = paste0("Watkins_1190", wilkins.data$`Watkins Number`)
wilkins.data$`Country of Origin` = gsub(" ", "_", wilkins.data$`Country of Origin`)

wilkins.data2 = wilkins.data[match(array.data.calls$Var_col, wilkins.data$`Watkins Number`), ]

array.pco$region = wilkins.data2$Region

ggplot(array.pco, aes(x = V1, y = V2, color = region, fill = region, shape = region)) + geom_point() + 
    scale_shape_manual(values = c(15, 16, 17, 18, 19, 4, 2))


#### PREPARE EXOME PCO ####

exome.pco = read.csv("rotation1scripts_v4/original_data/STRUCTURE/geno.pco.csv", stringsAsFactors = F)
names1 = multi.str.split(array.data.calls$Var_col, "_", 2)

exome.pco2 = exome.pco[, c(1, 2, unlist(sapply(names1, function(x) grep(x, colnames(exome.pco)))))]

exome.pco2 = exome.pco2[, -19]
exome.pco2 = as.data.frame(t(exome.pco2), stringsAsFactors = F)
exome.pco2 = exome.pco2[-(1:2), ]

s(exome.pco2, "rotation1scripts_v4/saved.objects/exome.pco2", "exome.array.comparison.R")

load("~/project.phd.main/rotation1scripts_v4/saved.objects/exome.pco2")
#server section
combinations1 = as.data.frame(combn(93, 2))
library(parallel)
count1 = 1



pairwise.dist = mclapply(combinations1[, ], function(x){
    a = as.character(exome.pco2[x[1], ])
    b = as.character(exome.pco2[x[2], ])
    
    g = length(which(a == b)) / ncol(exome.pco2)
    print(count1)
    count1 <<- count1 + 1
    g
}, mc.cores = 50)

library(dplyr)

p1 = bind_rows(combinations1, pairwise.dist)
p1 = as.data.frame(t(p1))

#add self comparisons
p1$V3 = p1$V3 * 100
p1$V3 = 100 - p1$V3

write.csv(p1.1, "rotation1scripts_v4/original_data/STRUCTURE/exome.pairwise.comparisons.csv", row.names = F)

p1.1 = bind_rows(p1, reset.colnames(data.frame(1:93, 1:93, 0)))

library(reshape2)
p2 = dcast(p1.1, V1 ~ V2, value.var = "V3")
exome.sim.matrix = p2
s(exome.sim.matrix, "rotation1scripts_v4/saved.objects/exome.sim.matrix", "exome.array.comparison.R")

#local code

load("rotation1scripts_v4/saved.objects/exome.sim.matrix")



exome.sim.matrix = exome.sim.matrix[, -1]

#make matrix symmetrical
exome.sim.matrix[lower.tri(exome.sim.matrix)] = t(exome.sim.matrix)[lower.tri(t(exome.sim.matrix))]
exome.sim.matrix = as.matrix(exome.sim.matrix)


as.df = as.data.frame(exome.sim.matrix)

distances1 = sapply(as.df, function(x) sum(x))
distances2 = sort(distances1, index.return = T)$ix
wilkins.data2[distances2[1], ]

a1 = as.dist(exome.sim.matrix)

plot(cmdscale(a1, k = 4)[, c(1, 2)])

cmdscale(a1)

exome.pco = cmdscale(exome.sim.matrix)
plot(exome.pco)

exome.pco = as.data.frame(exome.pco)

#### add geo data ####

library(readxl)
wilkins.data = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 2)
wilkins.data$`Watkins Number` = gsub("Watkins_", "", wilkins.data$`Watkins Number`)
wilkins.data$`Watkins Number` = paste0("Watkins_1190", wilkins.data$`Watkins Number`)
wilkins.data$`Country of Origin` = gsub(" ", "_", wilkins.data$`Country of Origin`)

wilkins.data2 = wilkins.data[match(array.data.calls$Var_col, wilkins.data$`Watkins Number`), ]

exome.pco$region = wilkins.data2$Region
wilkins.data2$short.an = multi.str.split(wilkins.data2$`Watkins Number`, "Watkins_1190", 2)

array.pco$an = wilkins.data2$short.an
exome.pco$an = wilkins.data2$short.an

plot1 = ggplot(array.pco, aes(x = V1, y = V2, color = region, fill = region, shape = region)) + geom_point() + 
    scale_shape_manual(values = c(15, 16, 17, 18, 19, 4, 2)) #+ geom_text_repel(aes(label = an), size = 4)

plot2 = ggplot(exome.pco, aes(x = V1, y = V2, color = region, fill = region, shape = region)) + geom_point() + 
    scale_shape_manual(values = c(15, 16, 17, 18, 19, 4, 2)) #+ geom_text_repel(aes(label = an), size = 4) + scale_y_reverse()

plot2

#main_plot
grid.arrange(plot1, plot2, ncol = 2)

pdf("rotation1scripts_v4/plots/exome.array.comparison/pco.plot1.pdf", 15, 10)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

png("rotation1scripts_v4/plots/exome.array.comparison/pco.plot1.png", 1000, 600, res = 100)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()


arrayp1 = read.csv("rotation1scripts_v4/original_data/STRUCTURE/array.pairwise.comparisons.csv")
exomep1 = read.csv("rotation1scripts_v4/original_data/STRUCTURE/exome.pairwise.comparisons.csv")



cor.test(arrayp1$V3, exomep1$V3)

arraylm = lm(exomep1$V3 ~ arrayp1$V3)

plot(arrayp1$V3, exomep1$V3)
abline(arraylm)

plot(arraylm$residuals)



biggest.deviations = sort(abs(arraylm$residuals), decreasing = T, index.return = T)$ix[1:10]

arrayp1[biggest.deviations, ]

g = c(arrayp1[biggest.deviations[2], ]$V1, arrayp1[biggest.deviations[2], ]$V2)
wilkins.data2[unique(g), ]


alm.df = data.frame(arrayp1$V3, exomep1$V3)
#main_plot
ggplot(alm.df, aes(x = arrayp1.V3, y = exomep1.V3)) + geom_point() + geom_abline(intercept = 2.7043, slope = 0.4788) + theme_bw()

png("rotation1scripts_v4/plots/exome.array.comparison/correlation.scatter.png", 600, 600)
plot(arrayp1$V3, exomep1$V3)
abline(arraylm)
dev.off()


#### PROBE DISTRIBUTION PLOTS #### 


exome.bases = exome.pco[, 1:2]
exome.bases2 = split(exome.bases, exome.bases$chr)

hist(exome.bases2[[2]]$pos, breaks = 1000)

allprobegenomeblast <<- read.csv("rotation1scripts_v4/processed_data/BLAST/820k/820k.full.unique.blast.ordered.csv", stringsAsFactors = F)
probes.1 = colnames(array.data.calls)
probes.1 = probes.1[-1]


probeblast = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/35k.probes.blast")

probeblast2 = grab.best.hits(probeblast)

probeblast2 = split(probeblast, probeblast$qseqid)
probeblast3 = lapply(probeblast2, function(x){
    x[which.min(x$evalue), ]
})

probeblast4 = bind_rows(probeblast3)

probeblast5 = split(probeblast4, probeblast4$sseqid)



par(mfrow = c(2, 1))
par(mar=c(1,1,1,1))
hist(exome.bases2[[1]]$pos, breaks = 1000, xlim = c(0, 750000000))
hist(probeblast5[[1]]$sstart, breaks = 1000, xlim = c(0, 750000000))

#probe distribution plot

probe.comp.plot.a = ggplot(exome.bases2[[1]], aes(x = pos)) + geom_histogram(bins = 1000) + xlab('Position (bp)') + ylab('Frequency') + ggtitle('a') + theme_gray(base_size = 22)
probe.comp.plot.b = ggplot(probeblast5[[1]], aes(x = sstart)) + geom_histogram(bins = 1000) + xlab('Position (bp)') + ylab('Frequency') + ggtitle('b') + theme_gray(base_size = 22)

grid.arrange(probe.comp.plot.a, probe.comp.plot.b)

#main_plot
grid.arrange(probe.comp.plot)

#### evanno plots for determining K #### 

arrayevanno = read.delim("software/strauto/array.results/harvester/evanno.txt", comment.char = "#", stringsAsFactors = F, header = F)

exomeevanno = read.delim("software/strauto/exome.results/harvester/evanno.txt", comment.char = "#", stringsAsFactors = F, header = F)

par(mfrow = c(1, 2))
plot(arrayevanno$V1, arrayevanno$V7)
plot(exomeevanno$V1, exomeevanno$V7)





##### PREPARE ARRAY PHYLOGENY ####


#load array tauschii

adat.phy = array.data

adat.phy2 = adat.phy
adat.phy2$base.comb = ""

# adat.phy2[which(adat.phy2$base1 == "C" & adat.phy2$base2 == "T"), ]$base.comb = "Y"
adat.phy2[which(adat.phy2$base1 == "A" & adat.phy2$base2 == "G"), ]$base.comb = "R"
adat.phy2[which(adat.phy2$base1 == "A" & adat.phy2$base2 == "T"), ]$base.comb = "W"
# adat.phy2[which(adat.phy2$base1 == "G" & adat.phy2$base2 == "C"), ]$base.comb = "S"
adat.phy2[which(adat.phy2$base1 == "T" & adat.phy2$base2 == "G"), ]$base.comb = "K"
# adat.phy2[which(adat.phy2$base1 == "C" & adat.phy2$base2 == "A"), ]$base.comb = "M"

adat.phy2[which(adat.phy2$base2 == "C" & adat.phy2$base1 == "T"), ]$base.comb = "Y"
# adat.phy2[which(adat.phy2$base2 == "A" & adat.phy2$base1 == "G"), ]$base.comb = "R"
# adat.phy2[which(adat.phy2$base2 == "A" & adat.phy2$base1 == "T"), ]$base.comb = "W"
adat.phy2[which(adat.phy2$base2 == "G" & adat.phy2$base1 == "C"), ]$base.comb = "S"
# adat.phy2[which(adat.phy2$base2 == "T" & adat.phy2$base1 == "G"), ]$base.comb = "K"
adat.phy2[which(adat.phy2$base2 == "C" & adat.phy2$base1 == "A"), ]$base.comb = "M"


adat.phy2$base.comb[which(adat.phy2$base.comb == "")] = adat.phy2$base1[which(adat.phy2$base.comb == "")]


array.phylo = adat.phy2
array.phylo$base.comb[which(array.phylo$base.comb == "10000")] = "-"

array.phylo2 = dcast(array.phylo, Var_col ~ Probe_row, value.var = "base.comb")

rownames(array.phylo2) = array.phylo2[, 1]
array.phylo2 = array.phylo2[, -1]
rownames(array.phylo2) = multi.str.split(rownames(array.phylo2), "Watkins_1190", 2)


#prepare exome data

load("rotation1scripts_v4/saved.objects/geno.df")

geno.df = geno.df[, c(1, 2, match(rownames(array.phylo2), colnames(geno.df)))]


codes1 = unique(unlist(lapply(geno.df[, 3:ncol(geno.df)], function(x){
    unique(x)
})))

codes1[which(nchar(codes1) == 3)]

codes2 = data.frame(c("G,T", "T,A", "A,T", "C,A", "G,A", "A,G", "C,G", "C,T", "T,C", "A,C", "G,C", "T,G"), 
                                        c("K", "W", "W", "M", "R", "R", "S", "Y", "Y", "M", "S", "K"), stringsAsFactors = F)


codes2 = reset.colnames(codes2)

geno.df[, 3:ncol(geno.df)] = lapply(geno.df[, 3:ncol(geno.df)], function(x){
    g = which(nchar(x) == 3)
    g2 = which(nchar(x) > 3)
    
    if(length(g) > 0){
        x[g] = codes2$V2[match(x[g], codes2$V1)]
    }
    
    if(length(g2) > 0) x[g2] = "-" #set triallelics to nocalls
    
    x
})

geno.df = geno.df[, -c(1, 2)]
geno.df = as.data.frame(t(geno.df), stringsAsFactors = F)

exome.headers = paste0(">", rownames(geno.df))

exome.seqs = apply(geno.df, 1, function(x){
    paste0(x, collapse = "")
})


interleave <- function(v1,v2)                                                                                                                                                                        
{                                                                                                                                                                                                                            
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}

exome.fa = interleave(exome.headers, exome.seqs)


array.headers = paste0(">", rownames(array.phylo2))

array.seqs = apply(array.phylo2, 1, function(x){
    paste0(x, collapse = "")
})


interleave <- function(v1,v2)                                                                                                                                                                        
{                                                                                                                                                                                                                            
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}

array.fa = interleave(array.headers, array.seqs)

writeLines(exome.fa, "rotation1scripts_v4/original_data/watkins.phylogenies/exome.array.comparison/exome.fa")
writeLines(array.fa, "rotation1scripts_v4/original_data/watkins.phylogenies/exome.array.comparison/array.fa")



#### plot compilation object

p1 = grid.arrange(plot1, plot2)
p2 = ggplot(alm.df, aes(x = arrayp1.V3, y = exomep1.V3)) + geom_point() + geom_abline(intercept = 2.7043, slope = 0.4788) + theme_bw()
p3 = grid.arrange(probe.comp.plot)
p4 = grid.arrange(splot1, splot2)

main_plots_exome_comp = list(
  p1,
  p2,
  p3,
  p4
)

s(main_plots_exome_comp, 'rotation1scripts_v4/saved.objects/main_plots_exome_comp', 'exome.array.comparison.R')


#### UPDATED PLOTTING 25/02/2021 ####
# previously most plots were saved as grobs, which are not able to be converted back to ggplot format
# and therefore cannot be modified. Here I've changed this so that the plots are saved in their original
# ggplot formats

p1 = list(plot1, plot2)

p2 = ggplot(alm.df, aes(x = arrayp1.V3, y = exomep1.V3)) +
    geom_point() +
    geom_abline(intercept = 2.7043, slope = 0.4788) +
    theme_bw(base_size = 22) + 
    xlab('Num. differences (array data)') +
    ylab('Num. differences (exome data)') 

p3 = list(probe.comp.plot.a, probe.comp.plot.b)

p4 = list(splot1, splot2)



main_plots_exome_comp = list(
  p1,
  p2,
  p3,
  p4
)

s(main_plots_exome_comp, 'rotation1scripts_v4/saved.objects/main_plots_exome_comp', 'exome.array.comparison.R')



