#### PARAGON X WATKINS49 ####
library(readxl)
pxw49.allen = read_excel("rotation1scripts_v4/original_data/luzie_genetic_mapping/Gardiner et al 2019 paragon x watkins axiom maps.xlsx", sheet = 3)



colnames(pxw49.allen)[1:3] = c("marker", "chr", "cM")
pxw49.allen = as.data.frame(pxw49.allen, stringsAsFactors = F)




pxw49.allen2 = split(pxw49.allen, pxw49.allen$chr)

#in this dataset from the Gardiner paper, there are multiple linkage groups listed as the same chromosome for several chromosomes
#this is obviously a mistake on their part. I'm selecting only the biggest linkage group in each case.
pxw49.allen2 = lapply(pxw49.allen2, function(x){
    x[longest_subseq.R(x$cM), ]
})

pxw49.allen = bind_rows(pxw49.allen2)

pxw49.allen[pxw49.allen == "X"] = "H"
pxw49.allen$chr[grep("[0-9][A-Z]1", pxw49.allen$chr)] = multi.str.split(pxw49.allen$chr[grep("[0-9][A-Z]1", pxw49.allen$chr)], "1$", 1)

#### PARAGON X WATKINS94 ####

pxw94.allen = read_excel("rotation1scripts_v4/original_data/luzie_genetic_mapping/Gardiner et al 2019 paragon x watkins axiom maps.xlsx", sheet = 4)
colnames(pxw94.allen)[1:3] = c("marker", "chr", "cM")

pxw94.allen2 = split(pxw94.allen, pxw94.allen$chr)

#in this dataset from the Gardiner paper, there are multiple linkage groups listed as the same chromosome for several chromosomes
#this is obviously a mistake on their part. I'm selecting only the biggest linkage group in each case.
pxw94.allen2 = lapply(pxw94.allen2, function(x){
    x[longest_subseq.R(x$cM), ]
})

pxw94.allen = bind_rows(pxw94.allen2)

pxw94.allen[pxw94.allen == "X"] = "H"
pxw94.allen$chr[grep("[0-9][A-Z]1", pxw94.allen$chr)] = multi.str.split(pxw94.allen$chr[grep("[0-9][A-Z]1", pxw94.allen$chr)], "1$", 1)



#------- analysis of tautology... 

# pxw942 = prepare.allen.map(pxw94.allen)

# pxw.6b = pxw942[, which(pxw942[1, ] == "6B")]
# g = apply(pxw.6b, 1, function(x){
#     which(x == "H")
# })
# 
# g2 = lapply(g, function(x){
#     length(x) = max(sapply(g, length))
#     x
# })
# 
# g3 = bind_rows(g2)
# g3 = g3[, -1]
# g4 = melt(g3)
# g4.y = lapply(1:108, function(x) rep(x, 198))
# g4.y = do.call(c, g4.y)
# 
# g4$y = g4.y
# ggplot(g4, aes(x = value, y = y)) + geom_point()
# 
# boxplot(sapply(g, mean))
# 
# pxw.recomb = detect.recombination(pxw942)
# pxw.recomb2 = pxw.recomb[which(pxw.recomb$chromo == "6B"), ]

# pxwmaps = qtlanalysis.all(pxw942, pxw94.map, "RIL7", geno.prob.step = 0, fgen = 4, skip.qtl.anal = T)

# pxwmaps2 = pxwmaps[[4]][which(pxwmaps[[4]]$chr == "6B"), ]


# --- end of tautology analysis



#### ALL SACHA POP. ANALYSIS V2 ####

s.x.r.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/sxr.allen.csv", stringsAsFactors = F)
a.x.p.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic map AxP.csv", stringsAsFactors = F)
oxs.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/o.x.s.allen.csv", stringsAsFactors = F)
csxp.map = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.cs.x.p.map.csv", stringsAsFactors = F)
axc.map = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.axc.csv")

all.allen.maps = list(s.x.r.allen, a.x.p.allen, oxs.allen, csxp.map, pxw49.allen, pxw94.allen, axc.map)

# all.allen.maps2 = lapply(all.allen.maps, prepare.allen.map, filter.chromosomes = F)
pop.type1 = list(7, 7, 7, 7, "dh", 4, 4)
pop.type2 = list("riself", "riself", "riself", "riself", "riself", "riself", "dh")

# all.allen.analysis.freq = Map(function(x, y){
#     ### Perform recombination frequency QTL analysis ####
#     geno1 = prepare.allen.map(x, F)
#     geno1[geno1 == "H"] = "-"
#     map1 = get.allen.genetic.map(x, geno1)
#     qtl1 = qtlanalysis.all(geno1, map1, "RIL6", geno.prob.step = 0, skip.liss = T, fgen = y)
#     return(qtl1)    
# }, all.allen.maps, pop.type2)


#multi core analysis for wilkins
# all.allen.analysis.freq = mclapply(all.allen.maps, function(x){
#     ### Perform recombination frequency QTL analysis ####
#     geno1 = prepare.allen.map(x, F)
#     geno1[geno1 == "H"] = "-"
#     map1 = get.allen.genetic.map(x, geno1)
#     
#     if(colnames(x[4]) == "AC1"){
#         poptype1 = "dh"
#     } else {
#         poptype1 = "riself"
#     }
#     
#     qtl1 = qtlanalysis.all(geno1, map1, "RIL6", geno.prob.step = 0, skip.liss = T, fgen = poptype1)
#     return(qtl1)    
# }, mc.cores = 8)

load("rotation1scripts_v4/saved.objects/all.allen.analysis.freq")

all.allen.analysis.freq[[7]][[1]]



# all.allen.analysis.dist = Map(function(x, y){
#     #### Perform recombination distribution QTL analysis #### 
#     geno2 = prepare.allen.map(x, T)
#     geno2[geno2 == "H"] = "-"
#     map2 = get.allen.genetic.map(x, geno2)
#     qtl2 = qtlanalysis.all(geno2, map2, "RIL6", geno.prob.step = 0, skip.liss = F, fgen = y)
#     
#     return(qtl2)    
# }, all.allen.maps, pop.type2)

#multi core analysis for wilkins
all.allen.analysis.dist = lapply(all.allen.maps, function(x){
    #### Perform recombination distribution QTL analysis ####
    geno2 = prepare.allen.map(x, T)
    geno2[geno2 == "H"] = "-"
    map2 = get.allen.genetic.map(x, geno2)
    
    if(colnames(x[4]) == "AC1"){
        poptype1 = "dh"
    } else {
        poptype1 = "riself"
    }
    
    qtl2 = qtlanalysis.all(geno2, map2, "RIL6", geno.prob.step = 0, skip.liss = F, fgen = poptype1, multi.core = T)
 
    return(qtl2)
 })

load("~/project.phd.main/rotation1scripts_v4/saved.objects/all.allen.analysis.dist")
all.allen.analysis.dist

ggplot(outmr, aes(x = pos, y = lod)) + geom_line(size = 1) + 
    facet_grid(cols = vars(chr), scales = "free_x", space = "free_x") + 
    xlab("Chromosome") + ylab("LOD") + geom_hline(yintercept = lod.threshold, linetype = "dashed") + 
    geom_hline(yintercept = lod.threshold2, linetype = "dotted") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                panel.spacing = unit(0, "lines"), strip.background = element_rect(),
                panel.grid = element_blank()) +
    ggtitle(phenotypes.to.test[x])



lapply(all.allen.analysis.freq, get.sig.qtl)



lapply(all.allen.analysis.dist, get.sig.qtl)


#distribution QTLs
qtl.d1 = all.allen.analysis.dist[[2]][[1]][[10]]

qtl.d3 = all.allen.analysis.dist[[3]][[1]][[11]]
qtl.d4 = all.allen.analysis.dist[[4]][[1]][[6]]
qtl.d5 = all.allen.analysis.dist[[7]][[1]][[6]]
qtl.d6 = all.allen.analysis.dist[[7]][[1]][[16]]


all.qtl1 = list(qtl.d1, qtl.d3, qtl.d4, qtl.d5, qtl.d6)
all.allen.maps = list(s.x.r.allen, a.x.p.allen, oxs.allen, csxp.map, pxw49.allen, pxw94.allen, axc.map)

qtl.pop = c("A x P F5", "O x S", "CS x P", "A x C DH", "A x C DH")

qtl.dat1 = bind_rows(Map(function(x, y){
    qtl.dat1 = parse.qtl.data(x)
    qtl.dat1$pop = y
    qtl.dat1
}, all.qtl1, qtl.pop))

colnames(qtl.dat1)
colnames(qtl.dat1) = c("Chromosome", "Position (cM)", "LOD", "LOD", "% Phenotypic Variance", "p-value (chi-square)", "p-value F", "AA Mean", "BB Mean", "AA SE",
                                             "BB SE", "ADD EST", "ADD SE", "ADD T", "Chromosome", "Position Previous Marker (cM)", "LOD previous marker", "Physical position of marker before (%)", "Physical position of marker before (bp)", "Chromosome after", "Position of marker after (cM)",
                                             "LOD of next marker", "Physical position of next marker (%)", "Physical position of next marker (bp)", "Phenotype", "Population")

qtl.dat1 = qtl.dat1[, -(c(4, 12, 13, 14, 15, 20))]
qtl.dat1 = qtl.dat1[, c(20, 19, 1:18)]

qtl.plots = Map(function(x, y){
    x[[2]] + ggtitle(paste0(y, ": ", x$phenotype)) + theme_bw(base_size = 22)
}, all.qtl1, qtl.pop)

qtl.plots = qtl.plots[-(4:5)]

do.call(grid.arrange, qtl.plots)

#### UPDATED PLOTTING ####


all.qtl1 = all.qtl1[1:3]
phenotypes1 = c('A x P F5: 6B', 'O x S: 5A', 'CS x P: 2D')

initial.qtl.plots = Map(function(x, y){
ggplot(x[[8]], aes(x = pos, y = lod)) + geom_line(size = 1) + 
    facet_grid(cols = vars(chr), scales = "free_x", space = "free_x") + 
    xlab("Chromosome") + ylab("LOD") + geom_hline(yintercept = x[[9]], linetype = "dashed") + 
    geom_hline(yintercept = x[[10]], linetype = "dotted") +
    theme_bw(base_size = 22) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                panel.spacing = unit(0, "lines"), strip.background = element_rect(),
                panel.grid = element_blank()) +
            ggtitle(y)
}, all.qtl1, phenotypes1)

save(initial.qtl.plots, file = '~/project.phd.main/rotation1scripts_v4/saved.objects/initial.qtl.plots')


##### F2 QTL ANALYSIS #####

#lgeno3comb2 = lgeno3comb

##can't remove all the heterozygotes for the F2 as this removes all of the data from the genetic map (50% hets)
## lgeno3comb2[lgeno3comb2 == "H"] = "-"
#lgeno3comb2.map = extract.centimorgan.from.genotypes(lgeno3comb2)
#lgeno3comb2.map = lgeno3comb2.map[, c(3, 2, 1)]
#colnames(lgeno3comb2.map) = c("marker", "chr", "cM")

#axcf2.qtl.analysis.step0.het = qtlanalysis.all(lgeno3comb2, lgeno3comb2.map, "RIL2", F, "-", geno.prob.step = 0, fgen = 2, qtl.algorithm = "ehk")
#axcf2.qtl.analysis.step0.het.noliss = qtlanalysis.all(lgeno3comb2, lgeno3comb2.map, "RIL2", F, ".", geno.prob.step = 0, fgen = 2, qtl.algorithm = "ehk", skip.liss = T)

#get.sig.qtl(axcf2.qtl.analysis.step0.het)
#get.sig.qtl(axcf2.qtl.analysis.step0.het.noliss)



