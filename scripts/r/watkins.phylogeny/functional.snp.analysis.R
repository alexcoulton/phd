setwd("/home/ac14037/project.phd.main/")
load("rotation1scripts_v4/saved.objects/integrated.snps5")
source("rotation1scripts_v4/scripts/r/functions.R")

read.vcf = function(x){
    read.table(x, comment.char = "#", sep = "\t", stringsAsFactors = F)
}

read.vep = function(x){
    g = read.table(x, comment.char = "#", sep = "\t", stringsAsFactors = F) 
    g$chromo = multi.str.split(g$V2, ":", 1)
    g$pos = multi.str.split(g$V2, ":", 2)
    g
}

output.paths1 = list.files("~/temp2/ensembl-vep/output.files", pattern = ".*output.sift.txt$", full.names = T)
input.paths1 = list.files("/local-scratch/alex/Watkins_exome_capture_data/samples", pattern = "Sample", full.names = T)
input.paths1 = paste0(input.paths1, "/all.filtered.coverage.vcf")
input.paths1 = c(input.paths1, "/local-scratch/alex/Watkins_exome_capture_data/tauschii.snps.vcf")

# vep.input = read.vcf("/local-scratch/alex/Watkins_exome_capture_data/Sample_96-116_1190788/all.filtered.coverage.vcf")
# vep.output = read.vep("~/temp2/ensembl-vep/output.files/Sample_96-116_1190788.output.txt")
# 
# 
# Map(function(x, y){
#     vep.input = read.vcf(x)
#     vep.output = read.vep(y)
#     
#     nunique(vep.input$V2) == nunique(vep.output$pos)
# }, input.paths1, output.paths1)


input.vcfs.all = lapply(input.paths1, function(x){
    read.vcf(x)
})

output.vcfs.all = lapply(output.paths1, function(x){
    read.vep(x)
})

unique(output.vcfs.all[[1]]$V7)
table(output.vcfs.all[[1]]$V7)


#update - don't want to limit the missense mutation analysis to the D genome.
output.vcfs.all.list.d = lapply(output.vcfs.all, function(x){
    #x = x[grep("D", x$V2), ]
    x[grep("missense", x$V7), ]
})

library(dplyr)
out.comb1 = bind_rows(output.vcfs.all.list.d)
#Number of genes with at least one missense mutation in a least one of the Watkins varieties
length(unique(out.comb1$V4)) #6942



missense.mutations1 = lapply(output.vcfs.all.list.d, function(x){
    x$V1
})

length(output.vcfs.all.list.d)
nrow(output.vcfs.all.list.d[[2]])

#total number of unique missense mutations among all watkins lines
nunique(unlist(missense.mutations1))

#mean number of mutations in all watkins lines
mean(sapply(missense.mutations1, function(x) length(x))) #1758.89
sd(sapply(missense.mutations1, function(x) length(x))) #240.7

#find the number of missense mutations that are shared by all Watkins varieties
length(Reduce(intersect, missense.mutations1)) #32

setdiff(missense.mutations1[[1]], missense.mutations1[[2]])

Reduce(setdiff, missense.mutations1)

pairwise.combiantions1 = as.data.frame((combn(104, 2)))

ggg = lapply(pairwise.combiantions1, function(x){
    intersect(missense.mutations1[[x[1]]], missense.mutations1[[x[2]]])
})

pairwise.shared.missense.mutations = unname(unlist(lapply(ggg, length)))

mean(pairwise.shared.missense.mutations) #640.75
sd(pairwise.shared.missense.mutations) #119.53

max(pairwise.shared.missense.mutations) #1485
min(pairwise.shared.missense.mutations) #328

output.paths1[pairwise.combiantions1[, which.max(pairwise.shared.missense.mutations)]]

output.paths1[pairwise.combiantions1[, which.min(pairwise.shared.missense.mutations)]]


length(intersect(missense.mutations1[[85]], missense.mutations1[[100]])) #434

#### analyse one output vep file ####

lapply(output.vcfs.all, function(x){
    
})

sift.scores.all = lapply(output.vcfs.all, function(g){
    g = g[grep("D", g$V2), ]
    g2 = g[grep("SIFT", g$V14), ]
    sift.scores = as.numeric(multi.str.split(multi.str.split(g2$V14, "\\(", 2), "\\)", 1))
    sift.scores
})

min(sapply(sift.scores.all, mean))
#examine SIFT scores within varieties
output.paths1[which.min(sapply(sift.scores.all, mean))]
output.paths1[which.max(sapply(sift.scores.all, mean))]
max(sapply(sift.scores.all, mean))
sapply(sift.scores.all, sd)

mean(sapply(sift.scores.all, length))
sd(sapply(sift.scores.all, length))

hist(sift.scores.all[[1]])

library(dplyr)
input.vcfs.all2 = bind_rows(input.vcfs.all)
output.vcfs.all2 = bind_rows(output.vcfs.all)

output.vcfs.all2.d = output.vcfs.all2[grep("D", output.vcfs.all2$chromo), ]

colnames(output.vcfs.all2.d) = c("Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
                                                                 "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
                                                                 "Extra", "chromo", "pos")

#grab only unique mutations
output.vcfs.all2.d = output.vcfs.all2.d[match(unique(output.vcfs.all2.d$Uploaded_variation), output.vcfs.all2.d$Uploaded_variation), ]

q1 = output.vcfs.all2.d$Extra[grep("SIFT", output.vcfs.all2.d$Extra)]
q2 = as.numeric(multi.str.split(multi.str.split(q1, "\\(", 2), "\\)", 1))
mean(q2)
sd(q2)

nunique(output.vcfs.all2.d$Gene)

#### EXAMINE EXOME COVERAGE ####

iwgsc.gff = read.delim("~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3", sep = "\t", stringsAsFactors = F, header = F)

iwgsc.d = iwgsc.gff[grep("D", iwgsc.gff$V1), ]
iwgsc.d = iwgsc.d[which(iwgsc.d$V3 == "exon"), ]
iwgsc.d2 = split(iwgsc.d, iwgsc.d$V1)


all.cov.d = all.coverage20x[grep("D", all.coverage20x$V1), ]
all.cov.d2 = split(all.cov.d, all.cov.d$V1)

library(parallel)
genes.in.coverage = Map(function(x, y){
    n = 50
    nr = nrow(x)
    x2 = split(x, rep(1:ceiling(nr/n), each=n, length.out=nr))
    
    genes2 = mclapply(x2, function(z){
        lapply(z$V2, function(z2){
            genes1 = y[which(y$V4 < z2 & y$V5 > z2), ]$V9
            genes1
        })
            
    }, mc.cores = 50)
    
    print("done")
    genes2
}, all.cov.d2, iwgsc.d2)

genes.table = table(unlist(genes.in.coverage))
genes.table = as.data.frame(genes.table)


iwgsc.d$length = iwgsc.d$V5 - iwgsc.d$V4

iwgsc.d3 = split(iwgsc.d, iwgsc.d$V9)
iwgsc.d4 = bind_rows(lapply(iwgsc.d3, function(x){
    x$total.exon.length = sum(x$length)
    x[1, ]
}))


iwgsc.d4$cov = ""
iwgsc.d4[match(genes.table$Var1, iwgsc.d4$V9), ]$cov = genes.table$Freq
iwgsc.d4$cov = as.numeric(iwgsc.d4$cov)
iwgsc.d4$cov[is.na(iwgsc.d4$cov)] = 0

sum(iwgsc.d4$total.exon.length)
sum(iwgsc.d4$cov) / 1000000

length(unique(unlist(genes.in.coverage)))


#### WHOLE GENOME EXOME COV #### 

iwgsc.all = iwgsc.gff
iwgsc.all = iwgsc.all[which(iwgsc.all$V3 == "exon"), ]
iwgsc.all2 = split(iwgsc.all, iwgsc.all$V1)


all.cov.all = all.coverage20x
all.cov.all2 = split(all.cov.all, all.cov.all$V1)

library(parallel)
genes.in.coverage.all = Map(function(x, y){
    n = 50
    nr = nrow(x)
    x2 = split(x, rep(1:ceiling(nr/n), each=n, length.out=nr))
    
    genes2 = mclapply(x2, function(z){
        lapply(z$V2, function(z2){
            genes1 = y[which(y$V4 < z2 & y$V5 > z2), ]$V9
            genes1
        })
        
    }, mc.cores = 50)
    
    print("done")
    genes2
}, all.cov.all2, iwgsc.all2)

genes.table.all = table(unlist(genes.in.coverage.all))
genes.table.all = as.data.frame(genes.table.all)


iwgsc.all$length = iwgsc.all$V5 - iwgsc.all$V4

iwgsc.all3 = split(iwgsc.all, iwgsc.all$V9)
iwgsc.all4 = bind_rows(lapply(iwgsc.all3, function(x){
    x$total.exon.length = sum(x$length)
    x[1, ]
}))


iwgsc.all4$cov = ""
iwgsc.all4[match(genes.table$Var1, iwgsc.all4$V9), ]$cov = genes.table$Freq
iwgsc.all4$cov = as.numeric(iwgsc.all4$cov)
iwgsc.all4$cov[is.na(iwgsc.all4$cov)] = 0

sum(iwgsc.all4$total.exon.length)
sum(iwgsc.all4$cov) / 1000000

(sum(iwgsc.all4$cov)/ sum(iwgsc.all4$total.exon.length)) * 100

length(unique(unlist(genes.in.coverage)))


#### PARSE FUNCTIONAL SNP INFORMATION #### 
load("rotation1scripts_v4/saved.objects/integrated.snps5")
integrated.snps6 = split(integrated.snps5, integrated.snps5$chr)

snp.anno1 = split(output.vcfs.all2.d, output.vcfs.all2.d$chromo)

snp.anno2 = Map(function(x, y){
    #get functional annotations with different consequences
    #not sure why more than one functional annotation is returned when the SNP is the same genotype in the same position...
    # g = lapply(x$pos, function(z){
    #     y2 = y[which(y$pos == z), ]
    #     y2[match(unique(y2$Consequence), y2$Consequence), ]
    # })
    
    #just going to take the first annotation for each SNP
    y2 = y[na.omit(match(x$pos, y$pos)), ]
    print(paste0("match result: ", nrow(x) == nrow(y2)))
    y2
    
}, integrated.snps6, snp.anno1)

#test
snp.anno3 = bind_rows(snp.anno2)

# all(snp.anno3$pos == integrated.snps5$pos)

barplot(table(snp.anno3$Consequence))
snp.anno.cons = as.data.frame(table(snp.anno3$Consequence))
library(ggplot2)
ggplot(snp.anno.cons, aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 60, hjust = 1))

#### PLOT COVERAGE AND HC GENE ANNOTATION ####

all.coverage20x = read.csv("/local-scratch/alex/Watkins_exome_capture_data/csv/all.coverage.20x.csv", stringsAsFactors = F)

hcgff = read.table("~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", stringsAsFactors = F, header = T)

ac2 = split(all.coverage20x, all.coverage20x$V1)
hcgff2 = split(hcgff, hcgff$V2)


par(mfrow = c(2, 1))
hist(ac2[[1]]$V2, breaks = 1000)
hist(hcgff2[[1]]$V5, breaks = 1000)

library(ggplot2)
formatter1000 = function(x) x / 1000
coverage.distribution = make.ggplot.histogram(ac2[[3]]$V2, breaks = seq(1, max(ac2[[3]]$V2), 100000), ylabel = 'Frequency', xlabel = 'Base Position (Mbp)')
coverage.distribution = coverage.distribution + ggtitle('Distribution of bases with > 20x coverage (Chr 1D)') + scale_x_continuous(labels = formatter1000)
gene.distribution = make.ggplot.histogram(hcgff2[[3]]$V5, breaks = seq(1, max(ac2[[3]]$V2), 100000), ylabel = 'Frequency', xlabel = 'Base Position (Mbp)')
gene.distribution = gene.distribution + ggtitle('Distribution of genes (Chr 1D)') + scale_x_continuous(labels = formatter1000)

library(gridExtra)
svg('~/project.phd.main/rotation1scripts_v4/plots/watkins.phylogeny/coverage.distribution.svg', 8, 8)
grid.arrange(coverage.distribution, gene.distribution)
dev.off()

png('~/project.phd.main/rotation1scripts_v4/plots/watkins.phylogeny/coverage.distribution.png', 1000, 800, res = 200)
grid.arrange(coverage.distribution, gene.distribution)
dev.off()


in.a.gene = lapply(ac2[[1]]$V2, function(x){
    which(hcgff2[[1]]$V5 < x & hcgff2[[1]]$V6 > x)
})

library(parallel)
in.a.gene.all = mclapply(1:22, function(z){
    in.a.gene = lapply(ac2[[z]]$V2, function(x){
        which(hcgff2[[z]]$V5 < x & hcgff2[[z]]$V6 > x)
    })
    length(unique(unlist(in.a.gene)))
}, mc.cores = 22)

#Number of genes at least partially covered by exome capture data
sum(unlist(in.a.gene.all)) #16867

#Number of gene in D genome at least partially covered by exome capture data
sum(unlist(in.a.gene.all[grep('D', unname(unlist(lapply(ac2, function(x) x$V1[[1]]))))]))


#number of positions of 20x coverage inside a gene
length(which(sapply(in.a.gene, length) == 1))
#number of positions of 20x coverage not inside a gene
length(which(sapply(in.a.gene, length) == 0))

# get distance of every 20x coverage SNP from the nearest HC gene

dist.from.gene = sapply(ac2[[1]]$V2, function(x){
    g = abs(x - hcgff2[[1]]$V6)
    g[which.min(abs(x - hcgff2[[1]]$V6))]
})

par(mfrow = c(1, 1))
hist(dist.from.gene, breaks = 1000)

dist.from.gene = Map(function(x, y){
    sapply(x$V2, function(z){
        g = abs(z - y$V6)
        g[which.min(abs(z - y$V6))]
    })
}, ac2, hcgff2)

