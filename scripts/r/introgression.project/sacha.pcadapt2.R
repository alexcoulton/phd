library(cowplot)
library(gridExtra)
setwd("/home/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
genotype.data = read.delim("rotation1scripts_v4/original_data/sacha.introgression.data/genotyping_data.txt", sep = "\t", stringsAsFactors = F)
genotype.data = t(genotype.data)
colnames(genotype.data) = genotype.data[1, ]
genotype.data = genotype.data[-1, ]

genotype.data[genotype.data == "AA"] = 0
genotype.data[genotype.data == "AB"] = 1
genotype.data[genotype.data == "BB"] = 2

class(genotype.data) = "numeric"

axiom.order = read.delim("rotation1scripts_v4/original_data/sacha.introgression.data/TABBED_INPUT_vcf_4-1_wheat_rel_sacha.vcf",
                                                 sep = "\t", comment.char = "#", stringsAsFactors = F, header = F)

axiom.order2 = multi.str.split(axiom.order$V3, ";", 2)

axiom.order3 = data.frame(axiom.order$V1, axiom.order$V2, axiom.order2, stringsAsFactors = F)
colnames(axiom.order3) = c("chr", "pos", "marker")



genotype.data2 = genotype.data[, na.omit(match(axiom.order2, colnames(genotype.data)))]
genotype.data3 = t(genotype.data2)

#### GET POPULATION DATA ####

pop1.vars = readLines("rotation1scripts_v4/original_data/sacha.introgression.data/pop1")
pop2.vars = readLines("rotation1scripts_v4/original_data/sacha.introgression.data/pop2")
pop3.vars = readLines("rotation1scripts_v4/original_data/sacha.introgression.data/pop3")
pop4.vars = readLines("rotation1scripts_v4/original_data/sacha.introgression.data/pop4")
genotype.data3.pop1 = genotype.data3[, which(colnames(genotype.data3) %in% pop1.vars[sample(149, 29)])]
genotype.data3.pop2 = genotype.data3[, which(colnames(genotype.data3) %in% pop2.vars)]
genotype.data3.pop3 = genotype.data3[, which(colnames(genotype.data3) %in% pop3.vars)]
genotype.data3.pop4 = genotype.data3[, which(colnames(genotype.data3) %in% pop4.vars)]



pcad1 = read.pcadapt(genotype.data3.pop4, type = "pcadapt")
pcad2 = pcadapt(pcad1, K = 20)

plot(pcad2, option = "screeplot")
plot(pcad2, option = "scores", i = 3, j = 4)

scree.values = c(5, 4, 3, 5)

pcad1.pop1 = read.pcadapt(genotype.data3.pop1, type = "pcadapt")
pcad1.pop2 = read.pcadapt(genotype.data3.pop2, type = "pcadapt")
pcad1.pop3 = read.pcadapt(genotype.data3.pop3, type = "pcadapt")
pcad1.pop4 = read.pcadapt(genotype.data3.pop4, type = "pcadapt")

pcad.list = list(pcad1.pop1, pcad1.pop2, pcad1.pop3, pcad1.pop4)

output.pcadapt = function(pcdat, scree.val){
    pcad2 = pcadapt(pcdat, K = scree.val)
    
    # plot(pcad2, option = "screeplot")
    # plot(pcad2, option = "scores", i = 3, j = 4)
    
    
    
    
    
    padj = p.adjust(pcad2$pvalues, method = "BH")
    
    # pc123 = get.pc(pcad2, which(padj < 0.0001))
    
    
    
    pvals1 = data.frame(padj)
    colnames(pvals1) = "pvalue"
    pvals1$padjlog = -log(pvals1$pvalue)
    pvals1$x = 1:nrow(pvals1)
    
    pvals1$sig = ""
    pvals1$sig[which(pvals1$pvalue < 0.001)] = T
    pvals1$chromo = axiom.order$V1
    pvals2 = split(pvals1, pvals1$chromo)
    
    pvals3 = lapply(pvals2, function(x){
        x$newx = 1:nrow(x)
        x
    })
    
    library(dplyr)
    
    pvals4 = bind_rows(pvals3)
    pvals4$marker = axiom.order2
    pvals4$position = axiom.order$V2
    
    plot1 = ggplot(pvals4, aes(x = newx, y = padjlog, color = sig)) + geom_point() + theme_bw() + 
        theme(panel.spacing = unit(0, "lines"), panel.background = element_blank(), panel.grid = element_blank(),
                    panel.border = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), 
                    axis.ticks.x = element_blank()) +
        facet_grid(cols = vars(chromo), scales = "free_x", space = "free_x") + ylab("-log(p)") + scale_color_brewer()
    
    pvals5 = split(pvals4, pvals4$chromo)
    count1 = 1
    plots.ind.chr = lapply(pvals5, function(z){
        chr.trunc.name = multi.str.split(names(pvals5)[count1], "chr", 2)
        plot2 = ggplot(z, aes(x = position, y = padjlog, color = sig)) + geom_point() + theme_bw() +
            ylab("-log(p)") + scale_color_brewer() + ggtitle(names(pvals5)[count1]) +
            coord_cartesian(xlim = c(0, 
                                                             chromosomecounts[which(chromosomecounts$chr == chr.trunc.name), ]$length.bp))
        count1 <<- count1 + 1
        plot2
    })
    
    
    list(plot1, pvals4, plots.ind.chr)
}

#### PERFORM PCADAPT ####

all.pcanal = Map(function(x, y){
    output.pcadapt(x, y)
}, pcad.list, scree.values)


pcplots1 = lapply(all.pcanal, function(x) x[[1]])

grid.arrange(pcplots1[[1]], pcplots1[[2]], pcplots1[[3]], pcplots1[[4]], ncol = 1)

sachapcadapt = read.delim("rotation1scripts_v4/original_data/sacha.introgression.data/pop1_pi.sites.pi", sep = "\t")
plot(sachapcadapt$POS, -log(sachapcadapt$PI))

nrow(all.pcanal[[4]][[2]])




# eight20 = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/820kprobes.vs.iwgsc.v1.blast")
# eight20.2 = split(eight20, eight20$qseqid)



# pvals1$chromo = genodat1[, 1]
pdf("rotation1scripts_v4/original_data/sacha.introgression.data/pcadapt.pop2.pdf", 10, 10)
ggplot(pvals1, aes(x = x, y = padjlog, color = sig)) + geom_point() + 
    theme_bw()
dev.off()

#### EXAMINE INTROGRESSIONS ####

introgressions = read.delim("rotation1scripts_v4/original_data/sacha.introgression.data/Introgression_table_280617.tdt", sep = "\t", header = F, stringsAsFactors = F)
colnames(introgressions) = c("elite", "chr", "start.bp", "length.bp", "cutoff", "source")
introgressions$end.bp = introgressions$start.bp + introgressions$length.bp

#filter introgressions file, remove introgressions of length 0
introgressions.large = introgressions[which(introgressions$length.bp > 0), ]

intr1 = split(introgressions.large, introgressions.large$elite)
intr2 = lapply(intr1, function(x){
    split(x, x$chr)
})

#filter introgressions file, remove non-unique introgressions (contained within other introgressions)
unique.intros2 = lapply(intr2, function(x){
    unique.intros = lapply(x, function(y){
        z = y[c("start.bp", "end.bp")]
        ranges1 = IRanges(z$start.bp, z$end.bp)
        ranges.overlap = as.data.frame(findOverlaps(ranges1, type = "within", drop.self = T, drop.redundant = T))
        y[-unique(ranges.overlap$subjectHits), ]
    })
    
    bind_rows(unique.intros)
})

unique.intros3 = bind_rows(unique.intros2)

introgressions.large = unique.intros3

pc.all = all.pcanal[[1]][[2]]
pc.all$sig[which(pc.all$pvalue < 0.05)] = T
pc.all.split = split(pc.all, pc.all$chromo)


pc.sig = all.pcanal[[1]][[2]][which(all.pcanal[[1]][[2]]$pvalue < 0.05), ]


pc.sig.split = split(pc.sig, pc.sig$chromo)

pc.sig.split.cons = lapply(pc.sig.split, function(x){
    #take only groups of significant markers; multiple significant markers in consecutive order
    x$x
    x$x[-which(diff(x$x) != 1)]
    
    cons1 = split(x$x, cumsum(c(1, diff(x$x) > 1)))
    
    cons2 = unlist(cons1[which(sapply(cons1, length) > 6)])
    x[which(x$x %in% cons2), ]
    # cons1[which.max(sapply(cons1, length))]
    
})

pc.sig.split.cons2 = lapply(pc.sig.split, function(x){
    #take only groups of significant markers; multiple significant markers in consecutive order
    x$x
    x$x[-which(diff(x$x) != 1)]
    
    cons1 = split(x$x, cumsum(c(1, diff(x$x) > 1)))
    
    cons2 = unlist(cons1[which(sapply(cons1, length) > 1)])
    x[which(x$x %in% cons2), ]
    cons1[which.max(sapply(cons1, length))]
    
})



pc.nosig = all.pcanal[[1]][[2]][which(all.pcanal[[1]][[2]]$pvalue > 0.05), ]

pc.nosig.split = split(pc.nosig, pc.nosig$chromo)


chromo1 = multi.str.split(names(pc.sig.split), "chr", 2)

introgressions.pop1 = introgressions.large[which(introgressions.large$elite %in% pop1.vars), ]

list.of.vars = unique(introgressions.pop1$elite)
list.of.vars = "ChineseSpring"

#### ANALYSE NUMBER OF PCADAPT SIG VALUES OVERLAPPING W INTROGRESSIONS #### 


find.overlapping.introgressions = function(list.of.vars, pcad.df1){
    pc.sig = bind_rows(pcad.df1)
    g1 = lapply(list.of.vars, function(z){
        #loop through varieties
        introgressions.large.var = introgressions.large[which(introgressions.large$elite == z), ]
        
        # chromo1 = sort(unique(introgressions.large.var$chr))
        chromo1 = listofwheatchromosomes
        
        
        percent.non.overlaps = Map(function(x, chromo1){
            #loop through chromosomes
            if(nrow(x) == 0){
                
            } else {
                x = as.data.frame(transpose(x))
                intro.chr = introgressions.large.var[which(introgressions.large.var$chr == chromo1), ]
                
                overlaps = lapply(x, function(y){
                    #loop through significant markers
                    y2 = as.numeric(as.character(y))
                    pos = y2[8]
                    
                    bigger1 = pos > intro.chr$start.bp
                    smaller1 = pos < intro.chr$end.bp
                    overlappingintro = intro.chr[which(bigger1 == T & smaller1 == T), ]
                    
                    # overlappingintro[which(overlappingintro$elite %in% pop1.vars), ]
                })
                
                overlaps2 = sapply(overlaps, nrow)
                sources1 = lapply(overlaps, function(x){
                    x$source[1]
                })
                
                data.frame(overlaps2, chromo1, unlist(sources1), stringsAsFactors = F)
                # length(which(overlaps2 > 0)) / (length(which(overlaps2 > 0)) + length(which(overlaps2 == 0)))
                # overlaps2
            }
            
        }, pcad.df1[paste0("chr", chromo1)], chromo1)
        
        overlap.all2 = bind_rows(percent.non.overlaps)
        colnames(overlap.all2) = c("overlap", "chr")
        overlap.all2$var = z
        overlap.all2
        # mean(unlist(percent.non.overlaps))
        
    })
    g2 = bind_cols(g1)
    
    g3 = g2[, grep("overlap", colnames(g2))]
    
    g.sources = g2[, grep("V", colnames(g2))]
    g.sources2 = as.data.frame(t(g.sources))
    
    o.allpc.sources3 = lapply(g.sources2, function(x){
        length(na.omit(unique(as.character(x))))
    })
    
    g4 = as.data.frame(t(g3))
    g4[g4 > 1] = 1
    
    par(mfrow = c(2, 1))
    num.pc.overlaps = unlist(lapply(g4, sum)) 
    barplot(unlist(lapply(g4, sum)))
    plot(pc.sig$padjlog)
    
    
    data.frame(num.pc.overlaps, unlist(o.allpc.sources3), pc.sig$padjlog, pc.sig$pvalue)
}

Z10 = find.overlapping.introgressions(list.of.vars, pc.sig.split.cons)

par(mfrow = c(3, 1))
barplot(Z10$num.pc.overlaps)
barplot(Z10$unlist.o.allpc.sources3.)
barplot(Z10$pc.sig.pvalue)

v(Z10)

overlaps.allpc1 = lapply(list.of.vars, function(z){
    #loop through varieties
    introgressions.large.var = introgressions.large[which(introgressions.large$elite == z), ]
    
    # chromo1 = sort(unique(introgressions.large.var$chr))
    chromo1 = listofwheatchromosomes
    
    
    percent.non.overlaps = Map(function(x, chromo1){
        #loop through chromosomes
        x = as.data.frame(transpose(x))
        intro.chr = introgressions.large.var[which(introgressions.large.var$chr == chromo1), ]
        
        overlaps = lapply(x, function(y){
            #loop through significant markers
            y2 = as.numeric(as.character(y))
            pos = y2[8]
            
            bigger1 = pos > intro.chr$start.bp
            smaller1 = pos < intro.chr$end.bp
            overlappingintro = intro.chr[which(bigger1 == T & smaller1 == T), ]
            
            # overlappingintro[which(overlappingintro$elite %in% pop1.vars), ]
        })
        
        overlaps2 = sapply(overlaps, nrow)
        sources1 = lapply(overlaps, function(x){
            x$source[1]
        })
        
        data.frame(overlaps2, chromo1, unlist(sources1), stringsAsFactors = F)
        # length(which(overlaps2 > 0)) / (length(which(overlaps2 > 0)) + length(which(overlaps2 == 0)))
        # overlaps2
    }, pc.sig.split[paste0("chr", chromo1)], chromo1)
    
    overlap.all2 = bind_rows(percent.non.overlaps)
    colnames(overlap.all2) = c("overlap", "chr")
    overlap.all2$var = z
    overlap.all2
    # mean(unlist(percent.non.overlaps))
    
})

unique(unlist(lapply(overlaps.allpc1, nrow)))




unlist(lapply(overlaps.allpc1, nrow))

overlaps.allpc2 = bind_cols(overlaps.allpc1)

overlaps.allpc3 = overlaps.allpc2[, grep("overlap", colnames(overlaps.allpc2))]

overlaps.allpc.sources = overlaps.allpc2[, grep("V", colnames(overlaps.allpc2))]
overlaps.allpc.sources2 = as.data.frame(t(overlaps.allpc.sources))

o.allpc.sources3 = lapply(overlaps.allpc.sources2, function(x){
    length(na.omit(unique(as.character(x))))
})

overlaps.allpc4 = as.data.frame(t(overlaps.allpc3))
overlaps.allpc4[overlaps.allpc4 > 1] = 1

par(mfrow = c(2, 1))
num.pc.overlaps = unlist(lapply(overlaps.allpc4, sum)) 
barplot(unlist(lapply(overlaps.allpc4, sum)))
plot(pc.sig$padjlog)


overlap.data1 = data.frame(num.pc.overlaps, unlist(o.allpc.sources3), pc.sig$padjlog)
colnames(overlap.data1) = c("overlaps", "num.sources", "pval")
overlap.data1$diversity = overlap.data1$overlaps / overlap.data1$num.sources

ggplot(overlap.data1, aes(x = pval, y = diversity)) + geom_point()

cor.test(overlap.data1$diversity, overlap.data1$pval)
summary(lm(overlap.data1$diversity ~ overlap.data1$pval))
summary(lm(overlap.data1$overlaps ~ overlap.data1$pval))

par(mfrow = c(2, 1))
barplot(overlap.data1$diversity)
barplot(overlap.data1$num.sources)
barplot(overlap.data1$pval)

par(mfrow = c(1, 1))
plot(overlap.data1$overlaps, overlap.data1$num.sources)

par(mfrow = c(1, 1))
plot(pc.sig$padjlog, num.pc.overlaps)

cor.test(pc.sig$padjlog, num.pc.overlaps)





overlaps.pc.cons1 = sapply(list.of.vars, function(z){
    introgressions.large.var = introgressions.large[which(introgressions.large$elite == z), ]
    # browser()
    
    chromo1 = sort(unique(introgressions.large.var$chr))
    
    
    
    
    percent.non.overlaps = Map(function(x, chromo1){
        x = as.data.frame(transpose(x))
        intro.chr = introgressions.large.var[which(introgressions.large.var$chr == chromo1), ]
        
        overlaps = lapply(x, function(y){
            y2 = as.numeric(as.character(y))
            pos = y2[8]
            
            bigger1 = pos > intro.chr$start.bp
            smaller1 = pos < intro.chr$end.bp
            overlappingintro = intro.chr[which(bigger1 == T & smaller1 == T), ]
            # overlappingintro[which(overlappingintro$elite %in% pop1.vars), ]
        })
        
        overlaps2 = sapply(overlaps, nrow)
        length(which(overlaps2 > 0)) / (length(which(overlaps2 > 0)) + length(which(overlaps2 == 0)))
        
    }, pc.sig.split.cons[paste0("chr", chromo1)], chromo1)
    
    mean(unlist(percent.non.overlaps))
    
})





chr2 = sort(unique(introgressions.large$chr))




chr2 = sort(unique(intro.var$chr))

#### PLOT INTROGRESSIONS #### 

plots.intro.lj = lapply(chr2, function(x){
    intro.var = introgressions.large[which(introgressions.large$elite == "Little_Joss_TurnerJIC"), ]
    intro1a = intro.var[which(intro.var$chr == x), ]
    # browser()
    # intro1a.2 = intro1a[which(intro1a$elite %in% pop1.vars), ]
    # intro1a.2 = intro1a[which(intro1a$elite %in% "ChineseSpring"), ]
    intro1a$intro.id = 1:nrow(intro1a)
    
    intro1a.3 = melt(intro1a, id.vars = "intro.id")
    intro1a.4 = intro1a.3[which(intro1a.3$variable %in% c("start.bp", "end.bp")), ]
    intro1a.4$value = as.numeric(intro1a.4$value)
    
    end.chrom = as.numeric(chromosomecounts[which(chromosomecounts$chr == x), ]$length.bp)
    
    ggplot(intro1a.4, aes(x = value, y = intro.id)) + geom_point() + geom_line(aes(group = intro.id)) +
        coord_cartesian(xlim = c(0, end.chrom)) + ggtitle(paste0(x, " Little Joss"))
    
    
})

plots.intro.cs = lapply(chr2, function(x){
    intro.var = introgressions.large[which(introgressions.large$elite == "ChineseSpring"), ]
    intro1a = intro.var[which(intro.var$chr == x), ]
    # browser()
    # intro1a.2 = intro1a[which(intro1a$elite %in% pop1.vars), ]
    # intro1a.2 = intro1a[which(intro1a$elite %in% "ChineseSpring"), ]
    intro1a$intro.id = 1:nrow(intro1a)
    
    intro1a.3 = melt(intro1a, id.vars = "intro.id")
    intro1a.4 = intro1a.3[which(intro1a.3$variable %in% c("start.bp", "end.bp")), ]
    intro1a.4$value = as.numeric(intro1a.4$value)
    
    end.chrom = as.numeric(chromosomecounts[which(chromosomecounts$chr == x), ]$length.bp)
    
    ggplot(intro1a.4, aes(x = value, y = intro.id)) + geom_point() + geom_line(aes(group = intro.id)) +
        coord_cartesian(xlim = c(0, end.chrom)) + ggtitle(paste0(x, " Chinese Spring"))
    
    
})

plots.intro.all = lapply(chr2, function(x){
    intro.var = introgressions.large[which(introgressions.large$elite %in% pop1.vars), ]
    intro1a = intro.var[which(intro.var$chr == x), ]
    # browser()
    # intro1a.2 = intro1a[which(intro1a$elite %in% pop1.vars), ]
    # intro1a.2 = intro1a[which(intro1a$elite %in% "ChineseSpring"), ]
    intro1a$intro.id = 1:nrow(intro1a)
    
    intro1a.3 = melt(intro1a, id.vars = "intro.id")
    intro1a.4 = intro1a.3[which(intro1a.3$variable %in% c("start.bp", "end.bp")), ]
    intro1a.4$value = as.numeric(intro1a.4$value)
    
    end.chrom = as.numeric(chromosomecounts[which(chromosomecounts$chr == x), ]$length.bp)
    
    ggplot(intro1a.4, aes(x = value, y = intro.id)) + geom_point() + geom_line(aes(group = intro.id)) +
        coord_cartesian(xlim = c(0, end.chrom)) + ggtitle(paste0(x, " All Varieties"))
    
    
})


plot_grid(plots.intro.lj[[5]], plots.intro.cs[[5]], plots.intro.all[[5]], all.pcanal[[1]][[3]][[20]], ncol = 1, align = "v", axis = "r")

all.pcanal[[4]][[1]]
all.pcanal[[1]][[1]]
all.pcanal[[4]][[3]][[1]]

chromosomecounts = read.delim("rotation1scripts_v4/original_data/chromosomecounts.txt", sep = " ",
                                                            header = F, stringsAsFactors = F)

chromosomecounts = chromosomecounts[, c(1, 3)]
chromosomecounts$V3 = multi.str.split(chromosomecounts$V3, "chr", 2)
colnames(chromosomecounts) = c("length.bp", "chr")



diff(which(all.pcanal[[1]][[2]]$pvalue < 0.05))
paul.csv = read.csv("rotation1scripts_v4/original_data/sacha.introgression.data/log10values_PCAdapt.csv",
                                        stringsAsFactors = F)

par(mfrow = c(1, 1))
plot(paul.csv$Position, paul.csv$log10.pop1.)
plot(paul.csv$Position, paul.csv$log10.pop2.)
plot(paul.csv$Position, paul.csv$log10.pop3.)
plot(paul.csv$Position, paul.csv$log10.pop4.)

#### EXAMINE GENES ####

process.blast.allsnps = function(blast.path, blasthit1){
    #returns the distance from the nearest Axiom SNP to the gene
    #we can then use this as a benchmark to determine how far significant selection signals (from pcadapt) are from the gene of interest
    if(missing(blasthit1)) blasthit1 = 1
    blast1 = read.blast(blast.path)
    chr1 = multi.str.split(as.character(blast1$sseqid[blasthit1]), ":", 1)
    chr2 = multi.str.split(chr1, "chr", 2)
    axiom.order4 = axiom.order3[which(axiom.order3$chr == chr1), ]
    # pos1 = multi.str.split(as.character(blast1$sseqid[blasthit1]), ":", 2)
    # pos2 = as.numeric(multi.str.split(pos1, "-", 1))
    pos2 = blast1$sstart[blasthit1]
    g = data.frame(abs(axiom.order4[which.min(abs(axiom.order4$pos - pos2)), ]$pos - pos2), 
    axiom.order4[which.min(abs(axiom.order4$pos - pos2)), ]$marker, stringsAsFactors = F)
    colnames(g) = c("pos", "marker")
    g
}

process.blast.sigsnps = function(blast.path, blasthit1){
    #returns the distance from the nearest Axiom SNP to the gene
    #we can then use this as a benchmark to determine how far significant selection signals (from pcadapt) are from the gene of interest
    # browser()
    if(missing(blasthit1)) blasthit1 = 1
    blast1 = read.blast(blast.path)
    chr1 = multi.str.split(as.character(blast1$sseqid[blasthit1]), ":", 1)
    chr2 = multi.str.split(chr1, "chr", 2)
    axiom.order4 = axiom.order3[which(axiom.order3$chr == chr1), ]
    # pos1 = multi.str.split(as.character(blast1$sseqid[blasthit1]), ":", 2)
    # pos2 = as.numeric(multi.str.split(pos1, "-", 1))
    pos2 = blast1$sstart[blasthit1]
    #rht gene
    analysis.pc1 = lapply(all.pcanal, function(x){
        # browser()
        x.chr = x[[2]][which(x[[2]]$chromo == chr1), ]
        x.sig = x.chr[which(x.chr$pvalue < 0.05), ]
        min.dist1 = min(na.omit(abs(x.sig$position - pos2)))
        g = data.frame(min.dist1, x.sig[which(na.omit(abs(x.sig$position - pos2)) == min.dist1), ]$marker, stringsAsFactors = F)
        colnames(g) = c("pos", "marker")
        g
        
    })
    
    bind_rows(analysis.pc1)
    
}


tags5 = read.blast("bioinf/blast/probe.vs.genes.blast/results.blast/TaGS5-3A.blast")

output.df1 = function(df1, df2, name1){
    df3 = bind_rows(df1, df2)
    df3$gene = name1
    df3 = as.data.frame(t(df3))
    write.csv(df3, paste0("rotation1scripts_v4/original_data/sacha.introgression.data/",
                                                name1, ".csv"))
}


#rht analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/rht-b1b.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/rht-b1b.blast")


output.df1(d1, d2, "rht")



read.blast("bioinf/blast/nr_gene/results.blast/rht-b1b.blast")

#vrn gene analysis

d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/vrn.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/vrn.blast")

output.df1(d1, d2, "vrn")

#ppd gene analysis 
# https://www.ncbi.nlm.nih.gov/nuccore/KJ147477.1
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/ppd.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/ppd.blast")

output.df1(d1, d2, "ppd")

ppdblast = read.blast("bioinf/blast/nr_gene/results.blast/ppd.blast")

all.pcanal[[2]][[1]]

#wfzp gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/wfzp.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/wfzp.blast")

output.df1(d1, d2, "wfzp")

#nfya gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/nfya.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/nfya.blast")

output.df1(d1, d2, "nfya")

nfya.blast = read.blast("bioinf/blast/nr_gene/results.blast/nfya.blast")
v(nfya.blast)

#nam gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/nam-b1.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/nam-b1.blast")

output.df1(d1, d2, "nam")

namb1.blast = read.blast("bioinf/blast/nr_gene/results.blast/nam-b1.blast")
v(namb1.blast)



#nac gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/nac2.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/nac2.blast")

output.df1(d1, d2, "nac")

nac2.blast = read.blast("bioinf/blast/nr_gene/results.blast/nac2.blast")


#tackx6-d1 gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/TaCKX6-D1.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/TaCKX6-D1.blast")

output.df1(d1, d2, "taxkx6-d1")

tackx6.blast = read.blast("bioinf/blast/nr_gene/results.blast/TaCKX6-D1.blast")

#ckx4 gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/ckx4.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/ckx4.blast")

output.df1(d1, d2, "ckx4")

ckx4.blast = read.blast("bioinf/blast/nr_gene/results.blast/ckx4.blast")

#tgw6 gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/tgw6.blast", 10)
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/tgw6.blast", 10)

output.df1(d1, d2, "tgw6")

tgw6.blast = read.blast("bioinf/blast/nr_gene/results.blast/tgw6.blast")


#gs-d1b gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/gs-d1b.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/gs-d1b.blast")

output.df1(d1, d2, "gs-d1b")

#cyp gene analysis
d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/cyp78a3-c.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/cyp78a3-c.blast")

output.df1(d1, d2, "cyp")




d1 = process.blast.allsnps("bioinf/blast/nr_gene/results.blast/TaGS5-3A.blast")
d2 = process.blast.sigsnps("bioinf/blast/nr_gene/results.blast/TaGS5-3A.blast")



##### MORE ANALYSIS ####


paul.log10 = read.csv("rotation1scripts_v4/original_data/sacha.introgression.data/log10values_PCAdapt.csv")

plot(1:nrow(paul.log10), paul.log10$log10.pop1.)

hist(as.numeric(na.omit(paul.log10$log10.pop1.)))

