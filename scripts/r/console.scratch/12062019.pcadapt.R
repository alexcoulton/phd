library(pcadapt)
library(dplyr)
library(gridExtra)
library(reshape2)
library(readr)
setwd("/home/ac14037/project.phd.main")
# setwd("E:/phd.project.main/")
cereals = read_csv("rotation1scripts_v4/original_data/genotype.data/cerealsdb.35k.genotype.data.csv", col_names = F)

cereals2 = dcast(cereals, X1 ~ X2)
source("rotation1scripts_v4/scripts/r/functions.R")

#get nucleotide bases for each SNP
mysql1 = setup.mysql.connection(T)
rs = dbSendQuery(mysql1, 'SELECT * FROM axiom35;')
rs2 = fetch(rs, n=-1)



luzie.order = read_csv("rotation1scripts_v4/original_data/IWGSC/Luzie.probes.JIC_RefSeq1_axiom_order_Nov18.csv")


cereals3 = cereals2[, c(1, na.omit(match(luzie.order$SNP.id, colnames(cereals2))))]

luzie.order2 = luzie.order[na.omit(match(colnames(cereals3), luzie.order$SNP.id)), ]


cereals3[cereals3 == "AA"] = "0/0"
cereals3[cereals3 == "AB"] = "0/1"
cereals3[cereals3 == "BB"] = "1/1"
cereals3[cereals3 == "NoCall"] = "./."

cereals3 = as.data.frame(t(cereals3))

cereals3 = as.data.frame(cbind(cereals3[, 1], cereals3), stringsAsFactors = F)
cereals3[, 1] = rownames(cereals3)

cereals3 = as.data.frame(cbind(cereals3[, 1], cereals3), stringsAsFactors = F)

# match(luzie.order$SNP.id, cereals3[, 1])

cereals3[, 1] = luzie.order[match(cereals3[, 1], luzie.order$SNP.id), 3]


add.row1 = function(x){
	x = as.data.frame(cbind(x[, 1], x), stringsAsFactors = F)
	x[, 1] = ""
	x
}

cereals3 = add.row1(cereals3)

cereals3[, 1] = luzie.order[match(rownames(cereals3), luzie.order$SNP.id), 2]
cereals3 = add.row1(cereals3)
cereals3 = add.row1(cereals3)
cereals3 = add.row1(cereals3)
cereals3 = add.row1(cereals3)
cereals3 = add.row1(cereals3)
cereals3 = add.row1(cereals3)





rs2 = rs2[match(rownames(cereals3), rs2$affycode35k), ]


#prepare nucleotide bases for vcf file
bases1 = regmatches(rs2$Sequence, regexpr("\\[[a-zA-Z\\/]*\\]", rs2$Sequence))
bases2 = substring(bases1, 2, 4)
bases3 = lapply(bases2, function(x){
	strsplit(x, "/")
	})


cereals3[, 1] = c("", sapply(bases3, function(x) x[[1]][1]))
cereals3[, 2] = c("", sapply(bases3, function(x) x[[1]][2]))
cereals3[, 3] = "."
cereals3[, 4] = "."
cereals3[, 5] = "."
cereals3[, 6] = "GT"

cereals4 = as.data.frame(cbind(cereals3[, 7:9], cereals3[, 1:6], cereals3[, 10:ncol(cereals3)]), stringsAsFactors = F)
cereals4[1, 1:9] = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
cereals4[] = lapply(cereals4, function(x) as.character(x))
colnames(cereals4) = cereals4[1, ]
cereals4 = cereals4[-1, ]

write.vcf = function(df1, name1){
    write.table(df1, paste0("rotation1scripts_v4/original_data/cerealsdb.data/", name1, ".vcf"), sep = "\t", quote = F, row.names = F)
    system(paste0("sed -i '1s/^/#/' rotation1scripts_v4/original_data/cerealsdb.data/", name1, ".vcf"))
    system(paste0("sed -i '1i ##fileformat=VCFv4' rotation1scripts_v4/original_data/cerealsdb.data/", name1, ".vcf"))    
    name1
}

write.table(cereals4, "rotation1scripts_v4/original_data/cerealsdb.data/35k.genotypes.vcf", sep = "\t", quote = F, row.names = F)
system("sed -i '1s/^/#/' rotation1scripts_v4/original_data/cerealsdb.data/35k.genotypes.vcf")
system("sed -i '1i ##fileformat=VCFv4' rotation1scripts_v4/original_data/cerealsdb.data/35k.genotypes.vcf")

#prepare files for analysis with iqtree

cereals5 = cereals4
cereals5[, 10:ncol(cereals5)] = lapply(cereals5[, 10:ncol(cereals5)], function(x){
    coord1 = which(x == "0/0")
    coord2 = which(x == "1/1")
    coord4 = 1:length(x)
    coord4 = coord4[-c(coord1, coord2)]
    
    
    
    x[coord1] = cereals5$REF[coord1]
    x[coord2] = cereals5$ALT[coord2]
    x[coord4] = "-"
    x
    
    
})

sequences = lapply(cereals5[, 10:ncol(cereals5)], function(x) paste0(x, collapse = ""))
watkins.coords = grep("Watkins", names(sequences))
sequences2 = sequences[watkins.coords]
sequences3 = do.call(c, sequences2)




interleave <- function(v1,v2)                                                                                                                                                                        
{                                                                                                                                                                                                                            
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}


sequences4 = interleave(gsub("'", "", paste0(">", names(sequences2))), sequences3)
# writeLines(sequences4, "/local-scratch/alex/Watkins_exome_capture_data/watkins.array.seq.fasta")

library(stringr)
watkins.tags = str_extract_all(names(sequences4), "_.*")
watkins.tags = gsub("_", "", watkins.tags)
watkins.tags = gsub("\\.1", "", watkins.tags) 
watkins.tags = gsub("\\.2", "", watkins.tags) 
watkins.tags = gsub("\\.3", "", watkins.tags) 
watkins.tags = gsub("\\.4", "", watkins.tags)
watkins.tags = gsub("'", "", watkins.tags)
watkins.tags[which(nchar(watkins.tags) == 1)] = paste0("119000", watkins.tags[which(nchar(watkins.tags) == 1)])
watkins.tags[which(nchar(watkins.tags) == 2)] = paste0("11900", watkins.tags[which(nchar(watkins.tags) == 2)])
watkins.tags[which(nchar(watkins.tags) == 3)] = paste0("1190", watkins.tags[which(nchar(watkins.tags) == 3)])
watkins.tags = watkins.tags[-which(watkins.tags == "character(0)")]

exome.watkins = list.files("/local-scratch/alex/Watkins_exome_capture_data", pattern = "Sample")
exome.watkins2 = str_extract(exome.watkins, "[0-9]*$")



sequences3.2 = sequences3[which(watkins.tags %in% exome.watkins2)]
watkins.tags2 = watkins.tags[which(watkins.tags %in% exome.watkins2)]

sequences3.2 = sequences3.2[-which(duplicated(watkins.tags2))]

watkins.tags2 = watkins.tags2[-which(duplicated(watkins.tags2))]

sequences5 = interleave(paste0(">", watkins.tags2), sequences3.2)
writeLines(sequences5, "/local-scratch/alex/Watkins_exome_capture_data/watkins.array.exome.subset.fa")


#### EXAMINE GEOGRAPHY ####

accessions = colnames(cereals4)
accessions = gsub("'", "", accessions)

library(readxl)
geography.df = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 3)

#### CODE FROM historical.recombination.mapping.R ##### 

library(readxl)
#geographical data for elite lines
geography.df = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 3)
#geo data for watkins lines
geo.df.watkins = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 2)

# s(first.cast, "rotation1scripts_v4/saved.objects/first.cast", "historical.recombination.mapping.R")
#genotyping data
load("rotation1scripts_v4/saved.objects/first.cast")
#accession names in genotyping data
g = multi.str.split(first.cast$V2, "'", 2)

coord.starting.w.digit = grep("^[[:digit:]]+\\_", g)

#only some of the accessions have names starting with digits, e.g. 93_robigus, I want to remove these digits. 
#other lines, such as the watkins lines, have digits that are useful and needed
g[coord.starting.w.digit] = multi.str.split(g[coord.starting.w.digit], "_", 2)

lines = data.frame(g)
watkins = data.frame(g[grep("Watkins", lines$g)])
row.names(watkins) = grep("Watkins", lines$g)
watkins$line = multi.str.split(as.character(watkins$g.grep..Watkins...lines.g..), "\\.", 1)

watkins$region = geo.df.watkins$Region[match(watkins$line, geo.df.watkins$`Watkins Number`)]

matches.coords = lapply(lines$g, function(x) grep(x, geography.df$Accession))
first.matches = unlist(lapply(matches.coords, function(x) x[1]))

lines$region = geography.df$Region[first.matches]
lines = lines[-grep("Watkins", lines$g), ] #remove watkins lines from elite lines

test.df.df = read_csv("rotation1scripts_v4/processed_data/genotypes/historical.recombination.mapping/1A.ske.incl.csv")
test.df.df = reset.colnames(test.df.df)
test.df.df = dcast(test.df.df, V2~V3)

table(lines$region)
elite.aus.row.coords = as.character(lines$g[which(lines$region == "Australia")])
elite.aus.row.coords2 = na.omit(match(elite.aus.row.coords, accessions))

elite.europe.row.coords = as.character(lines$g[which(lines$region == "Europe")])
elite.europe.row.coords2 = na.omit(match(elite.europe.row.coords, accessions))


table(watkins$region)
watkins.asia.row.coords = as.character(watkins$g.grep..Watkins...lines.g..[which(watkins$region == "Asia")])
watkins.asia.row.coords2 = na.omit(match(watkins.asia.row.coords, accessions))

watkins.eur.west.row.coords = as.character(watkins$g.grep..Watkins...lines.g..[which(watkins$region == "Europe (West)")])
watkins.eur.west.row.coords2 = na.omit(match(watkins.eur.west.row.coords, accessions))

watkins.eur.east.row.coords = as.character(watkins$g.grep..Watkins...lines.g..[which(watkins$region == "Europe (East)")])
watkins.eur.east.row.coords2 = na.omit(match(watkins.eur.east.row.coords, accessions))

watkins.mid.east.row.coords = as.character(watkins$g.grep..Watkins...lines.g..[which(watkins$region == "Middle East")])
watkins.mid.east.row.coords2 = na.omit(match(watkins.mid.east.row.coords, accessions))

all.vars1 = 10:ncol(cereals4)


all.row.coords = list(all.vars1, elite.aus.row.coords2, elite.europe.row.coords2,
                                                     watkins.asia.row.coords2, watkins.eur.east.row.coords2,
                                                     watkins.eur.west.row.coords2, watkins.mid.east.row.coords2)

cereals5 = lapply(all.row.coords, function(x){
    cereals4[, c(1:9, x)]
})


#### investigate missing snp ####

snps.to.rm = unlist(sapply(cereals5, function(z){
    coord1 = which(apply(z[, 10:ncol(z)], 1, function(x) length(unique(x))) == 1)
    
    which(z[coord1, 10] == "./.")
    coord1[which(z[coord1, 10] == "./.")]
}))

cereals6 = lapply(cereals5, function(x){
    x[-snps.to.rm, ]
})





count1 = 1
popnames1 = c("all.vars", "elite.aus", "elite.eur", "watkins.asia", "watkins.eur.east", "watkins.eur.west", "watkins.mid.east")
kvalues1 = c(9, 5, 8, 8, 9, 10, 7)
lapply(cereals6, function(x){
    write.vcf(x, popnames1[count1])
    count1 <<- count1 + 1
})

#### RUN PCADAPT ####

library(pcadapt)

run.pcadapt.anal = function(name1, kvalue1){
    pcfile1 = read.pcadapt(paste0("rotation1scripts_v4/original_data/cerealsdb.data/", name1, ".vcf"), type = "vcf")
    genodat1 = read.table(paste0("rotation1scripts_v4/original_data/cerealsdb.data/", name1, ".vcf"), sep = "\t", stringsAsFactors = F)
    pc2 = pcadapt(input = pcfile1, K = 20)
    
    dir.create(paste0("rotation1scripts_v4/plots/pcadapt/", name1))
    
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/screeplot.pdf"), 10, 10)
    plot(pc2, option = "screeplot")
    dev.off()
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/scoresplot.pdf"), 10, 10)
    plot(pc2, option = "scores")
    dev.off()
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/scoresplot2.pdf"), 10, 10)
    plot(pc2, option = "scores", i = 3, j = 4)
    dev.off()
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/scoresplot3.pdf"), 10, 10)
    plot(pc2, option = "scores", i = 5, j = 6)
    dev.off()
    
    pc3 = pcadapt(input = pcfile1, K = kvalue1)
    
    
    padj = p.adjust(pc3$pvalues, method = "BH")
    outliers = which(padj < 0.1)
    length(outliers)
    
    res <- pcadapt(pcfile1, K = 4)
    
    pvals1 = data.frame(padj)
    colnames(pvals1) = "pvalue"
    pvals1$padjlog = -log(pvals1$pvalue)
    pvals1$x = 1:nrow(pvals1)
    pvals1$sig = ""
    pvals1$sig[which(pvals1$pvalue < 0.001)] = T
    
    pvals1$chromo = genodat1[, 1]
    plot1 = ggplot(pvals1, aes(x = x, y = padjlog, color = sig)) + geom_point() + 
        facet_grid(cols = vars(chromo), scales = "free_x") + theme_bw() +
        theme(panel.spacing.x = unit(0, "lines"), axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(), panel.grid = element_blank()) + ggtitle(name1)
    
    #### OTHER PLOTTING ####
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/manhattanplot.pdf"), 10, 10)
    plot(pc3, option = "manhattan")
    dev.off()
    
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/qqplot.pdf"), 10, 10)
    plot(pc3, option = "qqplot")
    dev.off()
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/histplot.pdf"), 10, 10)
    hist(pc3$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
    dev.off()
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/histplot2.pdf"), 10, 10)
    plot(pc3, option = "stat.distribution")
    dev.off()
    
 
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/ldcheck.pdf"), 10, 10)
    par(mfrow = c(2, 2))
    for (i in 1:4)
        plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
    dev.off()
    
    res <- pcadapt(pcfile1, K = 20, LD.clumping = list(size = 200, thr = 0.1))
    # plot(res, option = "screeplot")
    
    res <- pcadapt(pcfile1, K = 4, LD.clumping = list(size = 10, thr = 0.1))
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/ldcheck3.pdf"), 10, 10)
    par(mfrow = c(2, 2))
    for (i in 1:4)
        plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
    dev.off()
    
    pdf(paste0("rotation1scripts_v4/plots/pcadapt/", name1, "/manhattan.after.thinning.pdf"), 10, 10)
    plot(res)
    dev.off()
    
    #### processing ####
    
    sig.markers1 = genodat1[which(pvals1$pvalue < 0.05), 1:3]
    
    return(list(plot1, sig.markers1))
    
    
    
}

library(parallel)
analyses2 = Map(function(x, y){
    run.pcadapt.anal(x, y)
}, popnames1, kvalues1)


load("rotation1scripts_v4/saved.objects/pcadapt/analyses2")
anal.plots = lapply(analyses2, function(x) x[[1]])

do.call(grid.arrange, anal.plots)

functional.anno = read.csv("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.csv", stringsAsFactors = F)

hcgff = read.delim("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", sep = "\t", stringsAsFactors = F)

#### GRAB GENES WITHIN 1MB ####

## run on wilkins (parallel)
#get dataframe of functional annotations nearest to markers exhibiting signicant selection signatures via pcadapt
genesfirst = lapply(analyses2, function(z){

    z2 = as.data.frame(t(z[[2]]))
    
    gene.func = mclapply(z2, function(x){
        chrom = x[1]
        pos = as.numeric(x[2])
        hcgff2 = hcgff[which(hcgff$V2 == paste0("chr", chrom)), ]
        nearby.genes = hcgff[which(hcgff2$V5 < (pos + 2500000) & hcgff2$V5 > (pos - 2500000)), ]
        nearby.genes2 = multi.str.split(nearby.genes$V10, "ID=", 2)
        nearby.genes3 = multi.str.split(nearby.genes2, ";", 1)
        func.coord = unlist(sapply(nearby.genes3, function(x){
            grep(x, functional.anno$Gene.ID)
        }))
        
        genes.funcs = functional.anno[func.coord, ]
        
        genes.funcs
    }, mc.cores = 30)
    

    
    gene.func2 = bind_rows(gene.func)
    print("done an interation")
    gene.func2
})


pont.genes = c("TraesCS2B01G069900", "TraesCS2B01G070100", "TraesCS2B01G070200", "TraesCS2B01G070300", "TraesCS2B01G070600", "TraesCS2D01G058000", "TraesCS2D01G058100", "TraesCS2D01G058200", "TraesCS3A01G101100", "TraesCS3A01G101200", "TraesCS3B01G118900", "TraesCS5A01G473800", "TraesCS5B01G256100", "TraesCS5B01G486900", "TraesCS5D01G486600", "TraesCS2A01G081900", "TraesCS2D01G079600", "TraesCS2A01G116900", "TraesCS2B01G350400", "TraesCS3A01G489000", "TraesCS4A01G271000", "TraesCS4B01G043100", "TraesCS5A01G391700", "TraesCS7B01G013100", "TraesCS7A01G070100", "TraesCS4A01G418200", "TraesCS1A01G008200", "TraesCS1A01G040200", "TraesCS1B01G010400", "TraesCS1B01G013500", "TraesCS1D01G007400", "TraesCS6A01G048900", "TraesCS6A01G049800", "TraesCS6B01G086000", "TraesCS6B01G086500")

lapply(pont.genes, function(x){
    grep(x, genesfirst[[6]]$Gene.ID)
})

functional.anno[grep("TraesCS5B01G256100", functional.anno$Gene.ID), ]



#### MANUAL ANALYSIS OF GENES ####


process.blast.allsnps = function(blast.path){
    #returns the distance from the nearest Axiom SNP to the gene
    #we can then use this as a benchmark to determine how far significant selection signals (from pcadapt) are from the gene of interest
    blast1 = read.blast(blast.path)
    chr1 = multi.str.split(as.character(blast1$sseqid[1]), ":", 1)
    chr2 = multi.str.split(chr1, "chr", 2)
    luzie.order3 = luzie.order2[which(luzie.order2$WGA.Chr == chr2), ]
    pos1 = multi.str.split(as.character(blast1$sseqid[1]), ":", 2)
    pos2 = as.numeric(multi.str.split(pos1, "-", 1))
    abs(luzie.order3[which.min(abs(luzie.order3$WGA.bp - pos2)), ]$WGA.bp - pos2)
}

process.blast.sigsnps = function(blast.path){
    #returns the distance from the nearest Axiom SNP to the gene
    #we can then use this as a benchmark to determine how far significant selection signals (from pcadapt) are from the gene of interest
    blast1 = read.blast(blast.path)
    chr1 = multi.str.split(as.character(blast1$sseqid[1]), ":", 1)
    chr2 = multi.str.split(chr1, "chr", 2)
    luzie.order3 = luzie.order2[which(luzie.order2$WGA.Chr == chr2), ]
    pos1 = multi.str.split(as.character(blast1$sseqid[1]), ":", 2)
    pos2 = as.numeric(multi.str.split(pos1, "-", 1))
    
    #rht gene
    lapply(analyses2, function(x){
        x.chr = x[[2]][which(x[[2]]$V1 == chr2), ]
        min(na.omit(abs(x.chr$V2 - pos2)))
    })
    
}


tags5 = read.blast("bioinf/blast/probe.vs.genes.blast/results.blast/TaGS5-3A.blast")


#rht analysis
process.blast.allsnps("bioinf/blast/probe.vs.genes.blast/results.blast/rht-b1b.blast")
process.blast.sigsnps("bioinf/blast/probe.vs.genes.blast/results.blast/rht-b1b.blast")

#vrn gene analysis

process.blast.allsnps("bioinf/blast/probe.vs.genes.blast/results.blast/vrn.blast")
process.blast.sigsnps("bioinf/blast/probe.vs.genes.blast/results.blast/vrn.blast")

#ppd gene analysis 
# https://www.ncbi.nlm.nih.gov/nuccore/KJ147477.1
process.blast.allsnps("bioinf/blast/probe.vs.genes.blast/results.blast/ppd.blast")
process.blast.sigsnps("bioinf/blast/probe.vs.genes.blast/results.blast/ppd.blast")






