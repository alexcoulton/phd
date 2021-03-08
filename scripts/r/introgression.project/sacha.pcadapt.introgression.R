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

genotype.data2 = genotype.data[, na.omit(match(colnames(genotype.data), axiom.order2))]


pcad1 = read.pcadapt(genotype.data2, type = "lfmm")
pcad2 = pcadapt(pcad1, K = 8)

padj = p.adjust(pcad2$pvalues, method = "BH")

pvals1 = data.frame(padj)
colnames(pvals1) = "pvalue"
pvals1$padjlog = -log(pvals1$pvalue)
pvals1$x = 1:nrow(pvals1)
pvals1$sig = ""
pvals1$sig[which(pvals1$pvalue < 0.001)] = T



eight20 = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/820kprobes.vs.iwgsc.v1.blast")


eight20.2 = split(eight20, eight20$qseqid)

lapply(eight20.2, function(x){
    browser()
})




# pvals1$chromo = genodat1[, 1]
pdf("rotation1scripts_v4/original_data/sacha.introgression.data/pcadapt.pdf", 10, 10)
ggplot(pvals1, aes(x = x, y = padjlog, color = sig)) + geom_point() + 
    theme_bw()
dev.off()
    
    







genotype.data2 = cbind(genotype.data[, 1], genotype.data) 
genotype.data2[, 1] = rownames(genotype.data2)
source("rotation1scripts_v4/scripts/r/functions.R")

genotype.data2 = reset.rownames(genotype.data2)

library(parallel)
genotype.data2[, ] = mclapply(genotype.data2[, ], as.character, mc.cores = 50)
colnames(genotype.data2) = genotype.data2[1, ]



introgressions = read.delim("rotation1scripts_v4/original_data/sacha.introgression.data/Introgression_table_280617.tdt", sep = "\t", header = F, stringsAsFactors = F)
colnames(introgressions) = c("elite", "chr", "start.bp", "length.bp", "cutoff", "source")

unique(introgressions$elite) %in% colnames(genotype.data) 

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'root', password = 'gamble', db = 'cerealsdb', host = '127.0.0.1')

rs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Avalon%'") #grab Avalon genotyping data for 35k array
rs2 = fetch(rs, n=-1)
rs2 = rs2[rs2$Var_col == "Avalon", ] 

pvals1 = read.delim("rotation1scripts_v4/original_data/sacha.introgression.data/pop1_pi.sites.pi", sep = "\t", stringsAsFactors = F)
pvals2 = pvals1[which(pvals1$CHROM == "chr1D" & pvals1$POS >= 484402678 & pvals1$POS <= 495229142), ]











