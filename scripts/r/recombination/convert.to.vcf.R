setwd("E:/phd.project.main/")
library(readr)
source("rotation1scripts_v4/scripts/r/functions.R")

axptemps1 = c("10(deg)C", "14(deg)C", "26(deg)C", "28(deg)C")

#apogee x paragon temp shift experiment - from recombination_analysisv3.R

load.data = function(x) read.csv(paste("rotation1scripts_v4/processed_data/genotypes/alexq", x, "_man_cur_sk_map_genogeno.csv", sep = ""), 
                                                                 header = T, stringsAsFactors = F)

#put in order of temperature, q1 = 10Â(deg)C, q3 = 14Â(deg)C, q2 = 26Â(deg)C, q4 = 28Â(deg)C
initial.data = lapply(c(1, 3, 2, 4), load.data)


convert.rqtl.to.vcf = function(rqtl.df){
    rqtl.df2 = rqtl.df[-1, -1]
    
    cereals3 = rqtl.df2
    colnames(cereals3)[2:ncol(cereals3)] = switch.affy.format(colnames(cereals3)[2:ncol(cereals3)])
    
    
    # get nucleotide bases for each SNP
    mysql1 = setup.mysql.connection(T)
    rs = dbSendQuery(mysql1, 'SELECT * FROM axiom35;')
    rs2 = fetch(rs, n=-1)
    dbDisconnect(mysql1)
    luzie.order = read_csv("rotation1scripts_v4/original_data/IWGSC/Luzie.probes.JIC_RefSeq1_axiom_order_Nov18.csv")
    # cereals3 = cereals2[, c(1, na.omit(match(luzie.order$SNP.id, colnames(cereals2))))]
    
    
    luzie.order2 = luzie.order[na.omit(match(colnames(cereals3), luzie.order$SNP.id)), ]
    
    cereals3[cereals3 == "A"] = "0/0"
    cereals3[cereals3 == "H"] = "0/1"
    cereals3[cereals3 == "B"] = "1/1"
    cereals3[cereals3 == "-"] = "./."
    
    cereals3 = as.data.frame(t(cereals3))
    
    cereals3 = as.data.frame(cbind(cereals3[, 1], cereals3), stringsAsFactors = F)
    cereals3[, 1] = rownames(cereals3)
    
    cereals3 = as.data.frame(cbind(cereals3[, 1], cereals3), stringsAsFactors = F)
    
    # match(luzie.order$SNP.id, cereals3[, 1])
    
    #add basepair positions of probes
    cereals3[, 1] = luzie.order[match(cereals3[, 1], luzie.order$SNP.id), 3]
    
    
    add.row1 = function(x){
        x = as.data.frame(cbind(x[, 1], x), stringsAsFactors = F)
        x[, 1] = ""
        x
    }
    
    cereals3 = add.row1(cereals3)
    
    #add chromosomes to dataframe
    cereals3[, 1] = luzie.order[match(rownames(cereals3), luzie.order$SNP.id), 2]
    
    cereals3 = add.row1(cereals3)
    cereals3 = add.row1(cereals3)
    cereals3 = add.row1(cereals3)
    cereals3 = add.row1(cereals3)
    cereals3 = add.row1(cereals3)
    cereals3 = add.row1(cereals3)
    
    
    
    
    #put nucleotide bases file in same order as VCF file
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
    cereals3[, 4] = "PASS"
    cereals3[, 5] = "ConversionType=PolyHighResolution"
    cereals3[, 6] = "GT"
    
    cereals4 = as.data.frame(cbind(cereals3[, 7:9], cereals3[, 1:6], cereals3[, 10:ncol(cereals3)]), stringsAsFactors = F)
    cereals4[1, 1:9] = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    cereals4[which(is.na(cereals4[, 2])), 2] = "."
    cereals4[] = lapply(cereals4, function(x) as.character(x))
    colnames(cereals4) = cereals4[1, ]
    cereals4 = cereals4[-1, ]    
    cereals4
}


initial.data.vcf = lapply(initial.data, convert.rqtl.to.vcf)

count1 = 1
lapply(initial.data.vcf, function(x){
    write.table(x, paste0("rotation1scripts_v4/original_data/genotype.data/vcf/axp", count1, ".vcf"), quote = F, sep = "\t")
    count1 <<- count1 + 1
})



g = data.frame(c(colnames(initial.data.vcf[[1]]), colnames(initial.data.vcf[[2]]), colnames(initial.data.vcf[[3]]), colnames(initial.data.vcf[[4]])))
write.csv(g, "rotation1scripts_v4/temp/axp.vcf.samp.names.csv")
