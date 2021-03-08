#qtl attempt
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/")
library(qtl)
library(dplyr)
library(tibble)

path = function(quadrant) paste("original_data/genotype.data/alexq", quadrant, "_man_cur_sk_map_genogeno.csv", sep = "")


#     ____________________________________________________________________________
#     CONVERT PHENO DATA INTO RQTL FORMAT                                                                         ####


genotypeinfo = lapply(1:4, function(x) read.csv(path(x), header = T, stringsAsFactors = F))

load.pheno = function(q) read.csv(paste("original_data/phenotype.data/awn.pheno.q", q, ".csv", sep = ""), header = T, stringsAsFactors = F)

pheno.data = lapply(1:4, load.pheno)

pheno.data = lapply(pheno.data, function(x){
    rownames(x) = x[, 1]
    x = x[, 2:13]
    return(x)
})

conv.q.to.temp = function(x){
    if(x == 1) y = 10
    if(x == 2) y = 26
    if(x == 3) y = 14
    if(x == 4) y = 28
    return(y)
}

addphenodata = function(quad, genotypedata){

    matches = regexpr(paste(conv.q.to.temp(quad), "_[[:upper:]][[:digit:]]{1,}\\.", sep = ""), genotypedata[grep(paste("_Q", quad, "_", sep = ""), genotypedata[,2]), 2])
    
    list.of.coords = substr(genotypedata[grep(paste("_Q", quad, "_", sep = ""), genotypedata[,2]), 2], matches + 3, matches + attributes(matches)$match.length - 2)
    
    conv.letter.to.coord = function(x){
        which(LETTERS == x)
    }
    
    pheno = lapply(list.of.coords, function(x){
        letter = substr(x, 1, 1)
        letter = conv.letter.to.coord(letter)
        pheno.data[[quad]][letter, as.numeric(substr(x, 2, nchar(x)))]
    })
    
    pheno = as.character(pheno)
    
    print(pheno)
    
    genotypedata$phenotype = ""
    
    genotypedata$phenotype[which(genotypedata[, 2] %in% genotypedata[grep(paste("_Q", quad, "_", sep = ""), genotypedata[,2]), 2])] = pheno
    
    pheno
    return(genotypedata)
}

genotypeinfo = Map(addphenodata, 1:4, genotypeinfo)

genotypeinfo = lapply(genotypeinfo, function(x){
    x = add_column(x, .before = "X", pheno = "")
    x$pheno = x$phenotype
    return(x)
})

code.pheno.values.integer = function(x){
    y = ""
    if(x == "Y") y = 1
    if(x == "N") y = 0 
    if(x == "-") y = "-"
    return(y)
}

genotypeinfo = lapply(genotypeinfo, function(x){
    x$pheno = lapply(x$pheno, code.pheno.values.integer)    
    return(x)
})

genotypeinfo = lapply(genotypeinfo, function(x){
    x = x[, -2]
    return(x)
})

genotypeinfo = lapply(genotypeinfo, function(x){
    x$pheno[which(x$pheno == "-")] = NA
    return(x)
})

genotypeinfo = lapply(genotypeinfo, function(x){
    rownames(x) = x[, 2]
    x = x[, -2]
    return(x)
})

combinedgenotypes = rbind(genotypeinfo[[1]], genotypeinfo[[2]][-(1:3), ], genotypeinfo[[3]][-(1:3), ], genotypeinfo[[4]][-(1:3), ])

combinedgenotypes$pheno[which(combinedgenotypes$pheno == "")] = NA
combinedgenotypes$pheno[1] = ""
combinedgenotypes = combinedgenotypes[, -ncol(combinedgenotypes)]


combinedgenotypes[] = lapply(combinedgenotypes, as.character)


write.csv(combinedgenotypes, "original_data/genotype.data/alex_combined_man_cur_sk_map_genogeno_w_phenotypes.csv", row.names = F)


#     ____________________________________________________________________________
#     BEGIN QTL ANALYSIS                                                                                                            ####

hyper = read.cross(format = "csv", dir = "./", file = "original_data/genotype.data/alex_combined_man_cur_sk_map_genogeno_w_phenotypes.csv", genotypes = c("A", "H", "B"), estimate.map = F)

a2 = sim.geno(hyper, step = 2.5, n.draws = 16, error.prob = 0.01)

out1 <- scanone(hyper, method="imp")
plot(out1)
max(out1)
g = pull.geno(hyper)[, "AX.94613491"]
mean(is.na(g))
g <- pull.geno(fill.geno(hyper))[,"AX.94613491"]
out1.c4 <- scanone(hyper, method="imp", addcovar=g)
plot(out1, out1.c4, col=c("blue", "red"))
plot(out1.c4 - out1, ylim=c(-3,3))
abline(h=0, lty=2, col="gray")
