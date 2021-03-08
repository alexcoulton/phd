#!/usr/bin/Rscript
#     ____________________________________________________________________________
#     HISTORICAL RECOMBINATION MAPPING                                                                                 ####
setwd("/home/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

args = commandArgs(trailingOnly = T)

if(!args %in% listofwheatchromosomes){
    print("Arguments: chromosome")
} else {
    chromosome.arg = args[[1]]
    library(dplyr)
    load("rotation1scripts_v4/saved.objects/sorted.chromosomes2")
    
    single.chromo = filter(sorted.chromosomes2, sseqid == p("chr", chromosome.arg))

    single.chromo.first.100 = single.chromo[1:100, ]    

    #extract genotyping information
    library(readr)
    all.geno.data = read_csv("rotation1scripts_v4/original_data/genotype.data/cerealsdb.35k.genotype.data.csv")

    print("all.geno.data")
    print(all.geno.data)

    extracted.geno.data = newdf(colnames(all.geno.data), no.rows = T)
    for(i in 1:nrow(single.chromo.first.100)){
        coords = which(all.geno.data[, 2] == single.chromo.first.100[i, 1])

        print("coords")
        print(coords)

        print(" ")
        print("all..geno.data[coords, ]")
        print(all.geno.data[coords, ])
        extracted.geno.data = rbind(extracted.geno.data, all.geno.data[coords, ])
    }

    write.csv(extracted.geno.data, p("rotation1scripts_v4/processed_data/genotypes/historical.recombination.mapping/zoom/", chromosome.arg, ".zoom.csv"))

}

# library(RMySQL)
# 
# #setup connection to mysql database
# mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')
# 
# ptm = proc.time()
# g2 = lapply(g1a$ï..marker, function(x){
#     
#     sql.query = p("SELECT * FROM 35k_array_oct2016_named WHERE Probe_row = '", x, "'")
#     
#     #pass mysql query to database
#     rs = dbSendQuery(mydb, sql.query)
#     rs2 = fetch(rs, n=-1)
#     
#     colnames(rs2)[3] = x
#     return(rs2)
# })
# 
# proc.time() - ptm

