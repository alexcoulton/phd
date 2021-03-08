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
        
    load("rotation1scripts_v4/saved.objects/cs.x.p.map.w.genocounts")

    library(dplyr)

    # -------- functions -
    library(compiler) #function for finding longest increasing subsequence
    #from https://www.r-bloggers.com/compute-longest-increasingdecreasing-subsequence-using-rcpp/
    longest_subseq.R = cmpfun(function(x) {
        P = integer(length(x))
        M = integer(length(x) + 1)
        L = newL = 0
        for (i in seq_along(x) - 1) {
            lo = 1
            hi = L
            while (lo <= hi) {
                mid = (lo + hi)%/%2
                if (x[M[mid + 1] + 1] < x[i + 1]) {
                    lo = mid + 1
                } else {
                    hi = mid - 1
                }
            }
            newL = lo
            P[i + 1] = M[newL]
            if (newL > L) {
                M[newL + 1] = i
                L = newL
            } else if (x[i + 1] < x[M[newL + 1] + 1]) {
                M[newL + 1] = i
            }
        }
        k = M[L + 1]
        re = integer(L)
        for (i in L:1) {
            re[i] = k + 1
            k = P[k + 1]
        }
        re
    })

    coords.of.skeleton.markers = function(cM.vec){
        match(unique(cM.vec), cM.vec)    
    }

    # --------- code

    #within bins, order cs.x.p.map.w.genocounts by physical sequence
    unique.chromo.cm.combos = unique(cs.x.p.map.w.genocounts[, 2:3])
    for(i in 1:nrow(unique.chromo.cm.combos)){
        g = cs.x.p.map.w.genocounts[which(cs.x.p.map.w.genocounts$chr == unique.chromo.cm.combos[i, 1] & cs.x.p.map.w.genocounts$cM == unique.chromo.cm.combos[i, 2]), ]
        sorted = sort(g$phys.dist.bp, na.last = F)
        g = g[match(sorted, g$phys.dist.bp), ]
        cs.x.p.map.w.genocounts[which(cs.x.p.map.w.genocounts$chr == unique.chromo.cm.combos[i, 1] & cs.x.p.map.w.genocounts$cM == unique.chromo.cm.combos[i, 2]), ] = g
        
    }

#     ____________________________________________________________________________
#     FIX FLIPPED CHROMOSOMES                                                                                                 ####


#check if any chromosomes are flipped 
flipped.chromosomes = lapply(listofwheatchromosomes, function(x){
    g = filter(cs.x.p.map.w.genocounts, chr == x)
    start.mean = mean(na.omit(g$phys.dist.bp[1:10]))
    end.mean = mean(na.omit(g$phys.dist.bp[(nrow(g)-10):nrow(g)]))
    # browser()
    if(start.mean > end.mean){
        return(x)
    } else{
        return()
    }
})

flipped.chromosomes = unlist(flipped.chromosomes)

#reverse selected chromosomes, reset centiMorgan values
flipped.chromosome.dfs = lapply(flipped.chromosomes, function(x){
    g = filter(cs.x.p.map.w.genocounts, chr == x)
    g = g[nrow(g):1, ]
    
    #reverse cm distance
    total.cm = max(as.numeric(g$cM))
    
    g$cM = abs(g$cM - total.cm)
    
    return(g)
})

names(flipped.chromosome.dfs) = flipped.chromosomes

#replace flipped chromosomes with newly inverted ones
for(i in names(flipped.chromosome.dfs)){
    cs.x.p.map.w.genocounts[which(cs.x.p.map.w.genocounts$chr == i), ] = flipped.chromosome.dfs[[i]]
}




    extract.monotonic.chromo = function(chromo){
        g = filter(cs.x.p.map.w.genocounts, chr == chromo)
             
        g[match(na.omit(g$phys.dist.bp)[longest_subseq.R(na.omit(g$phys.dist.bp))], g$phys.dist.bp), ]
    }

    g1a = extract.monotonic.chromo(chromosome.arg)
    s(g1a, p("rotation1scripts_v4/saved.objects/g1a.", chromosome.arg), "historical.recombination.mapping.server.script.unix.R")

    print("g1a")
    print(g1a)

    #extract genotyping information
    library(readr)
    all.geno.data = read_csv("rotation1scripts_v4/original_data/genotype.data/cerealsdb.35k.genotype.data.csv")

    print("all.geno.data")
    print(all.geno.data)

    extracted.geno.data = newdf(colnames(all.geno.data), no.rows = T)
    for(i in 1:nrow(g1a)){
        coords = which(all.geno.data[, 2] == g1a[i, 1])

        print("coords")
        print(coords)

        print(" ")
        print("all..geno.data[coords, ]")
        print(all.geno.data[coords, ])
        extracted.geno.data = rbind(extracted.geno.data, all.geno.data[coords, ])
    }

    write.csv(extracted.geno.data, p("rotation1scripts_v4/processed_data/genotypes/historical.recombination.mapping/", chromosome.arg, ".ske.incl.csv"))



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
}
