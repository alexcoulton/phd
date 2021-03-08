#load genetic map that was made using only cxa651 CEL files as input to axiom analysis suite
cxaf2.c651.asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa.axc.all.replicates.mapv3.newprobeset.csv", stringsAsFactors = F)
cxaf2.c651.asmap = cxaf2.c651.asmap[sort(cxaf2.c651.asmap$chromo, index.return = T)$ix, ]
# cxaf2.c651.asmap = read.csv("rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa651-redone-map.csv")

cxa.geno.all = read.csv("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/all.replicates.flipped.geno.data.newprobeset.csv", stringsAsFactors = F)

colnames(cxa.geno.all)[4:ncol(cxa.geno.all)] = switch.affy.format(colnames(cxa.geno.all)[4:ncol(cxa.geno.all)])

#update genotyping marker order to match genetic map marker order
cxa.geno.all2 = cxa.geno.all[, c(1, 3, na.omit(match(cxaf2.c651.asmap$marker, colnames(cxa.geno.all))))]
cxaf2.c651.asmap = cxaf2.c651.asmap[which(cxaf2.c651.asmap$marker %in% colnames(cxa.geno.all2)), ]


#add chromosome information
cxa.geno.all2[1, 3:ncol(cxa.geno.all2)] = cxaf2.c651.asmap$chromo

#remove parents
cxa.geno.all2 = cxa.geno.all2[-c(2, 3), ]

cxa.651 = cxa.geno.all2[c(1, grep("65-1", cxa.geno.all2$probeset_id)), ]
cxa.653 = cxa.geno.all2[c(1, grep("65-3", cxa.geno.all2$probeset_id)), ]
axc.512 = cxa.geno.all2[c(1, grep("51-2", cxa.geno.all2$probeset_id)), ]
axc.611 = cxa.geno.all2[c(1, grep("61-1", cxa.geno.all2$probeset_id)), ]



list.of.geno.dfs.each.pop3 = list(cxa.651, cxa.653, axc.512, axc.611)

recomb2 = lapply(list.of.geno.dfs.each.pop3, detect.recombination)

#get number of recombination events for each individual
individual.events = lapply(recomb2, function(x){
    sort(sapply(unique(x$individual), function(y){
        # browser()
        sum(x[which(x$individual == y), ]$num.events)
    }))
})

indivi.lar1 = names(which(unlist(individual.events) > 110))


#remove individuals with excessive (erroneous) number of recombination events
list.of.geno.dfs.each.pop3 = lapply(list.of.geno.dfs.each.pop3, function(x){
    to.rm1 = which(x$probeset_id %in% indivi.lar1)
    
    if(length(to.rm1) > 0){
        x = x[-to.rm1, ]
    }
    
    x
})


asmaplgs = unique(as.character(list.of.geno.dfs.each.pop3[[1]][1, ]))[-1]

axc.lg.lengths = table(as.character(list.of.geno.dfs.each.pop3[[1]][1, ]))
axc.lg.lengths = axc.lg.lengths[-c(1, length(axc.lg.lengths))]

