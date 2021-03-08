g10 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/all.snps.all.replicates.txt")

g10.axc611 = g10[grep("AxC 61-1", g10$probeset_id), ]
g10.axc512 = g10[grep("AxC 51-2", g10$probeset_id), ]
g10.cxa653 = g10[grep("CxA 65-3", g10$probeset_id), ]
g10.cxa651 = g10[grep("CxA 65-1", g10$probeset_id), ]

allg10 = list(g10.axc512, g10.axc611, g10.cxa653, g10.cxa651)
lapply(allg10, function(x){
    length(which(x == "-"))
})

allg10.count = lapply(allg10, function(x){
    lapply(x, geno.count)
})

allg10.mono = lapply(allg10, function(x){
    g = sapply(x, function(y){
        g = which(y == "-")
        if(length(g) > 0) y = y[-g]
        length(unique(y))
    })
    length(which(g == 1))
})



sample.g10 = read.table("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/all.replicates.sample.info.txt", sep = "\t", header = T)
sample.g10.axc611 = sample.g10[grep("AxC 61-1", sample.g10$Sample.Filename), ]
sample.g10.axc512 = sample.g10[grep("AxC 51-2", sample.g10$Sample.Filename), ]
sample.g10.cxa653 = sample.g10[grep("CxA 65-3", sample.g10$Sample.Filename), ]
sample.g10.cxa651 = sample.g10[grep("CxA 65-1", sample.g10$Sample.Filename), ]

all.sg10 = list(sample.g10.axc611, sample.g10.axc512, sample.g10.cxa653, sample.g10.cxa651)

lapply(all.sg10, function(x){
    mean(x$DQC)
})

lapply(all.sg10, function(x){
    mean(x$QC.call_rate)
})

lapply(all.sg10, function(x){
    mean(x$Average.call.rate.for.passing.samples)
})



#### AXP ####

q10 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.alex.f2/axp.all.quads.all.snps.txt")

q10.q1 = q10[grep("Q1", q10$probeset_id), ]
q10.q2 = q10[grep("Q2", q10$probeset_id), ]
q10.q3 = q10[grep("Q3", q10$probeset_id), ]
q10.q4 = q10[grep("Q4", q10$probeset_id), ]

allq = list(q10.q1, q10.q2, q10.q3, q10.q4)


lapply(allq, function(x){
    length(which(x == "-"))
})

allg10.count = lapply(allg10, function(x){
    lapply(x, geno.count)
})

allg10.mono = lapply(allg10, function(x){
    g = sapply(x, function(y){
        g = which(y == "-")
        if(length(g) > 0) y = y[-g]
        length(unique(y))
    })
    length(which(g == 1))
})




