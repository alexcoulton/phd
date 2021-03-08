source("rotation1scripts_v4/scripts/r/modelling/pedigree.sim.functions.R")

load("rotation1scripts_v4/saved.objects/cxa.6b.cm")



run.1000.sims("recombination.distribution.test.6b", 2, num.sims = 1000, w.selec = F, pop.size = 96, cm.pos = cxa.6b.cm)    

run.1000.sims("recombination.distribution.test.sparse", 2, num.sims = 100, w.selec = F, pop.size = 1000, cm.pos = seq(0, 100, 5))    


run.1000.sims("recombination.distribution.test.6b.500.pop", 2, num.sims = 100, w.selec = F, pop.size = 500, cm.pos = cxa.6b.cm)    




load("rotation1scripts_v4/saved.objects/axp.3b.cm")

run.1000.sims("recombination.distribution.test.axp3b.500.pop", 2, num.sims = 100, w.selec = F, pop.size = 500, cm.pos = axp.3b.cm)    
run.1000.sims("recombination.distribution.test.axp3b.500.pop.500.sims", 2, num.sims = 500, w.selec = F, pop.size = 500, cm.pos = axp.3b.cm)    



load("/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.rd.ind")

combinations1 = as.data.frame(combn(1:200, 4))

combinations1 = combinations1[, sample(ncol(combinations1), 100000)]
library(parallel)
library(dplyr)

# RUN ON WILKINS
num.sim.sig1 = mclapply(all.rd.ind, function(x){
    length(which(sapply(combinations1, function(y){
        
        
        m1 = x[[y[1]]]$phys.dist
        m2 = x[[y[2]]]$phys.dist
        m3 = x[[y[3]]]$phys.dist
        m4 = x[[y[4]]]$phys.dist
        
        m1 = data.frame(m1, "1")
        m2 = data.frame(m2, "2")
        m3 = data.frame(m3, "3")
        m4 = data.frame(m4, "4")
        
        colnames(m1) = c("phys.dist", "treatment")
        colnames(m2) = c("phys.dist", "treatment")
        colnames(m3) = c("phys.dist", "treatment")
        colnames(m4) = c("phys.dist", "treatment")
        
        m5 = bind_rows(m1, m2, m3, m4)
        # browser()
        # abs(mean(m1) - mean(m2))
        
        # browser()
        # hist(m1)
        # hist(m2)
        # wilcox.test(m1, m2)$p.value
        
        kruskal.test(m5$phys.dist ~ m5$treatment)$p.value
        
    }) < 0.05))
}, mc.cores = 50)