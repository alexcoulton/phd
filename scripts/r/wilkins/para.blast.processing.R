source("~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")
load("~/project.phd.main/rotation1scripts_v4/saved.objects/blast820.1")

library(dplyr)
library(parallel)
library(tibble)


count1 = make.counter()
blast820.1.list = mclapply(unique(blast820.1$qseqid), function(x){
	print(count1())
	g = filter(blast820.1, qseqid == x)
	g
}, mc.cores = 50)



# g2 = mclapply(blast820.1.list, function(x){
# 	remove.hits.with.same.best.e.value(x)
# }, mc.cores = 50)
