setwd('~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/')


list.files()
library(Biostrings)

all.genes = readDNAStringSet('genomic.genes.w.names.fa')

all.gc.content = sapply(all.genes, function(x){
	sum(letterFrequency(x, c('G', 'C'))) / length(x)
})


mean(all.gc.content)
sd(all.gc.content)

length(which(all.gc.content < 0.4))
