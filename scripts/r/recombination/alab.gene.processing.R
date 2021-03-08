library(Biostrings)

alab.genes1 = readDNAStringSet('/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/recombination.gene.analysis/alabdullah.genes1.fa')

alab.gene.names = multi.str.split(names(alab.genes1), ':', 1)

Map(function(x, y){
	writeXStringSet(DNAStringSet(x), paste0('~/project.phd.main/autoclonernew/primer/pipeline/sequences/alab.genes/', y, '.fa'))
}, alab.genes1, alab.gene.names)



