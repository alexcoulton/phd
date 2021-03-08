setwd("~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/")
library(Biostrings)

list.files()

g = read.delim('~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3', sep = "\t", header = F)


genomic = readDNAStringSet('~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/genomic.genes.w.names.fa')
genes = readDNAStringSet('~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa')

source('~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R')

names(genomic) = multi.str.split(names(genomic), "_", 2)
names(genes) = multi.str.split(names(genes), "\\.", 1)

genomic2 = genomic[na.omit(match(names(genes), names(genomic)))]

comparison.lengths = data.frame(width(genes), width(genomic2))

diffs1 = abs(comparison.lengths[, 1] - comparison.lengths[, 2])

diffs2 = sort(diffs1, decreasing = T)