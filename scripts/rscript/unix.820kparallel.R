#!/usr/bin/Rscript
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")

args = commandArgs(trailingOnly = T)

blast.df.name = args[1]


blast820 = read.csv(p("rotation1scripts_v4/processed_data/BLAST/820k/820.blast.subsets.4b.rev/", blast.df.name), header = T, stringsAsFactors = F)

# blast820 = blast820[, -1]

blast820.1 = remove.hits.with.same.best.e.value(blast820)
blast820.2 = grab.best.hits(blast820.1)

write.csv(blast820.2, p("rotation1scripts_v4/processed_data/BLAST/820k/only.hits.w.unique.best.e.value.4b.rev/", blast.df.name, ".csv"), row.names = F)
