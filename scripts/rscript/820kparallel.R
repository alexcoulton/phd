#!/usr/bin/Rscript
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")

args = commandArgs(trailingOnly = T)

blast.df.name = args[1]

blast820 = read.blast(p("rotation1scripts_v4/processed_data/BLAST/820k/820.blast.subsets.4b.rev/", blast.df.name))

blast820.1 = remove.hits.with.same.best.e.value(blast820)

write.csv(blast820.1, p("rotation1scripts_v4/processed_data/BLAST/820k/only.hits.w.unique.best.e.value.4b.rev/", blast.df.name, ".csv"), row.names = F)
