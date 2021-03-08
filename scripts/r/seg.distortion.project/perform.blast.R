#perform BLAST search of sequence against genomes for primer picking

library(tibble)
library(dplyr)
library(Biostrings)
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")
setwd("/home/ac14037/project.phd.main")

#provide path of file to be BLASTed against genomes as an argument
args = commandArgs(trailingOnly = T)
query.fa.path = args[1]
query.fa.name = strsplit(query.fa.path, "/")
query.fa.name = query.fa.name[[1]]
query.fa.name = query.fa.name[length(query.fa.name)]

#read the configuration file
config.file = readLines("rotation1scripts_v4/auto.primer.picker/config.txt")
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

lapply(1:number.genomes, function(x){
    #begin by parsing the config file for the genome name, fasta file path and blastdb path
    config.settings = config.file[grep(x, config.variables)]
    
    config.settings.temp = config.settings[[1]]
    genome.name = strsplit(config.settings.temp[1], "=")
    genome.name = genome.name[[1]][2]
    
    config.settings.temp = config.settings[[2]]
    fa.path = strsplit(config.settings.temp[1], "=")
    fa.path = fa.path[[1]][2]
    
    config.settings.temp = config.settings[[3]]
    blastdb.path = strsplit(config.settings.temp[1], "=")
    blastdb.path = blastdb.path[[1]][2]
    
    system(p("blastn -db ", blastdb.path, " -query ", query.fa.path,
                     " -outfmt 6 -out ./blast.results/", query.fa.name, ".vs.", genome.name, ".blast",
                     " -num_threads 40 -culling_limit 10"))
    
})
