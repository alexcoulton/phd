#!/usr/bin/Rscript
library(Biostrings)

args = commandArgs(trailingOnly = T)

fasta.file = args[1]
output.file = args[2]

#start.coord = args[3]
#end.coord = args[4]

sequence.name = args[3]
start.coord = args[4]
end.coord = args[5]


dna = readDNAStringSet(fasta.file)
dna = dna[sequence.name]


dna[[1]] = dna[[1]][start.coord:end.coord]

writeXStringSet(dna, output.file)