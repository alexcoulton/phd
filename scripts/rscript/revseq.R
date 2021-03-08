#!/usr/bin/Rscript
#Takes a fasta file as input, returns the same file but with the reverse complement
#of the specified chromosome

#args:
#script input.fa chromo.to.rev 2nd.chromo.to.rev ...
args = commandArgs(trailingOnly = T)

fasta.to.rev = args[1]

library(Biostrings)

dna = readDNAStringSet(fasta.to.rev)
dna = reverseComplement(dna)

writeXStringSet(dna, paste(fasta.to.rev, ".rev.fa", sep = ""))

