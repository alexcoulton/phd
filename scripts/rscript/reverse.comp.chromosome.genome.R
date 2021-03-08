#!/usr/bin/Rscript
#Takes a fasta file as input, returns the same file but with the reverse complement
#of the specified chromosome

#args:
#script input.fa chromo.to.rev 2nd.chromo.to.rev ...
args = commandArgs(trailingOnly = T)

base.path = args[1]
file.name = args[2]

if(length(args) < 2){
	print("args: script base.path input.fa chromo.to.rev 2nd.chromo.to.rev ...")
} else {
	library(Biostrings)

	genome = readDNAStringSet(paste(base.path, file.name, sep = ""))
	
	for(i in 3:length(args)){
		genome[[args[[i]]]] = reverseComplement(genome[[args[[i]]]])
	}

	writeXStringSet(genome, paste("./", file.name, ".rev.fa", sep = ""))

}
