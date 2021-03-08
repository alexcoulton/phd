#!/home/ac14037/bin/Rscript

setwd("/home/ac14037/project.phd.main/")
library(Biostrings)
library(dplyr)
library(tibble)
source("rotation1scripts_v4/scripts/r/functions.R")

#evaluate best primers

base.path = "rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/best.primers/set2/"
primer.files = list.files(base.path)

lapply(primer.files, function(x){

	parse.coords = strsplit(x, "\\.")
	parse.coords = parse.coords[[1]][2]
	parse.coords = strsplit(parse.coords, "-")

	print("parse.coords:")
	print(parse.coords)

	primer.file1 = readLines(p(base.path, x))
	primer.file1.seq = strsplit(primer.file1[grep("PRIMER_LEFT_0_SEQUENCE", primer.file1)], "=")
	primer.file1.seq = DNAString(primer.file1.seq[[1]][2])

	primer.file1.seq.rev = strsplit(primer.file1[grep("PRIMER_RIGHT_0_SEQUENCE", primer.file1)], "=")
	primer.file1.seq.rev = primer.file1.seq.rev[[1]][2]

	primer.file1.seq.rev = reverseComplement(DNAString(primer.file1.seq.rev))

	both.seqs = c(DNAStringSet(primer.file1.seq), DNAStringSet(primer.file1.seq.rev))
	names(both.seqs) = c(p("forward.primer.", parse.coords[[1]][1]), p("reverse.primer.", parse.coords[[1]][2]))

	writeXStringSet(both.seqs, p("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/extended/primers.set2/", x, ".seqs.fa"))	
})



