#evaluate best primers

primer.file1 = readLines("rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/primers/best.primers/set1/1primer3.186-3907.txt.output.txt")
primer.file1.seq = strsplit(primer.file1[grep("PRIMER_LEFT_0_SEQUENCE", primer.file1)], "=")
primer.file1.seq = DNAString(primer.file1.seq[[1]][2])

primer.file1.seq.rev = strsplit(primer.file1[grep("PRIMER_RIGHT_0_SEQUENCE", primer.file1)], "=")
primer.file1.seq.rev = primer.file1.seq.rev[[1]][2]

primer.file1.seq.rev = reverseComplement(DNAString(primer.file1.seq.rev))

both.seqs = c(DNAStringSet(primer.file1.seq), DNAStringSet(primer.file1.seq.rev))
names(both.seqs) = c("forward.primer", "reverse.primer")

writeXStringSet(both.seqs, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/extended/primer.seqs.fa")
