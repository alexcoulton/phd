library(Biostrings)
library(dplyr)
setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

g1 = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/TraesCS5A01G531300.genomic.vs.para.blast")

#remove low-value blast hits 
g1.1 = filter(g1, alignment.length > 50) # only four scaffolds remain
g1.2 = sort.blastdf(g1.1)


ms2.gene.cs = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/extended/seqs/TraesCS5A01G531300.genomic.fa")
ms2.gene.cs.cds = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300.fa")

load("rotation1scripts_v4/saved.objects/ms2.para.scaffolds")
q = pairwiseAlignment(ms2.gene.cs, paragon2[[1]])
q1.1 = pairwiseAlignment(ms2.gene.cs.cds, paragon2[[1]])



q2 = pairwiseAlignment(ms2.gene.cs, paragon2[[2]])
q3 = pairwiseAlignment(ms2.gene.cs, reverseComplement(paragon2[[3]]))
q4 = pairwiseAlignment(ms2.gene.cs, reverseComplement(paragon2[[4]]))

apo.vs.cs = pairwiseAlignment(apogee.consensus.trim, reverseComplement(ms2.gene.cs))
writePairwiseAlignments(apo.vs.cs, "rotation1scripts_v4/processed_data/TILLING/t300.apo.vs.cs.txt")

writePairwiseAlignments(q1.1, "rotation1scripts_v4/processed_data/TILLING/pairwise.ms2.align.1.1.txt")

apogee.consensus.trim = readDNAStringSet("rotation1scripts_v4/original_data/fasta/sanger.sequencing/consensus.rev.trim.fa")

testpair2 = pairwiseAlignment(apogee.consensus.trim, ms2.gene.cs.cds, type = "global")
writePairwiseAlignments(testpair, "rotation1scripts_v4/processed_data/TILLING/testalign2.txt")

testpair = pairwiseAlignment(reverseComplement(paragon2[[1]]), ms2.gene.cs.cds, type = "global")
writePairwiseAlignments(testpair, "rotation1scripts_v4/processed_data/TILLING/testalign.txt")

paragon.scaff = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/full.paragon.scaff.ms2.fa")

pa.1 = pairwiseAlignment(ms2.gene.cs.cds, reverseComplement(paragon.scaff[[1]]))
pa = pairwiseAlignment(ms2.gene.cs, paragon.scaff[[2]], type = "local")
pa2 = pairwiseAlignment(ms2.gene.cs, reverseComplement(paragon.scaff[[3]]), type = "local")
pa3 = pairwiseAlignment(ms2.gene.cs, reverseComplement(paragon.scaff[[4]]), type = "local")

writePairwiseAlignments(pa.1, "rotation1scripts_v4/processed_data/TILLING/para.ms2.cds.align.txt")


writePairwiseAlignments(pa3, "rotation1scripts_v4/processed_data/TILLING/pair.align.full.scaff.4.txt")

#extract only exons from paragon sequence
pattern(testpair)
mismatchTable(testpair)
ins = insertion(testpair)

DNAString(as.character(pattern(testpair)))[ins]

testpair[1]
aligned(testpair)
writeXStringSet(aligned(testpair), "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/t.300.apogee.consensus.exons.fa")

apo.vs.para = pairwiseAlignment(aligned(testpair), aligned(testpair2))
writePairwiseAlignments(apo.vs.para, "rotation1scripts_v4/processed_data/TILLING/apo.vs.para.exon.align.txt")


translate(aligned(testpair2)) == translate(ms2.gene.cs.cds)



protein.align = pairwiseAlignment(translate(aligned(testpair)), translate(ms2.gene.cs.cds))
protein.align2 = pairwiseAlignment(translate(aligned(testpair)), translate(aligned(testpair2)))
writePairwiseAlignments(protein.align, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/ms2.protein.align.txt")
writePairwiseAlignments(protein.align2, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/ms2.protein.align2.txt")



#     ____________________________________________________________________________
#     check male sterility genes                                                                                            ####

male.sterility.genes = functional.anno[grep("Male sterility", functional.anno$Pfam.IDs..Description.), ]

m.st = multi.str.split(as.character(male.sterility.genes$Gene.ID), "\\.", 1)

m.genes = genes[unlist(lapply(m.st, function(x) grep(x, genes$gene.names))), ]




#     ____________________________________________________________________________
#     Sequence alignment for second 5a candidate    TraesCS5A01G531700                 ####

g1 = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/second.5a.candidate.genomic.fa")
g1.1 = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/second.5a.candidate.fa")
g2 = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/Triticum_aestivum_Paragon_EIv1.1_scaffold_014159.fa")

g3 = readDNAStringSet("rotation1scripts_v4/original_data/fasta/sanger.sequencing/t.700.aligned.fa")


g2 = g2[[1]][50000:54000]

g.aligna.x.p = pairwiseAlignment(g2, g3)
writePairwiseAlignments(g.aligna.x.p, "rotation1scripts_v4/processed_data/TILLING/t.700.apo.para.align.txt")


g.align.cs.x.a = pairwiseAlignment(g1, g3)
writePairwiseAlignments(g.align.cs.x.a, "rotation1scripts_v4/processed_data/TILLING/t.700.apo.cs.txt")

g.align = pairwiseAlignment(g1, g2, type = "global")
g.align2 = pairwiseAlignment(g1.1, g2, type = "global")
writePairwiseAlignments(g.align, "rotation1scripts_v4/processed_data/TILLING/pair.second.5a.candidate.txt")
writePairwiseAlignments(g.align2, "rotation1scripts_v4/processed_data/TILLING/pair.second.5a.candidate.cds.txt")


g.align3 = pairwiseAlignment(g3, g1.1, type = "global")
writePairwiseAlignments(g.align3, "rotation1scripts_v4/processed_data/TILLING/pair.second.5a.candidate.cds.txt")

writeXStringSet(aligned(g.align3), "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/paragon.second.5a.candidate.fa")


apo.seq = aligned(g.align3)
para.seq = aligned(g.align3)

cs.vs.apo = pairwiseAlignment(translate(g1.1), translate(apo.seq))
apo.vs.para = pairwiseAlignment(translate(para.seq), translate(apo.seq))
cs.vs.para = pairwiseAlignment(translate(para.seq), translate(g1.1))
writePairwiseAlignments(g.align4, "rotation1scripts_v4/processed_data/TILLING/cs.x.apo.second.5a.candidate.aa.seq.txt")

cs.vs.apo.subs = as.data.frame(mismatchTable(cs.vs.apo))[, c(2, 4, 8)]
colnames(cs.vs.apo.subs) = c("Position", "Chinese Spring AA", "Apogee AA")
write.csv(cs.vs.apo.subs, "rotation1scripts_v4/processed_data/TILLING/cs.vs.apo.substitution.table.t700.csv")

apo.vs.para.subs = as.data.frame(mismatchTable(apo.vs.para))[, c(2, 4, 8)]
colnames(apo.vs.para.subs) = c("Position", "Paragon AA", "Apogee AA")
write.csv(apo.vs.para.subs, "rotation1scripts_v4/processed_data/TILLING/apo.vs.para.substitution.table.t700.csv")

cs.vs.para.subs = as.data.frame(mismatchTable(cs.vs.para))[, c(2, 4, 8)]
colnames(cs.vs.para.subs) = c("Position", "Paragon AA", "Chinese Spring AA")
write.csv(cs.vs.para.subs, "rotation1scripts_v4/processed_data/TILLING/cs.vs.para.substitution.table.t700.csv")



wider.blast = read.blast("bioinf/blast/inflo.transcriptome/results.blast/feng.genes.vs.refseq.5a.wider.blast")


#     ____________________________________________________________________________
#     GENE EXTRACTION FROM OTHER VARIETIES                                                                        ####

#TraesCS5A01G531300

cad.blast = read.blast("bioinf/blast/cadenza.genome/results.blast/TraesCS5A01G531300.genomic.vs.cadenza.blast")

robigus.blast = read.blast("bioinf/blast/robigus.genome/results.blast/TraesCS5A01G531300.genomic.vs.robigus.blast")

claire.blast = read.blast("bioinf/blast/claire.genome/results.blast/TraesCS5A01G531300.genomic.vs.claire.blast")

cad.scaff = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/Triticum_aestivum_Cadenza_EIv1.1_scaffold_003792.fa")

cad.scaff.align = pairwiseAlignment(reverseComplement(all.genomic.sequences3[[2]]), cad.scaff, robigus.scaff, type = "global")
writePairwiseAlignments(cad.scaff.align, "rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/cad.scaff.align.txt")

robigus.scaff = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/Triticum_aestivum_Robigus_EIv1.1_scaffold_003114.fa")

para.blast = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/TraesCS5A01G531300.genomic.vs.para.blast")


homeologue.5b = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/Traes300.5b.fa")

homo.align = pairwiseAlignment(homeologue.5b, all.genomic.sequences3[[2]])
writePairwiseAlignments(homo.align, "rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/5bhomoalign.txt")

#want to trim some of the sequences

Traes300.seqs = readDNAStringSet("processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/seq/all.fa")

apogee.consensus = readDNAStringSet("original_data/fasta/sanger.sequencing/consensus.fa")
apogee.rev = reverseComplement(apogee.consensus)
writeXStringSet(apogee.rev, "original_data/fasta/sanger.sequencing/consensus.rev.fa")

test.pair1 = pairwiseAlignment(Traes300.seqs[[7]], reverseComplement(apogee.consensus))
writePairwiseAlignments(test.pair1, "processed_data/TILLING/apo.traes300.cs.align.txt")

test.pair2 = pairwiseAlignment(Traes300.seqs[[3]], reverseComplement(apogee.consensus))
writePairwiseAlignments(test.pair2, "processed_data/TILLING/apo.traes300.paragon.align.txt")

traes.alignments = pairwiseAlignment(Traes300.seqs, Traes300.seqs[[7]], type = "global")

writePairwiseAlignments(traes.alignments, "rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531300/all.alignments.txt")

Traes300.seqs[[1]] = Traes300.seqs[[1]][800:nchar(Traes300.seqs[[1]])]
Traes300.seqs[[2]] = Traes300.seqs[[2]][0:(nchar(Traes300.seqs[[2]])-500)]
Traes300.seqs[[3]] = Traes300.seqs[[3]][900:(nchar(Traes300.seqs[[3]])-900)]
Traes300.seqs[[4]] = Traes300.seqs[[4]][1600:nchar(Traes300.seqs[[4]])]

writeXStringSet(Traes300.seqs, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531300/seq/all.trimmed.fa")



# --------------- TraesCS5A01G531700 --

cad1.blast = read.blast("bioinf/blast/cadenza.genome/results.blast/TraesCS5A01G531700.genomic.vs.cadenza.blast")
robigus1.blast = read.blast("bioinf/blast/robigus.genome/results.blast/TraesCS5A01G531700.genomic.vs.robigus.blast")
claire1.blast = read.blast("bioinf/blast/claire.genome/results.blast/TraesCS5A01G531700.genomic.vs.claire.blast")
para1.blast = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/TraesCS5A01G531700.genomic.vs.para.blast")


gen1 = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531700/TraesCS5A01G531700.genomic.fa")
gen2 = readDNAStringSet("rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/TraesCS5A01G531700/rob.cut.fa")

align1 = pairwiseAlignment(reverseComplement(gen1), gen2)
writePairwiseAlignments(align1, "rotation1scripts_v4/processed_data/TILLING/TraesCS5A01G531700/alignvscad.txt")








