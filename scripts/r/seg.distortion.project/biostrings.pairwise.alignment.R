library(Biostrings)

source("rotation1scripts_v4/scripts/r/generic.functions/fasta.R")


query = readDNAStringSet(filepath = "rotation1scripts_v4/processed_data/fasta/unique.seg.dist.genes.full.genonmic.sequence.v3.fa")
subject = readDNAStringSet(filepath = "rotation1scripts_v4/processed_data/fasta/Triticum_aestivum_Paragon_EIv1.1_scaffold_073295.fa")

coord = which(names(query) == "TraesCS6B01G444700.1")






align2 = pairwiseAlignment(subject, type = "local")

align3 = pairwiseAlignment(query[coord], rev(subject[[1]][10000:13000]))

writePairwiseAlignments(align3, "rotation1scripts_v4/processed_data/pairwise.alignments/temp.align2.txt")


summary(align)

writePairwiseAlignments(align2, "rotation1scripts_v4/processed_data/pairwise.alignments/tempalign.txt")

gff = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3")

gff2 = gff[gff$V3 == "gene", ]
gff3 = gff[gff$V3 == "mRNA", ]


write.table(gff2, "rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.gff")