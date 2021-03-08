source("rotation1scripts_v4/scripts/r/functions.R")
library(Biostrings)
probe.seq = readDNAStringSet("rotation1scripts_v4/original_data/fasta/35kprobes.fa")
load("rotation1scripts_v4/saved.objects/c5")

names(probe.seq) = switch.affy.format(names(probe.seq))

probe.seq2 = probe.seq[which(names(probe.seq) %in% c5$marker)]
probe.seq3 = probe.seq2[match(c5$marker, names(probe.seq2))]

snp.pos = sapply(probe.seq3, function(x){
	g = as.character(x)
	g2 = strsplit(g, "")[[1]]	
	which(!g2 %in% c("A", "G", "T", "C"))
	})

snp1 = sapply(probe.seq3, function(x){
	g = as.character(x)
	g2 = strsplit(g, "")[[1]]	
	g2[which(!g2 %in% c("A", "G", "T", "C"))]
	})


c5$snp.pos = snp.pos
c5$snp1 = snp1
writeXStringSet(probe.seq3, "bioinf/blast/probe.vs.genome.blast/fasta/axpf2.probes.fa")

axpblast = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/axpf2.probes.blast")

axpblast$sseqid = multi.str.split(as.character(axpblast$sseqid), "chr", 2)

axpblast2 = split(axpblast, axpblast$qseqid)

axpblast3 = lapply(axpblast2, function(x){
	chromo1 = c5[match(x$qseqid[[1]], c5$marker), ]$chr1
	x2 = x[which(x$sseqid == chromo1), ]
	x2
	})

library(dplyr)

axpblast4 = bind_rows(axpblast3)


iwgsc = readDNAStringSet("genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta")
names(iwgsc) = multi.str.split(names(iwgsc), "chr", 2)

c6 = c5[which(c5$marker %in% axpblast4$qseqid), ]

axpblast5 = axpblast4[match(c6$marker, axpblast4$qseqid), ]

c7 = bind_cols(c6, axpblast5)

c8 = split.df.into.list.of.dfs.by.column(c7, "marker", F)

c9 = lapply(c8, function(x){
	snp.pos2 = x$snp.pos - x$qstart	
	iwgsc.base1 = as.character(iwgsc[[x$chr1]][(x$sstart + snp.pos2)])
	x$iwgsc.snp.pos = (x$sstart + snp.pos2)
	x$iwgsc = iwgsc.base1
	x
	})

r.code = c("A", "G")
m.code = c("A", "C")
w.code = c("A", "T")
s.code = c("C", "G")
y.code = c("C", "T")
k.code = c("G", "T")

check1 = sapply(c9, function(x){
	if(x$snp1 == "M"){
		g = x$iwgsc %in% m.code
	}
	if(x$snp1 == "R"){
		g = x$iwgsc %in% r.code
	}
	if(x$snp1 == "W"){
		g = x$iwgsc %in% w.code
	}
	if(x$snp1 == "S"){
		g = x$iwgsc %in% s.code
	}
	if(x$snp1 == "Y"){
		g = x$iwgsc %in% y.code
	}
	if(x$snp1 == "K"){
		g = x$iwgsc %in% k.code
	}
	
	g

	})

c10 = c9[which(check1 == T)]
c11 = bind_rows(c10)

c12 = split.df.into.list.of.dfs.by.column(c11, "marker", F)

c13 = lapply(c12, function(x){
	if(x$snp1 == "M"){
		x$alt = m.code[-which(m.code == x$iwgsc)]
	}
	if(x$snp1 == "R"){
		x$alt = r.code[-which(r.code == x$iwgsc)]
	}
	if(x$snp1 == "W"){
		x$alt = w.code[-which(w.code == x$iwgsc)]
	}
	if(x$snp1 == "S"){
		x$alt = s.code[-which(s.code == x$iwgsc)]
	}
	if(x$snp1 == "Y"){
		x$alt = y.code[-which(y.code == x$iwgsc)]
	}
	if(x$snp1 == "K"){
		x$alt = k.code[-which(k.code == x$iwgsc)]
	}
	x
	})

c14 = bind_rows(c13)

library(readr)
vcf1 = read_delim("rotation1scripts_v4/original_data/genotype.data/vcf/axp10degreesv2.vcf", delim = "\t", comment = "#", col_names = F)

vcf1.1 = vcf1[na.omit(match(c14$marker, vcf1$X3)), ]

c15 = c14[na.omit(match(vcf1.1$X3, c14$marker)), ]

vcf1.1$X1 = c15$chr1
vcf1.1$X2 = c15$iwgsc.snp.pos
vcf1.1$X4 = c15$iwgsc
vcf1.1$X5 = c15$alt


write_delim(vcf1.1, "rotation1scripts_v4/original_data/genotype.data/vcf/axp10degreesv3.vcf", delim = "\t")

vcf1 = read_delim("rotation1scripts_v4/original_data/genotype.data/vcf/axp14degreesv2.vcf", delim = "\t", comment = "#", col_names = F)

vcf1.1 = vcf1[na.omit(match(c14$marker, vcf1$X3)), ]

c15 = c14[na.omit(match(vcf1.1$X3, c14$marker)), ]

vcf1.1$X1 = c15$chr1
vcf1.1$X2 = c15$iwgsc.snp.pos
vcf1.1$X4 = c15$iwgsc
vcf1.1$X5 = c15$alt

write_delim(vcf1.1, "rotation1scripts_v4/original_data/genotype.data/vcf/axp14degreesv3.vcf", delim = "\t")

vcf1 = read_delim("rotation1scripts_v4/original_data/genotype.data/vcf/axp26degreesv2.vcf", delim = "\t", comment = "#", col_names = F)

vcf1.1 = vcf1[na.omit(match(c14$marker, vcf1$X3)), ]

c15 = c14[na.omit(match(vcf1.1$X3, c14$marker)), ]

vcf1.1$X1 = c15$chr1
vcf1.1$X2 = c15$iwgsc.snp.pos
vcf1.1$X4 = c15$iwgsc
vcf1.1$X5 = c15$alt

write_delim(vcf1.1, "rotation1scripts_v4/original_data/genotype.data/vcf/axp26degreesv3.vcf", delim = "\t")


vcf1 = read_delim("rotation1scripts_v4/original_data/genotype.data/vcf/axp28degreesv2.vcf", delim = "\t", comment = "#", col_names = F)

vcf1.1 = vcf1[na.omit(match(c14$marker, vcf1$X3)), ]

c15 = c14[na.omit(match(vcf1.1$X3, c14$marker)), ]

vcf1.1$X1 = c15$chr1
vcf1.1$X2 = c15$iwgsc.snp.pos
vcf1.1$X4 = c15$iwgsc
vcf1.1$X5 = c15$alt

write_delim(vcf1.1, "rotation1scripts_v4/original_data/genotype.data/vcf/axp28degreesv3.vcf", delim = "\t")

