#here I am using transcriptome data from wheat-expression.com

source("rotation1scripts_v4/scripts/r/functions.R")
library(dplyr)

g = read.csv("rotation1scripts_v4/original_data/expression.data/bygene/custom5aregion/trimmed.data.csv", stringsAsFactors = F)

colnames(g)[2:ncol(g)] = paste(colnames(g)[2:ncol(g)], "_", g[1, 2:ncol(g)], sep = "")
g = g[2:nrow(g), ]

#convert columns to numeric
for(i in 2:ncol(g)){
    g[[i]] = as.numeric(g[[i]])
}

#grab only spike tissue
g.spike = g[grep("spike", g[, 1]), ]

#isolate genes with tpm higher than 10
good.hits = unlist(lapply(g.spike, function(x){
    length(which(x > 10))
}))

g.spike2 = g.spike[, which(good.hits > 0)]

genenames = colnames(g.spike2)
genenames2 = genenames[grep("_tpm", genenames)]
genenames3 = multi.str.split(genenames2, "_", 1)

genomic.sequences = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/genomic.genes.w.names.fa")

genomic.gene.names = multi.str.split(names(genomic.sequences), "_", 2)

gene.candidates = genomic.sequences[match(genenames3, genomic.gene.names)]

writeXStringSet(gene.candidates, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/christobal.spike.candidates.cs.genomic.fa")

chris.blast = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/christobal.candidates.blast")

chris.blast2 = sort.blastdf(chris.blast)

chris.blast3 = parse.scaffold.blast(chris.blast2, 2000)

chris.blast4 = chris.blast3[sort(chris.blast3$query, index.return = T)$ix, ]

#get best groups of hits
c.b5 = lapply(unique(chris.blast4$query), function(x){
    
    g = filter(chris.blast4, query == x)
    # browser()
    g[which(g$avg.bitscore == max(g$avg.bitscore)), ]
    
})

chris.blast6 = combine.list.of.data.frames(c.b5)

load("rotation1scripts_v4/saved.objects/paragon.candidates")

count1 = make.counter()
alignments = Map(function(cs.seq, para.seq){
    g = pairwiseAlignment(cs.seq, para.seq)
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/christo.candidate.align_", genenames3[count1()], "_.txt"))
    return(g)
}, gene.candidates, paragon.candidates)

cs.gene.cds = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")
cds.names = multi.str.split(names(cs.gene.cds), "=", 2)

gene.candidate.names = multi.str.split(names(gene.candidates), "_", 2)


gene.candidate.cds = cs.gene.cds[match(gene.candidate.names, cds.names)]

gene.candidate.cds[grep("TraesCS5A01G531300|TraesCS5A01G530800", names(gene.candidate.cds))] =
    reverseComplement(gene.candidate.cds[grep("TraesCS5A01G531300|TraesCS5A01G530800", names(gene.candidate.cds))])


count2 = make.counter()
paragon.cds = Map(function(paragon.seq, cs.cds.seq){
    g = pairwiseAlignment(paragon.seq, cs.cds.seq)
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/para.cds.align.txt_", genenames3[count2()], "_.txt"))
    g2 = aligned(g)
    return(g2)
}, paragon.candidates, gene.candidate.cds)

paragon.cds2 = list.of.stringset.to.stringset(paragon.cds)

paragon.cds3 = lapply(paragon.cds2, function(x){
    g = strsplit(as.character(x), "")
    coord.to.rm = which(unlist(lapply(g, function(x1) x1 == "-" | x1 == "N")))
    if(length(coord.to.rm) > 0) x = x[-coord.to.rm]
    DNAStringSet(x)
})


paragon.cds4 = list.of.stringset.to.stringset(paragon.cds3)

count2.1 = make.counter()
count2 = make.counter()
Map(function(paragon.cds, gene.cds){
    # print(count2.1())
    g = pairwiseAlignment(translate(paragon.cds), translate(gene.cds))
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/chris.align.aa.para.cs_", genenames3[count2()], "_.txt"))
    return(g)
}, paragon.cds4, gene.candidate.cds)


to.keep = grep("TraesCS5A01G530100|TraesCS5A01G530800|TraesCS5A01G531300|TraesCS5A01G532300|TraesCS5A01G532400", names(gene.candidates))
refined.gene.candidates = gene.candidates[to.keep]
refined.cs.cds = gene.candidate.cds[to.keep]

paragon.gene.candidates.refined = paragon.candidates[to.keep]
refined.paragon.cds = paragon.cds4[to.keep]

writeXStringSet(refined.gene.candidates, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/christobal.candidates/new.candidates.refined.fa")

genegff = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", stringsAsFactors = F)

refined.names = multi.str.split(names(refined.gene.candidates), "_", 2)

gff.coords = unlist(lapply(refined.names, function(x){
    grep(x, genegff$V10)    
}))

as.numeric(genegff[gff.coords, 5])


grep(refined.names, genegff$V10)

#mix2seq sequencing of apogee
apogee.genes = readDNAStringSet("rotation1scripts_v4/original_data/sanger.sequencing/mix2seq2/11104687993-1.fasta")

mix2seqblast = read.blast("bioinf/blast/probe.vs.genes.blast/results.blast/mix2seq2.vs.genes.blast")

#filter the blast dataframe
g2 = newdf(colnames(mix2seqblast), no.rows = T)
for(i in unique(mix2seqblast$qseqid)){
    g = filter(mix2seqblast, qseqid == i)
    g = g[which(g$bitscore == max(g$bitscore)), ]
    g2 = rbind(g2, g)
}

#do some further filtering
g2 = g2[which(g2$bitscore > 800), ]

apogee.genes.refined = apogee.genes[match(as.character(g2$qseqid), names(apogee.genes))]
apoids = multi.str.split(names(apogee.genes.refined), "\\.", 2)

matches = unlist(lapply(apoids, function(x){
    grep(x, names(refined.gene.candidates))    
}))


c2 = make.counter()
lapply(refined.gene.candidates, function(x){
    g = pairwiseAlignment(apogee.genes[1], x)
    
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/test/test", c2(), ".txt"))
})

apogee.genes.refined[6] = reverseComplement(apogee.genes.refined[6])

alignments1 = list()
for(i in 1:length(apogee.genes.refined)){
    g = pairwiseAlignment(apogee.genes.refined[i], refined.gene.candidates[matches[i]])
    alignments1 = c(alignments1, g)
}

c1 = make.counter()
lapply(alignments1, function(x){
    writePairwiseAlignments(x, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/apo.align", c1(), ".txt"))
})


alignments.para = list()
for(i in 1:length(apogee.genes.refined)){
    g = pairwiseAlignment(apogee.genes.refined[i], paragon.gene.candidates.refined[matches[i]])
    alignments.para = c(alignments.para, g)
}

c1 = make.counter()
lapply(alignments.para, function(x){
    writePairwiseAlignments(x, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/apo.align.para", c1(), ".txt"))
})

alignments.cs.cds = list()
for(i in 1:length(apogee.genes.refined)){
    g = pairwiseAlignment(apogee.genes.refined[i], refined.cs.cds[matches[i]])
    alignments.cs.cds = c(alignments.cs.cds, g)
}

c1 = make.counter()
lapply(alignments.cs.cds, function(x){
    writePairwiseAlignments(x, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/apo.align.cs.cds", c1(), ".txt"))
})


alignments.para.cds = list()
for(i in 1:length(apogee.genes.refined)){
    g = pairwiseAlignment(apogee.genes.refined[i], refined.paragon.cds[matches[i]])
    alignments.para.cds = c(alignments.para.cds, g)
}

c1 = make.counter()
lapply(alignments.para.cds, function(x){
    writePairwiseAlignments(x, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/apo.align.para.cds", c1(), ".txt"))
})


#have now combined the apogee sequences for each respective gene.
apo.t.31300 = readDNAStringSet("rotation1scripts_v4/original_data/sanger.sequencing/mix2seq2/t.31300.combined.fasta")
apo.t.32300 = readDNAStringSet("rotation1scripts_v4/original_data/sanger.sequencing/mix2seq2/t.32300.all.fasta")
apo.t.32400 = readDNAStringSet("rotation1scripts_v4/original_data/sanger.sequencing/mix2seq2/T.32400.combined.fasta")

apo.comb.genes = c(apo.t.31300, apo.t.32300, apo.t.32400)
names(paragon.gene.candidates.refined)
names(refined.gene.candidates)

refined.gene.candidates2 = refined.gene.candidates[c(3, 4, 5)]
para.ref2 = paragon.gene.candidates.refined[c(3, 4, 5)]

c1 = make.counter()
Map(function(x, y){
    g = pairwiseAlignment(x, y)
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/apovcs.combined", c1(), ".txt"))
}, apo.comb.genes, refined.gene.candidates2)

c1 = make.counter()
Map(function(x, y){
    g = pairwiseAlignment(x, y)
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/apovpara.combined", c1(), ".txt"))
}, apo.comb.genes, para.ref2)

c1 = make.counter()
Map(function(x, y){
    g = pairwiseAlignment(x, y)
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/csvpara.combined", c1(), ".txt"))
}, refined.gene.candidates2, para.ref2)

refined.cs.cds2 = refined.cs.cds[c(3, 4, 5)]
refined.para.cds2 = refined.paragon.cds[c(3, 4, 5)]

mult.align32400 = c(apo.comb.genes[3], refined.gene.candidates2[3], para.ref2[3], refined.cs.cds2[3], refined.para.cds2[3])
writeXStringSet(mult.align32400, "rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/sequences/32400.mult.align.fasta")

#30800 alignments

cs30800 = refined.gene.candidates[2]
p30800 = paragon.gene.candidates.refined[2]
a30800 = readDNAStringSet("rotation1scripts_v4/original_data/sanger.sequencing/mix2seq3/t.30800.combined.fa")

a30800 = reverseComplement(a30800)

p1 = pairwiseAlignment(cs30800, p30800)
writePairwiseAlignments(p1, "rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/csvsp30800.txt")

p1 = pairwiseAlignment(a30800, cs30800)
writePairwiseAlignments(p1, "rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/avscs30800.txt")

p1 = pairwiseAlignment(a30800, p30800)
writePairwiseAlignments(p1, "rotation1scripts_v4/processed_data/TILLING/christo.candidates.alignments/avsp30800.txt")

