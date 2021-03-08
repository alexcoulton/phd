#grab and process haplotype sequences for Sacha's segregation distortion markers

source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/generic.functions/fasta.R")


#     ____________________________________________________________________________
#     GRAB AFFYCODES FOR ALL HAPLOTYPES                                                                             ####

haplotypes = ReadFasta("rotation1scripts_v4/processed_data/fasta/extended.haplotypes.fa")


haplotypes$contigid = multi.str.split(as.character(haplotypes$name), "_", 1)


#setup connection to mysql database
mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

#pass mysql query to database
rs = dbSendQuery(mydb, "SELECT * FROM axiom820")
rs2 = fetch(rs, n=-1)

rs2$contigid = multi.str.split(as.character(rs2$id), "_", 1)
haplotypes$affycode = rs2$affycode[match(haplotypes$contigid, rs2$contigid)]

write.csv(haplotypes, "rotation1scripts_v4/processed_data/haplotypes/haplotypes.w.affycodes.csv", row.names = F)


#     ____________________________________________________________________________
#     GRAB HAPLOTYPES FOR SACHA SEG. DIST. MARKERS                                                        ####

load("rotation1scripts_v4/saved.objects/sacha.seg.dist.markers")

rownames(sacha.seg.dist.markers) = gsub("\\.", "-", rownames(sacha.seg.dist.markers))

haplotypes.seg.dist = na.omit(haplotypes[match(rownames(sacha.seg.dist.markers), haplotypes$affycode), ])

haplotypes.seg.dist.fasta = haplotypes.seg.dist[, c(4, 14)]

writefasta(haplotypes.seg.dist.fasta, "rotation1scripts_v4/processed_data/fasta/sacha.seg.dist.haplotypes.short.fa")

allprobe.sequences = ReadFasta("rotation1scripts_v4/original_data/fasta/35kprobes.fa")

sacha.seg.dist.probe.sequences = allprobe.sequences[match(gsub("-", ".", rownames(sacha.seg.dist.markers)), allprobe.sequences$name), ]

writefasta(sacha.seg.dist.probe.sequences, "rotation1scripts_v4/processed_data/fasta/sacha.seg.dist.probe.sequences.fa")


#     ____________________________________________________________________________
#     ANALYSE BLAST RESULTS                                                                                                     ####

probe.blast = read.table("bioinf/blast/probe.vs.genome.blast/results.blast/sacha.seg.dist.probes.blast")
haplotype.blast = read.table("bioinf/blast/probe.vs.genome.blast/results.blast/sacha.seg.dist.haplotypes.blast")

probe.blast2 = clean.blast(probe.blast)


haplotype.blast = remove.hits.with.same.best.e.value(haplotype.blast)

dif1 = get.list.of.differences.between.best.and.second.best.e.values(haplotype.blast)

names(which(dif1 >= 2.141728e-64))

#sorting of probe.blast df 

probe.blast2 = sort.blastdf(probe.blast2)



#     ____________________________________________________________________________
#     FURTHER HAPLOTYPE BLAST ANALYSIS                                                                                ####

#want to compare the best blast hits for haplotypes with the best blast hits for probe sequences. how often are these
#concordant?

hap.blast = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/extended.haplotypes.vs.genome.blast")
hap.blast$hap.id = multi.str.split(as.character(hap.blast$qseqid), character.split = "_", split.id = 1)
hap.blast2 = split.df.into.list.of.dfs.by.column(hap.blast, "hap.id")

count1 = make.counter()
hap.blast3 = lapply(hap.blast2, function(x){
    print(count1())
    min.e = min(x$evalue)
    x[which(x$evalue == min.e), ]
})

g = lapply(hap.blast3, function(x){
    nunique(x$sseqid)
})

#grab haplotype blast hits where best hits are concordant on the same chromosome sequence
hap.blast4 = hap.blast3[which(g == 1)]

hap.blast5 = lapply(hap.blast4, function(x){
    x[1, ]
})

hap.blast6 = do.call(rbind, hap.blast5)

#want to match contig ids from hap.blast to marker IDs

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

#pass mysql query to database
rs = dbSendQuery(mydb, "SELECT * FROM axiom820")
rs2 = fetch(rs, n=-1)

rs2$id.1 = multi.str.split(rs2$id, "_", 1)

rs3 = rs2[match(hap.blast6$hap.id, rs2$id.1), ]
hap.blast6$affycode = rs3$affycode

blast820 = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/820kprobes.vs.iwgsc.4b.rev.blast")

blast820.1 = blast820[which(blast820$qseqid %in% rs3$affycode), ]

s(blast820.1, "rotation1scripts_v4/saved.objects/blast820.1", "sacha.markers.haplotype.processing.R")

#grabbed the best hits on a Wilkins session, here is the object with the results
load("rotation1scripts_v4/saved.objects/blast820.5")

blast820.5 = blast820.5[match(hap.blast6$affycode, blast820.5$qseqid), ]

#percentage of concordant hits between haplotype blast results and probe blast results
table(hap.blast6$sseqid == blast820.5$sseqid)[["TRUE"]] / nrow(blast820.5) * 100

