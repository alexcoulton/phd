#want to match contig ids from hap.blast to marker IDs

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

#pass mysql query to database
rs = dbSendQuery(mydb, "SELECT * FROM axiom820")
rs2 = fetch(rs, n=-1)

rs2$id.1 = multi.str.split(rs2$id, "_", 1)

rs3 = rs2[match(hap.blast6$hap.id, rs2$id.1), ]

blast820 = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/820kprobes.vs.iwgsc.4b.rev.blast")

blast820.1 = blast820[which(blast820$qseqid %in% rs3$affycode), ]

s(blast820.1, "rotation1scripts_v4/saved.objects/blast820.1", "sacha.markers.haplotype.processing.R")

blast820.2 = remove.hits.with.same.best.e.value(blast820.1)

#experimenting with code profiling
source("rotation1scripts_v4/scripts/r/generic.functions/rm.best.e.value.R")
blast820.3 = blast820.1[1:10, ]
l = lineprof(remove.hits.with.same.best.e.value(blast820.3))
shine(l)

ptm = proc.time()
remove.hits.with.same.best.e.value(blast820.3)
proc.time() - ptm

bt.dt = as.data.table(blast820)
setkey(bt.dt, qseqid)
bt.dt[.("AX-94381226"), nomatch = 0L]



microbenchmark(bt.dt[.("AX-94381226"), nomatch = 0L])

microbenchmark(filter(bt.dt, qseqid == "AX-94381226"))

microbenchmark(bt.dt["AX-94381226", ])

microbenchmark(filter(blast820, qseqid == "AX-94381226"))


blast820.test = add_column(blast820, index = 1:nrow(blast820), .before = "qseqid")

microbenchmark(filter(blast820.test, index == 488))




blast820.unique.qseq = blast820[match(unique(blast820$qsedid), blast820$qseqid), ]
blast820.unique.qseq = as.data.table(blast820.unique.qseq)
setkey(blast820.unique.qseq, qseqid)
microbenchmark(blast820.unique.qseq[.("AX-94381226"), ])

blast820.test2 = blast820.test[1:10, ]
blast820.test.dt = as.data.table(blast820.test2)
setkey(blast820.test.dt, qseqid)
microbenchmark(blast820.test.dt[.("AX-94381226"), ])




g = data.frame(1:1000000, sort(rep(seq(100001, 200000, 1), 10)), sort(rep(seq(100001, 200000, 1), 10)))
g = as.data.table(g)
colnames(g) = c("V1", "V2", "V3")
# setkey(g, V1)
setkey(g, V2)
microbenchmark(g[.(100000), ])
microbenchmark(filter(g, V1 == 10))

