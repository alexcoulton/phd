library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')
apogee = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Apogee'")
apogee2 = fetch(apogee, n=-1)

chinesespring = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Chinese_Spring'")
chinesespring2 = fetch(chinesespring, n=-1)

paragon = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE 'Paragon'")
paragon2 = fetch(paragon, n=-1)

(length(which(apogee2$Matrix_value == chinesespring2$Matrix_value)) / nrow(chinesespring2)) * 100
(length(which(apogee2$Matrix_value == paragon2$Matrix_value)) / nrow(chinesespring2)) * 100
(length(which(chinesespring2$Matrix_value == paragon2$Matrix_value)) / nrow(chinesespring2)) * 100

cs.x.p5a = all.m.8.iwgsc.4b.rev[[4]][["5A"]]
marker5a = cs.x.p5a$marker

a5a = apogee2[match(marker5a, apogee2$Probe_row), ]
cs5a = chinesespring2[match(marker5a, chinesespring2$Probe_row), ]
p5a = paragon2[match(marker5a, paragon2$Probe_row), ]

length(which(a5a$Matrix_value == cs5a$Matrix_value))
length(which(a5a$Matrix_value == p5a$Matrix_value))
length(which(cs5a$Matrix_value == p5a$Matrix_value))

#can't use 5a SNPs from the genetic map of P x CS as these have been specifically selected to be polymorphic between the two varieties. 

snp.blast = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/allprobegenome.blast")

snp.blast2 = grab.best.hits(snp.blast)

snp.blast3 = grab.best.hits(snp.blast, T)





