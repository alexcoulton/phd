
#     ____________________________________________________________________________
#     CHECK MISSING VARS FROM CEREALSDB                                                                             ####


library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

q1 = dbSendQuery(mydb, "SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='cerealsdb' AND `TABLE_NAME`='axiom820dataab';")
q2 = fetch(q1, n=-1)

list.of.w.relatives = as.character(unique(introgression$wild.relative))
list.of.elites = as.character(unique(introgression$elite))
w.rel2 = list.of.w.relatives[-which(list.of.w.relatives %in% q2[, 1])]
eli2 = list.of.elites[-which(list.of.elites %in% q2[, 1])]

fg = function(x) grep(x, q2[, 1])

write(w.rel2, "rotation1scripts_v4/original_data/sacha.introgression.data/listof.w.relatives.needed.txt")
write(eli2, "rotation1scripts_v4/original_data/sacha.introgression.data/list.of.missing.hexa.txt")


#     ____________________________________________________________________________
#     PROCESS PROBES                                                                                                                    ####


mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

sql.query = "SELECT affycode, Sequence FROM axiom820"

#pass mysql query to database
rs = dbSendQuery(mydb, sql.query)
rs2 = fetch(rs, n=-1)


non.match.coords = grep("\\[.*\\]", rs2[, 2], invert = T)
rs3 = rs2[non.match.coords, ]
#grab anything that isn't A, G, C or T
unique(unlist(regmatches(rs3[, 2], gregexpr("[^AGCT]", rs3[, 2]))))

rs2 = rs2[-non.match.coords, ]

uni.seq = unlist(unique(regmatches(rs2[, 2], gregexpr("\\[.*\\]", rs2[, 2]))))

rs2[, 2] = gsub("\\[R\\]", "R", rs2[, 2])
rs2[, 2] = gsub("\\[Y\\]", "Y", rs2[, 2])
rs2[, 2] = gsub("\\[M\\]", "M", rs2[, 2])
rs2[, 2] = gsub("\\[K\\]", "K", rs2[, 2])
rs2[, 2] = gsub("\\[S\\]", "S", rs2[, 2])
rs2[, 2] = gsub("\\[W\\]", "W", rs2[, 2])

check.worked = grep("\\[.*\\]", rs2[, 2])

rs4 = rbind(rs2, rs3)

writefasta(rs4, "rotation1scripts_v4/original_data/fasta/820kprobes.fa")


#     ____________________________________________________________________________
#     BLAST PROCESSING                                                                                                                 ####

blast820 = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/820kprobes.vs.iwgsc.4b.rev.blast")

#splitting the dataset into multiple parts
list.of.coords.to.split.by = as.numeric()
counter = 32357 #divide 820 dataset into ~ 50 sets
for(i in 1:55){
    if(counter > 1617840){ #THIS NUMBER SHOULD BE 1 LESS THAN THE NUMBER OF ROWS IN blast820
        
    } else {
        while(blast820[(counter+1), 1] == blast820[counter, 1]){
            counter = counter + 1    
        }
        if(counter > 1617840){
            list.of.coords.to.split.by = c(list.of.coords.to.split.by, nrow(blast820))
        } else {
            list.of.coords.to.split.by = c(list.of.coords.to.split.by, counter)    
        }
        counter = counter + 32357
    }
}
list.of.coords.to.split.by = c(list.of.coords.to.split.by, nrow(blast820))
start.ranges = c(1, (list.of.coords.to.split.by + 1)[1:(length(list.of.coords.to.split.by)-1)])
end.ranges = list.of.coords.to.split.by

blast.subsets = Map(function(s.ranges, e.ranges){
    blast820[s.ranges:e.ranges, ]
}, start.ranges, end.ranges)

count2 = make.counter()
lapply(blast.subsets, function(x){
    write.csv(x, p("rotation1scripts_v4/processed_data/BLAST/820k/820.blast.subsets.4b.rev/subset.", count2(), ".csv"), row.names = F)
})

#     ____________________________________________________________________________
#     PROCESS UNIQUE                                                                                                                    ####


mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

sql.query3 = "SELECT * FROM axiom820"

#pass mysql query to database
consensus = dbSendQuery(mydb, sql.query3)
consensus2 = fetch(consensus, n=-1)
consensus3 = consensus2[-which(consensus2$consensus == ""), ]
consensus4 = consensus3[-which(consensus3$consensus == "none"), ]

blast820.unique = read.csv("rotation1scripts_v4/processed_data/BLAST/820k/only.hits.w.unique.best.e.value.4b.rev/all.comb.4b.rev.csv", header = T, stringsAsFactors = F)

unuseful.cons = unique(unlist(regmatches(consensus4$consensus, gregexpr("[^DBA]", consensus4$consensus))))
consensus5 = consensus4[-which(consensus4$consensus %in% unuseful.cons), ]
consensus5$con.trunc = substr(consensus5$consensus, 1, 2)


blast820.non.unique = blast820[which(!blast820$qseqid %in% blast820.unique$qseqid), ]
blast820.non.unique.w.consensus = blast820.non.unique[which(blast820.non.unique$qseqid %in% consensus5$affycode), ]

blast820.non.unique.w.consensus$cons.chromo = ""

blast820.non.unique.w.consensus$cons.chromo = consensus5$con.trunc[match(blast820.non.unique.w.consensus$qseqid, consensus5$affycode)]
b.n.uniq.w.2 = blast820.non.unique.w.consensus[which(substr(blast820.non.unique.w.consensus$sseqid, 4, 5) == blast820.non.unique.w.consensus$cons.chromo), ]
more.to.del = names(which(table(as.character(b.n.uniq.w.2$qseqid)) > 1))

b.n.uniq.w.3 = b.n.uniq.w.2[-which(as.character(b.n.uniq.w.2$qseqid) %in% more.to.del), ]
b.n.4 = b.n.uniq.w.3[, -13]

final.blast = rbind(blast820.unique, b.n.4)
final.blast2 = split.df.into.list.of.dfs.by.column(final.blast, "sseqid")
final.blast2.1 = lapply(final.blast2, function(x) arrange(x, sstart))
final.blast3 = combine.list.of.data.frames(final.blast2.1)

# write.csv(final.blast3, "rotation1scripts_v4/processed_data/BLAST/820k/820k.full.unique.blast.ordered.blast", row.names = F)
write.csv(final.blast3, "rotation1scripts_v4/processed_data/BLAST/820k/820k.full.unique.blast.4brev.ordered.csv", row.names = F)

probes820fa = ReadFasta("rotation1scripts_v4/original_data/fasta/820kprobes.fa")

tauachii.blast = read.blast("bioinf/blast/probe.vs.tauschii.blast/results.blast/820kprobes.d.genome.vs.tauschii.blast")

d.genome.blast = final.blast3[grep("D", final.blast3$sseqid), ] 

a.genome.blast = final.blast3[grep("A", final.blast3$sseqid), ] 

b.genome.blast = final.blast3[grep("B", final.blast3$sseqid), ]

urartu.blast = read.blast("bioinf/blast/urartu.blast/results.blast/820kprobes.vs.urartu4rev.blast")

urartu.blast2 = urartu.blast[match(a.genome.blast$qseqid, urartu.blast$qseqid), ]
uru.b3 = split.df.into.list.of.dfs.by.column(urartu.blast2, "sseqid")
u.b4 = lapply(uru.b3, function(x) arrange(x, sstart))
u.b5 = combine.list.of.data.frames(u.b4)
m.coord2 = match(a.genome.blast$qseqid, u.b5$qseqid)

ur.df1 = data.frame(1:length(m.coord2), m.coord2)
colnames(ur.df1) = c("n", "matches")

ur.df1.2 = na.omit(ur.df1)
qual.probes.uru = u.b5$qseqid[ur.df1.2$matches[longest_subseq.R(ur.df1.2$matches)]]


ur.df1[which(is.na(ur.df1$matches)), 2] = -100000
t.df[which(is.na(t.df$matches)), 2] = -100000


barley.blast.1 = read.blast("bioinf/blast/barley.blast/results.blast/820kvsbarley.blast")
# barley.blast.1 = barley.blast[-which(barley.blast$sseqid == "chrUn"), ]
make.barley.order.df = function(wheat.genome.blast.queries){
    barley.blast2 = barley.blast.1[match(wheat.genome.blast.queries, barley.blast.1$qseqid), ]
    barley.3 = split.df.into.list.of.dfs.by.column(barley.blast2, "sseqid")
    bar.4 = lapply(barley.3, function(x) arrange(x, sstart))
    bar.5 = combine.list.of.data.frames(bar.4)
    bar.coord2 = match(wheat.genome.blast.queries, bar.5$qseqid)
    
    bar.df1 = data.frame(1:length(bar.coord2), bar.coord2)
    colnames(bar.df1) = c("n", "matches")
    
    bar.df.return = na.omit(bar.df1)
    
    list(bar.df.return, bar.5)
}

bar.df.a = make.barley.order.df(a.genome.blast$qseqid)
bar.df.b = make.barley.order.df(b.genome.blast$qseqid)
bar.df.d = make.barley.order.df(d.genome.blast$qseqid)

#     ____________________________________________________________________________
#     GET QUALIFYING PROBES FROM LONGEST INCR. SUBSEQ.                                                ####

bar.df1.2$dist = unlist(Map(function(x, y){
    point.dist(x, y, 0, 0.48)
}, bar.df1.2$n, bar.df1.2$matches))
bar.df1.2$dist = as.numeric(bar.df1.2$dist)

bar.df1.2$group = ""
bar.df1.2$group[which(bar.df1.2$dist < 3000)] = "< 100"
bar.df1.2$group[which(bar.df1.2$dist > 3000)] = "> 100"

bar.df1.2$group.subseq = ""
bar.df1.2$group.subseq[longest_subseq.R(bar.df1.2$matches)] = "T"

ur.df1.2$group.subseq = ""
ur.df1.2$group.subseq[longest_subseq.R(ur.df1.2$matches)] = "T"

t.df.1$group.subseq = ""
t.df.1$group.subseq[longest_subseq.R(t.df.1$matches)] = "T"

qualifying.probes.bar.d = bar.df.d[[2]]$qseqid[bar.df.d[[1]]$matches[longest_subseq.R(bar.df.d[[1]]$matches)]]
qualifying.probes.bar.b = bar.df.b[[2]]$qseqid[bar.df.b[[1]]$matches[longest_subseq.R(bar.df.b[[1]]$matches)]]
qualifying.probes.bar.a = bar.df.a[[2]]$qseqid[bar.df.a[[1]]$matches[longest_subseq.R(bar.df.a[[1]]$matches)]]

qual.probes.uru = u.b5$qseqid[ur.df1.2$matches[longest_subseq.R(ur.df1.2$matches)]]

qual.probes.tau = t.bla4$qseqid[t.df.1$matches[longest_subseq.R(t.df.1$matches)]]

qual.probes.all.d = qual.probes.tau[which(qual.probes.tau %in% qualifying.probes.bar.d)]
qual.probes.all.a = qual.probes.uru[which(qual.probes.uru %in% qualifying.probes.bar.a)]
qual.probes.all.b = qualifying.probes.bar.b

all.qual.probes = c(as.character(qual.probes.all.a), as.character(qual.probes.all.b), as.character(qual.probes.all.d))

cons.order.blast = final.blast3[match(all.qual.probes, final.blast3$qseqid), ]
# write.csv(cons.order.blast, "rotation1scripts_v4/processed_data/BLAST/820k/820k.cons.order.w.barley.csv", row.names = F)


#qualifying.probes.bar.b = bar.5$qseqid[bar.df1.2$matches[which(bar.df1.2$dist < 2000)]]


library(gridExtra)

plot0 = ggplot(bar.df1.2, aes(n, matches, color = group.subseq)) + geom_point(alpha = 1/20) + 
    xlab("IWGSC Assembly Position") + ylab("H. vulgare Position") + theme_bw() +
    theme(axis.title = element_text(size = 15), legend.position = "none", plot.title = element_text(size = 15)) + ggtitle("(a)") +
    coord_cartesian(ylim = c(0, 101085), xlim = c(0, 203852)) + scale_color_manual(values = c("#000000", "#eef442"))

plot1 = ggplot(ur.df1.2, aes(n, matches, color = group.subseq)) + geom_point(alpha = 1/20) + 
    xlab("IWGSC Assembly Position") + ylab("T. urartu Position") + theme_bw() +
    theme(axis.title = element_text(size = 15), legend.position = "none", plot.title = element_text(size = 15)) + ggtitle("(b)") +
    coord_cartesian(ylim = c(0, 171310), xlim = c(0, 203852)) + scale_color_manual(values = c("#000000", "#eef442"))

plot2 = ggplot(t.df.1, aes(n, matches, color = group.subseq)) + geom_point(alpha = 1/20) + 
    xlab("IWGSC Assembly Position") + ylab("A. tauschii Position") + theme_bw() +
    theme(axis.title = element_text(size = 15), legend.position = "none", plot.title = element_text(size = 15)) + ggtitle("(c)") +
    coord_cartesian(ylim = c(0, 233759)) + scale_color_manual(values = c("#000000", "#eef442"))


# pdf("rotation1scripts_v4/plots/iwgsc.vs.a.tauschii.probe.blast.scatterplot/urartu.and.tauschii.pdf", 10, 10)
# png("rotation1scripts_v4/plots/iwgsc.vs.a.tauschii.probe.blast.scatterplot/urartu.and.tauschii2.png", 1800, 900)
png("rotation1scripts_v4/plots/iwgsc.vs.a.tauschii.probe.blast.scatterplot/urartu.and.tauschiiv5.png", 1000, 900)
grid.arrange(plot0, plot1, plot2, ncol = 2)
dev.off()


tauschii.blast2 = split.df.into.list.of.dfs.by.column(tauachii.blast, "sseqid")
t.bla3 = lapply(tauschii.blast2, function(x) arrange(x, sstart))
t.bla4 = combine.list.of.data.frames(t.bla3)
t.matches = match(d.genome.blast$qseqid, t.bla4$qseqid)
t.df = data.frame(1:length(t.matches), t.matches)
colnames(t.df) = c("n", "matches")

t.df.1 = na.omit(t.df)

qual.probes.tau = t.bla4$qseqid[t.df.1$matches[longest_subseq.R(t.df.1$matches)]]

png("rotation1scripts_v4/plots/iwgsc.vs.a.tauschii.probe.blast.scatterplot/scat.png", 400, 400)

dev.off()

plot(1:length(t.matches), t.matches, alpha = 0.2)

