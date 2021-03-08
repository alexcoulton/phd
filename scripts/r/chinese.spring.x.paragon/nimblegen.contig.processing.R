#     ____________________________________________________________________________
#     INTERFACE WITH MYSQL DATABASE                                                                                     ####

setwd("/Users/ac14037/project.phd.main/")

library(RMySQL)

#setup connection to mysql database
mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

#make mysql query for passing to database
seg.markers.sql = paste("affycode = '", seg.dist.markers, "'", " OR ", sep = "", collapse = "") #note seg.dist.markers is from seg.dist.marker.analysis.R
seg.markers.sql = substr(seg.markers.sql, 1, nchar(seg.markers.sql)-4)
seg.markers.sql = paste("SELECT * FROM axiom820 WHERE ", seg.markers.sql, sep = "")
seg.markers.sql = gsub("\\.", "-", seg.markers.sql)

#pass mysql query to database
rs = dbSendQuery(mydb, seg.markers.sql)
rs2 = fetch(rs, n=-1)
rs2 = cereals.axiom.table


#     ____________________________________________________________________________
#     GRAB HAPLOTYPE SEQUENCES FOR SEG. DIST. MARKERS                                                 ####

haplotype.seqs = ReadFasta("rotation1scripts_v4/original_data/fasta/170317_haplotype_sequences.fa")

q = lapply(cereals.axiom.table$id, function(x){
    g = strsplit(x, "_")[[1]]
    grep(g, haplotype.seqs$name)
})

length(which(unlist(lapply(q, function(x) length(x) == 0)))) #number of probes w/ no representatives in haplotype_sequences.fa

#grab only the first haplotype for each probe
q.coord = lapply(q, function(x) x[1])
#add to original dataframe
cereals.axiom.table$haplotype.seq = haplotype.seqs$sequence[unlist(q.coord)]

g2 = haplotype.seqs[unlist(q.coord), ]
g3 = g2[!is.na(g2$name), ]

nimblegen.contigs = ReadFasta("rotation1scripts_v4/original_data/fasta/NimblegenRef.fa")

q2 = lapply(cereals.axiom.table$id, function(x){
    g = strsplit(x, "_")[[1]]
    grep(g, nimblegen.contigs$name)
})

cereals.axiom.table$nimblecontig = nimblegen.contigs[unlist(q2), 2]

save(cereals.axiom.table, file = "rotation1scripts_v4/saved.objects/cereals.axiom.table")


#     ____________________________________________________________________________
#     ALIGN HAPLOTYPES TO NIMBLEGEN CONTIGS                                                                     ####

counter = 0

extend.haplotype = function(haplotype.sequence, nimblegen.sequence){
    g = pairwiseAlignment(DNAString(haplotype.sequence), DNAString(nimblegen.sequence))
    pos.mis = mismatchTable(g)
    nimblecontig = as.character(nimblegen.sequence)
    
    #split string up into vector then replace vector positions according to pos.mis
    nimblecontig = strsplit(nimblecontig, "")
    nimblecontig = nimblecontig[[1]]
    #print(nrow(pos.mis))
    
    if(nrow(pos.mis) != 0){
        for(i in 1:nrow(pos.mis)){
            nimblecontig[pos.mis$SubjectStart[[i]]] = as.character(pos.mis$PatternSubstring[[i]])
        }
    }
    
    nimblecontig = paste(nimblecontig, collapse = "")
    
    counter <<- counter + 1
    print(counter)
    
    return(nimblecontig)
}

cereals.axiom.table$extended.haplotype = ""

cereals.table.only.hap = cereals.axiom.table[!is.na(cereals.axiom.table$haplotype.seq), ]

cereals.table.only.hap$extended.haplotype = Map(extend.haplotype, cereals.table.only.hap$haplotype.seq, cereals.table.only.hap$nimblecontig)


for(i in 1:nrow(cereals.table.only.hap)){
    g = pairwiseAlignment(cereals.table.only.hap$nimblecontig[[i]], cereals.table.only.hap$extended.haplotype[[i]])    
    writePairwiseAlignments(g, file = p("rotation1scripts_v4/processed_data/pairwise.alignments/test", i, ".txt"))
}



#     ____________________________________________________________________________
#     MASKING WITH Ns                                                                                                                 ####

#grab entire axiom820 table
all.axiom820 = dbSendQuery(mydb, "SELECT * FROM axiom820")
all.axiom820.2 = fetch(all.axiom820, n=-1)

write.csv(all.axiom820.2, "rotation1scripts_v4/original_data/cerealsdb.data/all.axiom820.csv", row.names = F)

all.axiom820.2$contigid = unlist(lapply(all.axiom820.2$id, function(x) strsplit(x, "_")[[1]][1]))

haplotype.seqs = convert.to.character.data.frame(haplotype.seqs)
haplotype.seqs$contigid = unlist(lapply(haplotype.seqs$name, function(x) strsplit(x, "_")[[1]][1]))

haplotype.seqs$nimble.seq = nimblegen.contigs[match(haplotype.seqs$contigid, nimblegen.contigs$name), 2]

haplotype.seqs$contigid[[1]]



multi.str.split = function(x, character.split, split.id){
    #An implementation of strsplit() that allows selection of a particular element of the output value of strsplit(), 
    #for all elements of the input vector to which the split is applied.
    #--------------
    #args:
    #x = character vector
    #character.split = character specifying what character to split with 
    #split.id = integer specifying which element of the split you want
    
    unlist(lapply(x, function(q) strsplit(q, character.split)[[1]][split.id]))
}




all.axiom820.2$snp.location = multi.str.split(all.axiom820.2$id, "_", 2)

haplotype.seqs$snps = multi.str.split(haplotype.seqs$name, "_", 2)

nrow(haplotype.seqs) == length(grep("[[:upper:]]\\d*[[:upper:]]\\d*[[:upper:]]\\d*", haplotype.seqs$snps)) #check whether all of the haplotypes contain three snps

haplotype.seqs$snp.pos.1 = multi.str.split(haplotype.seqs$snps, "[[:upper:]]", 2)
haplotype.seqs$snp.pos.2 = multi.str.split(haplotype.seqs$snps, "[[:upper:]]", 3)
haplotype.seqs$snp.pos.3 = multi.str.split(haplotype.seqs$snps, "[[:upper:]]", 4)

haplotype.seqs$snp.1 = multi.str.split(haplotype.seqs$snps, "\\d*", 1)
haplotype.seqs$snp.2 = multi.str.split(haplotype.seqs$snps, "\\d*", 3)
haplotype.seqs$snp.3 = multi.str.split(haplotype.seqs$snps, "\\d*", 5)

haplotype.seqs$extendedhaplotype = ""


#loop for producing extended haplotypes
for(i in 16199:nrow(haplotype.seqs)){
    q = strsplit(as.character(haplotype.seqs$nimble.seq[[i]]), "")[[1]]
    
    q[as.numeric(haplotype.seqs$snp.pos.1[[i]])] = haplotype.seqs$snp.1[[i]]
    
    q[as.numeric(haplotype.seqs$snp.pos.2[[i]])] = haplotype.seqs$snp.2[[i]]
    
    q[as.numeric(haplotype.seqs$snp.pos.3[[i]])] = haplotype.seqs$snp.3[[i]]
    
    
    g = all.axiom820.2[which(all.axiom820.2$contigid == haplotype.seqs$contigid[[i]]), ] #grab rows from all.axiom820.2 that match the haplotype
    pos.to.mask = as.numeric(g$snp.location[which(!as.numeric(g$snp.location) %in% as.numeric(haplotype.seqs[i, 6:8]))]) #grab the SNP positions that are not in the haplotype
    
    q[pos.to.mask] = "N"
    
    
    
    q1 = paste(q, collapse = "")
    
    haplotype.seqs$extendedhaplotype[[i]] = q1
    print(i)
}
#done to 16200

#check if above loop is working (it seems to be!)
for(i in 1:30){
    g = pairwiseAlignment(haplotype.seqs$nimble.seq[[i]], haplotype.seqs$extendedhaplotype[[i]])
    writePairwiseAlignments(g, p("rotation1scripts_v4/processed_data/pairwise.alignments/newtest", i, ".txt")) 
}

fa.ext.haplotypes = haplotype.seqs[, c(1, 12)]

writefasta(fa.ext.haplotypes, "rotation1scripts_v4/processed_data/fasta/extended.haplotypes.fa")




#write alignments for inspection
for(i in 1:nrow(haplotype.seqs)){
    g = pairwiseAlignment(haplotype.seqs$sequence[[i]], haplotype.seqs$nimble.seq[[i]])    
    writePairwiseAlignments(g, file = p("rotation1scripts_v4/processed_data/pairwise.alignments/test2", i, ".txt"))
}


#cut extended haplotypes down

haplotype.seqs.seg.dist.genes.only = haplotype.seqs[match(cereals.table.only.hap$contigid, haplotype.seqs$contigid), ]


haplotype.seqs.seg.dist.genes.only$shorter.ext.haplotype = ""

for(i in 1:nrow(haplotype.seqs.seg.dist.genes.only)){
    if((as.numeric(haplotype.seqs.seg.dist.genes.only$snp.pos.1[[i]]) - 50) < 0){ 
        sub.start = 0
    } else { 
        sub.start = as.numeric(haplotype.seqs.seg.dist.genes.only$snp.pos.1[[i]]) - 50
    }
    
    if((as.numeric(haplotype.seqs.seg.dist.genes.only$snp.pos.3[[i]]) + 50) > nchar(haplotype.seqs.seg.dist.genes.only$extendedhaplotype[[i]])){
        sub.end = nchar(haplotype.seqs.seg.dist.genes.only$extendedhaplotype[[i]])
    } else { 
        sub.end = as.numeric(haplotype.seqs.seg.dist.genes.only$snp.pos.3[[i]]) + 50
    }
    
    
    g = substr(haplotype.seqs.seg.dist.genes.only$extendedhaplotype[[i]], sub.start, sub.end)
    haplotype.seqs.seg.dist.genes.only$shorter.ext.haplotype[[i]] = g
    print(i)
}


writefasta(haplotype.seqs.seg.dist.genes.only[, c(3, 13)], "rotation1scripts_v4/processed_data/fasta/shorter.extended.haplotypes.fa")


#     ____________________________________________________________________________
#     ADDING EXTENDED HAPLOTYPES TO SEG. DIST. GENES                                                    ####

cereals.table.only.hap$contigid = multi.str.split(cereals.table.only.hap$id, "_", 1)

cereals.table.only.hap$extended.haplotype = haplotype.seqs$extendedhaplotype[match(cereals.table.only.hap$contigid, haplotype.seqs$contigid)]

only.hap.fa = cereals.table.only.hap[, c(1, 28)]

writefasta(only.hap.fa, "rotation1scripts_v4/processed_data/fasta/cs.x.p.seg.dist.genes.extended.haplotypes.fa")



#     ____________________________________________________________________________
#     ANALYSE BLAST                                                                                                                     ####

hap.blast = read.table("bioinf/blast/probe.vs.genome.blast/results.blast/cs.x.para.ext.haplotypes.vs.genome.blast")
hap.blast2 = clean.blast(hap.blast, 5.0e-510)

hap.shorter.blast = read.table("bioinf/blast/probe.vs.genome.blast/results.blast/shorter.ext.haplotypes.vs.genome.blast")
#hap.shorter.blast2 = clean.blast(hap.shorter.blast, 2.0e-69)



hap.shorter.blast2.2 = remove.hits.with.same.best.e.value(hap.shorter.blast)



differences = get.list.of.differences.between.best.and.second.best.e.values(hap.shorter.blast2.2)
hits.to.remove = names(differences[which(differences < 1.679223e-70)])

hap.shorter.blast2.2 = hap.shorter.blast2.2[-which(hap.shorter.blast2.2$V1 %in% hits.to.remove), ]

hap.2.3 = grab.best.hits(hap.shorter.blast2.2)

hap.2.3$contigid = cereals.table.only.hap[match(hap.2.3$V1, cereals.table.only.hap$contigid), ]$affycode

load(file = "rotation1scripts_v4/saved.objects/newblast")

#check for concordance between probe BLAST and haplotype BLAST 

probes.also.in.hap = newblast[which(newblast$V1 %in% hap.2.3$contigid), ]

hap.corresponding.to.unique.probes = hap.2.3[match(probes.also.in.hap$V1, hap.2.3$contigid), ]

probes.also.in.hap = cbind(probes.also.in.hap, hap.corresponding.to.unique.probes[, 9:10])

non.concordant = probes.also.in.hap[which(abs(as.numeric(probes.also.in.hap$V9) - as.numeric(probes.also.in.hap[, 13])) > 1000), ] #grab the haplotype BLAST results that are not concordant

hap.2.3 = hap.2.3[-which(hap.2.3$contigid %in% non.concordant$V1), ] #remove non-concordant hits between haplotype sequences and probe sequences from the haplotype blast

save(hap.2.3, file = "rotation1scripts_v4/saved.objects/hap.2.3")

newblast = newblast[-which(newblast$V1 %in% non.concordant$V1), ] #remove non-concordant hits from probe blast

newblast$contigid = ""

newblast = newblast[which(!newblast$V1 %in% hap.2.3$contigid), ] # remove entries from probe blast that are also in the haplotype blast

combined = rbind(newblast, hap.2.3)

combined = combined[sort(combined$V2, index.return = T)$ix, ]


comb2 = newdf(colnames(combined), no.rows = T)

#do some sorting
for(i in unique(combined$V2)){
    temp = combined[combined$V2 == i, ]
    temp = temp[sort(as.numeric(temp$V9), index.return = T)$ix, ]
    comb2 = rbind(comb2, temp)
}

comb2$contigid[grep("AX", comb2$V1)] = comb2$V1[grep("AX", comb2$V1)] #populate the empty values in comb2$contigid with Axiom IDs

save(comb2, file = "rotation1scripts_v4/saved.objects/comb2")
load(file = "rotation1scripts_v4/saved.objects/comb2")

plot.blastdf(comb2)

