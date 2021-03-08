#     ____________________________________________________________________________
#     PROCESSING BEDTOOLS EXTRACTION FILE                                                                         ####


fafile = ReadFasta("rotation1scripts_v4/original_data/IWGSC/bed.tools.genes.only.extraction.fa")

fafile[] = lapply(fafile, as.character)

fafile$coords = lapply(strsplit(fafile[, 1], ":"), function(x) x[[2]])

gfffile = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.genes.only.gff3")


fafile[] = lapply(fafile, as.character)

fafile$startcoord = lapply(strsplit(fafile$coords, "-"), function(x) x[[1]])
fafile$endcoord = lapply(strsplit(fafile$coords, "-"), function(x) x[[2]])

fafile$genename = ""

fafile$genename = gfffile$V9

fafile$chr = lapply(strsplit(fafile$name, ":"), function(x) x[[1]])

fafile$name = paste(fafile$chr, fafile$genename)

fafilenew = fafile[, 1:2]

writefasta(fafilenew, "rotation1scripts_v4/original_data/IWGSC/iwgsc.full.gene.sequences.fa")

fafilenew2 = ReadFasta("rotation1scripts_v4/processed_data/fasta/unique.seg.dist.genes.v2.fasta")

fafilenew2 = convert.to.character.data.frame(fafilenew2)

list.of.seg.genes = unlist(lapply(lapply(strsplit(fafilenew2$name, " "), function(x) x[[1]]), function(x) substr(x, 1, (nchar(x)-2))))

fafile.seg.genes.extracted = newdf(colnames(fafilenew))

for(i in 1:length(list.of.seg.genes)){
    print(grep(list.of.seg.genes[i], fafilenew$name))
    fafile.seg.genes.extracted = rbind(fafile.seg.genes.extracted, fafilenew[grep(list.of.seg.genes[i], fafilenew$name), ])
}

fafile.seg.genes.extracted = fafile.seg.genes.extracted[-1, ]
fafile.seg.genes.extracted$name = unlist(lapply(lapply(strsplit(fafile.seg.genes.extracted$name, " "), function(x) x[[2]]), function(x) substr(x, 4, nchar(x))))

fafile.seg.genes.extracted$name = paste(unlist(lapply(fafile.seg.genes.extracted$name, function(x) substr(x, 1, (nchar(x)-12)))), ".1", sep = "")


writefasta(fafile.seg.genes.extracted, "rotation1scripts_v4/processed_data/fasta/unique.seg.dist.genes.full.genonmic.sequence.v3.fa")



blastfilenew = read.table("bioinf/blast/genes.vs.paragon.genome/results.blast/unique.seg.dist.genes.cs.x.para.vs.para.genome.v3.outputfmt6.word_sizev3.blast")
blastfilenew2 = grab.best.groups.of.hits(blastfilenew)




#     ____________________________________________________________________________
#     NEW EXTRACTION 08032018                                                                                                 ####

bedfa = ReadFasta("rotation1scripts_v4/processed_data/fasta/bedtools.genomic.extraction08032018.fa")
bedfa = convert.to.character.data.frame(bedfa)

max.orf.len = function(seq){
    max.orf = max_orf(seq, reverse.strand = T)
    g = c(max.orf$ORF.Forward$ORF.Max.Len, max.orf$ORF.Reverse$ORF.Max.Len)
    g2 = max(as.numeric(g))
    if(which(g == g2) == 1) for.or.rev = "forward"
    if(which(g == g2) == 2) for.or.rev = "reverse"
    
    if(for.or.rev == "forward"){
        attr(g2, "max.orf.seq") = max.orf$ORF.Forward$ORF.Max.Seq
    } else { 
        attr(g2, "max.orf.seq") = max.orf$ORF.Reverse$ORF.Max.Seq
    }
    attr(g2, "for.or.rev") = for.or.rev
    
    
    return(g2)
}

orf.lengths = unlist(lapply(bedfa$sequence, max.orf.len))

add.probe.names.to.bedtools.ex.file = function(bed.tools.df, blastdf){
    fafile = bed.tools.df
    
    fafile[] = lapply(fafile, as.character)
    
    fafile$coords = lapply(strsplit(fafile[, 1], ":"), function(x) x[[2]])    
    fafile[] = lapply(fafile, as.character)
    
    fafile$startcoord = lapply(strsplit(fafile$coords, "-"), function(x) x[[1]])
    fafile$endcoord = lapply(strsplit(fafile$coords, "-"), function(x) x[[2]])
    
    fafile$genename = ""
    fafile$startcoord = as.numeric(fafile$startcoord)
    fafile$endcoord = as.numeric(fafile$endcoord)
    
    fafile$probenames = as.character(blastdf$V1[unlist(lapply(fafile$startcoord, function(x){
        which(blastdf$V9 == (x + 2001) | blastdf$V10 == (x + 2001))
    }))])
    
    return(fafile)
}

bedfa2 = add.probe.names.to.bedtools.ex.file(bedfa, newblast)

bedfa2 = bedfa2[which(orf.lengths > 1000), ] #NEED TO ADD PROBE NAMES TO bedfa; SEE CODE ABOVE

bedfa2$max.orf.seq = unlist(lapply(bedfa2$sequence, function(x) attr(max.orf.len(x), "max.orf.seq")))

fastadf = bedfa2[, 7:8]
colnames(fastadf) = c("header", "sequence")
fastadf$header = unlist(lapply(fastadf$header, function(x)    p(">", x)))

writefasta(fastadf, "rotation1scripts_v4/processed_data/fasta/bedtools.genomic.orf.bigger.1000.fa")


#     ____________________________________________________________________________
#     bedtools blast                                                                                                                    ####

b = read.table("bioinf/blast/genes.vs.paragon.genome/results.blast/bedtools.orf.1000.cs.x.para.vs.para.genome.v3.outputfmt6.word_sizev3.blast")
b = read.table("bioinf/blast/probe.vs.genes.blast/results.blast/bedtools.orfs.refined.cs.x.para.vs.genes.blast")




