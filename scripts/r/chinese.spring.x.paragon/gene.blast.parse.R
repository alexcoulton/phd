#Parse the BLAST results of the CS segregation distortion genes against the Paragon genome assembly.
setwd("C:/Users/ac14037/Google Drive/University/PhD/project.phd.main/")
library(Rsamtools)
library(Biostrings)
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/generic.functions/fasta.R")

#     ____________________________________________________________________________
 #     DEFINE FUNCTIONS                                                                                                                ####

extract.sequences.from.genome.assembly = function(assembly.path, list.of.sequences.to.extract){
    if(exists("fa") == F){
        fa <<- FaFile(assembly.path)    
    }
    
    gr = as(seqinfo(fa), "GRanges")
    
    #extract listed sequences from assembly
    #NB this might take a long time (~ 15 mins to extract ~ 100 sequences from the Paragon genome assembly)
    sequences = lapply(unique(list.of.sequences.to.extract), function(x){
        idx = which(names(gr) == x)
        seq = getSeq(fa, gr[idx])
        return(seq)
    }
    )
    
    
    if(length(sequences) == 1) sequences = sequences[[1]][[1]]
    
    return(sequences)
    
}

extract.seq.from.paragon.assembly = function(list.of.seq.to.extract){
    extract.sequences.from.genome.assembly("rotation1scripts_v4/original_data/genome.assemblies/Triticum_aestivum_Paragon_EIv1.1.fa", list.of.seq.to.extract)
}

extract.seq.from.seg.dist.genes = function(name.of.seq, path.to.seg.dist.sequences){
    #GET QUERY SEQUENCE TO PASS TO perform.pairwise.alignment()
    if(missing(path.to.seg.dist.sequences)) path.to.seg.dist.sequences = "rotation1scripts_v4/processed_data/fasta/sacha.seg.dist.genes.full.genomic.sequence.2a.fa"
    if(exists("query.set") == FALSE){
        query.set = readDNAStringSet(filepath = path.to.seg.dist.sequences)
    }
    
    coord = grep(name.of.seq, names(query.set))
    query = query.set[coord]
    #query = query[[1]]
    return(query)
}

extract.seq.from.seg.dist.genes.cds = function(name.of.seq){
    #GET QUERY SEQUENCE TO PASS TO perform.pairwise.alignment()
    if(exists("query.set2") == FALSE){
        .GlobalEnv$query.set2 = readDNAStringSet(filepath = "rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")
    }
    
    coord = grep(name.of.seq, names(query.set2))
    query = query.set2[coord]
    #query = query[[1]]
    return(query)
}

perform.pairwise.alignment = function(query, subject, subject.start, subject.end, output.filepath, write.to.file){
    if(missing(write.to.file)) write.to.file = T
    subject = subject[subject.start:subject.end]
    align3 = pairwiseAlignment(query, subject)
    if(write.to.file == T) writePairwiseAlignments(align3, output.filepath)
    return(align3)
}

extract.relevant.sequences = function(query.name, scaffold.name){
    query = extract.seq.from.seg.dist.genes(query.name)
    scaffold = extract.seq.from.paragon.assembly(scaffold.name)
    
}

run.align = function(query.name, subject.name, blast.df, rev.subject, subset.start.buffer,
                                         subset.end.buffer, write.to.file){
    #query.name, subject.name: single character strings 
    #blast.df: the dataframe containing tabular blast output; required for subsetting
    #rev.subject: Boolean value specifying whether to reverse the subject sequence
    #subset.start.buffer: an integer specifying how many nt before the beginning of the BLAST
    # hit to subset the subject sequence; defaults to 20 if not specified.
    #subset.end.buffer: an integer specifying how many nt after the end of the BLAST
    # hit to subset the subject sequence; defaults to 20 if not specified.
    
    if(missing(rev.subject)){
        rev.subject = F
    }
    if(missing(subset.start.buffer)){
        subset.start.buffer = 20
    }
    if(missing(subset.end.buffer)){
        subset.end.buffer = 20
    }
    if(missing(write.to.file)){
        write.to.file = T
    }
    
    #This function originally ran a loop. I was too lazy to rewrite the whole thing so the following
    #lines simply convert the looping code structure to a single run.
    unique.blast.hit.pairs = newdf(c("V1", "V2"))
    unique.blast.hit.pairs[1, 1] = query.name
    unique.blast.hit.pairs[1, 2] = subject.name
    i = 1
    
    
    #GET QUERY SEQUENCE TO PASS TO perform.pairwise.alignment()
    query = all.sequences[[unique.blast.hit.pairs[i, 1]]] #see all.sequences object below -- saves computation as sequences are already extracted
    
    #subset main dataframe for ease of processing
    sub = blast.df[blast.df$V1 == unique.blast.hit.pairs[i, 1] & blast.df$V2 == unique.blast.hit.pairs[i, 2], ]
    sub[, 3:ncol(sub)][] = lapply(sub[, 3:ncol(sub)], as.numeric) #convert some columns to numeric
    
    #GET SUBJECT SEQUENCE TO PASS TO perform.pairwise.alignment()
    scaffold.name = unique.blast.hit.pairs[i, 2]
    scaffold = all.sequences[[unique.blast.hit.pairs[i, 2]]]
    
    # browser()
    
    print(head(as.character(scaffold[1:100])))
    print("--------------")
    if(rev.subject == T){
        scaffold = reverseComplement(scaffold)
    }
    print(head(as.character(scaffold[1:100])))
    scaffold.length = nchar(scaffold)
    
    
    print(" ")
    print(p(c("______ Hit # ", i, " ______")))
    # get start and end ranges for subsetting the subject based on the length of the query sequence
    if(sub$V9[[1]] < sub$V10[[1]]){ #if start position of blast hit against subject is SMALLER than end position of blast hit against subject
        # if((sub$V9[[1]] - sub$V14[[1]]) < 0 ){ #subtract the length of the gene from the start of the subject subset
        s.start = sub$V9[[1]] - subset.start.buffer
        s.end = (s.start + sub$V14[[1]] + subset.end.buffer) #sub$V14[[1]] is the length of the gene / query
        print("s.start = 0")
        # } else {
        #     s.start = (sub$V9[[1]] - 100)
        #     s.end = (s.start + sub$V14[[1]] + sub$V14[[1]] + sub$V14[[1]])
        #     print("s.start = (sub$V9[[1]] - sub$V14[[1]])")
        # }
    } else { #if start position of blast hit against subject is BIGGER than end position of blast hit against subject
        s.start = scaffold.length - sub$V9[[1]] - subset.start.buffer
        s.end = scaffold.length - sub$V10[[nrow(sub)]] + subset.end.buffer
        print("s.start = (sub$V9[[1]] + sub$V14[[1]])")
    }
    
    if(s.start > s.end & s.start > scaffold.length){
        s.start = scaffold.length
        print("s.start = scaffold.length")
    } else if(s.start < s.end & s.end > scaffold.length){
        s.end = scaffold.length
        print("s.end = scaffold.length")
    }
    
    if(s.start < 0) s.start = 1
    
    
    print(p(c(query.name, " vs ", subject.name)))
    
    print(p(c("Scaffold length: ", scaffold.length)))
    print(p(c("Query length: ", nchar(query))))
    print(p(c("# subhits: ", nrow(sub))))
    for(q in 1:nrow(sub)){
        print(p("        ---"))
        print(p(c("        Start of query subhit ", q, ": ", sub$V7[[q]])))
        print(p(c("        End of query subhit ", q, ": ", sub$V8[[q]])))
        print(p("        ------"))
        if(rev.subject == F){
            print(p(c("        Start of subject subhit ", q, ": ", sub$V9[[q]])))
            print(p(c("        End of subject subhit ", q, ": ", sub$V10[[q]])))
        } else {
            print(p(c("        Start of subject subhit ", q, ": ", scaffold.length - sub$V9[[1]])))
            print(p(c("        End of subject subhit ", q, ": ", scaffold.length - sub$V10[[1]])))
        }
        print("        ---")
    }
    print(p(c("Start of subject subset: ", s.start)))
    print(p(c("End of subject subset: ", s.end)))
    print(p(c("Length of subset: ", abs(s.end - s.start))))
    
    if(write.to.file == T){
        b = perform.pairwise.alignment(query, scaffold, s.start, s.end, p(c("rotation1scripts_v4/processed_data/pairwise.alignments/alignment.", query.name, "-", subject.name, ".txt")), write.to.file = T)
    }
    
    if(write.to.file == F){
        b = perform.pairwise.alignment(query, scaffold, s.start, s.end, p(c("rotation1scripts_v4/processed_data/pairwise.alignments/alignment.", query.name, "-", subject.name, ".txt")), write.to.file = F)
    }
    
    return(b)
    
}

get.orf.diff = function(list.of.alignments){ 
    #get the difference in maximum orf length between query and subject sequences in a list of pairwise alignments
    sapply(list.of.alignments, function(x){
        max.orf.pat = max_orf(gsub("-*", "", as.character(pattern(x))), reverse.strand = T)
        g = c(max.orf.pat$ORF.Forward$ORF.Max.Len, max.orf.pat$ORF.Reverse$ORF.Max.Len)
        g = max(as.numeric(g))
        
        max.orf.sub = max_orf(gsub("-*", "", as.character(subject(x))), reverse.strand = T)
        g2 = c(max.orf.sub$ORF.Forward$ORF.Max.Len, max.orf.sub$ORF.Reverse$ORF.Max.Len)
        g2 = max(as.numeric(g2))
        
        q = abs(g - g2)
        
        return(q)
    })
}

grab.query.seq = function(pairwise.alignment){
    #pairwise.alignment: a PairwiseAlignmentsSingleSubject object from bioconductor
    g = gsub("-*", "", as.character(pattern(pairwise.alignment)))
    return(g)
}

grab.subject.seq = function(pairwise.alignment){
    #pairwise.alignment: a PairwiseAlignmentsSingleSubject object from bioconductor
    g = gsub("-*", "", as.character(subject(pairwise.alignment)))
    return(g)
}

return.only.gene.titles.from.orf.output = function(char.vec){
    unlist(lapply(strsplit(char.vec, "\\."), function(x) paste(x[[1]], x[[2]], sep = ".")))
}

grab.and.write.sequences.to.file = function(seq.name){
    ex = extract.seq.from.seg.dist.genes(seq.name)
    ex2 = extract.seq.from.seg.dist.genes.cds(seq.name)    
    write(as.character(ex), p(c("rotation1scripts_v4/processed_data/fasta/individual.gene.sequences/", seq.name, ".genomic.fa")))
    write(as.character(ex2), p(c("rotation1scripts_v4/processed_data/fasta/individual.gene.sequences/", seq.name, ".cds.fa")))
}


#     ____________________________________________________________________________
#     BIOSTRINGS PAIRWISE ALIGNMENT                                                                                     ####

#some initial processing -- load seg. dist. genes from Sacha's genotype data set of CS X P

load(file = "rotation1scripts_v4/saved.objects/seg.unique.v3")

seg.unique.v3 = reset.colnames(seg.unique.v3, 1:14)

seg.unique.v3$V1 = multi.str.split(seg.unique.v3$V1, "=", 2)
seg.unique.v3$V1 = multi.str.split(seg.unique.v3$V1, ";", 1)

newg = seg.unique.v3

unique.blast.hit.pairs = unique(seg.unique.v3[, 1:2])

newg.high.coverage = newg[newg$per.cov > 90, ]
newg.low.coverage = newg[newg$per.cov < 90, ]

unique.blast.hit.pairs = unique(newg.high.coverage[, 1:2])

newg.high.coverage.rev.hits = newg.high.coverage[which(newg.high.coverage$V9 > newg.high.coverage$V10), ]
newg.high.coverage.forward.hits = newg.high.coverage[which(newg.high.coverage$V9 < newg.high.coverage$V10), ]

newg.high.coverage = newg.high.coverage.forward.hits
unique.blast.hit.pairs = unique(newg.high.coverage[, 1:2])

unique.blast.hit.pairs.rev = unique(newg.high.coverage.rev.hits[, 1:2])

unique.blast.hit.pairs = reset.rownames(unique.blast.hit.pairs)
newg.high.coverage = reset.rownames(newg.high.coverage)


#put sequences into RAM 
for.seq.genes = sapply(unique.blast.hit.pairs[, 1], extract.seq.from.seg.dist.genes)
for.seq.scaffold = sapply(unique.blast.hit.pairs[, 2], extract.seq.from.paragon.assembly)
forward.sequences = c(for.seq.genes, for.seq.scaffold)

s(forward.sequences, "rotation1scripts_v4/saved.objects/forward.sequences", "gene.blast.parse.R")


rev.seq.genes = sapply(unique.blast.hit.pairs.rev[, 1], extract.seq.from.seg.dist.genes)
rev.seq.scaffold = sapply(unique.blast.hit.pairs.rev[, 2], extract.seq.from.paragon.assembly)
rev.sequences = c(rev.seq.genes, rev.seq.scaffold)

all.sequences = c(forward.sequences, rev.sequences)
s(all.sequences, "rotation1scripts_v4/saved.objects/all.sequences", "gene.blast.parse.R")




#     ____________________________________________________________________________
#     DO PAIRWISE ALIGNMENTS AGAIN FOR ELIGIBLE                                                             ####


for.alignments = Map(function(x, y) run.align(x, y, newg.high.coverage, subset.start.buffer = 50, subset.end.buffer = 200, write.to.file = T), unique.blast.hit.pairs[, 1], unique.blast.hit.pairs[, 2])

for(i in 1:length(for.alignments)){
    writePairwiseAlignments(for.alignments[[i]], p("rotation1scripts_v4/processed_data/pairwise.alignments/latest", i, ".txt"))
}

rev.alignments = Map(function(x, y) run.align(x, y, newg.high.coverage.rev.hits, T, write.to.file = T, subset.start.buffer = 25, subset.end.buffer = 50), unique.blast.hit.pairs.rev[, 1], unique.blast.hit.pairs.rev[, 2])


#double check: alignment.TraesCS2A01G575200-Triticum_aestivum_Paragon_EIv1.1_scaffold_067961
#going to run this one manually:
manual.rev.alignment = run.align("TraesCS2A01G575200", "Triticum_aestivum_Paragon_EIv1.1_scaffold_067961", newg.high.coverage.rev.hits, T, write.to.file = T, subset.start.buffer = 50, subset.end.buffer = 60)

#grab the subject (Paragon) sequences from the alignments 
subject.seqs = sapply(for.alignments, subject)
forward.cds = sapply(unique.blast.hit.pairs[, 1], extract.seq.from.seg.dist.genes.cds)

#aligning chinese spring CDS sequences to Paragon sequences from the previous alignments
cds.alignments = lapply(names(subject.seqs), function(x){
    rev.cds = reverseComplement(forward.cds[[x]])
    align.pattern.forward = pairwiseAlignment(forward.cds[[x]], subject.seqs[[x]])
    align.pattern.rev = pairwiseAlignment(rev.cds, subject.seqs[[x]])
    
    #do some automatic parsing to check whether the CDS should be forward or reverse complement based on the number of rows in mismatchTable()
    if(nrow(mismatchTable(align.pattern.forward)) > nrow(mismatchTable(align.pattern.rev))) {
        writePairwiseAlignments(align.pattern.rev, p("rotation1scripts_v4/processed_data/pairwise.alignments/cds.to.paragon/", x, ".rev.txt"))    
        return(align.pattern.rev)
    } else {
        writePairwiseAlignments(align.pattern.forward, p("rotation1scripts_v4/processed_data/pairwise.alignments/cds.to.paragon/", x, ".txt"))    
        return(align.pattern.forward)
    }
})


as.character(pattern(cds.alignments[[2]]))


raw.subject.string = as.character(subject(cds.alignments[[2]]))
raw.subject.string = strsplit(raw.subject.string, "")[[1]]

subject.string = subject(cds.alignments[[2]])
del.table = deletion(cds.alignments[[2]])
ins.table = insertion(cds.alignments[[2]])

g = as.numeric()
starts = start(del.table)[[1]]
starts = as.list(starts)
ends = end(del.table)[[1]]
ends = as.list(ends)
for(i in 1:length(starts)){
    g = c(g, starts[[i]]:ends[[i]])
}

raw.subject.string.cut = raw.subject.string[-g]
raw.subject.string.cut = paste(raw.subject.string.cut, collapse = "")

raw.cut = DNAString(raw.subject.string.cut)

pair.test = pairwiseAlignment(DNAString(as.character(pattern(cds.alignments[[2]]))), raw.cut)
writePairwiseAlignments(pair.test, "rotation1scripts_v4/processed_data/pairwise.alignments/testertt.txt")
writePairwiseAlignments(cds.alignments[[2]], "rotation1scripts_v4/processed_data/pairwise.alignments/sdlkfj.txt")



library(LncFinder) #need this for max_orf function



orf.differences = get.orf.diff(all.align)

#The following code is used in conjunction with geneblast2 object in seg.dist.marker.analysis.R
#here we are simply finding the corresponding probe names for our genes of interest in orf.differences
probes = geneblast[unlist(lapply(names(orf.differences), function(x) which(x == geneblast$V2))), 1:2]
#now we want to grab only the probes that fall under the stricter criteria for no homeologous interference
probes.refined = probes[which(probes$V1 %in% newblast$V1), ]

#SAVE probes.refined to file
save(probes.refined, file = "rotation1scripts_v4/saved.objects/gene.blast.parse.R_probes.refined")

#the following g is defined on line 34 of seg.dist.marker.analysis.R
seg.dist.data = g[match(probes.refined$V1, names(g))]

orf.diff.refined = orf.differences[which(names(orf.differences) %in% probes.refined$V2)]

all.align.refined = all.align[which(names(all.align) %in% probes.refined$V2)]



forward.orfs$genetitles = return.only.gene.titles.from.orf.output(forward.orfs$cdsorf)
rev.orfs$genetitles = return.only.gene.titles.from.orf.output(rev.orfs$cdsorf)


forward.orfs[na.omit(match(names(all.align.refined), forward.orfs$genetitles)), ]
rev.orfs[na.omit(match(names(all.align.refined), rev.orfs$genetitles)), ]

all.orfs = rbind(forward.orfs, rev.orfs)

#some code to compare orf output produced with max_orf and getorf command line application
# all.orfs$patternorf = pattern.orfs[match(all.orfs$genetitles, names(pattern.orfs))]
# all.orfs$subject.orfs = subject.orfs[match(all.orfs$genetitles, names(subject.orfs))]
# 
# all.orfs.refined = all.orfs[-which(unname(unlist(lapply(all.orfs$patternorf, is.null)))), ]




#the alignments in all.align.refined are all forward alignments -- I had to convert reverse alignments when using run.align()


#remove some alignments that have big gaps between cds and genomic sequence


#grab both the biggest orfs and their dna sequences for pattern in alignments

alignment.seq.analysis = function(query.or.subject, alignment){
    if(query.or.subject != "q" & query.or.subject != "s"){ 
        print("Please specify either query (q) or subject (s)")
        } else {
            if(query.or.subject == "q") query.seq = grab.query.seq(alignment)
            if(query.or.subject == "s") query.seq = grab.subject.seq(alignment)
            
            length.of.q = nchar(query.seq)
            
            g = max_orf(query.seq, reverse.strand = T)
            biggest.orf = max(c(g$ORF.Forward$ORF.Max.Len, g$ORF.Reverse$ORF.Max.Len))
            diff.between.len.and.max.orf = abs(length.of.q - biggest.orf)
            
            for.len = nchar(g$ORF.Forward$ORF.Max.Seq)
            rev.len = nchar(g$ORF.Reverse$ORF.Max.Seq)
            
            if(for.len > rev.len){
                q = list(biggest.orf, g$ORF.Forward$ORF.Max.Seq, diff.between.len.and.max.orf)
                names(q) = c("Max.orf.len", "orf.seq", "orf.full.len.diff")
                return(q)
            } else {
                q = list(biggest.orf, g$ORF.Reverse$ORF.Max.Seq, diff.between.len.and.max.orf)
                names(q) = c("Max.orf.len", "orf.seq", "orf.full.len.diff")
                return(q)
            }
        }
}


pattern.orfs = lapply(all.align.refined, function(x) alignment.seq.analysis("q", x))
pattern.orfs2 = pattern.orfs[-which(names(pattern.orfs) %in% c("TraesCS2D01G123200.1", "TraesCS5B01G010200.1", "TraesCS5B01G295500.1"))] #remove genes with discrepancies between genomic and cds sequences
pattern.orfs2 = pattern.orfs2[-which(names(pattern.orfs2) == "TraesCSU01G090700.1" | names(pattern.orfs2) == "TraesCS3A01G272400.1")] #remove genes with discrepancies between probe vs genome blast and gff coordinates

#grab both the biggest orfs and their dna sequences for subject in alignments
subject.orfs = lapply(all.align.refined, function(x) alignment.seq.analysis("s", x))
subject.orfs2 = subject.orfs[-which(names(subject.orfs) %in% c("TraesCS2D01G123200.1", "TraesCS5B01G010200.1", "TraesCS5B01G295500.1"))] #remove genes with discrepancies between genomic and cds sequences
subject.orfs2 = subject.orfs2[-which(names(subject.orfs2) == "TraesCSU01G090700.1" | names(subject.orfs2) == "TraesCS3A01G272400.1")] #remove genes with discrepancies between probe vs genome blast and gff coordinates


pattern.aa.seqs = lapply(pattern.orfs2, function(x){
    Biostrings::translate(DNAString(x$orf.seq))
})

subject.aa.seqs = lapply(subject.orfs2, function(x){
    Biostrings::translate(DNAString(x$orf.seq))
})


Map(function(x, y, z){
    g = pairwiseAlignment(x, y)
    writePairwiseAlignments(g, p(c("rotation1scripts_v4/processed_data/pairwise.alignments/aa.alignments/", z, ".txt")))
}, pattern.aa.seqs, subject.aa.seqs, names(pattern.aa.seqs))

probes.refined2 = probes.refined[match(names(subject.aa.seqs), probes.refined$V2), ]

seg.dist.data[match(as.character(probes.refined[match(names(subject.aa.seqs), probes.refined$V2), ]$V1), names(seg.dist.data))]

functional.anno = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", sep = "\t", header = T, fill = T)
functional.anno$gene2 = multi.str.split(as.character(functional.anno$Gene.ID), "\\.", 1)



which(functional.anno$Gene.ID %in% as.character(probes.refined2$V2))

which(length(unlist(lapply(strsplit(as.character(functional.anno$Gene.ID), "\\."), function(x) x[1])), function(x) grep(x, probes.refined$V2))    > 0)

q1 = unlist(lapply(strsplit(as.character(functional.anno$Gene.ID), "\\."), function(x) x[1]))

q2 = unlist(lapply(strsplit(as.character(probes.refined2$V2), "\\."), function(x) x[1]))

functional.anno[which(q1 %in% q2), ]

dnaset = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")

select.probe.vs.genome.blast = newblast[match(as.character(probes.refined2$V1), newblast$V1), ]

select.probe.vs.genome.blast$genename = probes.refined2$V2

genome.gff.file = read.delim("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3", header=F, comment.char="#")

genome.gff.file.subset = genome.gff.file[unlist(lapply(select.probe.vs.genome.blast$genename, function(x) grep(x, genome.gff.file$V9))), ]
gen.sub.2 = genome.gff.file.subset[which(genome.gff.file.subset$V3 == "mRNA"),]

gen.sub.2$diff.prob.gff = as.numeric(gen.sub.2$V4) - as.numeric(select.probe.vs.genome.blast$V9)

save(gen.sub.2, file = "rotation1scripts_v4/saved.objects/gene.blast.parse.R_gen.sub.2")

genes.w.big.diff.between.blast.and.gff = c("TraesCSU01G090700.1", "TraesCS3A01G272400.1")





#     ____________________________________________________________________________
#     some investigation                                                                                                            ####

#want to check concordance between CDS sequences and genomic sequences for all of the genes in the list
lapply(names(pattern.aa.seqs), function(x){
    genomic = extract.seq.from.seg.dist.genes(x)
    cds = extract.seq.from.seg.dist.genes.cds(x)
    pair = pairwiseAlignment(genomic, cds)
    writePairwiseAlignments(pair, p(c("rotation1scripts_v4/processed_data/pairwise.alignments/cds.to.genomic/", x, ".txt")))
})

#some of the previous alignments look like they need to be using the reverseComplement of the genomic sequence
genes.to.rev = c("TraesCS5B01G010200.1", "TraesCS6A01G101200.1", "TraesCSU01G070900.1", "TraesCSU01G090700.1")

lapply(genes.to.rev, function(x){
    genomic = reverseComplement(extract.seq.from.seg.dist.genes(x))
    cds = extract.seq.from.seg.dist.genes.cds(x)
    pair = pairwiseAlignment(genomic, cds)
    writePairwiseAlignments(pair, p(c("rotation1scripts_v4/processed_data/pairwise.alignments/cds.to.genomic/", x, ".rev.txt")))
})



#     ____________________________________________________________________________
#     SCRATCH CODE                                                                                                                        ####



lapply(unique.blast.hit.pairs[, 1], grab.and.write.sequences.to.file)

###

grab.and.write.sequences.to.file(unique.blast.hit.pairs[3, 1])

unique.blast.hit.pairs[1, 1]

ex = extract.seq.from.seg.dist.genes("TraesCS1B01G480100.1")

ex2 = extract.seq.from.seg.dist.genes.cds("TraesCS1B01G480100.1")

ex3 = pairwiseAlignment(ex, ex2, type = "global")

writePairwiseAlignments(ex3, "rotation1scripts_v4/processed_data/pairwise.alignments/blablabalaa.txt")

###

lapply(unique.blast.hit.pairs.rev[, 1], grab.and.write.sequences.to.file)
