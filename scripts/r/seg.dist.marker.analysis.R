#segregation distortion marker analysis
#processes axiom analysis suite output and grabs probes showing significant segregation distortion (chi.sq p < 0.05)

setwd("C:/Users/ac14037/project.phd.main/")

library(qtl)
library(dplyr)
library(tibble)

source("rotation1scripts_v4/scripts/r/functions.R")

# g = convert.aas.to.rqtl("original_data/axiom analysis suite output/all.polyhighres.snps.txt")
# write.csv(g, "processed_data/genotypes/allpolyhighresmarkers.csv", row.names = F)


#     ____________________________________________________________________________
#     Grab seg.dist markers with deviation from 1:1 ratio (custom function)     ####

cs.x.p.path = "rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/cs.x.para.flipped.csv"
a.x.c.path = "rotation1scripts_v4/processed_data/genotypes/avalon.x.cadenza/a.x.c.flipped.csv"
o.x.s.path = "rotation1scripts_v4/processed_data/genotypes/opata.x.synthetic/o.x.s.flipped.csv"
r.x.s.path = "rotation1scripts_v4/processed_data/genotypes/rialto.x.savannah/rialto.x.savannah.flipped.csv"
s.x.s.path = "rotation1scripts_v4/processed_data/genotypes/shamrock.x.shango/shamrock.x.shango.flipped.csv"


#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####


geno.count.chi = function(genovector){
    #genovector - a character vector of genotypes "A", "B" and "H"
    A = length(which(genovector == "A"))
    B = length(which(genovector == "B"))
    H = length(which(genovector == "H"))
    all = c(A, B)
    c = chisq.test(all)
    res = c(c[[1]], c[[2]], c[[3]], A, B, H)
    names(res) = c("X-squared", "df", "p-value", "A", "B", "H")
    return(res)
}

geno.count.chi.f2 = function(genovector){
    #genovector - a character vector of genotypes "A", "B" and "H"
    total.genotypes = length(genovector)
    
    A = length(which(genovector == "A"))
    B = length(which(genovector == "B"))
    H = length(which(genovector == "H"))
    
    all.obs = c(A, B, H)
    
    
    exp.A = (total.genotypes / 4)
    exp.B = (total.genotypes / 4)
    exp.H = (total.genotypes / 2)
    
    all.exp = c(exp.A, exp.B, exp.H)
    
    sum = 0
    for(i in 1:3){
        sum = sum + (((all.obs[[i]] - all.exp[[i]])^2) / all.exp[[i]])
    }
    
    p.value = pchisq(sum, df = 2, lower.tail = F) #calculate p-value from chi-square value
    
    res = c(sum, 2, p.value, A, B, H)
    names(res) = c("X-squared", "df", "p-value", "A", "B", "H")
    return(res)
}

grab.seg.dist.markers = function(path.to.genotype.data, is.f2, p.threshold){
    if(missing(is.f2)) is.f2 = F
    if(missing(p.threshold)) p.threshold = 0.005
    
    if(is.f2) p.threshold = 0.0005 #my custom weighted chi-square function seems to produce much smaller p values than the inbuilt one, so a lower threshold is needed
    
    
    geno.data = read.csv(path.to.genotype.data, header = T, stringsAsFactors = F)
    
    pop.size = nrow(geno.data)
    
    if(is.f2 == F){
        g = lapply(geno.data[, 3:ncol(geno.data)], function(x){
            geno.count.chi(x[2:length(x)])
        })
    } else {
        g = lapply(geno.data[, 3:ncol(geno.data)], function(x){
            geno.count.chi.f2(x[2:length(x)])
        })
    }
    
    
    if(is.f2 == F){
        seg.dist.markers = g[which(unlist(lapply(g, function(x){
            x[[3]] < p.threshold
        })))]
    } else {
        seg.dist.markers = g[which(unlist(lapply(g, function(x){
            x[[3]] < p.threshold
        })))]
    }
    
    seg.dist.markers = t(as.data.frame(seg.dist.markers))
    
    if(nrow(seg.dist.markers) == 0){
        print("No seg. dist. markers found for this threshold")
        return()
    } else {
        if(is.f2 == F & length(which(seg.dist.markers[, 6] > 20)) != 0) seg.dist.markers = seg.dist.markers[-which(seg.dist.markers[, 6] > 20), ]
        
        #seg.dist.markers = names(seg.dist.markers)
        return(seg.dist.markers)
    }
    
}


#     ____________________________________________________________________________
#     GRAB SEG. DIST. MARKERS                                                                                                 ####


#load comb2, a df from nimblegen.contig.processing.R that contains BLAST results from haplotype data and probes combined
load(file = "rotation1scripts_v4/saved.objects/comb2")

cs.x.p.seg.dist.markers = grab.seg.dist.markers(cs.x.p.path) #chinese spring x paragon

#since removing false seg.dist. markers (containing > 20 heterozygotes), comb2 now contains some invalid markers
#fixing this:
comb2 = comb2[which(comb2$contigid %in% rownames(cs.x.p.seg.dist.markers)), ]
comb2.1 = comb2
save(comb2.1, file = "rotation1scripts_v4/saved.objects/comb2.1")

#order the markers 
cs.x.p.seg.dist.markers2 = cs.x.p.seg.dist.markers[na.omit(match(comb2$contigid, rownames(cs.x.p.seg.dist.markers))), ]
#add BLAST information to df
cs.x.p.seg.dist.markers2 = cbind(cs.x.p.seg.dist.markers2, comb2$V2, comb2$V9, comb2$V1)

attr(cs.x.p.seg.dist.markers2, "original.script") = "seg.dist.marker.analysis.R"
save(cs.x.p.seg.dist.markers2, filepath = "rotation1scripts_v4/saved.objects/cs.x.p.seg.dist.markers2")


#grab seg.dist markers for other crosses
a.x.c.seg.dist.markers = grab.seg.dist.markers(a.x.c.path) #avalon x cadenza
o.x.s.seg.dist.markers = grab.seg.dist.markers(o.x.s.path) #opata x synthetic
r.x.s.seg.dist.markers = grab.seg.dist.markers(r.x.s.path) #rialto x savannah
s.x.s.seg.dist.markers = grab.seg.dist.markers(s.x.s.path) #shamrock x shango
a.x.p.s.d.2 = grab.seg.dist.markers(a.x.p.path, T)


cross.list = list(cs.x.p.seg.dist.markers, a.x.c.seg.dist.markers, o.x.s.seg.dist.markers, r.x.s.seg.dist.markers, s.x.s.seg.dist.markers, a.x.p.s.d.2)
names(cross.list) = c("cs.x.p", "ava.x.cad", "opa.x.syn", "rial.x.sav", "sham.x.shango", "apo.x.para")

cross.list = lapply(cross.list, as.data.frame)



#     ____________________________________________________________________________
#     grab seg.dist markers using rqtl (deviation from 1:2:1)                                 ####

#not needed now that I've made my own function for F2 seg.dist. identification

# library(qtl)
# a.x.p.path = "rotation1scripts_v4/original_data/genotype.data/manuallycuratedsnps_flipped.csv"
# geno.data = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/original_data/genotype.data/manuallycuratedsnps_flipped.csv", genotypes = c("A", "H", "B"), estimate.map = F)
# # 
# gt = geno.table(geno.data)
# #check deviation from 1:2:1 (note this is for an F2)
# a.x.p.seg.dist.markers = gt[gt$P.value < 0.05/totmar(geno.data), ][, -c(1, 2, 6, 7)]
# 





#     ____________________________________________________________________________
#     EXAMINE SHARED MARKERS                                                                                                    ####

grab.shared.markers = function(cross1, cross2){
    shared = cross1[which(rownames(cross1) %in% rownames(cross2)), ]
    shared.markers.reciprocal = cross2[which(rownames(cross2) %in% rownames(cross1)), ]
    
    shared.markers.reciprocal = shared.markers.reciprocal[match(rownames(shared), rownames(shared.markers.reciprocal)), ]
    
    g = list(shared, shared.markers.reciprocal)
    names(g) = c("cross1", "cross2")
    return(g)
}

combinations = t(combn(1:6, 2)) #grab all unique comparisons of crosses to index cross.list

comparisons = Map(function(x, y){
    cross = c(names(cross.list)[x], names(cross.list)[y])
    g = grab.shared.markers(cross.list[[x]], cross.list[[y]])
    
    attr(g, "cross") = cross
    return(g)
    }, combinations[, 1], combinations[, 2]
)

comparisons =    comparisons[-which(unlist(lapply(comparisons, function(x) nrow(x[[1]]))) == 0)]

lapply(1:5, function(x) lapply(rownames(comparisons[[x]]$cross1), grab.individual.hit)) #grab physical locations for shared seg. dist. markers 


#     ____________________________________________________________________________
#     grab probe sequences for seg.dist.markers                                                             ####


convert.probes.to.fasta = function(listofprobes, output.filepath){
    #take a character vector of probe names and output a fasta file containing the sequences for those probes
    allprobeseq = read.table("rotation1scripts_v4/original_data/fasta/35kprobes.fa")
    
    nrow(allprobeseq)
    g = allprobeseq[seq(1, nrow(allprobeseq), 2),]
    g1 = allprobeseq[seq(2, nrow(allprobeseq), 2),]
    
    allprobeseq2 = data.frame(g, g1)
    
    listofprobes = allprobeseq2[which(allprobeseq2$g %in% paste(">", listofprobes, sep = "")), ]
    listofprobes[] = lapply(listofprobes, as.character)
    
    g = character()
    
    for(i in 1:nrow(listofprobes)){
        g = c(g, listofprobes[i, 1], listofprobes[i, 2])    
    }
    
    fileconn = file(output.filepath)
    writeLines(g, fileconn)
    close(fileconn)
}

# convert.probes.to.fasta(seg.dist.probe.seqs, "processed_data/fasta/seg.dist.genes.cs.x.para.fasta")

#     ____________________________________________________________________________
#     ANALYSE CS.X.PARA PROBES VS GENOME BLAST                                                                ####
#             (identify probes which are least likely to suffer from homeologous interference)
csblast = read.table("rotation1scripts_v4/pull/cs.x.para.probes.vs.genome.blast")

bedconversion = read.table("bioinf/blast/probe.vs.genome.blast/results.blast/cs.x.para.probes.vs.genome.cleaned.blast.bed")
as.numeric(bedconversion$V2 + 1) %in% as.numeric(csblast$V9) | as.numeric(bedconversion$V2 + 1) %in% as.numeric(csblast$V10) #check bed conversion worked properly (start coords same as in BLAST file [except for a difference of + 1])

#do some editing of bed file
all(bedconversion$V2 < bedconversion$V3) #check all start smaller than all end
bedconversion$V2 = bedconversion$V2 - 2000 #4000 bp to coords for genomic sequence extraction
bedconversion$V3 = bedconversion$V3 + 2000

write.blast(bedconversion, "bioinf/blast/probe.vs.genome.blast/results.blast/cs.x.para.probes.vs.genome.cleaned.blast.+4000.bed")

newblast = clean.blast(csblast)
save(newblast, file = "rotation1scripts_v4/saved.objects/newblast")
write.blast(newblast, "bioinf/blast/probe.vs.genome.blast/results.blast/cs.x.para.probes.vs.genome.cleaned.blast")


convert.probes.to.fasta(newblast$V1, "rotation1scripts_v4/processed_data/fasta/unique.seg.dist.probes.refined.cs.x.para.fa")

#     ____________________________________________________________________________
#     ANALYSE BLAST VS GENES                                                                                                    ####

source("scripts/r/generic.functions/fasta.R")

#note: this "refined" blast search has a stricter set of criteria with which probes were selected from the probe vs genome blast.
#namely, instead of taking the hit with the minimum e-value from a group of hits with the same query-subject pair,
#the e-value had to be smaller than the next best hit by an empirically determined margin that amounted to a difference of more than 3 query-subject mismatches between
#the best and the second best hit.
geneblast.upd = read.table("bioinf/blast/probe.vs.genes.blast/results.blast/unique.seg.dist.probes.refined.cs.x.para.vs.genes.blast")
geneblast.upd2 = clean.blast(geneblast.upd)


geneblast = read.table("rotation1scripts_v4/processed_data/BLAST/unique.seg.dist.probes.cs.x.para.vs.genes.blast")

geneblast2 = clean.blast(geneblast)

duplicate.elements = find.dupes(geneblast2$V2)
geneblast2 = reset.rownames(geneblast2)

functional.anno = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", sep = "\t", header = T, fill = T)

functional.anno2 = functional.anno[which(functional.anno$Gene.ID %in% geneblast2$V2), ]

# g = ReadFasta("original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")
# 
# g = convert.to.character.data.frame(g)
# 
# g2 = g[which(unlist(lapply(strsplit(g$name, split = " "), function(x) x[[1]])) %in% geneblast2$V2), ]

#length(which(geneblast2$V2 %in% lapply(strsplit(g2[, 1], split = " "), function(x) x[[1]]))) #check all genes are present in fasta

# writefasta(g2, "processed_data/fasta/unique.seg.dist.genes.v2.fasta")


#     ____________________________________________________________________________
#     ANALYSE BLAST GENES VS GENOME                                                                                     ####
#seg. dist. genes cs x para vs Paragon genome assembly

g = read.table("processed_data/BLAST/unique.seg.dist.genes.cs.x.para.vs.para.genome.v2.outputfmt6.word_size.blast")




