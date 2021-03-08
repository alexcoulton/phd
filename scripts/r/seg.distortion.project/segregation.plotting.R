setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
load("rotation1scripts_v4/saved.objects/all.m.8.iwgsc.4b.rev")
load("rotation1scripts_v4/saved.objects/all.maps.geno.data")

other.p.cross.genotypes = list.files("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/")[-2]

other.p.cross.genotypes2 = lapply(other.p.cross.genotypes, function(x){
    read.csv(p("rotation1scripts_v4/processed_data/genotypes/other.paragon.mapping.populations/", x), stringsAsFactors = F, header = T)
})

other.p.cross.maps = list.files("rotation1scripts_v4/processed_data/genetic_maps/other.paragon.mapping.populations/")
other.p.cross.maps2 = lapply(other.p.cross.maps, function(x){
    read.csv(p("rotation1scripts_v4/processed_data/genetic_maps/other.paragon.mapping.populations/", x), stringsAsFactors = F, header = T)
})

#order genotypes according to genetic map and add chromosomes to linkage groups
other.p.cross.genotypes3 = Map(function(geno, map){
    # browser()
    geno = convert.to.character.data.frame(geno)
    map = convert.to.character.data.frame(map)
    match.coords = match(as.character(switch.affy.format(map$X)), as.character(colnames(geno)))
    geno = geno[, c(1, 2, match.coords)]
    geno[1, 3:ncol(geno)] = map$chromo
    geno
}, other.p.cross.genotypes2, other.p.cross.maps2)

other.p.cross.maps3 = lapply(other.p.cross.maps2, function(x){
    g = split.df.into.list.of.dfs.by.column(x, "chromo")
    names(g) = unlist(lapply(g, function(x) unique(x$chromo)))
    g
})

names(other.p.cross.maps3) = other.p.cross.maps

all.maps.geno.data = c(all.maps.geno.data, other.p.cross.genotypes3)
all.m.8.iwgsc.4b.rev = c(all.m.8.iwgsc.4b.rev, other.p.cross.maps3)


#grab markers that are present in the 5a map from cs x p and also the genotype data from the a x p cross
coords = na.omit(unlist(lapply(all.m.8.iwgsc.4b.rev$cs.x.p[["5A"]]$marker, function(x){
    match(x, switch.affy.format(colnames(other.p.cross.genotypes2[[1]])))
})))

shared.genotypes = other.p.cross.genotypes2[[1]][, c(1, 2, coords)]

a.x.p.shared.seg = unlist(lapply(shared.genotypes, function(x){
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    h = length(which(x == "H"))
    a / b
}))

a.x.p.shared.phys = genetictophysical("5A", names(a.x.p.shared.seg), ".", return.bp = T)

plot(attr(a.x.p.shared.phys, "bp"), a.x.p.shared.seg, ylim = c(0.5, 1.5)) + abline(1, 0) + axis(1, at = seq(0, 8e+08, 2e+07))

shared.data = as.data.frame(list(attr(a.x.p.shared.phys, "bp"), a.x.p.shared.seg))
colnames(shared.data) = c("phys", "seg")

pdf("rotation1scripts_v4/plots/seg.dist.5a/a.x.p.wide.pdf", 100, 10)
ggplot(shared.data, aes(x = phys, y = seg)) + geom_point() + coord_cartesian(ylim = c(0.5, 1.5)) + geom_hline(yintercept = 1) + scale_x_continuous(breaks = seq(0, 8e+08, 2e+07))
dev.off()

#     ____________________________________________________________________________
#     GENE DENSITY FUNCTIONS                                                                                                    ####

#examine gene density
genes = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", header = T, stringsAsFactors = F)
genes$gene.names = multi.str.split(as.character(genes$V10), "=", 2)
genes$gene.names = multi.str.split(genes$gene.names, ";", 1)
genes = genes[, -1]
genes = reset.colnames(genes)
colnames(genes) = c(colnames(genes)[1:9], "gene.names")
genes$V4 = as.numeric(genes$V4)

check.gene.density = function(chromosome, start.region.bp, end.region.bp, return.df){
    if(missing(return.df)) return.df = F
    if(("genes" %in% ls(.GlobalEnv)) == F) genes = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.gff")
    genes2d = genes[genes$V1 == p("chr", chromosome), ]
    
    if(return.df == F){
        #how many genes are in this particular region from cs.x.p.seg.dist.markers2    
        return(nrow(genes2d[genes2d$V4 > start.region.bp & genes2d$V4 < end.region.bp, ]))
    } else {
        return(genes2d[genes2d$V4 > start.region.bp & genes2d$V4 < end.region.bp, ])
    }
}



library(readr)
# functional.anno = read_delim("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", delim = "\t")

functional.anno = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", sep = "\t", header = T, fill = T)

functional.anno$Gene.ID = as.character(functional.anno$Gene.ID)

grab.genes.and.functions = function(chromosome, start.bp, end.bp){
    func.temp = lapply(check.gene.density(chromosome, start.bp, end.bp, T)$gene.names, function(x) grep(x, functional.anno$Gene.ID))
    genes.func = functional.anno[na.omit(unlist(lapply(func.temp, function(x) x[1]))), ]
    genes.func = convert.to.character.data.frame(genes.func)
    return(genes.func)
}


#     ____________________________________________________________________________
#     SEGREGATION PLOTTING                                                                                                        ####

all.maps.geno.data.axp.rm = all.maps.geno.data[-4]


#setup a closure / funcitonal that applies new.func() to counts of A and B genotypes
grab.seg.info.closure = function(new.func){
    grab.seg.info.base = function(map.number, chromo.name){
        lapply(all.maps.geno.data.axp.rm[[map.number]][, which(all.maps.geno.data.axp.rm[[map.number]][1, ] == chromo.name)], function(x){
            x = x[2:length(x)]
            x.a = length(which(x == "A"))
            x.b = length(which(x == "B"))
            new.func(x.a, x.b)
        })
    }    
}

#setup functions to grab segregation ratios and chi-square values
grab.seg.ratios = grab.seg.info.closure(function(x, y) x / y)
grab.chi.vals = grab.seg.info.closure(function(x, y) chisq.test(c(x, y))$p.value)

#setup function to plot the data
library(ggplot2)
gen.plot = function(map.number, chromo.name, return.df){
    if(missing(return.df)) return.df = F
    # browser()
    g1 = data.frame(all.m.8.iwgsc.4b.rev[[map.number]][[chromo.name]]$phys.dist, unlist(grab.seg.ratios(map.number, chromo.name)),
                                    unlist(grab.chi.vals(map.number, chromo.name)),
                                    all.m.8.iwgsc.4b.rev[[map.number]][[chromo.name]]$phys.dist.bp)
    g1$sig = ""
    colnames(g1) = c("phys.dist", "seg.ratios", "chi", "phys.dist.bp", "sig")
    g1$sig[which(g1$chi < 0.05)] = "T"
    g1$sig[which(!g1$chi < 0.05)] = "F"
    # browser()
    if(return.df == T){
        return(g1)
    } else {
        return(ggplot(data = g1, aes(x = phys.dist.bp, y = seg.ratios, group = sig, color = sig)) + geom_point() + geom_hline(yintercept = 1) + geom_vline(xintercept = 689955376) +
                         labs(title = p(names(all.m.8.iwgsc.4b.rev)[map.number], ": ", chromo.name),
                                    x = "Physical distance (bp)", y = "Segregation ratio (A/B)") + scale_x_continuous(breaks = seq(0, 8e+08, 2e+07))) #+
                         # coord_cartesian(xlim = c(0, 100), ylim = c(0.5, 1.5)))
    }
}

pdf(p("rotation1scripts_v4/plots/seg.dist.5a/cs.x.p/", "5a.w.vline", "cs.x.p.wide2.pdf"), 100, 10)
print(gen.plot(4, "5A"))
dev.off()


for(x in names(all.m.8.iwgsc.4b.rev[[4]])){
    
    pdf(p("rotation1scripts_v4/plots/seg.dist.5a/cs.x.p/", x, "cs.x.p.wide2.pdf"), 100, 10)
    print(gen.plot(4, x))
    dev.off()
}

v(gen.plot(4, "5A", T))

gen.plot(4, "5A")


gen.plot.wo.colour.groups = function(map.number, chromo.name, return.df){
    if(missing(return.df)) return.df = F
    g1 = data.frame(all.m.8.iwgsc.4b.rev[[map.number]][[chromo.name]]$phys.dist, unlist(grab.seg.ratios(map.number, chromo.name)),
                                    unlist(grab.chi.vals(map.number, chromo.name)),
                                    all.m.8.iwgsc.4b.rev[[map.number]][[chromo.name]]$phys.dist.bp)
    g1$sig = ""
    colnames(g1) = c("phys.dist", "seg.ratios", "chi", "phys.dist.bp", "sig")
    g1$sig[which(g1$chi < 0.05)] = "T"
    g1$sig[which(!g1$chi < 0.05)] = "F"
    g1$seg.ratios = ((g1$seg.ratios - 1) * 1.7) + 1
    # browser()
    if(return.df == T){
        return(g1)
    } else {
        return(ggplot(data = g1, aes(x = phys.dist, y = seg.ratios )) + geom_point(size = 2.5, aes(color = seg.ratios)) + 
                         scale_colour_gradientn(limits = c(0.5, 1.5), breaks = c(0.5, 1.5), 
                                                                        labels = c(expression(paste("P"[2], " Allele")), expression(paste("P"[1], " Allele"))), 
                                                                        colours = c("#0015ff", "#3a4aff", "#7f89ff", "#ff7f7f", "#ff4444", "#ff0000"),
                                                                        name = "Ratio") +
                         geom_hline(yintercept = 1) + geom_vline(xintercept = 97.16) +
                         labs(title = "Empirical F5 cross with lots of distortion",
                                    x = "Physical distance along chromosome (%)", y = expression(paste("Segregation ratio (P"[1], "/ P"[2], ")"))) +
                         coord_cartesian(xlim = c(0, 100), ylim = c(0.5, 1.5)) + theme_bw(base_size = 15)) 
    }
}

gen.plot.wo.colour.groups(map.number = 4, chromo.name = "5A", return.df = F)

no.dist.df = data.frame(seq(1, 100, 0.2), rnorm(496, mean = 1, sd = 0.01))
no.dist.df = no.dist.df[-sample(496, 300), ]
colnames(no.dist.df) = c("phys.dist", "seg.ratios")

no.seg.dist.plot = ggplot(data = no.dist.df, aes(x = phys.dist, y = seg.ratios)) + geom_point(size = 2.5, aes(color = seg.ratios)) + 
    scale_colour_gradientn(limits = c(0.5, 1.5), breaks = c(0.5, 1.5), labels = c(expression(paste("P"[2], " Allele")), expression(paste("P"[1], " Allele"))), 
                                                 colours = c("#0015ff", "#3a4aff", "#7f89ff", "#ff7f7f", "#ff4444", "#ff0000"),
                                                 name = "Ratio") + 
    geom_hline(yintercept = 1) + #geom_vline(xintercept = 97.16) +
    labs(title = "Hypothetical cross with no distortion",
             x = "Physical distance along chromosome (%)", 
             y = expression(paste("Segregation ratio (P"[1], "/ P"[2], ")"))) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0.5, 1.5)) + theme_bw(base_size = 15)



pdf("plots/swbioconference/seg.dist14.pdf", 15, 7)
grid.arrange(no.seg.dist.plot, gen.plot.wo.colour.groups(map.number = 4, chromo.name = "5A", return.df = F), ncol = 1)
dev.off()

#find out what chromosomes are in which maps
map.chromos = lapply(all.m.8.iwgsc.4b.rev, names)

library(gridExtra)


gen.plot2 = function(chromo.name){
    #only make plots for those maps that contain that particular chromosome
    g = which(unlist(lapply(map.chromos, function(x) chromo.name %in% x)))
    g1 = lapply(g, function(x){
        gen.plot(x, chromo.name)
    })
    
    #render the plots 
    do.call(grid.arrange, g1)
} 

lapply(listofwheatchromosomes, function(x){
    pdf(p("rotation1scripts_v4/plots/seg.dist.all.maps.comp/", x, "comp.pdf"), 10, 10)
    gen.plot2(x)
    dev.off()    
})



#     ____________________________________________________________________________
#     FILTER GENES FROM FENG ET AL                                                                                        ####

#MAKE SURE FUNCTIONS IN DEFINE FUNCTIONS SECTION ARE LOADED BEFORE RUNNING THIS CODE

#load gene datasets for both IWGSC assemblies
iwgsc.css.genes = readDNAStringSet("genome.assemblies.v2/iwgsc.2.25/survey_sequence_gene_models_MIPS_v2.2_Jul2014/ta_IWGSC_MIPSv2.2_HighConf_CDS_2014Jul18.fa")

iwgsc.ref.genes = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")

#grab the gene annotations from the region on 6A that shows common distortion between a.x.c, o.x.s and cs.x.p
v(gen.plot(4, "6A", T))
v(gen.plot(4, "5A", T))
gen.plot(4, "5B", F)

#grab.genes.and.functions() should not be used as it filters the functional annotations table, which does not
#contain annotations for all genes, and therefore could miss some genes out

# genes1 = grab.genes.and.functions("6A", 14354024, 61036697) #from sacha.consensus.map.processing.R
# genes1.2 = grab.genes.and.functions("6A", 30032847, 40686978) #from sacha.consensus.map.processing.R
# genes2 = grab.genes.and.functions("5A", 688783153, 691034248) #from sacha.consensus.map.processing.R

genes.fil.5a = filter(genes, V1 == "chr5A" & V4 > 688783153 & V4 < 691034248)

genes.fil.6a = filter(genes, V1 == "chr6A" & V4 > 34032847 & V4 < 36686978)





feng.6a.genes = grab.feng.genes("6A")
css.genes = filter.iwgsc.genes(iwgsc.css.genes, feng.6a.genes)

ref.genes = filter.iwgsc.genes(iwgsc.ref.genes, genes.fil.6a$gene.names)


writeXStringSet(iwgsc.css.genes2, "rotation1scripts_v4/processed_data/fasta/inflorescense.transcriptome.genes/feng.genes.6a.iwgsc.css.fa")

writeXStringSet(genes.6a, "rotation1scripts_v4/processed_data/fasta/inflorescense.transcriptome.genes/6a.region.of.interest.narrow.genes.refseq.fa")


#     ____________________________________________________________________________
#     PROCESS BLAST                                                                                                                     ####
#did a blast search of genes found in inflorescense transcriptome vs genes in my region
#of interest on 6A

inflo.blast = read.blast("bioinf/blast/inflo.transcriptome/results.blast/feng.genes.vs.refseq.blast")

#check 5a BLAST 
inflo.blast2 = read.blast("bioinf/blast/inflo.transcriptome/results.blast/feng.genes.vs.refseq.5a.blast")

#second 5a candidate BLAST
inflo.blast3 = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/TraesCS5A01G531700.genomic.vs.para.blast")

inflo.blast4 = read.blast("bioinf/blast/inflo.transcriptome/results.blast/feng.genes.vs.refseq.6a.narrow.blast")


#     ____________________________________________________________________________
#     5A seg genes                                                                                                                        ####

fivea.seg.genes = iwgsc.ref.genes[match(genes2$Gene.ID, multi.str.split(names(iwgsc.ref.genes), " ", 1))]

writeXStringSet(fivea.seg.genes, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/5a.seg.dist.genes.fa")

f.blast = read.blast("bioinf/blast/refseq.genes.vs.css.genes/results.blast/5a.seg.refseq.vs.css.blast")
f.blast5a = f.blast[which(multi.str.split(as.character(f.blast$sseqid), "\\_", 2) == "5AL"), ]


#grab genomic sequences for 5a seg dist genes
all.genomic.sequences = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/iwgsc.full.gene.sequences.fa")

genomic.seq.names = multi.str.split(multi.str.split(names(all.genomic.sequences), "=", 2), ";", 1)

fivea.seg.gene.names = multi.str.split(multi.str.split(names(fivea.seg.genes), " ", 1), "\\.", 1)

all.genomic.sequences2 = all.genomic.sequences[match(fivea.seg.gene.names, genomic.seq.names)]

all.genomic.sequences3 = all.genomic.sequences2[c(25, 21)]

writeXStringSet(all.genomic.sequences2, "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/5a.seg.dist.genes.genomic.fa")

#do some automation for searching the TILLING website
#1. write sequences to a text file one at a time for easy copy and paste to the site

g1 = newdf("seq.name", "seq", "notes", no.rows = T)
write_excel_csv(g1, "rotation1scripts_v4/processed_data/TILLING/tilling5asearch.csv")


for(i in 1:length(all.genomic.sequences2)){
    
    g1 = read_csv("rotation1scripts_v4/processed_data/TILLING/tilling5asearch.csv")
    
    writeXStringSet(all.genomic.sequences2[i], "rotation1scripts_v4/processed_data/fasta/tilling.seg.dist.gene.candidates/5a.seg.dist.genes.genomic.single.fa")
    
    g1 = add_row(g1)
    g1$seq.name[i] = names(all.genomic.sequences2)[i]
    g1$seq[i] = as.character(all.genomic.sequences2[i])
    write_excel_csv(g1, "rotation1scripts_v4/processed_data/TILLING/tilling5asearch.csv")
    
    
    
    
    
}




#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####

grab.feng.genes = function(chromo1){
    #returns character vector of gene names for a particular chromosome from feng et al., study; in IWGSC CSS format
    #args:
    #chromo1: character string of chromosome to search for
    chromo2 = p(chromo1, "L|", chromo1, "S")
    
    library(readxl)
    
    feng.data = lapply(1:4, function(x){
        as.data.frame(read_excel("rotation1scripts_v4/original_data/inflorescence.transcriptome/feng.et.al.2017.supp.data/PP2017-RA-00310DR1_Supplemental_Dataset_1.xlsx", x))
    })
    
    feng.data2 = lapply(feng.data, function(x){
        x1 = x[grep(chromo2, x[, 1]), ]
        x1
    })
    
    feng.genes = lapply(feng.data2, function(x){
        x[, 1]
    })
    
    feng.genes2 = do.call(c, feng.genes)
    feng.genes2
}

filter.iwgsc.genes = function(gene.data.set, list.of.gene.names){
    #returns a subset of a DNAStringSet object by a list of gene names
    #args:
    #gene.data.set: a DNAStringSet object
    #list.of.gene.names: a character vector of gene names to filter the DNAStringSet object by
    
    library(Biostrings)
    
    
    gene.coords = unlist(lapply(list.of.gene.names, function(x){
        grep(x, names(gene.data.set))
    }))
    
    filtered.gene.data.set = gene.data.set[gene.coords]
    filtered.gene.data.set
}


#     ____________________________________________________________________________
#     HOMEOLOGUE BLASTS                                                                                                             ####


g = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/TraesCS5A01G531300.vs.genome.blast")
g1 = read.blast("bioinf/blast/probe.vs.genome.blast/results.blast/TraesCS5A01G531700.vs.genome.blast")

