
#     ____________________________________________________________________________
#     INITIAL PROCESSING                                                                                                            ####


setwd("C:/Users/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

# #load sacha maps
# cs.x.p.map = read.csv("C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/Allen Characterization of Wheat breeders array suppl. docs/allen.cs.x.p.map.csv", header = T, stringsAsFactors = F)
# #remove unwanted columns
# cs.x.p.map = cs.x.p.map[, 1:3]
# 
# #going to use Sacha's cs.x.p map due to problems with consensus
# 
# phys.distances.cs.x.p = lapply(listofwheatchromosomes, function(x){
#     g = cs.x.p.map[cs.x.p.map$chr == x, 1]
#     genetictophysical(x, g, "-", 30, F)
# })
# 
# phys.distances.cs.x.p.bp = lapply(listofwheatchromosomes, function(x){
#     g = cs.x.p.map[cs.x.p.map$chr == x, 1]
#     genetictophysical(x, g, "-", 30, T)
# })
# 
# names(phys.distances.cs.x.p) = listofwheatchromosomes
# 
# cs.x.p.map$phys.dist.bp = unlist(phys.distances.cs.x.p.bp)
# 
# #load seg. dist. markers of cs.x.p cross from seg.dist.marker.analysis.R
# load(file = "rotation1scripts_v4/saved.objects/cs.x.p.seg.dist.markers2")
# 
# #change chromosome format from "chr1A" to "1A"
# cs.x.p.seg.dist.markers2[, 7] = substr(cs.x.p.seg.dist.markers2[, 7], 4, 5)

# s(cs.x.p.map, "rotation1scripts_v4/saved.objects/cs.x.p.map", "sacha.consensus.map.processing.R")

#     ____________________________________________________________________________
#     GENE DENSITY                                                                                                                        ####


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
#     FURTHER PROCESSING                                                                                                            ####

# cs.x.p.map.full = read.csv("C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/Allen Characterization of Wheat breeders array suppl. docs/allen.cs.x.p.map.csv", header = T, stringsAsFactors = F)
# 
# cs.x.p.map.full = t(cs.x.p.map.full)
# 
# cs.x.p.map.full.correct.format = read.csv("rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/sacha.genotypes.csv", header = T, stringsAsFactors = F)
# 
# sacha.seg.dist.markers2 = grab.seg.dist.markers("rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/sacha.genotypes.csv", p.threshold = 0.05)
# 
# sacha.genotype.data = read.csv("rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/sacha.genotypes.csv", header = T, stringsAsFactors = F)
# 
# all.sacha.counts = lapply(sacha.genotype.data[, 3:ncol(sacha.genotype.data)], geno.count.chi)
# all.sacha.counts = t(as.data.frame(all.sacha.counts))
# 
# rownames(all.sacha.counts) = switch.affy.format(rownames(all.sacha.counts))
# 
# cs.x.p.map = cs.x.p.map[-1, ]
# 
# all.sacha.counts = all.sacha.counts[match(rownames(all.sacha.counts), cs.x.p.map$ï..marker), ]
# 
# cs.x.p.map.w.genocounts = cbind(cs.x.p.map, all.sacha.counts)

s(cs.x.p.map.w.genocounts, "rotation1scripts_v4/saved.objects/cs.x.p.map.w.genocounts", "sacha.consensus.map.processing.R")

#     ____________________________________________________________________________
#     EXAMINE ENTIRE CHROMOSOME COUNTS                                                                                ####

ind.chromo.stats = lapply(listofwheatchromosomes, function(x){
    g = filter(cs.x.p.map.w.genocounts, chr == x)
    n.A.bigger = length(which(g$A > g$B))
    total = nrow(g)
    all = c(x, total, n.A.bigger)
    names(all) = c("chromosome", "total.markers", "n.A.bigger")
    return(all)
})

ind.chromo.stats = as.data.frame(t(as.data.frame(ind.chromo.stats)))
ind.chromo.stats = reset.rownames(ind.chromo.stats)
ind.chromo.stats = add_column(ind.chromo.stats, per.A.big = "")
ind.chromo.stats$per.A.big = (as.numeric(as.character(ind.chromo.stats$n.A.bigger)) / as.numeric(as.character(ind.chromo.stats$total.markers)))*100

ind.chromo.stats = ind.chromo.stats[sort(as.numeric(ind.chromo.stats$per.A.big), index.return = T)$ix, ]

v(ind.chromo.stats)


cs.x.p.map.w.genocounts$seg.dist = ""
#mark segregation distortion markers in map dataframe with new column "seg.dist"
for(x in listofwheatchromosomes){
    cs.x.p.map.w.genocounts$seg.dist[match(filter(probe.blast3, V2 == p("chr", x))$V1, cs.x.p.map.w.genocounts$ï..marker)] = "T"    
}

cs.x.p.map.w.genocounts$linear.count = ""
for(i in unique(cs.x.p.map.w.genocounts$chr)){
    cs.x.p.map.w.genocounts[which(cs.x.p.map.w.genocounts$chr == i), ]$linear.count = 1:nrow(filter(cs.x.p.map.w.genocounts, chr == i))
}
cs.x.p.map.w.genocounts$linear.count = as.numeric(cs.x.p.map.w.genocounts$linear.count)
cs.x.p.map.w.genocounts$linear.count = as.factor(cs.x.p.map.w.genocounts$linear.count)

cs.x.p.map.w.genocounts$seg.dist[which(cs.x.p.map.w.genocounts$`p-value` < 0.05)] = "T" #tag markers exhibiting significant segregation distortion

library(ggplot2)

#plot segregation of markers as a function of physical distance along chromosome
#setup a closure so the dataframe supplied to plot.chromo can be customized
plot.chromo.df = function(geno.dataframe, physical){
    if(missing(physical)) physical = F
    geno.df = geno.dataframe
    plot.chromo = function(chromosome){
        g = filter(geno.df, chr == chromosome)
        g$ratio = ((g$A/(g$A + g$B))*100)
        if(physical == F){
        ggplot(g, aes(x = linear.count, y = ratio, color = seg.dist, shape = seg.dist)) + 
            geom_point() + geom_hline(yintercept = 50) + scale_y_continuous(breaks = seq(0, 100, 5)) +
            scale_x_discrete(breaks = seq(-50, (nrow(g)+50), 50), labels = seq(-50, (nrow(g)+50), 50)) +
            coord_cartesian(xlim = c(-50, (nrow(g)+50)), ylim = c(40, 60)) +
            scale_color_manual(values = c("#000000", "#ff0000"))
        } else {
            ggplot(g, aes(x = phys.dist.bp, y = ratio, color = seg.dist, shape = seg.dist)) + 
                geom_point() + geom_hline(yintercept = 50) + scale_y_continuous(breaks = seq(0, 100, 5)) +
                scale_color_manual(values = c("#000000", "#ff0000"))
        }
        # plot(g$phys.dist.bp, (g$A/(g$A + g$B))*100)
        # abline(h = 50)
    }
}

plot.chromo = plot.chromo.df(cs.x.p.map.w.genocounts)

cs.x.p.map.physically.ordered = cs.x.p.map.w.genocounts

cs.x.p.map.phys.ord.2 = newdf(colnames(cs.x.p.map.physically.ordered), no.rows = T)
for(i in unique(cs.x.p.map.physically.ordered$chr)){
    g = cs.x.p.map.physically.ordered[cs.x.p.map.physically.ordered$chr == i, ]
    g = g[-which(is.na(g$phys.dist)), ]
    g = g[sort(g$phys.dist, index.return = T)$ix, ]
    g$linear.count = 1:nrow(g)
    cs.x.p.map.phys.ord.2 = rbind(cs.x.p.map.phys.ord.2, g)    
    
}

plot.chromo.phys = plot.chromo.df(cs.x.p.map.phys.ord.2, T)

#     ____________________________________________________________________________
#     EXAMINE RECOMBINATION EVENTS                                                                                        ####

#Here I'm making a new genotyping dataframe with the marker order specified by physical order 
#This will allow me to match up the recombination histograms with my physically ordered seg. dist. plots
cs.x.p.geno.phys.ordered = cs.x.p.map.full.correct.format[, c(1, 2, (match(switch.affy.format(rownames(cs.x.p.map.phys.ord.2)), colnames(cs.x.p.map.full.correct.format)[3:ncol(cs.x.p.map.full.correct.format)]) + 2))]

cs.x.p.geno.phys.ordered[1, 1] = ""
write.csv(cs.x.p.geno.phys.ordered, "rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/sacha.genotypes.phys.ordered.csv", row.names = F)

cs.x.p.geno.list = makegenotypelist(cs.x.p.map.full.correct.format)
cs.x.p.geno.list.phys = makegenotypelist(cs.x.p.geno.phys.ordered)

s(cs.x.p.geno.list, "rotation1scripts_v4/saved.objects/cs.x.p.geno.list", "sacha.consensus.map.processing.R")
s(cs.x.p.geno.list.phys, "rotation1scripts_v4/saved.objects/cs.x.p.geno.list.phys", "sacha.consensus.map.processing.R")

library(gridExtra)

histlists = lapply(listofwheatchromosomes, function(x){
    g = makehistlist(x, cs.x.p.geno.list, return_recombinationcounts = T, homozygous = T)
    return(g)
})

names(histlists) = listofwheatchromosomes

histlists.phys = lapply(listofwheatchromosomes, function(x){
    g = makehistlist(x, cs.x.p.geno.list.phys, return_recombinationcounts = T, homozygous = T)
    return(g)
})

names(histlists.phys) = listofwheatchromosomes

for(i in listofwheatchromosomes){
    g = filter(cs.x.p.map.phys.ord.2, chr == i)
    h1 = histlists.phys[[i]]
    linear.to.phys = g$phys.dist.bp[match(h1$position, g$linear.count)]
    h1$phys.pos = linear.to.phys
    histlists.phys[[i]] = h1
}



plot.side.by.side = function(chromo, phys){
    #phys: boolean - if T, use physical ordering of markers
    if(missing(phys)) phys = F
    if(phys == T){
        hist = histlists.phys
        plot1 = plot.chromo.phys(chromo)
        plot2 = qplot(sort(as.numeric(unique(cs.x.p.histlist2$phys.pos))), geom = "histogram", binwidth = 1000)
    } else {
        hist = histlists
        plot1 = plot.chromo(chromo)
        plot2 = qplot(sort(as.numeric(unique(cs.x.p.histlist2$position))), geom = "histogram", binwidth = 4)
        plot4 = plot.chromo.phys(chromo)
    }
    cs.x.p.histlist2 = hist[[chromo]]
    
    plot3 = plot.introgression(chromo)
    
    
    grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2, top = p("Segregation distortion: ", chromo))
}


plot.side.by.side("1A", phys = F)

#print all chromosomes
lapply(listofwheatchromosomes, function(x){
    ggsave(p("rotation1scripts_v4/plots/4x4seg.dist.plots/", x, ".pdf"), plot.side.by.side(x, phys = F), width = 10, height = 10)    
})


#     ____________________________________________________________________________
#     ASMAP                                                                                                                                     ####
#need to reestimate genetic distances for physical ordering of markers in order to get cM/Mb graph

library(ASMap)
library(qtl)
sacha.cross = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/sacha.genotypes.phys.ordered.csv", genotypes = c("A", "B"), estimate.map = F)
sacha.cross = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/processed_data/genotypes/chinese.spring.x.paragon/sacha.genotypes.csv", genotypes = c("A", "B"), estimate.map = F)
sacha.map.reesti = quickEst(sacha.cross, map.function = "kosambi")
sacha.map.reesti2 = pull.map(sacha.map.reesti, as.table = T)


#the above approach doesn't seem to be working, genetic distances heavily inflated with physical order

#     ____________________________________________________________________________
#     ANALYSIS OF CHROMOSOME 2A                                                                                             ####

chr2a = filter(cs.x.p.map.w.genocounts, chr == "2A")
chr2a.no.na = chr2a[-which(is.na(chr2a$phys.dist.bp)), ]
chr2a.phys.sort = chr2a.no.na[sort(chr2a.no.na$phys.dist.bp, index.return = T)$ix, ]
chr2a.phys.sort$linear.count = as.factor(1:nrow(chr2a.phys.sort))

plot.chromo.phys = plot.chromo.df(chr2a.phys.sort)
plot.chromo.phys("2A")

genes2a.seg.dist.region = check.gene.density("2A", 767792030, 793811875, T)
genes2a.seg.dist.region$V9 = as.character(genes2a.seg.dist.region$V9)
genes2a.seg.dist.region$gene = multi.str.split(genes2a.seg.dist.region$V9, "\\;", 1)
genes2a.seg.dist.region$gene = multi.str.split(genes2a.seg.dist.region$gene, "\\=", 2)

all.full.genes = ReadFasta("rotation1scripts_v4/original_data/IWGSC/iwgsc.full.gene.sequences.fa")


#     ____________________________________________________________________________
#     PROCESS BLAST                                                                                                                     ####

gen.blast.script("sacha.seg.dist.genes.full.genomic.sequence.2a.v2.fa", "genes.vs.paragon.genome", culling_limit = 1)

seg.2a.blast = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/sacha.seg.dist.genes.full.genomic.sequence.2a.v2.fa.blast")

seg.2a.blast.v2 = grab.best.groups.of.hits(seg.2a.blast)

seg.2a.blast.v3 = select.best.groups.of.hits(seg.2a.blast.v2)
seg.2a.blast.v3$coverage = seg.2a.blast.v3$coverage + 1 #I think bedtools or samtools listed the genelengths as one shorter than they actually are
seg.2a.blast.v3$per.cov = (seg.2a.blast.v3$coverage / seg.2a.blast.v3$gene.length)*100 #update per.cov to reflect new gene lengths

seg.unique.v3 = remove.100.percent.cov.and.ident.hits(seg.2a.blast.v3) # this isn't valid as there could be mismatches

s(seg.unique.v3, "rotation1scripts_v4/saved.objects/seg.unique.v3", "sacha.consensus.map.processing.R")

#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####

remove.100.percent.cov.and.ident.hits = function(df){
    g = unique(df[, 1:2])
    g2 = newdf(colnames(df), no.rows = T)
    for(i in 1:nrow(g)){
        q = filter(df, qseqid == g[i, 1] & sseqid == g[i, 2])
        #complicated logic gate follows:
        #if both the coverage and the identity are 100%, for the same hit,
        #exclude that hit from the final dataframe
        per.cov.check = which(100 %in% as.numeric(q$per.cov))
        per.ident.check = which(100 %in% as.numeric(q$percentage.identical))
        
        if(length(per.cov.check) == 0 | length(per.ident.check) == 0){
            both.check = F
        } else {
            both.check = all(per.cov.check == per.ident.check)
        }
        
        if(both.check == F) {
            g2 = rbind(g2, q)
        }
    }
    return(g2)
}

#grab best groups
select.best.groups.of.hits = function(df){
    #to be used after grab.best.groups.of.hits();
    #returns best hits by analysing the % coverage of the query sequence
    g = df
    
    newg = newdf(colnames(g))
    
    #grab top two groups of hits by % coverage of the query
    for(i in unique(g$qseqid)){
        q = g[g$qseqid == i, ]
        q$range = abs(q$qend - q$qstart)
        
        best.per.cov = q[which(q$per.cov %in% unique(sort(q$per.cov, decreasing = T))[1]), ]
        
        
        
        
        newg = rbind(newg, best.per.cov)
    }
    
    newg = newg[-1, ]
    
    newg[, 3:ncol(newg)][] = lapply(newg[, 3:ncol(newg)], as.numeric) #convert some columns to numeric
    return(newg)
}

clean.within.q.s.pairs = function(df){
    un = unique(df[, 1:2])
    for(i in 1:nrow(un)){
        g = filter(df, qseqid == un[1, 1] & sseqid == un[1, 2]) #grab query-subject pair from main df
        if(nrow(g) > 1){ #don't proceed if only 1 row
            g$range = abs(g$qend - g$qstart)
            biggest.range.pos = which.max(g$range)
            if(length(biggest.range.pos) == 1){ #don't proceed if two hits have the same range
                
            }
            
            
        }
    }
}


