#     ____________________________________________________________________________
#     HISTORICAL RECOMBINATION MAPPING                                                                                 ####

setwd("E:/phd.project.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
load("rotation1scripts_v4/saved.objects/cs.x.p.map.w.genocounts")

cs.x.p.map.w.genocounts2 = split(cs.x.p.map.w.genocounts, f = cs.x.p.map.w.genocounts$chr)

attributes(cs.x.p.map.w.genocounts)



# --- functions

coords.of.skeleton.markers = function(cM.vec){
    match(unique(cM.vec), cM.vec)    
}

# --------- code

#within bins, order cs.x.p.map.w.genocounts by physical sequence
unique.chromo.cm.combos = unique(cs.x.p.map.w.genocounts[, 2:3])
for(i in 1:nrow(unique.chromo.cm.combos)){
    g = cs.x.p.map.w.genocounts[which(cs.x.p.map.w.genocounts$chr == unique.chromo.cm.combos[i, 1] & cs.x.p.map.w.genocounts$cM == unique.chromo.cm.combos[i, 2]), ]
    sorted = sort(g$phys.dist.bp, na.last = F)
    g = g[match(sorted, g$phys.dist.bp), ]
    cs.x.p.map.w.genocounts[which(cs.x.p.map.w.genocounts$chr == unique.chromo.cm.combos[i, 1] & cs.x.p.map.w.genocounts$cM == unique.chromo.cm.combos[i, 2]), ] = g
    
}


#     ____________________________________________________________________________
#     FIX FLIPPED CHROMOSOMES                                                                                                 ####


#check if any chromosomes are flipped

fix.flipped.chromo = function(df.to.check){
    flipped.chromosomes = lapply(unique(df.to.check$chr), function(x){
        print(x)
        g = filter(df.to.check, chr == x)
        start.mean = mean(na.omit(g$phys.dist.bp[1:20]))
        end.mean = mean(na.omit(g$phys.dist.bp[(nrow(g)-20):nrow(g)]))
        if(is.na(start.mean)) browser()
        if(is.na(end.mean)) browser()
    
        if(start.mean > end.mean){
            return(x)
        } else{
            return()
        }
    })
    
    flipped.chromosomes = unlist(flipped.chromosomes)
    print(p(flipped.chromosomes, " are flipped"))
    
    #reverse selected chromosomes, reset centiMorgan values
    flipped.chromosome.dfs = lapply(flipped.chromosomes, function(x){
        g = filter(df.to.check, chr == x)
        g = g[nrow(g):1, ]
        
        #reverse cm distance
        total.cm = max(as.numeric(g$cM))
        
        g$cM = abs(g$cM - total.cm)
        
        return(g)
    })
    
    names(flipped.chromosome.dfs) = flipped.chromosomes
    
    #replace flipped chromosomes with newly inverted ones
    for(i in names(flipped.chromosome.dfs)){
        df.to.check[which(df.to.check$chr == i), ] = flipped.chromosome.dfs[[i]]
    }
    
    return(df.to.check)

}

cs.x.p.map.w.genocounts = fix.flipped.chromo(cs.x.p.map.w.genocounts)

extract.monotonic.chromo = function(chromo, fil.df){
    if(missing(fil.df)) fil.df = cs.x.p.map.w.genocounts
    g = filter(fil.df, chr == chromo)
    
    g[match(na.omit(g$phys.dist.bp)[longest_subseq.R(na.omit(g$phys.dist.bp))], g$phys.dist.bp), ]
}

g2a = extract.monotonic.chromo("4B")

#genotyping information should now be extracted with historical.recombination.mapping.server.script.unix.R; supply chromosome as argument

#     ____________________________________________________________________________
#     LOAD GENOTYPING DATA / PREPARE FOR PHASE 2.1.1                                                    ####

prepare.for.phase = function(chromosome.name, row.numbers, descriptor, zoom, output.data.frame){
    #args:
    #chromosome.name: string indicating chromosome name e.g. 1A
    #row.numbers: numeric vector of row coordinates with which to subset entire dataframe for varieties
    #descriptor: string indicating e.g. region and line type (watkins or elite)
    #zoom: boolean indicating whether or not this is the "zoom" analysis
    #output.data.frame: boolean, if T, don't write anything and return dataframe in non-PHASE format
    
    if(missing(zoom)) zoom = F
    if(missing(output.data.frame)) output.data.frame = F
    
    if(zoom == F){
        first = read_csv(p("rotation1scripts_v4/processed_data/genotypes/historical.recombination.mapping/", chromosome.name, ".ske.incl.csv"))
        g1a = extract.monotonic.chromo(chromosome.name)
    } else {
        first = read_csv(p("rotation1scripts_v4/processed_data/genotypes/historical.recombination.mapping/zoom/", chromosome.name, ".zoom.csv"))
        load("rotation1scripts_v4/saved.objects/sorted.chromosomes2")
        g1a = filter(sorted.chromosomes2, sseqid == p("chr", chromosome.name))[1:100, ]
    }
    first = reset.colnames(first)
    
    library(reshape2)
    # browser()
    
    # first2 = read_csv(p("rotation1scripts_v4/original_data/historical.recombination.analysis/", chromosome.name, ".csv"))
    # 
    # first2.cast = dcast(first, V2~V3)
    # first2.cast2 = first2.cast[row.numbers, ]
    # 
    first.cast = dcast(first, V2~V3)
    first.cast = first.cast[row.numbers, ]
    
    # browser()
    
    first.cast = add_row(first.cast, .before = 1)
    # browser()
    if(zoom == F) first.cast[1, 2:ncol(first.cast)] = g1a$phys.dist.bp
    if(zoom == T) first.cast[1, 2:ncol(first.cast)] = g1a$sstart
    
    if(output.data.frame == T){
        return(first.cast)
    } else {
        #assign variables for PHASE 2.1.1 software format
        no.ind = (nrow(first.cast) - 1)
        no.loci = (ncol(first.cast) - 1)
        
        positions = paste(unname(unlist(c("P", first.cast[1, 2:ncol(first.cast)]))), collapse = " ")
        locustype = paste(rep("S", (ncol(first.cast) - 1)), collapse = "")
        
        test = lapply(first.cast[2:nrow(first.cast), 1], function(x){
            g = unlist(unname(first.cast[which(first.cast[, 1] == x), 2:ncol(first.cast)]))
            row1 = as.numeric()
            row2 = as.numeric()
            
            lapply(g, function(y){
                if(y == "AB"){
                    row1 = c(row1, 0)
                    row2 = c(row2, 1)
                }
                
                if(y == "AA"){
                    row1 = c(row1, 0)
                    row2 = c(row2, 0)
                }
                
                if(y == "BB"){
                    row1 = c(row1, 1)
                    row2 = c(row2, 1)
                }
                
                if(y == "NoCall"){
                    row1 = c(row1, "?")
                    row2 = c(row2, "?")
                }
                
                list(row1, row2)
                
            })
        }
        )
        
        test2 = lapply(test, function(q){
            row1 = unlist(lapply(q, function(x) x[[1]]))
            row2 = unlist(lapply(q, function(x) x[[2]]))
            list(row1, row2)
        })
        
        make.counter = function(){
            counter = 0
            ret.count = function(){
                counter <<- counter + 1
                counter
            }
        }
        
        first.counter = make.counter()
        
        test3 = lapply(test2, function(x) list(p("#", first.counter()), x[[1]], x[[2]]))
        
        all.test = list()
        for(i in 1:length(test3)){
            all.test = c(all.test, test3[[i]])
        }
        
        all.test2 = lapply(all.test, paste, collapse = " ")
        
        
        to.write = c(no.ind, no.loci, positions, locustype, all.test2)
        
        
        if(zoom == F){
            lapply(to.write, function(x){
                write(x, p("rotation1scripts_v4/processed_data/PHASE/phaseinput.", chromosome.name, ".", descriptor, ".100.ind.txt"), append = T)
            }
            )
        } else {
            lapply(to.write, function(x){
                write(x, p("rotation1scripts_v4/processed_data/PHASE/zoom/", chromosome.name, ".", descriptor, ".zoom.100.ind.txt"), append = T)
            }
            )
        }
    }
}

#     ____________________________________________________________________________
#     PARSE PHASE OUTPUT / ANALYSE ZOOM                                                                             ####

parse.phase.output = function(phase.output.path){
    g = readLines(phase.output.path)
    g1 = lapply(g, strsplit, split = " ")
    g2 = lapply(g1, function(x) as.numeric(unlist(x)))
    
    g3 = reset.rownames(as.data.frame(t(as.data.frame(g2))))
    
    g4 = g3[2:nrow(g3), ]
    g4 = reset.rownames(g4)
    
    g5 = unname(unlist(lapply(g4, mean)))
    
    phys.pos = as.numeric(g3[1, ])
    
    historical.recomb.df = data.frame(phys.pos, g5)
    return(historical.recomb.df)
}

phase.path = "rotation1scripts_v4/processed_data/PHASE/zoom/"



recomb.files = list.files(phase.path)[grep("recom", list.files(phase.path))]

pull.chromo.data = function(chromo){
    onea.files = recomb.files[grep(chromo, recomb.files)]
    
    onea.data = lapply(onea.files, function(x){
        hist.recomb.df = parse.phase.output(p(phase.path, x))
    })
    
    onea.data2 = lapply(onea.data, function(x){
        x$phys.pos.diff = calculate.cm.diff(x$phys.pos)
        x$lambda.diff = calculate.cm.diff(x$g5)
        return(x)
    })
    
    names(onea.data2) = onea.files
    
    return(onea.data2)
}


plot.zoom = function(data.to.plot){
    count1 = make.counter()
    plots = lapply(data.to.plot, function(q){
        ggplot(data = q, aes(x = phys.pos, y = g5)) + geom_point() + geom_line() + ggtitle(names(data.to.plot)[count1()])
    })
    
    library(gridExtra)
    grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]])
}

oneb.data = pull.chromo.data("1A")

all.hist.data.zoom = lapply(listofwheatchromosomes, pull.chromo.data)

#     ____________________________________________________________________________
#     PCA PLOT                                                                                                                                ####

stack.zoom.data = function(list.of.zoom.dfs){
    list.of.zoom.dfs.stack = newdf(colnames(as.data.frame(t(list.of.zoom.dfs[[1]]))), no.rows = T)
    for(i in 1:length(list.of.zoom.dfs)){
        temp.one = as.data.frame(t(list.of.zoom.dfs[[i]]))
        temp.one = temp.one[2, ]
        list.of.zoom.dfs.stack = rbind(list.of.zoom.dfs.stack, temp.one)
    }
    return(list.of.zoom.dfs.stack)
}

all.hist.data.zoom.stack = lapply(all.hist.data.zoom, stack.zoom.data)
all.hist.data.zoom.stack2 = lapply(all.hist.data.zoom.stack, function(x){
        rownames(x) = wheat.groups
        return(x)
    })
all.stack2 = all.hist.data.zoom.stack[[1]]
for(i in 2:length(all.hist.data.zoom.stack)){
    all.stack2 = cbind(all.stack2, all.hist.data.zoom.stack[[i]])
}

all.stack.var = unlist(lapply(all.stack2, function(x){ 
    sd(x)
}))
all.stack.var = sort(all.stack.var, decreasing = T)[10:length(all.stack.var)]


ggplot(as.data.frame(all.stack.var), aes(x = all.stack.var)) + geom_histogram(binwidth = 0.05)





all.stack2.log = log(all.stack2)
all.stack2.log = reset.colnames(all.stack2.log)
all.stack2.pca = prcomp(all.stack2.log, center = T, scale. = T)

library(devtools)
install_github("ggbiplot", "vqv")

substr(names(oneb.data), 1, 11)

wheat.groups = c("e.aus", "e.eur", "wat.asia", "wat.eu.e", "wat.eu.w", "wat.mid.e")

library(ggbiplot)
g <- ggbiplot(all.stack2.pca, obs.scale = 1, var.scale = 1,
                            circle = TRUE, var.axes = F, groups = c("e.aus", "e.eur", "wat.as", "wat.eu.e", "wat.eu.w", "wat.mid.e"))
g <- g + scale_shape_discrete(name = '') + geom_point(aes(shape = wheat.groups), size = 4)
g <- g + theme(legend.direction = 'horizontal', 
                             legend.position = 'top')
print(g)


plot.zoom(oneb.data)

lambda.diffs = data.frame(lapply(oneb.data, function(x) x[, 4]))

g = apply(lambda.diffs[, ], 1, function(x){
    all(x > 0) | all(x < 0) | all(x == 0)
})

which(lambda.diffs$X3A.watkins.eur.west.row.coords.zoom.100.ind.output_recom == max(lambda.diffs$X3A.watkins.eur.west.row.coords.zoom.100.ind.output_recom))

grab.orig.data = function(chromosome, dataset.num, start.marker, end.marker){
    chromo.coord = which(chromosome == listofwheatchromosomes)
    all.dat[[chromo.coord]][[dataset.num]][, start.marker:end.marker]
}

grab.orig.data("3A", 1, 21, 23)


q = apply(all.dat[[7]][[1]][, 23:25], 1, function(x){
    p(x[1], x[2], x[3])
})

q2.1 = apply(all.dat[[7]][[5]][, 23:24], 1, function(x){
    p(x[1], x[2])
})

q2.2 = apply(all.dat[[7]][[5]][, 24:25], 1, function(x){
    p(x[1], x[2])
})

q2.3 = apply(all.dat[[7]][[5]][, 25:26], 1, function(x){
    p(x[1], x[2])
})


q3 = apply(all.dat[[7]][[4]][, 23:25], 1, function(x){
    p(x[1], x[2], x[3])
})

smallerg5 = lapply(oneb.data, function(x) sort(x$g5)[2])

which(oneb.data[[1]]$g5 == smallerg5[1])










plot(historical.recomb.df)

generate.gff3.for.sequence.extraction("1A", 1159750, 1174830, "high.recomb.region.1a", "rotation1scripts_v4/processed_data/sequence.extraction.gff.files/high1a.gff3")

g = ReadFasta("rotation1scripts_v4/processed_data/sequence.extraction.gff.files/extracted.sequences/high5a.fa")
g = convert.to.character.data.frame(g)
max_orf(g$sequence)


#     ____________________________________________________________________________
#     EXAMINE COs IN MAPPING POP                                                                                            ####

markerlist = colnames(first.cast)[2:ncol(first.cast)]

geno.1a = cs.x.p.map.full.correct.format[, c(1, which(cs.x.p.map.full.correct.format[1, ] == "1A"))]
geno.1a = geno.1a[-2, ]

geno.1a = geno.1a[, c(1, which(switch.affy.format(colnames(geno.1a)) %in% markerlist))]

genotypelist1a = makegenotypelist(geno.1a)
geno.histlist = makehistlist("1A", genotypelist1a, T, homozygous = T)

geno.histlist$phys = ""

for(i in 1:nrow(geno.histlist)){
    geno.histlist$phys[i] = phys.pos[as.numeric(geno.histlist$position[i])]
}


pdf("rotation1scripts_v4/plots/PHASE/1a.both.pdf")
par(mfrow=c(1, 2))
plot(phys.pos, g5)
hist(as.numeric(geno.histlist$phys), breaks = 100)
dev.off()

plot1.data = data.frame(phys.pos, g5)
plot2.data = as.data.frame(table(geno.histlist$phys))
plot2.data$Var1 = as.numeric(as.character(plot2.data$Var1))

library(ggplot2)
library(gridExtra)
plot1 = ggplot(plot1.data, aes(x = phys.pos, y = g5)) + geom_line()
plot2 = ggplot(plot2.data, aes(x = Var1, y = Freq)) + geom_line()
grid.arrange(plot1, plot2)



#     ____________________________________________________________________________
#     SQL TESTING                                                                                                                         ####



# mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

# sql.query = "SELECT * FROM 35k_array_oct2016_named WHERE Var_col = '117_Moulin'"

#pass mysql query to database
# rs = dbSendQuery(mydb, sql.query)
# rs2 = fetch(rs, n=-1)
# rs2 = cereals.axiom.table


#     ____________________________________________________________________________
#     GEOGRAPHICAL DATA                                                                                                             ####

library(readxl)
geography.df = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 3)
geo.df.watkins = read_excel("rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 2)

# s(first.cast, "rotation1scripts_v4/saved.objects/first.cast", "historical.recombination.mapping.R")
load("rotation1scripts_v4/saved.objects/first.cast")
g = multi.str.split(first.cast$V2, "'", 2)

coord.starting.w.digit = grep("^[[:digit:]]+\\_", g)

#only some of the accessions have names starting with digits, e.g. 93_robigus, I want to remove these digits. 
#other lines, such as the watkins lines, have digits that are useful and needed
g[coord.starting.w.digit] = multi.str.split(g[coord.starting.w.digit], "_", 2)

lines = data.frame(g)
watkins = data.frame(g[grep("Watkins", lines$g)])
row.names(watkins) = grep("Watkins", lines$g)
watkins$line = multi.str.split(as.character(watkins$g.grep..Watkins...lines.g..), "\\.", 1)

watkins$region = geo.df.watkins$Region[match(watkins$line, geo.df.watkins$`Watkins Number`)]

matches.coords = lapply(lines$g, function(x) grep(x, geography.df$Accession))
first.matches = unlist(lapply(matches.coords, function(x) x[1]))

lines$region = geography.df$Region[first.matches]
lines = lines[-grep("Watkins", lines$g), ] #remove watkins lines from elite lines

test.df.df = read_csv("rotation1scripts_v4/processed_data/genotypes/historical.recombination.mapping/1A.ske.incl.csv")
test.df.df = reset.colnames(test.df.df)
test.df.df = dcast(test.df.df, V2~V3)

table(lines$region)
elite.aus.row.coords = as.numeric(rownames(lines)[which(lines$region == "Australia")])[1:100]
elite.europe.row.coords = as.numeric(rownames(lines)[which(lines$region == "Europe")])[1:100]

table(watkins$region)
watkins.asia.row.coords = as.numeric(rownames(watkins)[which(watkins$region == "Asia")])[1:100]
watkins.eur.west.row.coords = as.numeric(rownames(watkins)[which(watkins$region == "Europe (West)")])[1:100]
watkins.eur.east.row.coords = as.numeric(rownames(watkins)[which(watkins$region == "Europe (East)")])[1:100]
watkins.mid.east.row.coords = as.numeric(rownames(watkins)[which(watkins$region == "Middle East")])[1:100]

row.coords.df = data.frame(elite.aus.row.coords, elite.europe.row.coords,
                                                     watkins.asia.row.coords, watkins.eur.east.row.coords,
                                                     watkins.eur.west.row.coords, watkins.mid.east.row.coords)


v(test.df.df[watkins.mid.east.row.coords, ])
v(test.df.df[watkins.asia.row.coords, ])

# s(row.coords.df, "rotation1scripts_v4/saved.objects/row.coords.df", "historical.recombination.mapping.R")

#generate phase files for all chromosomes for all regions etc
all.dat = lapply(listofwheatchromosomes, function(x){
    counter2 = make.counter()
    lapply(row.coords.df, function(y){
        # browser()
        prepare.for.phase(x, na.omit(y), colnames(row.coords.df[counter2()]), F, zoom = F) #set last argument here to F if want to generate PHASE files, otherwise output df
    })
})


#generate phase files for all chromosomes for all regions etc
all.dat = lapply(listofwheatchromosomes[11:21], function(x){
    counter2 = make.counter()
    lapply(row.coords.df, function(y){
        # browser()
        prepare.for.phase(x, na.omit(y), colnames(row.coords.df[counter2()]), F, zoom = F) #set last argument here to F if want to generate PHASE files, otherwise output df
    })
})


files1 = list.files("rotation1scripts_v4/processed_data/PHASE/", pattern = "100.ind.txt", full.names = T)
files1.1 = list.files("rotation1scripts_v4/processed_data/PHASE/", pattern = "100.ind.txt", full.names = F)
files2 = list.files("rotation1scripts_v4/processed_data/PHASE/watkins.and.elite.input/", pattern = "100.ind.txt", full.names = T)
files2.1 = list.files("rotation1scripts_v4/processed_data/PHASE/watkins.and.elite.input/", pattern = "100.ind.txt", full.names = F)

files1 = multi.str.split(files1, "phaseinput.", 2)
files1.1 = multi.str.split(files1.1, "phaseinput.", 2)

files1[which(!unlist(lapply(files1, file.size)) == unlist(lapply(files2, file.size)))]

files3 = data.frame(files1.1, files2.1)

# s(all.dat, "rotation1scripts_v4/saved.objects/all.dat", "historical.recombination.mapping.R")

#     ____________________________________________________________________________
#     EXAMINE ALL CHROMOSOMAL MAPS                                                                                        ####

obs.to.load = paste("rotation1scripts_v4/saved.objects/", list.files("rotation1scripts_v4/saved.objects/")[grep("g1a...", list.files("rotation1scripts_v4/saved.objects/"))], sep = "")

load(obs.to.load[1])

functional.anno = read_delim("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.TAB", "\t")
genomegff = read_delim("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3", "\t")
genomegff = reset.colnames(genomegff)
colnames(genomegff) = c("chr", "assembly", "ele.type", "bp.start", "bp.end", "V6", "V7", "V8", "gene.name")
genomegff$gene.name.trunc = genomegff$gene.name
genomegff$gene.name.trunc = multi.str.split(genomegff$gene.name.trunc, "=", 2)

genomegff$gene.name.trunc[grep("\\;", genomegff$gene.name.trunc)] = multi.str.split(genomegff$gene.name.trunc[grep("\\;", genomegff$gene.name.trunc)], "\\;", 1)

start = g1a$phys.dist.bp[1]
end = g1a$phys.dist.bp[2]
diff = abs(start - end)
g = filter(genomegff, chr == "chr1B" & bp.start > start & bp.end < end)

cum(g1a$phys.dist.bp)
?cummean


#     ____________________________________________________________________________
#     ZOOM IN - EXAMINE ALL PROBE BLAST                                                                             ####

all.probe.blast = read.blast("rotation1scripts_v4/original_data/allprobegenome.blast")
all.probe.blast2 = remove.hits.with.same.best.e.value(all.probe.blast)
all.probe.blast3 = grab.best.hits(all.probe.blast2)
all.probe.blast3$sstart = as.numeric(all.probe.blast3$sstart)


sorted.chromosomes = lapply(sort(unique(all.probe.blast3$sseqid)), function(x){
    chr.df = filter(all.probe.blast3, sseqid == x)
    chr.df = arrange(chr.df, sstart)
    chr.df$interval = c(0, diff(chr.df$sstart))
    chr.df$interval.kb = chr.df$interval / 1000
    return(chr.df)
})

sorted.chromosomes2 = newdf(colnames(sorted.chromosomes[[1]]), no.rows = T)
lapply(sorted.chromosomes, function(x){
    sorted.chromosomes2 <<- rbind(sorted.chromosomes2, x)
    return()
})

#save the object (processing steps take a long time)
# s(sorted.chromosomes2, "rotation1scripts_v4/saved.objects/sorted.chromosomes2", "historical.recombination.mapping.R")

hist(sorted.chromosomes2$interval.kb, breaks = 100000, xlim = c(0, 10)) #quite a lot of smaller intervals







#### 19/08/2019 NEW ANALYSIS ####

source("rotation1scripts_v4/scripts/r/recombination/recombination.distribution.functions.R")
csxp.map.orig = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.cs.x.p.map.csv", stringsAsFactors = F)

csxp.geno = prepare.allen.map(csxp.map.orig)
csxp.geno.nofilter = prepare.allen.map(csxp.map.orig, F)

csxp.map1.nofilter = get.allen.genetic.map(csxp.map.orig, csxp.geno.nofilter)

csxp.phys = get.phys.for.chromos(csxp.geno.nofilter, F, T, "-")

csxp.map = csxp.map1.nofilter
csxp.map$phys = csxp.phys$phys.pos
csxp.map$phys.bp = csxp.phys$phys.bp

csxp.map = order.by.physical.dist.within.bins(csxp.map, "chr", "cM", "phys")
csxp.map2 = split.df.into.list.of.dfs.by.column(csxp.map, "chr", F)
lapply(csxp.map2, check.if.inverted2, cm.colname = "cM", phys.pos.colname = "phys") #all chromosomes are in proper order compared to physical assembly

s(csxp.map, "rotation1scripts_v4/original_data/historical.recombination.analysis/csxp.map", "historical.recombination.mapping.R")

#### SERVER CODE ####
#prepare genotyping data for all varieties in cerealsDB
load("rotation1scripts_v4/original_data/historical.recombination.analysis/csxp.map")
geno = read.csv("rotation1scripts_v4/original_data/genotype.data/cerealsdb.35k.genotype.data.csv", stringsAsFactors = F, header = F)



geno2 = dcast(geno, V1 ~ V2, value.var = "V3")
geno3 = geno2[, c(1, match(csxp.map$marker, colnames(geno2)))]

geno4 = melt(geno3, "V1")


geno4$chr = unlist(lapply(csxp.map$chr, function(x) rep(x, 2119)))

colnames(geno4) = c("var", "marker", "geno", "chr")

geno5 = split(geno4, f = geno4$chr)

lapply(unique(geno4$chr), function(x){
    write.csv(geno5[[x]], paste0("rotation1scripts_v4/original_data/historical.recombination.analysis/", x, ".csv"), row.names = F)
})














