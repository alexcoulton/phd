#install.packages("ASMap")
setwd("C:/Users/ac14037/project.phd.main/rotation1scripts_v4/")
source("C:/Users/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")
library(ASMap)
library(qtl)
library(dplyr)
library(tibble)

list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)
path = function(quadrant, lg) paste("original_data/genotype.data/seperated.lgs.w.individuals.over.threshold.recombination.count.removed/q",
                                                                        quadrant, "lg", lg, ".csv", sep = "")

repten = function(x) rep(x, 10)
genotype.data.paths = Map(path, unlist(lapply(1:4, repten)), rep(list.of.full.lgs.to.test, 4))

#make rQTL cross object for all four quadrants
axiom.data = lapply(genotype.data.paths, function(x) read.cross(format = "csv", dir = "./", file = x, genotypes = c("A", "H", "B"), 
                                        estimate.map = F))


#     ____________________________________________________________________________
#     ESTIMATE MAPS (PRESERVE ORDER)                                                                                    ####

#estimate the maps
estimated.maps.rqtl = lapply(axiom.data, est.map, map.function = "kosambi") #want to compare est.map() to quickest()
estimated.maps.rqtl.quickest = lapply(axiom.data, quickEst, map.function = "kosambi")
estimated.maps = lapply(estimated.maps.rqtl.quickest, pull.map, as.table = T)

estimated.maps.adj.dist = lapply(estimated.maps, function(x){
    g = x$pos
    for(i in 2:nrow(x)){
        g[i] = x$pos[i] - x$pos[i-1]
    }
    return(g)
})

#grab mean and sd centimorgan distances for each quadrant 
grab.mean = function(lgs) mean(unlist(lapply(estimated.maps[lgs], function(x) x$pos[nrow(x)])))
grab.sd = function(lgs) sd(unlist(lapply(estimated.maps[lgs], function(x) x$pos[nrow(x)])))

lapply(list(1:10, 11:20, 21:30, 31:40), grab.mean)
lapply(list(1:10, 11:20, 21:30, 31:40), grab.sd)


#round cM to 2 dec. places
estimated.maps = lapply(estimated.maps, function(x){
    x$pos = as.numeric(x$pos)
    x$pos = round(x$pos, 2)
    return(x)
})

#do some formatting
estimated.maps = lapply(estimated.maps, function(x){
    x$ske = ""
    x$marker = row.names(x)
    x$chr = paste("LG", x$chr, sep = "")
    colnames(x) = c("lg", "pos", "ske", "marker")
    return(x)
})


for(i in seq_along(estimated.maps)){
    write.csv(estimated.maps[i], paste("processed_data/genetic_maps/q", i ,"_multipoint_ordercons_asmap_est.csv", sep = ""),
                        row.names = F)
}

# plotMap(estimated.maps.rqtl[[1]], show.marker.names = T, chr = 1)

#     ____________________________________________________________________________
#     PRODUCE MAPS (NEW ORDERING)                                                                                         ####

#make genetic map with all quadrants included    
axiom.data = read.cross(format = "csv", dir = "./", file = "original_data/genotype.data/manuallycuratedsnps_flipped.csv", genotypes = c("A", "H", "B"), 
                                                estimate.map = F)

axiomdata=convert2bcsft(axiom.data,BC.gen=0,F.gen=2)
axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1e-20, anchor=F)
q1_asmap_map=pull.map(axiomdata_map, as.table=T)
r=rownames(q1_asmap_map)
q1_asmap_map=add_column(q1_asmap_map, ske="")    
q1_asmap_map$marker=r
q1_asmap_map[,1]=as.character(q1_asmap_map[,1])

for(i in 1:nrow(q1_asmap_map)){
    q1_asmap_map[i,1]=paste("LG",q1_asmap_map[i,1],sep="")    
}

colnames(q1_asmap_map)=c("lg","pos","ske","marker")


remove.small.lgs = function(geneticmap, lgsize){
    linkage.groups = unique(geneticmap$chr)
    
    g = unlist(lapply(linkage.groups, function(x){
        g = geneticmap[which(geneticmap$chr == x), ]
        nrow(g)
        if(nrow(g) < lgsize){
            return(F)
        } else {
            return(T)
        }
    }))
    
    lg.to.keep = linkage.groups[which(g)]
    
    geneticmap[which(geneticmap$chr %in% lg.to.keep), ]
    
    
}

q1map.wo.small.lgs = remove.small.lgs(q1_asmap_map, 30)

lg.chromo.as = lapply(unique(q1map.wo.small.lgs$chr), function(x){
    g = q1map.wo.small.lgs[which(q1map.wo.small.lgs$chr == x), ]
    grab.consensus.chromosome(switch.affy.format(rownames(g)))
})

separate.chromo = function(map){
    lapply(unique(map$chr), function(x){
        g = map[which(map$chr == x), ]
        return(g)
    })
}

q = separate.chromo(q1map.wo.small.lgs)
q = q[-which(unlist(lapply(q, nrow)) < 50)]

chromo.candidates = lapply(q, function(x){
    g = switch.affy.format(rownames(x))
    grab.consensus.chromosome(g)
})

chromo.assignments = names(unlist(lapply(chromo.candidates, function(x){
    if(length(x) > 0){
        return(x[which(x == max(x))])
    } else {
        return("")
    }
})))

q2 = Map(function(x, y){
    x$chromo = y
    return(x)
}, q, as.list(chromo.assignments))

#find chromosome assignments for those linkage groups with no information from nullisomic consensus lines
phys.lists.9 = lapply(listofwheatchromosomes, function(x){
    genetictophysical(x, rownames(q2[[9]]), ".")
})

q2[[9]]$chromo = "5B"


phys.lists.10 = lapply(listofwheatchromosomes, function(x){
    genetictophysical(x, rownames(q2[[10]]), ".")
})

q2[[10]]$chromo = "1B"
q2[[11]]$chromo = "3B"


#continue with analysis


q3 = lapply(q2, function(x){
    phys = genetictophysical(x$chromo[[1]], row.names(x), markerformat = ".", threshold = 30)
    
    x$phys.dist = phys
    x$phys.dist.bp = attr(phys, "bp")
    return(x)
})

q4 = lapply(q3, function(x){
    x$cm.diff = c(0, diff(x$pos))
    return(x)
})

names(q4) = unlist(lapply(q4, function(x) x$chromo[[1]]))


q.comb = combine.list.of.data.frames(q4)


#     ____________________________________________________________________________
#     INTROGRESSION ANALYSIS                                                                                                    ####

introgression = read.table("original_data/sacha.introgression.data/Introgression_table_280617.tdt")
colnames(introgression) = c("elite", "chr", "bp.start", "bp.length", "snp.score", "wild.relative")

grab.only.large.intro = function(df){
    df[df$bp.length != 0, ]
}

introgression = grab.only.large.intro(introgression)

intro.cs = introgression[grep("ChineseSpring", introgression$elite), ] # grab chinese spring introgression data
intro.para = filter(introgression, elite == "Paragon")
intro.para2 = split.df.into.list.of.dfs.by.column(intro.para, "chr")
intro.para3 = lapply(intro.para2, function(x){
    arrange(x, bp.start)
})


intro.apo = introgression[grep("Apogee", introgression$elite), ] # grab apogee intro. data
split.and.arrange = function(intro.df1){
    intro.apo2 = split.df.into.list.of.dfs.by.column(intro.df1, "chr")
    intro.apo3 = lapply(intro.apo2, function(x){
        arrange(x, bp.start)
    })
    return(intro.apo3)
}
intro.apo3 = split.and.arrange(intro.apo)

gen.intr = function(list.of.df1, list.of.intro.df1, list.of.intro.df2, hotspot){
    #args:
    #hotspot: boolean, if T, grabs introgression associated with hotspots instead of coldspots
    if(missing(hotspot)) hotspot = F
    Map(function(gen.map, chromo.name){
        if(hotspot == T){
            phys = gen.map[which(gen.map$cm.diff != 0), ]$phys.dist.bp
        } else {
            phys = gen.map[which(gen.map$cm.diff == 0), ]$phys.dist.bp
        }
        check.elite.variety.for.matching.introgressions = function(list.of.chromo.dfs){
            if(hotspot == T){
                phys = gen.map[which(gen.map$cm.diff != 0), ]$phys.dist.bp
                phys.per = gen.map[which(gen.map$cm.diff != 0), ]$phys.dist
            } else {
                phys = gen.map[which(gen.map$cm.diff == 0), ]$phys.dist.bp
                phys.per = gen.map[which(gen.map$cm.diff == 0), ]$phys.dist
            }
            
            intro.candidates = lapply(phys, function(physical){
                which(list.of.chromo.dfs[[chromo.name]]$bp.start < (physical + 100000) & 
                                list.of.chromo.dfs[[chromo.name]]$bp.start > (physical - 100000))
            })
            
            intro.cand2 = unlist(intro.candidates)
            
            #grab introgressions that match up to non-recombinant markers     
            matching.introgressions = list.of.chromo.dfs[[chromo.name]][intro.cand2, ]    
            #remove non-recombinant markers that do not have any matching introgressions
            
            non.matching.markers = which(unlist(lapply(intro.candidates, length)) == 0)
            if(hotspot == T){
                non.skele.markers = gen.map[which(gen.map$cm.diff != 0), ][-non.matching.markers, ] #grab markers that do have introgressions
            } else {
                non.skele.markers = gen.map[which(gen.map$cm.diff == 0), ][-non.matching.markers, ] #grab markers that do have introgressions
            }
            phys.pos.of.non.matching.markers = phys.per[non.matching.markers]
            
            if(length(matching.introgressions) == 0 | length(non.skele.markers) == 0){
                return(list("no matching intro", "no.non.skele.markers", 
                                        p("Number of non-recombinant markers wo/ matching introgressions: ", 
                                            length(non.matching.markers)), non.matching.markers))
                
            } else if(nrow(matching.introgressions) == 0 | nrow(non.skele.markers) == 0) {
                
                return(list("no matching intro", "no.non.skele.markers", 
                                        p("Number of non-recombinant markers wo/ matching introgressions: ", 
                                            length(non.matching.markers)), non.matching.markers))
                
            } else {
                return(list(matching.introgressions, non.skele.markers, length(non.matching.markers), non.matching.markers))
            }
        }
        
        list(check.elite.variety.for.matching.introgressions(list.of.intro.df1),
                 check.elite.variety.for.matching.introgressions(list.of.intro.df2))
        
    }, list.of.df1, names(list.of.df1))
} 

intr = gen.intr(q8, intro.apo3, intro.para3)
# s(intr, "saved.objects/intr", "asmapv2.R")
# s(q10, "saved.objects/q10", "asmapv2.R")

#find regions with no introgressions in either paragon or apogee
no.intro.regions = lapply(intr, function(x){
    x[[1]][[4]][which(x[[1]][[4]] %in% x[[2]][[4]])]
})


to.rev = which(unlist(lapply(q4, function(x){
    check.if.inverted(x$phys.dist.bp)
})))

#reverse inverted chromosomes
q4.1 = lapply(q4, function(x){
    g = check.if.inverted(x$phys.dist.bp)
    if(g == T){
        return(reverse.genetic.map(x, "pos"))
    } else {
        return(x)
    }
})

#order the chromosomes by physical position within bins
q4.2 = lapply(q4.1, order.by.physical.dist.within.bins, cm.colname = "pos", phys.pos.colname = "phys.dist.bp")

#add information from sacha data 
q5 = Map(function(sacha, mine){
    coords = match(rownames(sacha), rownames(mine))
    mine$sacha.pos = ""
    mine$sacha.cm.diff = ""
    mine$sacha.pos[coords] = sacha$pos
    mine$sacha.cm.diff[coords] = sacha$cm.diff
    return(mine)
    
}, sacha.map.split3, q4.2)

q5.1 = lapply(q5, function(x){
    x$cm.diff = calculate.cm.diff(as.numeric(x$pos))
    x$sacha.cm.diff = calculate.cm.diff(as.numeric(x$sacha.pos))
    return(x)
})

intro.positions = lapply(intr, function(x){
    if(x[[1]][2] == "no.non.skele.markers" & x[[2]][2] == "no.non.skele.markers"){
        return()
    } else if(x[[1]][2] == "no.non.skele.markers") {
        return(x[[2]][2])
    } else if(x[[2]][2] == "no.non.skele.markers"){
        return(x[[1]][2])
    } else {
        return(c(x[[1]][2], x[[2]][2]))
    }
})

intro.positions2 = intro.positions[-which(unlist(lapply(intro.positions, is.null)))]
intro.positions3 = lapply(intro.positions2, function(x) x[[1]])

q7 = q6

#mark markers that do have introgressions
q7[match(names(intro.positions3), names(q7))] = Map(function(main.maps, intro.posi){
    main.maps$contains.intro = ""
    main.maps$contains.intro[match(intro.posi$phys.dist.bp, main.maps$phys.dist.bp)] = T
    return(main.maps)
}, q7[match(names(intro.positions3), names(q7))], intro.positions3)

q8 = lapply(q7, function(x){
    x$phys.dist.diff.kb = c(0, diff(x$phys.dist.bp)/1000)
    return(x)
})

#mark regions with no introgressions on the genetic map dataframe
q9 = Map(function(my.map, no.intro){
    my.map$no.introg = ""
    my.map[my.map$cm.diff == 0, ]$no.introg[no.intro] = "N"
    return(my.map)
}, q8, no.intro.regions)

#mark the markers that have no recombination events in my map but do have recombination events in sacha's map
q10 = lapply(q9, function(x){
    x$rec.w.h.freq = ""
    coords = which(x$cm.diff == 0)
    inspect.coords = which(x$cm.diff[coords] != x$sacha.cm.diff[coords])
    x$rec.w.h.freq[coords][inspect.coords] = "T"
    return(x)
})





#some separate analysis starts here... 

g = q4[[1]][which(q4[[1]]$cm.diff == 0), ]$phys.dist.bp
intro.candidates = lapply(g, function(x){
    which(intro.para3[["2A"]]$bp.start < (x + 100000) & intro.para3[["2A"]]$bp.start > (x - 100000))
})

intro.cand2 = unlist(intro.candidates)

intro.para3[["2A"]][intro.cand2, ]
q4[[1]][which(q4[[1]]$cm.diff == 0), ][-which(unlist(lapply(intro.candidates, length)) == 0), ]


plot.intro.closure = function(first.var, second.var, first.var.name, second.var.name){
    intro.cs = first.var
    intro.para = second.var
    plot.introgression = function(chromo){
        intro.cs.chromo = filter(intro.cs, chr == chromo)
        intro.cs.chromo = intro.cs.chromo[sort(intro.cs.chromo$bp.start, index.return = T)$ix, ]
        
        intro.para.chromo = filter(intro.para, chr == chromo)
        intro.para.chromo = intro.para.chromo[sort(intro.para.chromo$bp.start, index.return = T)$ix, ]
        
        intro.comb = rbind(intro.para.chromo, intro.cs.chromo)
        intro.comb = arrange(intro.comb, bp.start)
        
        intro.comb = intro.comb[as.numeric(rownames(unique(intro.comb[, c(1, 3, 6)]))), ]
        
        ggplot(intro.comb, aes(x = bp.start, y = bp.length, color = elite, shape = elite, group = elite)) + 
            geom_point() +
            scale_color_discrete(name = "", breaks = c(first.var.name, second.var.name), labels = c(substr(first.var.name, 1, 2), substr(second.var.name, 1, 2))) +
            scale_shape_discrete(name = "", breaks = c(first.var.name, second.var.name), labels = c(substr(first.var.name, 1, 2), substr(second.var.name, 1, 2)))
    }
}

plot.apo.para = plot.intro.closure(intro.apo, intro.para, "Apogee", "Paragon")
plot.cs.para = plot.intro.closure(intro.cs, intro.para, "Chinese Spring", "Paragon")



#organise genotype data and split into quadrants (this isn't actually needed for this analysis -- only realised this after doing)

man.geno.data = read.csv("original_data/genotype.data/manuallycuratedsnps_flipped.csv", header = T, stringsAsFactors = F)

man.geno.data.mst.map = man.geno.data[, c(1, 2, match(row.names(q.comb), colnames(man.geno.data)))]

man.geno.data.mst.map[1, 3:ncol(man.geno.data.mst.map)] = q.comb$chromo
man.geno.data.mst.map[1, 1] = ""

quads = c("Q1", "Q2", "Q3", "Q4")
quad.geno.data = lapply(quads, function(x){
    man.geno.data.mst.map[c(1, grep(x, man.geno.data.mst.map$probeset_id)), ]
})

Map(function(x, y){
    write.csv(x, p("original_data/genotype.data/rotation1.quads.mstmap.w.ske.incl/", y, ".csv"), row.names = F)
}, quad.geno.data, quads)





# genetictophysical("2D", row.names(q1map.wo.small.lgs[which(q1map.wo.small.lgs$chr == 1.1), ]), ".")


#     ____________________________________________________________________________
#     CHECK SACHA APOGEE X PARAGON COLD SPOTS                                                                 ####

sacha.a.x.p = read.csv("original_data/genotype.data/apogee.x.paragon.sacha/a.x.p.map.correct.format.csv")

v(q.comb)

aas = read_delim("original_data/genotype.data/apogee.x.paragon.sacha/rec.snps.txt", "\t")
aas2 = convert.aas.to.rqtl("original_data/genotype.data/apogee.x.paragon.sacha/rec.snps.txt", c("Apogee", "Paragon"))
#following functions from a.x.c.initial
aas3 = cleanup(aas2, c(2, 5))
aas3 = convert.to.character.data.frame(aas3)
aas4 = assign.parental.genotypes(aas3, c(2, 5))
aas5 = remove.excess.parents(aas4, c("Paragon"))
aas5 = remove.excess.parents(aas5, c("Apogee"))

aas6 = aas5[, c(1, 2, na.omit((match(rownames(q.comb), switch.affy.format(colnames(aas5))))))]
aas6 = convert.to.character.data.frame(aas6)
aas6[1, 3:ncol(aas6)] = q.comb$chromo[na.omit(match(switch.affy.format(colnames(aas6)), rownames(q.comb)))]

write.csv(aas6, "processed_data/genotypes/apogee.x.paragon.sacha/genotypes.matching.alex.f2.map.csv", row.names = F)

sacha.qtl.class = read.cross(format = "csv", dir = "./", file = "processed_data/genotypes/apogee.x.paragon.sacha/genotypes.matching.alex.f2.map.csv", genotypes = c("A", "H", "B"), 
                     estimate.map = F)

sacha.map.est = quickEst(sacha.qtl.class, map.function = "kosambi")
sacha.map = pull.map(sacha.map.est, as.table = T)
rownames(sacha.map) = switch.affy.format(rownames(sacha.map))

sacha.map.split = split.df.into.list.of.dfs.by.column(sacha.map, "chr")

#check which of sacha's chromosomes are missing markers from my map
Map(function(sacha, mine){
    length(which(!rownames(mine) %in% rownames(sacha)))
}, sacha.map.split, q4)

calculate.cm.diff = function(list.of.centimorgan.dist) c(0, diff(list.of.centimorgan.dist))
sacha.map.split2 = lapply(sacha.map.split, function(x){
    x$cm.diff = calculate.cm.diff(x$pos)
    return(x)
})

sacha.map.split3 = sacha.map.split2

sacha.map.split3[to.rev] = lapply(sacha.map.split3[to.rev], function(x){
    reverse.genetic.map(x, "pos")
})


#     ____________________________________________________________________________
#     OTHER PROCESSING                                                                                                                ####


for(i in seq_along(axiom.data)){
    write.csv(axiom.data[[i]], paste("processed_data/genetic_maps/q", i, "_asmap_map_nonmonotonic_removed.csv", sep = ""), row.names=F)
}


summary(axiomdata)




#     ____________________________________________________________________________
#     REDUCED DATASET ANALYSIS                                                                                                ####


q1_asmap_map=read.csv("reduced_datasets/q1_asmap_map.csv", header=T, stringsAsFactors = F)
q2_asmap_map=read.csv("reduced_datasets/q2_asmap_map.csv", header=T, stringsAsFactors = F)
q3_asmap_map=read.csv("reduced_datasets/q3_asmap_map.csv", header=T, stringsAsFactors = F)
q4_asmap_map=read.csv("reduced_datasets/q4_asmap_map.csv", header=T, stringsAsFactors = F)

listoflgs=unique(q1_asmap_map$lg)
listoflgs

#reverse a particular linkage group --- see lab notebook page 16/12/2017
#reverse all LGs in a map
q2=q4_asmap_map
for(i in listoflgs){
    lgtoreverse=i
    q2[q2$lg==lgtoreverse,]=q2[q2$lg==lgtoreverse,][sort(which(q2[q2$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
    q2[q2$lg==lgtoreverse,2]=-(q2[q2$lg==lgtoreverse,2]-max(q2[q2$lg==lgtoreverse,2]))
    q2[q2$lg==lgtoreverse,]
}
q4_asmap_map=q2

#reverse a single LG
q2=q4_asmap_map
lgtoreverse="LG3"
q2[q2$lg==lgtoreverse,]=q2[q2$lg==lgtoreverse,][sort(which(q2[q2$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
q2[q2$lg==lgtoreverse,2]=-(q2[q2$lg==lgtoreverse,2]-max(q2[q2$lg==lgtoreverse,2]))
q2[q2$lg==lgtoreverse,]
q4_asmap_map=q2



for(i in unique(q1_asmap_map$lg)){
    filmap=filter(q1_asmap_map, lg==i)
    filmap2=filter(q4_asmap_map, lg==i)
    print(i)
    print(filmap$marker==filmap2$marker)
}

lgtotest="LG3"
filmap=filter(q1_asmap_map, lg==lgtotest)
filmap2=filter(q2_asmap_map, lg==lgtotest)
filmap3=filter(q3_asmap_map, lg==lgtotest)
filmap4=filter(q4_asmap_map, lg==lgtotest)
filmap$marker==filmap2$marker
filmap$marker==filmap3$marker
filmap$marker==filmap4$marker

?barplot
windows()
par(mfrow=c(2,2))
barplot(filmap$pos, ylim=c(0,160))
barplot(filmap2$pos, ylim=c(0,160))
barplot(filmap3$pos, ylim=c(0,160))
barplot(filmap4$pos, ylim=c(0,160))

list(filmap$pos)
list(filmap2$pos)
windows()
plotMap(list(filmap$pos))
?plotMap
