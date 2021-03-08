#install.packages("ASMap")
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/")
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


setwd("C:/Users/ac14037/project.phd.main/")
data.ian = read.csv("rotation1scripts_v4/original_data/ian.henderson.wei.f2/HEI10_PolyHighResolution_F2s.csv", header = T, stringsAsFactors = F)

g = quickEst(axiom.data, map.function = "kosambi")
g = mstmap(axiom.data, bychr = F, dist.fun = "kosambi", trace = T)


library(qtl)
library(ASMap)
axiom.data = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/original_data/ian.henderson.wei.f2/HEI10_PolyHighResolution_F2s.csv", genotypes = c("A", "H", "B"))
axiom.data = convert2bcsft(axiom.data, BC.gen=0, F.gen=2)
map.ian = mstmap(axiom.data, trace=T, dist.fun="kosambi", id="probeset_id")
map.ian2 = pull.map(map.ian, as.table=T)
    

mstmap.data.frame(data.ian, pop.type = "F2", dist.fun = "kosambi")

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
    

#make all four genetic maps
axiom.data = lapply(axiom.data, function(x){
    axiomdata=convert2bcsft(x,BC.gen=0,F.gen=2)
    axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1, anchor=T)
    q1_asmap_map=pull.map(axiomdata_map, as.table=T)
    r=rownames(q1_asmap_map)
    q1_asmap_map=add_column(q1_asmap_map, ske="")    
    q1_asmap_map$marker=r
    q1_asmap_map[,1]=as.character(q1_asmap_map[,1])
    
    for(i in 1:nrow(q1_asmap_map)){
        q1_asmap_map[i,1]=paste("LG",q1_asmap_map[i,1],sep="")    
    }
    
    colnames(q1_asmap_map)=c("lg","pos","ske","marker")
    
    return(q1_asmap_map)
})

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
