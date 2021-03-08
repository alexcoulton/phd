#generic ASMap genetic map estimation

library(ASMap)
library(qtl)
library(dplyr)
library(tibble)
library(parallel)

# setwd("E:/phd.project.main/")

#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####

remove.small.lgs = function(geneticmap, lgsize){
    #Removes linkage groups that have less than lgsize markers
    #args:
    # geneticmap - Dataframe of the genetic map, with columns chr, pos, ske and marker
    #lgsize - Integer, remove linkage groups that are less than this
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

separate.chromo = function(map){
    #Seperates dataframe into list of dataframes by chromosome
    #args:
    # map - Dataframe containing genetic map
    lapply(unique(map$chr), function(x){
        g = map[which(map$chr == x), ]
        return(g)
    })
}

grab.consensus.chromosome = function(list.of.markers){
    #args: list.of.markers - character vector of marker names in "." format 
    #examines cerealsdb nullisomic line consensus data and returns a table of the most common 
    #chromosome assignment values for this group of markers 
    if(("rs2" %in% ls(.GlobalEnv)) == F){
        library(RMySQL)    
        mydb = setup.mysql.connection()
        
        rs = dbSendQuery(mydb, "SELECT affycode, consensus FROM axiom820") #grab Avalon genotyping data for 35k array
        rs2 <<- fetch(rs, n=-1)
    }
    
    g = rs2$consensus[match(list.of.markers, rs2$affycode)]
    g2 = na.omit(g[-c(which(g == ""), which(g == "none"))])
    g3 = gsub("L", "", g2)
    g4 = gsub("S", "", g3)
    
    g4 = g4[-grep("^\\d$", g4)]
    
    table(g4)
}

#     ____________________________________________________________________________
#     START ANALYSIS                                                                                                                    ####

#first read in data using rQTL function
axiom.data = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/all.replicates.flipped.geno.data.newprobeset.csv", genotypes = c("A", "H", "B"), 
                                                estimate.map = F)

#Specify the Filial generation of the cross with the F.gen parameter
axiomdata = convert2bcsft(axiom.data, BC.gen = 0, F.gen = 2, estimate.map = F)

#calculate genetic map
# axiomdata_map = mstmap(axiomdata, trace = T, dist.fun = "kosambi", id = "probeset_id", p.value = 1e-50, anchor = F)
higher.clust = c(1e-50, 1e-51, 1e-52, 1e-53, 1e-54, 1e-55, 1e-56, 1e-57, 1e-58, 1e-59, 1e-60)

high.clust = c(1e-40, 1e-41, 1e-42, 1e-43, 1e-44, 1e-45, 1e-46, 1e-47, 1e-48, 1e-49, 1e-50)
med.clust = c(1e-30, 1e-31, 1e-32, 1e-33, 1e-33, 1e-35, 1e-36, 1e-37, 1e-38, 1e-39, 1e-40)
low.clust = c(1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15)
#Run mstmap with multiple clustering parameter settings in parallel. (must be done on Wilkins)
maps.parallel_high = mclapply(high.clust, function(x){
    axiomdata_map = mstmap(axiomdata, trace = T, dist.fun = "kosambi", id = "probeset_id", p.value = x, anchor = F)
}, mc.cores = 20)

#if only a clustering parameter has been used to generate the map:
# maps.parallel = list(axiomdata_map)


maps.parallel2 = lapply(maps.parallel_high, function(x){
    q1_asmap_map = pull.map(x, as.table=T)
    
    #do some formatting
    r = rownames(q1_asmap_map)
    q1_asmap_map = add_column(q1_asmap_map, ske="")    
    q1_asmap_map$marker = r
    q1_asmap_map[,1] = as.character(q1_asmap_map[,1])
    
    for(i in 1:nrow(q1_asmap_map)){
        q1_asmap_map[i,1] = paste("LG", q1_asmap_map[i,1], sep = "")    
    }
    
    colnames(q1_asmap_map) = c("chr", "pos", "ske", "marker")
    q1_asmap_map
})


grab.lg.max.cm = function(genetic.map.df){
    unlist(lapply(unique(genetic.map.df$chr), function(x){ 
        max(filter(genetic.map.df, chr == x)$pos    )
    }))    
}





#what is the highest cM value for each of the clustering parameter settings
assess.maps = lapply(maps.parallel2, function(x) table(cut(grab.lg.max.cm(x), c(0, 20, 80, 100, 200, 400, 600, 1000, 10000))))

#how many of the linkage groups are over 50 cM
assess.maps2 = lapply(maps.parallel2, function(x) length(which(grab.lg.max.cm(x) > 50)))

#I will pick map 11 from the various clustering parameter settings
q1_asmap_map = maps.parallel2[[7]]


#remove small linkage groups from genetic map
q1map.wo.small.lgs = remove.small.lgs(q1_asmap_map, 30)

q = separate.chromo(q1map.wo.small.lgs)
q = q[-which(unlist(lapply(q, nrow)) < 50)]

#query cerealsDB mySQL database for nullisomic line information for chromosome identification
chromo.candidates = lapply(q, function(x){
    # g = switch.affy.format(rownames(x))
    grab.consensus.chromosome(rownames(x))
})

#analyse chromo.candidates and pick the most common chromosome 
chromo.assignments = names(unlist(lapply(chromo.candidates, function(x){
    if(length(x) > 0){
        return(x[which(x == max(x))][1]) #if two equally selected chromosomes, pick the first
    } else {
        return("")
    }
})))

names(chromo.assignments) = 1:length(chromo.assignments)

#add chromosome assignments to the genetic map
q2 = Map(function(x, y){
    x$chromo = y
    return(x)
}, q, as.list(chromo.assignments))

#find chromosome assignments for those linkage groups with no information from nullisomic consensus lines
# s(q2, "rotation1scripts_v4/saved.objects/q2", "generic.asmap.R")

#the following code should be executed on Wilkins

# load("~/project.phd.main/rotation1scripts_v4/saved.objects/q2")
# source("~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")
g.phys = mclapply(q2, function(q2.chromo){

    phys.lists.9 = lapply(listofwheatchromosomes, function(x){
        genetictophysical(x, rownames(q2.chromo), "-")
    })

    #the output of genetictophysical() that has the least amount of NAs is probably the
    #correct chromosome assignment for this particular linkage group
    process.g2p.output = unlist(lapply(phys.lists.9, function(x) length(which(is.na(x)))))
    listofwheatchromosomes[which(process.g2p.output == min(process.g2p.output))]

}, mc.cores = 50)

physical.chromosome.assignments = unlist(g.phys)
names(physical.chromosome.assignments) = 1:length(physical.chromosome.assignments)

#check mismatches between chromosome assignments from nullisomics and physical BLAST search
chromo.assignments[which(chromo.assignments != physical.chromosome.assignments)]
physical.chromosome.assignments[which(chromo.assignments != physical.chromosome.assignments)]

for(i in which(chromo.assignments == "")){
    q2[[i]]$chromo = physical.chromosome.assignments[[i]]
}


#some manual chromosome assignments if nullisomic cerealsdb data doesn't work
q2[[5]]$chromo = "3B"
q2[[21]]$chromo = "5B"
# q2[[16]]$chromo = "3B"
# q2[[20]]$chromo = "1A"
# q2[[25]]$chromo = "3D"


#Evaluate which of the chromosomes a particular linkage group has the most BLAST hits to
# phys.lists.10 = lapply(listofwheatchromosomes, function(x){
#     genetictophysical(x, rownames(q2[[10]]), "-")
# })

#Add physical locations of markers to genetic map
q3 = mclapply(q2, function(x){
    phys = genetictophysical(x$chromo[[1]], row.names(x), markerformat = "-", threshold = 30)
    
    x$phys.dist = phys
    x$phys.dist.bp = attr(phys, "bp")
    return(x)
}, mc.cores = 27)

inverted.lgs = unlist(lapply(q3, function(x){
    # print(na.omit(x$phys.dist.bp))
    check.if.inverted(na.omit(x$phys.dist.bp))
    
}))

for(i in which(inverted.lgs)){
    q3[[i]] = reverse.genetic.map(q3[[i]], "pos")
}

#Add difference in centiMorgans between adjacent markers column (cm.diff) to genetic map
q4 = lapply(q3, function(x){
    x$cm.diff = c(0, diff(x$pos))
    return(x)
})

names(q4) = unlist(lapply(q4, function(x) x$chromo[[1]]))

#remove excessively large linakge groups
# large.lgs = which(unlist(lapply(q4, function(x){
#     max(x$pos)
# })) > 250)

# if(length(large.lgs) > 0){
#     q4 = q4[-large.lgs]    
# }


min.max.phys = lapply(q4, function(x){
    min.bp = min(na.omit(x$phys.dist))
    max.bp = max(na.omit(x$phys.dist))
    list.bp = c(min.bp, max.bp)
    names(list.bp) = c("min.bp", "max.bp")
    list.bp
})

#for which linkage groups do we have markers that span more than 80% of the physical distance of the chromosome?
good.lg.coords = which(unlist(lapply(min.max.phys, function(x){
    abs(x[1] - x[2]) > 30
})))

#sort by physical distance within bins (for more markers after LISS)
q5 = lapply(q4, function(x){
    x2 = split(x, x$pos)
    bind_rows(lapply(x2, function(y){
        y[sort(y$phys.dist, index.return = T, na.last = T)$ix, ]
    }))
})

q.comb = combine.list.of.data.frames(q5[good.lg.coords])

write.csv(q.comb, "rotation1scripts_v4/processed_data/genetic_maps/cxaf2/cxa.axc.all.replicates.mapv3.newprobeset.csv")



