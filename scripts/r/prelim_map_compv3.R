#setwd("~/Google Drive/R/scripts/rotation1scripts_v2/") #macbook settings
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/")
source("scripts/r/functions.R")
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)

#load up each quadrant and assign bins (sequentially)
list.of.quadrant.names=c("q1","q2","q3","q4")
map.path="processed_data/genetic_maps/"

# ---------------------------------------------#
#DEFINE FUNCTIONS---------------------------####
read.quadrant.genetic.map=function(quadrant){
    #args:
    #quadrant --- a string, the name of the quadrant to read data for (q1, q2 etc.)
        read.csv(paste(map.path, quadrant, "_multipoint_ordercons_asmap_est_w_chr.csv", sep=""), header=T, stringsAsFactors = F)
}
    
#assign bins to each LG
assignbins=function(quadrant){
    #args:
    #quadrant --- a genetic map dataframe
    for(g in unique(quadrant$lg)){
        counter=1
        for(i in unique(quadrant[quadrant$lg==g,2])){
            quadrant[quadrant$pos==i,3]=counter
            counter=counter+1
        }
    }
    return(quadrant)
}

reverse.lg=function(lgtoreverse, quadrant){
    #args:
    #lgtoreverse --- string, the name of the LG (as specified in genetic map) to reverse
    #quadrant --- genetic map dataframe
    quadrant.fil=filter(quadrant, lg==lgtoreverse)
    quadrant.fil=quadrant.fil[sort(which(quadrant.fil[quadrant.fil$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
    quadrant.fil[quadrant.fil$lg==lgtoreverse,2]=-(quadrant.fil[quadrant.fil$lg==lgtoreverse,2]-max(quadrant.fil[quadrant.fil$lg==lgtoreverse,2]))
    quadrant.fil[quadrant.fil$lg==lgtoreverse,]
    quadrant[quadrant$lg==lgtoreverse,]=quadrant.fil
    return(quadrant)
}

#---------------#
#LOAD DATA---------------------####
list.of.quadrants=lapply(list.of.quadrant.names, read.quadrant.genetic.map)

#assign markers
list.of.quadrants.w.bins=lapply(list.of.quadrants, assignbins)

#NEED TO PERFORM CHROMOSOME ASSIGNMENT WITH cerealscomp_multipointv2.R BEFORE USING THIS SCRIPT

#     ____________________________________________________________________________
#     BEGINING OF PAIRWISE LG COMPARISON FOR Q3 AND Q4                                                ####

listselect1=1
listselect2=2
nonmonotonic.markers=newdf(colnames = c("thesemarkers","thesebins","thosebins","thosemarkers","monotony.indicator"))
alldatabins=newdf(colnames = c("thesemarkers","thesebins","thosebins","thosemarkers","monotony.indicator"))
q1=list.of.quadrants.w.bins[[listselect1]]
q2=list.of.quadrants.w.bins[[listselect2]]
#begin main comparisons of bin and marker order
listof.lg=unique(q1$lg)
#listof.lg=listof.lg[-5]
first.quad.to.compare=list.of.quadrants.w.bins[[listselect1]]; first.quad.to.compare.name=paste("q",listselect1,sep="")
second.quad.to.compare=list.of.quadrants.w.bins[[listselect2]]; second.quad.to.compare.name=paste("q",listselect2,sep="")
flipped.lgs=list()
for(b in listof.lg){
#DEBUGGING: GOING THROUGH CODE WITHIN LOOP MANUALLY
onebq1=first.quad.to.compare[first.quad.to.compare$lg==b,] #loop through list
onebq2=second.quad.to.compare[second.quad.to.compare$lg==b,]



# onebq1=q1[q1$lg=="LG1",] #choose chromosome manually
# onebq2=q2[q2$lg=="LG1",]

nrow(onebq1)
nrow(onebq2)

#check common markers between LGs
length(which(onebq1[,4] %in% onebq2[,4]))
length(which(onebq2[,4] %in% onebq1[,4]))

onebq1_comb=cbind(onebq1,onebq2)
onebq1_comb_chr_trunc=onebq1_comb[,c(6,12)]
onebq1_comb=onebq1_comb[,-c(6,12)]

onebq1_comb=add_column(onebq1_comb, .after=5, q1markerposinq2="" )
onebq1_comb=add_column(onebq1_comb, .after=6, q2markerposinq1="" )

for(i in 1:nrow(onebq1_comb)){
    onebq1_comb[i,7]=match(onebq1_comb[i,11], onebq1_comb[,4])
}
for(i in 1:nrow(onebq1_comb)){
    onebq1_comb[i,6]=match(onebq1_comb[i,4], onebq1_comb[,11])
}



onebq1_comb=onebq1_comb[,c(1,2,5,4,3,6,7,10,11,12,9,8)] #reorganize columns for ease of comparison
onebq1_comb=add_column(onebq1_comb, .after=6, bins="") 
onebq1_comb=add_column(onebq1_comb, .after=7, bins2="")
colnames(onebq1_comb)=c("lg1","cM1","chr1","marker1","skeleton1","q1markerposinq2","bins","bins2","q2markerposinq1","skeleton2","marker2","chr2","cM2","lg2")
onebq1_comb_backup=onebq1_comb
onebq1_comb=onebq1_comb_backup

onebq1_comb$bins=onebq1_comb$skeleton1
onebq1_comb$bins2=onebq1_comb$skeleton2

#which bins in the other quadrant contain markers that are in this bin?

thesebins=list()
thosebins=list()
thesemarkers=list()
thosemarkers=list()
thismarker_colnum=4
thatmarker_colnum=11
thisbin_colnum=7
thatbin_colnum=8
list_bins2=unique(onebq1_comb[,thisbin_colnum])
for(i in list_bins2){
    for(g in 1:length(onebq1_comb[onebq1_comb[,thisbin_colnum]==i, thismarker_colnum])){
        rowpos=match(onebq1_comb[onebq1_comb[,thisbin_colnum]==i,thismarker_colnum][g], onebq1_comb[,thatmarker_colnum])
        thesemarkers=c(thesemarkers, onebq1_comb[onebq1_comb[,thisbin_colnum]==i,thismarker_colnum][g])
        thosemarkers=c(thosemarkers, onebq1_comb[rowpos, thatmarker_colnum])
        thesebins=c(thesebins, i)
        thosebins=c(thosebins, onebq1_comb[rowpos,thatbin_colnum])
    }
}

#store this information in dataframe databins
databins=do.call(rbind, Map(data.frame, thesemarkers=thesemarkers, thesebins=thesebins, thosebins=thosebins, thosemarkers=thosemarkers)) #convert two lists into a dataframe https://stackoverflow.com/questions/28630024/combine-two-lists-in-a-dataframe-in-r
databins$thesebins=as.character(databins$thesebins)
databins$thosebins=as.character(databins$thosebins)

databins_sorted=databins[1:2,]
databins$thesebins=as.numeric(as.character(databins$thesebins))
databins$thosebins=as.numeric(as.character(databins$thosebins))

databins_sorted$thosebins


#check if "thosebins" is monotonic after sorting within each unique element of thesebins
checkformonotony=list()
for(i in unique(databins$thesebins)){
    sort_index=unname(unlist(sort(databins[databins$thesebins==i, 3], index.return=T)[2]))
    databins_sorted=rbind(databins_sorted, databins[databins$thesebins==i,][sort_index,])
    checkformonotony=c(checkformonotony, sort(databins[databins$thesebins==i, 3])) #sort "thosebins" from databins (monotony doesn't have to be retained within thosebins, it only needs to be retained for each block of thesebins)
}
databins_sorted=databins_sorted[-c(1,2),]
checkformonotony

#label non-monotonic bins in databins_sorted df
databins_sorted$monotony.indicator=""
databins_sorted$thosebins=as.numeric(databins_sorted$thosebins)
for(i in 1:(nrow(databins_sorted)-1)){
    if(databins_sorted$thosebins[i] <= databins_sorted$thosebins[i+1]) {
        databins_sorted$monotony.indicator[i+1] = "S"
    } else {
        databins_sorted$monotony.indicator[i+1] = "N"
    }
}
nonmonotonic.markers=rbind(nonmonotonic.markers, filter(databins_sorted, monotony.indicator=="N"))
alldatabins=rbind(alldatabins, databins_sorted)

onebq1_comb$dist_adjacent1=""
for(i in 1:(nrow(onebq1_comb)-1)){
    onebq1_comb[i,15]=(onebq1_comb[(i+1),2]-onebq1_comb[i,2])
}

onebq1_comb$dist_adjacent2=""
for(i in 1:(nrow(onebq1_comb)-1)){
    onebq1_comb[i,16]=(onebq1_comb[(i+1),13]-onebq1_comb[i,13])
}

onebq1_comb=cbind(onebq1_comb, onebq1_comb_chr_trunc)

#process chromosome vs genome blasts


listof.blastfiles=list.files("original_data/blast-files-LG-vs-genome/")





auto.chromosome.assignment=q1$chr[which(q1$lg==b)[1]]
chromosome_truncated=q1$chr_truncated[which(q1$lg==b)[1]]
#print(c("chromosome assignment", auto.chromosome.assignment))

if(file.exists(paste("original_data/blast-files-LG-vs-genome/listofprobesandsequences",auto.chromosome.assignment,".fa.blast",sep=""))==T){
    #print(c(auto.chromosome.assignment, "exists"))
    #chromo.blast.results=read.table(paste("original_data/blast-files-LG-vs-genome/listofprobesandsequences",b,".fa.blast",sep=""))
    chromo.blast.results=read.table(paste("original_data/blast-files-LG-vs-genome/listofprobesandsequences",auto.chromosome.assignment,".fa.blast",sep=""))
    nrow(chromo.blast.results[chromo.blast.results[,2]==paste("chr",chromosome_truncated,sep=""),])
    chromo.blast.results.filtered=chromo.blast.results[chromo.blast.results[,2]==paste("chr",chromosome_truncated,sep=""),]
    
    chromo.blast.results.filtered=filter(chromo.blast.results.filtered, V11<1e-24) #remove hits with low e-values
    chromo.blast.results.filtered=distinct(chromo.blast.results.filtered, V1, .keep_all=T)    
    
    #grab and process individual chromosome character lengths
    chromo.character.counts=read.table("original_data/fullchromosomesequences/chromosome_character_counts.txt",sep="\n", stringsAsFactors = F)
    chromo.character.counts=add_column(chromo.character.counts, one="", two="")
    for(i in 1:nrow(chromo.character.counts)){
        splitted=strsplit(chromo.character.counts[i,1], split=" ")
        chromo.character.counts[i,2]=splitted[[1]][1]
        chromo.character.counts[i,3]=splitted[[1]][3]
    }
    
    chromo.character.counts[,3]=substr(chromo.character.counts[,3],4,5)
    
    
    
    #add blast information to onebq1_comb dataframe
    matches=match(chromo.blast.results.filtered[,1], onebq1_comb$marker1)
    count=1
    onebq1_comb$blast_start=""
    for(i in matches){
        onebq1_comb$blast_start[i]=chromo.blast.results.filtered[count,9]
        count=count+1
    }
    
    #convert column to numeric and remove NA values for subsequent section to work properly
    onebq1_comb$blast_start=as.numeric(onebq1_comb$blast_start)
    onebq1_comb$blast_start[which(is.na(onebq1_comb$blast_start))]=""
    
    #add the percentage of the physical sequence at which the marker blast hit occurs
    onebq1_comb$blast_percentage=""
    for(i in 1:nrow(onebq1_comb)){
        if(onebq1_comb$blast_start[i] != ""){
            chromo.count.index=match(chromosome_truncated, chromo.character.counts$two)
            onebq1_comb$blast_percentage[i]=(as.numeric(onebq1_comb$blast_start[i])/as.numeric(chromo.character.counts[chromo.count.index,2]))*100
        }
    }
    print(as.numeric(onebq1_comb$q1markerposinq2[1:10]))
    if(!(1 %in% as.numeric(onebq1_comb$q1markerposinq2[1:10]))){
        print("The following LG appears to be flipped")
        print(onebq1_comb$lg1[1])
        print(as.numeric(onebq1_comb$q1markerposinq2[1:10]))
        flipped.lgs=c(flipped.lgs, onebq1_comb$lg1[1])
    }
    write.csv(onebq1_comb, paste("processed_data/pairwiseLGcomparisons/", first.quad.to.compare.name, "_", second.quad.to.compare.name, "/onebq1_comb", first.quad.to.compare.name, "_", second.quad.to.compare.name, "_",b,"_",auto.chromosome.assignment,".csv",sep=""), row.names = F)
    
}

onebq1_comb_q1_q2=onebq1_comb



} #temporary end of loop

#setup multiple line plot comparison
q1list=list.files("processed_data/pairwiseLGcomparisons/q1_q2/")
q3list=list.files("processed_data/pairwiseLGcomparisons/q3_q4/")
for(i in 1:length(q1list)){
    print(q1list[i])
    print(q3list[i])
    one=read.csv(paste("processed_data/pairwiseLGcomparisons/q1_q2/",q1list[i],sep=""), header=T)
    two=read.csv(paste("processed_data/pairwiseLGcomparisons/q3_q4/", q3list[i],sep=""), header=T)
    multi.line.df=one[,c(2,13)]
    multi.line.df=cbind(multi.line.df, two[,c(2,13)])
    colnames(multi.line.df)=c("10","26","14","28")
    multi.line.df$id=seq(1,nrow(multi.line.df),1)
    multi.line.df.melted=melt(multi.line.df, id.vars="id")
    
    graph = ggplot(data=multi.line.df.melted, aes(x=id, y=value)) + geom_line(aes(colour=variable, linetype=variable), size=1.2)
    png(filename=paste("multilineplots/",q1list[i],".png",sep=""))
    print(graph + theme_classic() + labs(x="Marker ID", y="Distance along LG (cM)")) #need to use print() function to output ggplot graphs to png
    dev.off()
}

?aes
?geom_line
?ggplot
??melt

#     ____________________________________________________________________________
#     END OF LOOP                                                                                                                         ####




#some graphs

#attempting to plot only unique values... not working... need to investigate further
# png(filename=paste("scatterplotsq1vsq2_adja/q1vsq1-",b,".png",sep=""))
# par(mfrow=c(1,2))
# plot(unique(onebq1_comb$cM1), unique(onebq1_comb$dist_adjacent1), xlim=c(0,160), ylim=c(0,12))
# plot(unique(onebq1_comb$cM2), unique(onebq1_comb$dist_adjacent2), xlim=c(0,160), ylim=c(0,12))
# dev.off()


#now do the whole process again but in reverse (comparing thosebins to thesebins)
#if checkformonotony is monotonic for both, then bin order is conserved.


# par(mfrow=c(1,1))
# png(filename=paste("scatterplotsq1vsq2/q1vsq1-",b,".png",sep=""))
# plot(onebq1_comb$cM1, onebq1_comb$cM2, xlim=c(0,120), ylim=c(0,120))
# abline(0,1) #add y=x for illustration
# dev.off()




#     ____________________________________________________________________________
#     probe sequence extraction                                                                                             ####

list.of.probes.in.chromosome=onebq1_comb$marker1

all.axiom.probes=read.table("35kprobes.fa",sep="\n", stringsAsFactors = F)

listofprobesandsequences=list()
for(i in list.of.probes.in.chromosome){
    indexofprobe=match(paste(">",i,sep=""), all.axiom.probes[,1])
    #print(indexofprobe)
    seqs=all.axiom.probes[(indexofprobe+1),1]
    #print(seqs)
    listofprobesandsequences=c(listofprobesandsequences, paste(">",i,sep=""), seqs)
}

listofprobesandsequences

d<-lapply(listofprobesandsequences, write, file=paste("fasta_files_chromosome_probes/listofprobesandsequences",b,".fa",sep=""), append=T);


#will eventually need to check if order is conserved between all quadrants in a pairwise manner
#all unique subsets of 2 elements from a set of 4 elements: 6 total comparisons



#     ____________________________________________________________________________
#     OTHER STUFF                                                                                                                         ####



length(list_bin_comparison)
nrow(onebq1_comb)

onebq1_comb[match(onebq1_comb[onebq1_comb$bins==list_bins[4],4][5], onebq1_comb[,11]), 9]
onebq1_comb[onebq1_comb$bins==list_bins[4],4]

?add_column
nrow(onebq1)
nrow(onebq2)




#make dataframe of the index at which marker from 1 dataset occurs in the other dataset
orderlist=list()
for(i in 1:nrow(onebq1)){
    numsame=as.numeric(which(onebq1[,4][i] == onebq2[,4]))
    orderlist=c(orderlist,numsame)
}
ba=data.frame(orderlist)
ba=rbind(ba,ba)
ba[2,]=""


#populate dataframe with ordinal indicator values - S if marker order is preserved between datasets, N if it is not.
counter=1
while(counter<ncol(ba)){
    if(as.numeric(ba[1,counter]) < as.numeric(ba[1,(counter+1)])){
        ba[2,counter]="S"
    } else {
        ba[2,counter]="N"
    }
    counter=counter+1
}


colnames(ba)=1:ncol(ba)
ba=t(ba)

table(ba[,2])






#     ____________________________________________________________________________
#     dist vs gap size graph                                                                                                    ####

path = function(x) read.csv(paste("processed_data/genetic_maps/q", x, "_multipoint_ordercons_asmap_est.csv", sep = ""),
                                                        header = T, stringsAsFactors = F)
genetic.maps = list(path(1), path(2), path(3), path(4))

#add dist_adj col
genetic.maps = lapply(genetic.maps, function(x){
    x$dist_adj = ""
    for(i in 2:(nrow(x))){
        x$dist_adj[i] = (x$pos[i] - x$pos[(i-1)])
    }
    list.of.lg = unique(genetic.maps[[1]]$lg)
    start.of.lg.coords = match(list.of.lg, genetic.maps[[1]]$lg)
    x$dist_adj[start.of.lg.coords] = 0
    x$dist_adj = as.numeric(x$dist_adj)
    return(x)
})



plot.chr = function(map, chromo) filter(genetic.maps[[map]], lg == chromo)

windows(2,2)
par(mfrow=c(2,2))
chr.to.test = "LG15" #note, this variable is sensitive to the linkage group format in the genetic map (e.g. either "LGX" or simply "X")
lapply(1:4, function(x){
    barplot(plot.chr(x, chr.to.test)$dist_adj, names.arg = plot.chr(x, chr.to.test)$marker, las = 2, cex.names = 0.5, ylim = c(0, 35))
})

#add quad column
allgeneticmaps = newdf(c(colnames(genetic.maps[[1]]), "quad"))
for(i in 1:4){
    genetic.maps[[i]]$quad = i
    allgeneticmaps = rbind(allgeneticmaps, genetic.maps[[i]])
}

#remove blank first row
allgeneticmaps = allgeneticmaps[-1,]

ggplot(data = allgeneticmaps, aes(x = marker, y = dist_adj, group = quad, colour = quad)) + geom_line()
