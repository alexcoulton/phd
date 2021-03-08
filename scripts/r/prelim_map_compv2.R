# q1=read.csv("man_cur_skq1_ord_split_all.csv", header=T, stringsAsFactors=F) #read csv for q1
# q2=read.csv("man_cur_skq2_ord_split_all.csv", header=T, stringsAsFactors=F) #read csv for q2
# q3=read.csv("man_cur_skq3_ord_split_all.csv", header=T, stringsAsFactors=F) #read csv for q3
# q4=read.csv("man_cur_skq4_ord_split_all.csv", header=T, stringsAsFactors=F) #read csv for q4

# install.packages("devtools")
# devtools::install_github("lorenzwalthert/strcode")
setwd("~/Google Drive/R/scripts/rotation1scripts/") #macbook settings
library(dplyr)
library(tibble)
library(ggplot2)

#load up each quadrant and assign bins (sequentially)

q1=read.csv("q1_asmap_map_w_chr.csv", header=T, stringsAsFactors = F)

#assign bins to each LG
for(g in unique(q1$lg)){
    counter=1
    for(i in unique(q1[q1$lg==g,2])){
        q1[q1$pos==i,3]=counter
        counter=counter+1
    }
}

q2=read.csv("q2_asmap_map_w_chr.csv", header=T, stringsAsFactors = F)

#reverse a particular linkage group --- see lab notebook page 16/12/2017
lgtoreverse="LG5"
q2[q2$lg==lgtoreverse,]=q2[q2$lg==lgtoreverse,][sort(which(q2[q2$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
q2[q2$lg==lgtoreverse,2]=-(q2[q2$lg==lgtoreverse,2]-max(q2[q2$lg==lgtoreverse,2]))
q2[q2$lg==lgtoreverse,]

#assign bins to each LG
for(g in unique(q2$lg)){
    counter=1
    for(i in unique(q2[q2$lg==g,2])){
        q2[q2$pos==i,3]=counter
        counter=counter+1
    }
}

q3=read.csv("q3_asmap_map_w_chr.csv", header=T, stringsAsFactors = F)

for(i in unique(q3$lg)){
    lgtoreverse=i
    q3[q3$lg==lgtoreverse,]=q3[q3$lg==lgtoreverse,][sort(which(q3[q3$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
    q3[q3$lg==lgtoreverse,2]=-(q3[q3$lg==lgtoreverse,2]-max(q3[q3$lg==lgtoreverse,2]))
    q3[q3$lg==lgtoreverse,]
}


#assign bins to each LG
for(g in unique(q3$lg)){
    counter=1
    for(i in unique(q3[q3$lg==g,2])){
        q3[q3$pos==i,3]=counter
        counter=counter+1
    }
}

q4=read.csv("q4_asmap_map_w_chr.csv", header=T, stringsAsFactors = F)

# reverse a particular linkage group
lgtoreverse="LG15"
q4[q4$lg==lgtoreverse,]=q4[q4$lg==lgtoreverse,][sort(which(q4[q4$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
q4[q4$lg==lgtoreverse,2]=-(q4[q4$lg==lgtoreverse,2]-max(q4[q4$lg==lgtoreverse,2]))
q4[q4$lg==lgtoreverse,]

for(i in unique(q4$lg)){
    lgtoreverse=i
    q4[q4$lg==lgtoreverse,]=q4[q4$lg==lgtoreverse,][sort(which(q4[q4$lg==lgtoreverse,1]==lgtoreverse),decreasing=T),]
    q4[q4$lg==lgtoreverse,2]=-(q4[q4$lg==lgtoreverse,2]-max(q4[q4$lg==lgtoreverse,2]))
    q4[q4$lg==lgtoreverse,]
}

#assign bins to each LG
for(g in unique(q4$lg)){
    counter=1
    for(i in unique(q4[q4$lg==g,2])){
        q4[q4$pos==i,3]=counter
        counter=counter+1
    }
}




#NEED TO PERFORM CHROMOSOME ASSIGNMENT WITH cerealscomp_multipointv2.R BEFORE USING THIS SCRIPT


#Initial analysis: compare the linkage groups between quadrants
#this redundant in the case of using ASMap, as all LGs are the same

##quadrant1
# s=as.data.frame(matrix(ncol=3,nrow=2000))
# s[is.na(s)==T]=""
# colnames(s)=c("quadrant","lg","nomarkers")
# 
# listlg=unique(q1$chr)
# for(i in 1:length(listlg)){
#     p=q1[q1[,5]==listlg[i],]
#     s[i,1]="1"
#     s[i,2]=listlg[i]
#     s[i,3]=nrow(p)
#     
# }
# 
# s=s[1:length(listlg),]
# s[,3]=as.numeric(s[,3])
# s=s[sort(s[,3],index.return=T, decreasing=T)[[2]],]
# rownames(s)=1:nrow(s)
# 
# q1s=s
# 
# ##quadrant2
# s=as.data.frame(matrix(ncol=3,nrow=2000))
# s[is.na(s)==T]=""
# colnames(s)=c("quadrant","lg","nomarkers")
# 
# listlg=unique(q2$chr)
# for(i in 1:length(listlg)){
#     p=q2[q2[,5]==listlg[i],]
#     s[i,1]="2"
#     s[i,2]=listlg[i]
#     s[i,3]=nrow(p)
#     
# }
# 
# s=s[1:length(listlg),]
# s[,3]=as.numeric(s[,3])
# s=s[sort(s[,3],index.return=T, decreasing=T)[[2]],]
# rownames(s)=1:nrow(s)
# 
# q2s=s

##quadrant3
# s=as.data.frame(matrix(ncol=3,nrow=2000))
# s[is.na(s)==T]=""
# colnames(s)=c("quadrant","lg","nomarkers")
# 
# listlg=unique(q3$chr)
# for(i in 1:length(listlg)){
#     p=q3[q3[,5]==listlg[i],]
#     s[i,1]="3"
#     s[i,2]=listlg[i]
#     s[i,3]=nrow(p)
#     
# }
# 
# s=s[1:length(listlg),]
# s[,3]=as.numeric(s[,3])
# s=s[sort(s[,3],index.return=T, decreasing=T)[[2]],]
# rownames(s)=1:nrow(s)
# 
# q3s=s

##quadrant4
# s=as.data.frame(matrix(ncol=3,nrow=2000))
# s[is.na(s)==T]=""
# colnames(s)=c("quadrant","lg","nomarkers")
# 
# listlg=unique(q4$chr)
# for(i in 1:length(listlg)){
#     p=q4[q4[,5]==listlg[i],]
#     s[i,1]="4"
#     s[i,2]=listlg[i]
#     s[i,3]=nrow(p)
#     
# }
# 
# s=s[1:length(listlg),]
# s[,3]=as.numeric(s[,3])
# s=s[sort(s[,3],index.return=T, decreasing=T)[[2]],]
# rownames(s)=1:nrow(s)
# 
# q4s=s

#sort by chromosome
# q1s=q1s[sort(q1s[,2],index.return=T)[[2]],]
# q2s=q2s[sort(q2s[,2],index.return=T)[[2]],]
# q3s=q3s[sort(q3s[,2],index.return=T)[[2]],]
# q4s=q4s[sort(q4s[,2],index.return=T)[[2]],]

#allq=rbind(q1s,q2s,q3s,q4s)
# allq=rbind(q1s,q2s)

# write.csv(allq, "allq_cur_ske_asmap_test.csv", row.names=F)

#check size of LGs for each population
#nrow(q1[q1$chr=="1AL",])
#nrow(q2[q2$chr=="1AL",])
# nrow(q3[q3$chr=="1AL",])
# nrow(q4[q4$chr=="1AL",])

#make new dataframes for comparison #IMPORTANT NB. FILTERING BY "chr" COLUMN CAN LEAD TO AN ERROR;
#TWO SEPERATE LGs WITH THE SAME CHROMOSOME DESIGNATION CAN BE INCLUDED
#THIS DOES HOWEVER YIELD CANDIDATES FOR MERGING


#     ____________________________________________________________________________
#     BEGINING OF PAIRWISE LG COMPARISON FOR Q1 AND Q2                                                ####


#begin main comparisons of bin and marker order
listof.lg=unique(q1$lg)
listof.lg[1]
#listof.selected.chr=c("6B","3D","5A","1AL","6A","7 7BL","4AL none","2","1B","5BL","3A","2D","7A","2B","4 4BL","5BS","6DL","1DL","2BS","BLANK 2DS")
for(b in listof.lg){
onebq1=q1[q1$lg==b,] #loop through list
onebq2=q2[q2$lg==b,]
    
# onebq1=q1[q1$lg=="LG1",] #choose chromosome manually
# onebq2=q2[q2$lg=="LG1",]

nrow(onebq1)
nrow(onebq2)

#check common markers between LGs
length(which(onebq1[,4] %in% onebq2[,4]))
length(which(onebq2[,4] %in% onebq1[,4]))

#check non-common markers between LGs
# which(!onebq1[,4] %in% onebq2[,4])


#assess where non-common markers are in the population (not finished)
# onebq1_noncommon=onebq1[which(!onebq1[,4] %in% onebq2[,4]),]
# onebq1_noncommon[,4] %in% q2[,4]
# q2[na.omit(match(onebq1_noncommon[,4], q2[,4])),] #check this rest of the dataset for the noncommon markers - here they are in another cluster which identifies to the partner arm


#grab only common markers between LGs
# onebq1=onebq1[which(onebq1[,4] %in% onebq2[,4]),]
# onebq2=onebq2[which(onebq2[,4] %in% onebq1[,4]),]

onebq1_comb=cbind(onebq1,onebq2)
onebq1_comb_chr_trunc=onebq1_comb[,c(6,12)]
onebq1_comb=onebq1_comb[,-c(6,12)]
#sorting rows by marker within bins
#first marker col index is 4
#first ske col index is 3
for(i in unique(onebq1_comb[,3])){
    sort.by.marker=as.vector(sort(onebq1_comb[onebq1_comb[,3]==i,4], index.return=T)[[2]])
    onebq1_comb[onebq1_comb[,3]==i,1:5]=onebq1_comb[onebq1_comb[,3]==i,1:5][sort.by.marker,]
    
}

#second marker col index is 11
#second ske col index is 8
for(i in unique(onebq1_comb[,8])){
    sort.by.marker=as.vector(sort(onebq1_comb[onebq1_comb[,8]==i,9], index.return=T)[[2]])
    onebq1_comb[onebq1_comb[,8]==i,6:10]=onebq1_comb[onebq1_comb[,8]==i,6:10][sort.by.marker,]
    
}


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

##bin assignment -- redundant for the asmap data
#
# skeleton_markercolumn_no=5
# newbin_column_no=7
# count=1 #give each bin a unique identifier
# for(i in 1:nrow(onebq1_comb)){
#     if(onebq1_comb[i,skeleton_markercolumn_no]=="S"){
#         cur_s=count
#         onebq1_comb[i,newbin_column_no]=cur_s
#         count=count+1
#     } else if(onebq1_comb[i,skeleton_markercolumn_no]=="B"){
#         onebq1_comb[i,newbin_column_no]=cur_s
#     }
# }
# 
# skeleton_markercolumn_no=10
# newbin_column_no=8
# count=1 #give each bin a unique identifier
# for(i in 1:nrow(onebq1_comb)){
#     if(onebq1_comb[i,skeleton_markercolumn_no]=="S"){
#         cur_s=count
#         onebq1_comb[i,newbin_column_no]=cur_s
#         count=count+1
#     } else if(onebq1_comb[i,skeleton_markercolumn_no]=="B"){
#         onebq1_comb[i,newbin_column_no]=cur_s
#     }
# }

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


listof.blastfiles=list.files("blast-files-LG-vs-genome/")



    

auto.chromosome.assignment=q1$chr[which(q1$lg==b)[1]]
chromosome_truncated=q1$chr_truncated[which(q1$lg==b)[1]]
print(c("chromosome assignment", auto.chromosome.assignment))

if(file.exists(paste("blast-files-LG-vs-genome/listofprobesandsequences",auto.chromosome.assignment,".fa.blast",sep=""))==T){
print(c(auto.chromosome.assignment, "exists"))
#chromo.blast.results=read.table(paste("blast-files-LG-vs-genome/listofprobesandsequences",b,".fa.blast",sep=""))
chromo.blast.results=read.table(paste("blast-files-LG-vs-genome/listofprobesandsequences",auto.chromosome.assignment,".fa.blast",sep=""))
nrow(chromo.blast.results[chromo.blast.results[,2]==paste("chr",chromosome_truncated,sep=""),])
chromo.blast.results.filtered=chromo.blast.results[chromo.blast.results[,2]==paste("chr",chromosome_truncated,sep=""),]

chromo.blast.results.filtered=filter(chromo.blast.results.filtered, V11<1e-24) #remove hits with low e-values
chromo.blast.results.filtered=distinct(chromo.blast.results.filtered, V1, .keep_all=T)    

#grab and process individual chromosome character lengths
chromo.character.counts=read.table("fullchromosomesequences/chromosome_character_counts.txt",sep="\n", stringsAsFactors = F)
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

write.csv(onebq1_comb, paste("pairwiseLGcomparisons/q1_q2/onebq1_comb_q1_q2","_",b,"_",auto.chromosome.assignment,".csv",sep=""), row.names = F)

}

onebq1_comb_q1_q2=onebq1_comb



} #temporary end of loop

#     ____________________________________________________________________________
#     BEGINING OF PAIRWISE LG COMPARISON FOR Q3 AND Q4                                                ####


#begin main comparisons of bin and marker order
listof.lg=unique(q1$lg)
listof.lg[1]
#listof.selected.chr=c("6B","3D","5A","1AL","6A","7 7BL","4AL none","2","1B","5BL","3A","2D","7A","2B","4 4BL","5BS","6DL","1DL","2BS","BLANK 2DS")
first.quad.to.compare=q1
first.quad.to.compare.name="q1"
second.quad.to.compare=q4
second.quad.to.compare.name="q4"
for(b in listof.lg){
    onebq1=first.quad.to.compare[first.quad.to.compare$lg==b,] #loop through list
    onebq2=second.quad.to.compare[second.quad.to.compare$lg==b,]
    
    # onebq1=q1[q1$lg=="LG1",] #choose chromosome manually
    # onebq2=q2[q2$lg=="LG1",]
    
    nrow(onebq1)
    nrow(onebq2)
    
    #check common markers between LGs
    length(which(onebq1[,4] %in% onebq2[,4]))
    length(which(onebq2[,4] %in% onebq1[,4]))
    
    #check non-common markers between LGs
    # which(!onebq1[,4] %in% onebq2[,4])
    
    
    #assess where non-common markers are in the population (not finished)
    # onebq1_noncommon=onebq1[which(!onebq1[,4] %in% onebq2[,4]),]
    # onebq1_noncommon[,4] %in% q2[,4]
    # q2[na.omit(match(onebq1_noncommon[,4], q2[,4])),] #check this rest of the dataset for the noncommon markers - here they are in another cluster which identifies to the partner arm
    
    
    #grab only common markers between LGs
    # onebq1=onebq1[which(onebq1[,4] %in% onebq2[,4]),]
    # onebq2=onebq2[which(onebq2[,4] %in% onebq1[,4]),]
    
    onebq1_comb=cbind(onebq1,onebq2)
    onebq1_comb_chr_trunc=onebq1_comb[,c(6,12)]
    onebq1_comb=onebq1_comb[,-c(6,12)]
    #sorting rows by marker within bins
    #first marker col index is 4
    #first ske col index is 3
    for(i in unique(onebq1_comb[,3])){
        sort.by.marker=as.vector(sort(onebq1_comb[onebq1_comb[,3]==i,4], index.return=T)[[2]])
        onebq1_comb[onebq1_comb[,3]==i,1:5]=onebq1_comb[onebq1_comb[,3]==i,1:5][sort.by.marker,]
        
    }
    
    #second marker col index is 11
    #second ske col index is 8
    for(i in unique(onebq1_comb[,8])){
        sort.by.marker=as.vector(sort(onebq1_comb[onebq1_comb[,8]==i,9], index.return=T)[[2]])
        onebq1_comb[onebq1_comb[,8]==i,6:10]=onebq1_comb[onebq1_comb[,8]==i,6:10][sort.by.marker,]
        
    }
    
    
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
    
    ##bin assignment -- redundant for the asmap data
    #
    # skeleton_markercolumn_no=5
    # newbin_column_no=7
    # count=1 #give each bin a unique identifier
    # for(i in 1:nrow(onebq1_comb)){
    #     if(onebq1_comb[i,skeleton_markercolumn_no]=="S"){
    #         cur_s=count
    #         onebq1_comb[i,newbin_column_no]=cur_s
    #         count=count+1
    #     } else if(onebq1_comb[i,skeleton_markercolumn_no]=="B"){
    #         onebq1_comb[i,newbin_column_no]=cur_s
    #     }
    # }
    # 
    # skeleton_markercolumn_no=10
    # newbin_column_no=8
    # count=1 #give each bin a unique identifier
    # for(i in 1:nrow(onebq1_comb)){
    #     if(onebq1_comb[i,skeleton_markercolumn_no]=="S"){
    #         cur_s=count
    #         onebq1_comb[i,newbin_column_no]=cur_s
    #         count=count+1
    #     } else if(onebq1_comb[i,skeleton_markercolumn_no]=="B"){
    #         onebq1_comb[i,newbin_column_no]=cur_s
    #     }
    # }
    
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
    
    
    listof.blastfiles=list.files("blast-files-LG-vs-genome/")
    
    
    
    
    
    auto.chromosome.assignment=q1$chr[which(q1$lg==b)[1]]
    chromosome_truncated=q1$chr_truncated[which(q1$lg==b)[1]]
    print(c("chromosome assignment", auto.chromosome.assignment))
    
    if(file.exists(paste("blast-files-LG-vs-genome/listofprobesandsequences",auto.chromosome.assignment,".fa.blast",sep=""))==T){
        print(c(auto.chromosome.assignment, "exists"))
        #chromo.blast.results=read.table(paste("blast-files-LG-vs-genome/listofprobesandsequences",b,".fa.blast",sep=""))
        chromo.blast.results=read.table(paste("blast-files-LG-vs-genome/listofprobesandsequences",auto.chromosome.assignment,".fa.blast",sep=""))
        nrow(chromo.blast.results[chromo.blast.results[,2]==paste("chr",chromosome_truncated,sep=""),])
        chromo.blast.results.filtered=chromo.blast.results[chromo.blast.results[,2]==paste("chr",chromosome_truncated,sep=""),]
        
        chromo.blast.results.filtered=filter(chromo.blast.results.filtered, V11<1e-24) #remove hits with low e-values
        chromo.blast.results.filtered=distinct(chromo.blast.results.filtered, V1, .keep_all=T)    
        
        #grab and process individual chromosome character lengths
        chromo.character.counts=read.table("fullchromosomesequences/chromosome_character_counts.txt",sep="\n", stringsAsFactors = F)
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
        
        write.csv(onebq1_comb, paste("pairwiseLGcomparisons/", first.quad.to.compare.name, "_", second.quad.to.compare.name, "/onebq1_comb", first.quad.to.compare.name, "_", second.quad.to.compare.name, "_",b,"_",auto.chromosome.assignment,".csv",sep=""), row.names = F)
        
    }
    
    onebq1_comb_q1_q2=onebq1_comb
    
    
    
} #temporary end of loop


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

setchrom="7A"

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### quad 1                                                                                                                                    ####


oneborigq1=q1[q1$chr==setchrom,]

#grab only skeleton markers and put in "new" dataframe
new=oneborigq1[1:2,]
for(i in 1:(nrow(oneborigq1)-1)){
    if(oneborigq1[i,2]!=oneborigq1[(i+1),2]){
        new=rbind(new, oneborigq1[i,])
    }
}

new=add_column(new, dist="")
new=new[-c(1,2),]

#populate distance column of new 
for(i in 1:(nrow(new)-1)){
    new$dist[i]=-(new[i,2]-new[i+1,2])
}

newq1=new
#plot centimorgan (x) vs dist (y)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### quad 2                                                                                                                                    ####


oneborigq2=q2[q2$chr==setchrom,]

#grab only skeleton markers and put in "new" dataframe
new=oneborigq2[1:2,]
for(i in 1:(nrow(oneborigq2)-1)){
    if(oneborigq2[i,2]!=oneborigq2[(i+1),2]){
        new=rbind(new, oneborigq2[i,])
    }
}

new=add_column(new, dist="")
new=new[-c(1,2),]

#populate distance column of new 
for(i in 1:(nrow(new)-1)){
    new$dist[i]=-(new[i,2]-new[i+1,2])
}

newq2=new
#plot centimorgan (x) vs dist (y)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### quad 3                                                                                                                                    ####


oneborigq3=q3[q3$chr==setchrom,]

#grab only skeleton markers and put in "new" dataframe
new=oneborigq3[1:2,]
for(i in 1:(nrow(oneborigq3)-1)){
    if(oneborigq3[i,2]!=oneborigq3[(i+1),2]){
        new=rbind(new, oneborigq3[i,])
    }
}

new=add_column(new, dist="")
new=new[-c(1,2),]

#populate distance column of new 
for(i in 1:(nrow(new)-1)){
    new$dist[i]=-(new[i,2]-new[i+1,2])
}

newq3=new
#plot centimorgan (x) vs dist (y)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### quad 4                                                                                                                                    ####


oneborigq4=q4[q4$chr==setchrom,]

#grab only skeleton markers and put in "new" dataframe
new=oneborigq4[1:2,]
for(i in 1:(nrow(oneborigq4)-1)){
    if(oneborigq4[i,2]!=oneborigq4[(i+1),2]){
        new=rbind(new, oneborigq4[i,])
    }
}

new=add_column(new, dist="")
new=new[-c(1,2),]

#populate distance column of new 
for(i in 1:(nrow(new)-1)){
    new$dist[i]=-(new[i,2]-new[i+1,2])
}

newq4=new
#plot centimorgan (x) vs dist (y)


#do some side by side plots
windows(2,2)
par(mfrow=c(2,2))
plot(newq1[,2],newq1[,6], xlim=c(0,250), ylim=c(0,20), main=paste("q1", setchrom, "(10 deg)"), xlab="Distance along chromosome (cM)", ylab="Distance between markers (cM)")
plot(newq3[,2],newq3[,6], xlim=c(0,250), ylim=c(0,20), main=paste("q3", setchrom, "(14 deg)"), xlab="Distance along chromosome (cM)", ylab="Distance between markers (cM)")
plot(newq2[,2],newq2[,6], xlim=c(0,250), ylim=c(0,20), main=paste("q2", setchrom, "(26 deg)"), xlab="Distance along chromosome (cM)", ylab="Distance between markers (cM)")
plot(newq4[,2],newq4[,6], xlim=c(0,250), ylim=c(0,20), main=paste("q4", setchrom, "(28 deg)"), xlab="Distance along chromosome (cM)", ylab="Distance between markers (cM)")
