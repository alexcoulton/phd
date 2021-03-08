setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#setwd("~/Google Drive/R/scripts/rotation1scripts/")
library(qtl)
library(dplyr)
library(tibble)
load("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts/axiomrqtl5.RData") #load rqtl data
cereals=read.csv("C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/cerealsdb/axiom820_2.csv", header=T, stringsAsFactors = F)
cereals$affycode=gsub('-','.',cereals$affycode) #fix affycode field so matches with flipped_data

multirf="q4_0.2_postsplit" #used to name files

#######################################################################
#########PROCESS MULTIPOINT OUTPUT FILE INTO USABLE DATAFRAME##########
#######################################################################
#aka. "Result file_Sk&Ext.txt"
mult_dat=read.table("mp_results_q4_rf0.2_presplit.txt",sep="\n", stringsAsFactors = F)

mult_dat2=cbind(mult_dat,mult_dat,mult_dat,mult_dat)

length(strsplit(mult_dat[10,], split="\\s+"))

#put data into multiple columns
for(i in 1:nrow(mult_dat)){
    b=strsplit(mult_dat[i,], split="\\s+")
        for(l in 1:length(b[[1]])){
            mult_dat2[i,l]=b[[1]][l]
    }
}

#remove blank rows
rowtodel=vector()
for(i in 1:nrow(mult_dat)){
    if(((grepl("", mult_dat2[i,1])==T)==T & (grepl("\\s+",mult_dat2[i,2])==T)==T & (grepl("\\s+",mult_dat2[i,3])==T)==T & (grepl("\\s+",mult_dat2[i,4])==T)==T)==T){
        rowtodel=c(rowtodel, i)
    }
}
rowtodel
mult_dat3=mult_dat2[-c(rowtodel),]

#populate LG column
for(i in 2:nrow(mult_dat3)){
    if((grepl("LG\\d{1,2}",mult_dat3[(i-1),1]) == T) & mult_dat3[i,1]==""){
        mult_dat3[i,1]=mult_dat3[(i-1),1]
    }
}

rowstodel=vector()
for(i in 1:nrow(mult_dat3)){
    if((mult_dat3[i,1] == mult_dat3[i,2]) & (mult_dat3[i,1] == mult_dat3[i,3])){
        rowstodel=c(rowstodel,i)
    }
}
mult_dat3=mult_dat3[-c(rowstodel),]

##########################################################
#GRAB LIST OF MARKERS IN MULTIPOINT LINKAGE GROUP mlg#####
##########################################################
######## HAVING PROBLEMS WITH THIS CODE SO COMMENTED OUT (droplevels no longer required, probably because I added stringsAsFactors=F to cereals data import##############
# mlg="LG4"
# listofmlg=unique(mult_dat3[,1])
# #COMPARE TO CEREALSDB DATA
# for(i in listofmlg){
# bn=mult_dat3[mult_dat3=="LG1",][,4]
# t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group
# print(table(droplevels(t$consensus)))
# }

##############################################################
#######STAGE 2 - INVESTIGATE A PARTICULAR LINKAGE GROUP#######
##############################################################
#compare a single LG to cerealsdb
listofmlg=unique(mult_dat3[,1])
listofmlg
mlg="LG31"
bn=mult_dat3[mult_dat3==mlg,][,c(2,4)]
bnrow=nrow(mult_dat3[mult_dat3==mlg,])
t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group

#not all of the markers in linkage group data from multiqtl are in cerealsdb...
#grab those that are
bn=bn[c(which(bn[,2] %in% cereals[,1])),]

#COMPARE MARKERS TO CEREALSDB, OUTPUT MATRIX WITH INDIVIDUAL MARKER COMPARISONS IN ORDER OF LG. USE TO SEE WHERE TO SPLIT BIG CHROMOSOMES
matofcons=matrix(ncol=3,nrow=bnrow)
matofcons=data.frame(matofcons)
matofcons <- data.frame(lapply(matofcons, as.character), stringsAsFactors=FALSE)

for(i in 1:nrow(bn)){
    g=which(bn[i,2] == cereals[,1]) #NB. this works but could of used the match() function instead
    matofcons[i,1]=bn[i,2]
    matofcons[i,2]=bn[i,1]
    matofcons[i,3]=cereals$consensus[g]
}

barplot(table(matofcons[,3])) #check barplot of the comparison


##############################################################################################################
###CAUTION: NEED TO SUPPLY COORD OF SPLIT AFTER MANUAL INSPECTION OF matofcons, DON'T JUST RUN THIS CODE######
##############################################################################################################
splitindex=251 #split is initiated before this marker (index of matofcons)

#AFTER INSPECTION OF matofcons, SPLIT LINKAGE GROUPS IN mult_dat3
#give marker name to split before
mtosplit=matofcons[splitindex,1]
mcoord=which(mult_dat3[mult_dat3==mlg,][,4]==mtosplit) #grab index of that marker in mult_dat3
endcoord=nrow(mult_dat3[mult_dat3==mlg,]) #grab end of linkage group in mult_dat3
mult_dat3[mult_dat3==mlg,][,1][mcoord:endcoord] = paste("LG", (length(listofmlg)+1),sep="") #split linkage group into two (second half is assigned
#an LG number one more than the biggest in the whole dataset)
newlgname=paste("LG", (length(listofmlg)+1),sep="")
listofmlg=c(listofmlg, newlgname) #append new LG to listofmlg
#also need to repair distances...
#first marker in LG is 0 cM.
#simply subtract the value of the first marker from all markers in the new LG.
mult_dat3[,2]=as.vector(sapply(mult_dat3[,2], as.numeric)) #first, convert column to numeric
mult_dat3[mult_dat3==newlgname,][,2]=(mult_dat3[mult_dat3==newlgname,][,2]-mult_dat3[mult_dat3==newlgname,][,2][1])


#################################################################################################################
###########TO PERFORM A SECOND SPLIT, GO UP TO STAGE 2 AND CHANGE mlg TO DESIRED LINKAGE GROUP###################
#################################################################################################################
quadrant="q4"

write.csv(mult_dat3, paste(multirf,"_ord_split.csv",sep=""), row.names = F)

#SETUP GENOTYPE INFORMATION WITH CLUSTERING AND ORDERING PROVIDED BY MULTIPOINT. I.E. USE mult_dat3 to organise alexq1.csv
#ignore variable names here --- only quadrant variable specified above matters 
alexq1=read.csv(paste("alex",quadrant,".csv",sep=""), header=T, stringsAsFactors = F)
alexq1=alexq1[,-c((which(!(colnames(alexq1)[3:ncol(alexq1)] %in% mult_dat3[,4]))+2))] #delete markers in alexq1 that arn't in mult_dat3

#match returns position of the first argument in the second argument
alexq1org=alexq1
for(i in 1:nrow(mult_dat3)){
    ind=(i+2) #don't alter first two columns of alexq1org
    pos=match(mult_dat3[,4][i],colnames(alexq1))
    alexq1org[,ind]=alexq1[,pos]
    colnames(alexq1org)[ind]=colnames(alexq1)[pos]
    alexq1org[1,ind]=mult_dat3[i,1]
}

alexq1org[1,]=gsub("LG(\\d)","\\1",alexq1org[1,]) #rename linkage group row in format suitable for rQTL \\d seems to match multiple digits... why?
alexq1org[1,1]=""

#write and load into rQTL, plot rf
write.csv(alexq1org, paste("alex",quadrant,"_",multirf,"geno.csv",sep=""), row.names=F)
alexq1org_rqtl=read.cross("csv", "./", "alexq1org_multi3.33_split.csv", estimate.map=F)
summary(alexq1org_rqtl)
plotRF(alexq1org_rqtl, alternate.chrid = T)


#     ____________________________________________________________________________
#     ADD CEREALSDB MARKER INFO TO MULT_DAT3                                                                    ####

mult_dat4=add_column(mult_dat3, cerealsdb="")
bn=mult_dat4[,4]
t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group
for(i in 1:length(bn)){
    pos=match(bn[i],t$affycode)
    mult_dat4[i,5]=t$consensus[pos]
}
write.csv(mult_dat4, "mult_dat4.csv",row.names = F)


###################################################################
###SEARCH WITHIN EACH LG IN mult_dat3 FOR cM GAPS BIGGER THAN X####
###################################################################

###INDICATIVE OF LGs THAT SHOULD BE SPLIT. OUTPUT = newmat2
X = 20 #size of cM gap to search for
newmat=as.data.frame(matrix(nrow=1,ncol=3))
colnames(newmat)=c("lg","pos_in_lg","centi")
newmat2=newmat
for(i in listofmlg){
    searchframe=mult_dat3[mult_dat3==i,]
    for(p in 1:(nrow(searchframe)-1)){
            if(searchframe[p,2] < (searchframe[(p+1),2]-X)){
            newmat[1,1]=i; newmat[1,2]=p; newmat[1,3]=searchframe[p,2]
            newmat2=rbind(newmat2,newmat)
        }
    }
}

############################################################################################################################################
######ITERATIVE SCRIPT --- OUTPUTS CHROMOSOME HISTOGRAM PLOTS AS WELL AS EXCEL SPREADSHEET CONTAINING CHROMOSOME DIST. INFO FOR EACH LG#####
############################################################################################################################################
#POPULATE MYMAT
dformat="multipoint"

    #SETUP DATAFRAME FOR POPULATION WITH CHROMOSOME COUNT DATA FOR ITERATIONS ON LINKAGE GROUP, LOD AND RF
    #####
    mymat=matrix(data=NA,nrow=2000,ncol=52)
    mymat=data.frame(mymat)
    mymat[1,1]="linkage group"
    cerealstitles=data.frame(table(cereals$consensus))
    cerealstitles <- rapply(cerealstitles, as.character, classes="factor", how="replace")
    cerealstitles[1,1]="BLANK"
    cerealstitles[,1]
    length(cerealstitles[,1])
    mymat[1,2:(length(cerealstitles[,1])+1)]=cerealstitles[,1]
    mymat=mymat[,-c(ncol(mymat),(ncol(mymat)-1))]
    colnames(mymat)=mymat[1,]
    mymat[1,]=NA
    mymat[is.na(mymat)]=""
    mymat=sapply(mymat, as.numeric)
    #supply the function with a character vector of marker names as well as the data source (either rQTL or multipoint)
    if(dformat=="rqtl"){
        numchromosomes=(length(table(lg[,2]))-10)
    } else if(dformat=="multipoint"){
        numchromosomes=length(listofmlg)
    }
    



    for(p in 1:length(listofmlg)){
            bn=mult_dat3[mult_dat3==listofmlg[p],][,4]
            t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group

        
        #SAVE CHROMOSOME COUNT BARPLOTS FOR EACH LINKAGE GROUP
        print("creating directories")
        if(dir.exists(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, "_rf_", multirf, sep=""))==F){
            dir.create(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, "_rf_", multirf, sep=""))
        }
        if(dir.exists(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, "_rf_", multirf, sep=""))==F){
            dir.create(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, "_rf_", multirf, sep=""))
        }
        png(filename=paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, "_rf_", multirf, "/","chromo_", listofmlg[p], ".png", sep=""))
        barplot(table(t$consensus))
        dev.off()
        
        #make dataframe containing table output
        mydat=data.frame(table(t$consensus), stringsAsFactors = F) #tabulate counts of consensus. NB: droplevels() removes zero counts
        mydat <- rapply(mydat, as.character, classes="factor", how="replace")
        
        if(mydat[1,1]==""){ mydat[1,1]="BLANK"}
        #populate mymat with chromosome frequency data for linkage group p
        for(x in 1:nrow(mydat)){
            coord=which(colnames(mymat)==mydat[x,1])
            mymat[p,coord]=mydat[x,2]
        }
    }
    
    mymat=cbind(mymat[,1],mymat)
    colnames(mymat)[1]="linkage group"
    colnames(mymat)[2]="most frequent"
    
    
    for(i in 1:numchromosomes){
        gogo=data.frame(na.exclude(mymat[i,4:(ncol(mymat)-1)]))
        biggest=max(as.numeric(as.character(gogo[,1])))
        mymat[i,2]=paste(rownames(data.frame(which(mymat[i,]==biggest))), collapse=" ")
    }
    
    mymat=cbind(mymat[,1],mymat)
    colnames(mymat)[1]="linkage group"
    colnames(mymat)[2]="totalinlg"
    
    for(i in 1:numchromosomes){
        gogo=data.frame(na.exclude(mymat[i,4:(ncol(mymat)-1)]))
        total=sum(as.numeric(as.character(gogo[,1])))
        mymat[i,2]=total
    }    
    
    for(i in 1:numchromosomes){
        mymat[1:numchromosomes,1]=listofmlg
    }

    if(dir.exists(paste("linkage_grouping_it_quad",quadrant, "_form_", dformat, sep=""))==F){
        dir.create(paste("linkage_grouping_it_quad",quadrant, "_form_", dformat, sep=""))
    }
    write.csv(mymat,paste("linkage_grouping_it_quad",quadrant, "_form_", dformat,"/lgroups_","_multirf_",multirf,".csv",sep=""), row.names=F)

 

#     ____________________________________________________________________________
#     ASSIGN CHROMOSOMES TO EXISTING PRELIMINARY MAPS                                                 ####
#     Repeat of above code but wo/ iterative saving of plots etc.

    
    
#mult_dat3=read.csv("q1_rf0.33postmulti_ord_split.csv", header=T,stringsAsFactors = F) #load first quadrant
#mult_dat3=read.csv("q2_0.31_postsplit_ord_split.csv", header=T,stringsAsFactors = F) #load second quadrant
mult_dat3=read.csv("q3_0.2_postsplit_ord_split.csv", header=T,stringsAsFactors = F) #load third quadrant
    
listofmlg=unique(mult_dat3[,1])
    
dformat="multipoint"

#SETUP DATAFRAME FOR POPULATION WITH CHROMOSOME COUNT DATA FOR ITERATIONS ON LINKAGE GROUP, LOD AND RF
####
mymat=matrix(data=NA,nrow=2000,ncol=52)
mymat=data.frame(mymat)
mymat[1,1]="linkage group"
cerealstitles=data.frame(table(cereals$consensus))
cerealstitles <- rapply(cerealstitles, as.character, classes="factor", how="replace")
cerealstitles[1,1]="BLANK"
cerealstitles[,1]
length(cerealstitles[,1])
mymat[1,2:(length(cerealstitles[,1])+1)]=cerealstitles[,1]
mymat=mymat[,-c(ncol(mymat),(ncol(mymat)-1))]
colnames(mymat)=mymat[1,]
mymat[1,]=NA
mymat[is.na(mymat)]=""
mymat=sapply(mymat, as.numeric)
#supply the function with a character vector of marker names as well as the data source (either rQTL or multipoint)
if(dformat=="rqtl"){
    numchromosomes=(length(table(lg[,2]))-10)
} else if(dformat=="multipoint"){
    numchromosomes=length(listofmlg)
}




for(p in 1:length(listofmlg)){
    bn=mult_dat3[mult_dat3==listofmlg[p],][,4]
    t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group
    
    #make dataframe containing table output
    mydat=data.frame(table(t$consensus), stringsAsFactors = F) #tabulate counts of consensus. NB: droplevels() removes zero counts
    mydat <- rapply(mydat, as.character, classes="factor", how="replace")
    
    if(mydat[1,1]==""){ mydat[1,1]="BLANK"}
    #populate mymat with chromosome frequency data for linkage group p
    for(x in 1:nrow(mydat)){
        coord=which(colnames(mymat)==mydat[x,1])
        mymat[p,coord]=mydat[x,2]
    }
}

mymat=cbind(mymat[,1],mymat)
colnames(mymat)[1]="linkage group"
colnames(mymat)[2]="most frequent"


for(i in 1:numchromosomes){
    gogo=data.frame(na.exclude(mymat[i,4:(ncol(mymat)-1)]))
    biggest=max(as.numeric(as.character(gogo[,1])))
    mymat[i,2]=paste(rownames(data.frame(which(mymat[i,]==biggest))), collapse=" ")
}

mymat=cbind(mymat[,1],mymat)
colnames(mymat)[1]="linkage group"
colnames(mymat)[2]="totalinlg"

for(i in 1:numchromosomes){
    gogo=data.frame(na.exclude(mymat[i,4:(ncol(mymat)-1)]))
    total=sum(as.numeric(as.character(gogo[,1])))
    mymat[i,2]=total
}    

for(i in 1:numchromosomes){
    mymat[1:numchromosomes,1]=listofmlg
}


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### automatic assignment of chromosomes                                                                         ####

mymat[is.na(mymat)==T]=""
mult_dat3=add_column(mult_dat3, chr="")

#this loop examines the most frequent chromosome assignment, and if that happens to be an individual chromosome arm, examines the corresponding matching arm (i.e.
#if most freq. is 6BL, it checks 6BS). 
#If the corresponding matching arm has a count higher than 0 in the LG,
#the loop truncates the most frequent chromosome assignment to just the chromosome number and subgenome, then assigns it to the LG
#NB. the loop ignores cases in which there are two chromosomes in the most frequent column of mymat. These should be inspected manually.

#NB this loop was developed after q1 and q2 were processed manually, so was first used on q3.

for(i in unique(mymat[,1])){
    b=0
    part_no=0
    mf=mymat[mymat[,1]==i,][3] #grab most freq. LG from mymat

    if(grepl("S",mf)==T){
        b=nchar(as.character(grep("S",mf,value=T)[[1]]))
        tok="S"
        most=grep("S",mf,value=T)[[1]]
        partner=paste(substr(most,1,2),"L",sep="")
        if(b==3){
            if(mymat[mymat[,1]==i,colnames(mymat)==partner][[1]] != "") part_no=as.numeric(mymat[mymat[,1]==i,colnames(mymat)==partner][[1]])
        }
    }
    
    if(grepl("L",mf)==T){
        b=nchar(as.character(grep("L",mf,value=T)[[1]]))
        
        tok="L"
        most2=grep("L",mf,value=T)[[1]]
        partner=paste(substr(most2,1,2),"S",sep="")
        
        if(b==3){
            if(mymat[mymat[,1]==i,colnames(mymat)==partner][[1]] != "") part_no=as.numeric(mymat[mymat[,1]==i,colnames(mymat)==partner][[1]])
        }
    }
    
    if(mf != ""){
        if(part_no > 0){
            mult_dat3[mult_dat3[,1]==i,][,5]=substr(mf,1,2)
        } else {
            mult_dat3[mult_dat3[,1]==i,][,5]=mf
        }
    } 
}





### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### manual assignment of chromosomes (q1)                                                                     ####

mult_dat3[mult_dat3[,1]=="LG1",][,5]="2A"
mult_dat3[mult_dat3[,1]=="LG2",][,5]="3D"
mult_dat3[mult_dat3[,1]=="LG5",][,5]="7B"
mult_dat3[mult_dat3[,1]=="LG25",][,5]="5D"

#replacing arm notation when LGs clearly represent a whole chromosome
mult_dat3[mult_dat3[,1]=="LG3",][,5]="3A"
mult_dat3[mult_dat3[,1]=="LG34",][,5]="5D"
mult_dat3[mult_dat3[,1]=="LG4",][,5]="5B"
mult_dat3[mult_dat3[,1]=="LG6",][,5]="7A"
mult_dat3[mult_dat3[,1]=="LG7",][,5]="4A"
mult_dat3[mult_dat3[,1]=="LG8",][,5]="1A"
mult_dat3[mult_dat3[,1]=="LG9",][,5]="1D"
mult_dat3[mult_dat3[,1]=="LG10",][,5]="2D"
mult_dat3[mult_dat3[,1]=="LG11",][,5]="6B"
mult_dat3[mult_dat3[,1]=="LG12",][,5]="6A"
mult_dat3[mult_dat3[,1]=="LG13",][,5]="4D"
mult_dat3[mult_dat3[,1]=="LG12",][,5]="6A"
mult_dat3[mult_dat3[,1]=="LG20",][,5]="1B"
mult_dat3[mult_dat3[,1]=="LG21",][,5]="2B"
mult_dat3[mult_dat3[,1]=="LG22",][,5]="6D"


#     ____________________________________________________________________________
#     manual assignment of chromosomes (q2)                                                                     ####

mult_dat3[mult_dat3[,1]=="LG1",][,5]="4A"
mult_dat3[mult_dat3[,1]=="LG2",][,5]="2D"
mult_dat3[mult_dat3[,1]=="LG3",][,5]="5A"
mult_dat3[mult_dat3[,1]=="LG4",][,5]="5B"
mult_dat3[mult_dat3[,1]=="LG6",][,5]="2B"
mult_dat3[mult_dat3[,1]=="LG9",][,5]="3A"
mult_dat3[mult_dat3[,1]=="LG13",][,5]="3A"
mult_dat3[mult_dat3[,1]=="LG14",][,5]="2D"
mult_dat3[mult_dat3[,1]=="LG15",][,5]="6B"
mult_dat3[mult_dat3[,1]=="LG16",][,5]="7A"
mult_dat3[mult_dat3[,1]=="LG17",][,5]="7A"
mult_dat3[mult_dat3[,1]=="LG18",][,5]="1A"
mult_dat3[mult_dat3[,1]=="LG19",][,5]="6A"
mult_dat3[mult_dat3[,1]=="LG21",][,5]="1B"
mult_dat3[mult_dat3[,1]=="LG22",][,5]="3D"
mult_dat3[mult_dat3[,1]=="LG26",][,5]="7D"
mult_dat3[mult_dat3[,1]=="LG27",][,5]="4D"
mult_dat3[mult_dat3[,1]=="LG31",][,5]="3D"
mult_dat3[mult_dat3[,1]=="LG32",][,5]="7B"
mult_dat3[mult_dat3[,1]=="LG33",][,5]="2A"
mult_dat3[mult_dat3[,1]=="LG36",][,5]="5D"


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### manual assignment of chromosomes (q3)                                                                     ####


#check which LGs need to be tidied up manually
lgs=unique(mult_dat3[,1])
lc=list()
for(i in lgs){
    print(c(i,unique(mult_dat3[mult_dat3[,1]==i,][,5])))
    lc=c(unique(mult_dat3[mult_dat3[,1]==i,][,5]), lc)
}










lc=as.data.frame(lc)
lc=t(lc)
lc

#make a mymat temp file
write.csv(mymat,"mymattemp.csv",row.names=F)

#write.csv(mult_dat3, "q1_rf0.33postmulti_ord_split_chromosomes.csv", row.names=F) #write csv for q1
write.csv(mult_dat3, "q2_0.31_postsplit_ord_split_chromosomes.csv", row.names=F) #write csv for q2

