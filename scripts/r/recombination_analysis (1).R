setwd("~/Google Drive/R/scripts/rotation1scripts/")
library(dplyr)
library(tibble)

alexq1=read.csv("alexq1_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)
alexq2=read.csv("alexq2_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)
alexq3=read.csv("alexq3_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)
alexq4=read.csv("alexq4_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)


# all(colnames(alexq1)==colnames(alexq2)) #check these genotype files are based on the same genetic map

quadrant.to.analyse=alexq4

#setup dataframe - for counting the number of consecutive genotypes that are the same
genotypelist=as.data.frame(matrix(nrow=1, ncol=5))
genotypelist[is.na(genotypelist)]=""
colnames(genotypelist)=c("genotype","number_consec","lg", "individual", "individual_id")

geno_count=1
for(p in 4:nrow(quadrant.to.analyse)){
    
    #setup init. variables
    first=quadrant.to.analyse[p,3]
    genotypelist$genotype[geno_count]=first
    genotypelist$number_consec[geno_count]=1
    genotypelist$lg[geno_count]=quadrant.to.analyse[1,3]
    genotypelist$individual[geno_count]=quadrant.to.analyse[p,2]
    genotypelist$individual_id[geno_count]=p
    firstlg=quadrant.to.analyse[1,3]
    
    for(i in 4:ncol(quadrant.to.analyse)){
        if(quadrant.to.analyse[p,i] == first & quadrant.to.analyse[1,i] == firstlg){
            genotypelist$number_consec[geno_count]=(as.numeric(genotypelist$number_consec[geno_count])+1)
        } else {
        firstlg=quadrant.to.analyse[1,i]
        genotypelist=add_row(genotypelist)
        genotypelist[nrow(genotypelist),]=""
        
        geno_count=geno_count+1
        first=quadrant.to.analyse[p,i]
        genotypelist$genotype[geno_count]=first
        genotypelist$number_consec[geno_count]=1
        genotypelist$lg[geno_count]=quadrant.to.analyse[1,i]
        genotypelist$individual[geno_count]=quadrant.to.analyse[p,2]
        genotypelist$individual_id[geno_count]=p
        }
    }
    genotypelist=add_row(genotypelist)
    geno_count=nrow(genotypelist)
}

#remove NoCalls
genolist2=genotypelist[-which(genotypelist$genotype=="-"),]
genolg1=filter(genolist2, lg==1)

#concatenate consecutive genotypes within LGs within individuals that are the same (i.e. previously separated by NoCalls / "-")
#NOTE: using only first LG at the moment to test
for(i in 1:(nrow(genolg1)-1)){
    if(genolg1[i,1] == genolg1[(i+1),1] 
         & genolg1$lg[i] == genolg1$lg[(i+1)] 
         & genolg1$individual[i] == genolg1$individual[(i+1)]){
        genolg1[i,2]=(as.numeric(genolg1[i,2])+as.numeric(genolg1[(i+1),2]))
        genolg1[(i+1),]="del"
        }
}

genolg1=genolg1[-which(genolg1[,1]=="del"),]

#count recombination events and note where they occur
#first, setup a cumulative count column
genolg1=add_column(genolg1, .after="number_consec", cumulative_count="")
for(i in unique(as.numeric(genolg1$individual_id))){
    tempgeno=filter(genolg1, individual_id==i)
    tempgeno$cumulative_count=cumsum(tempgeno$number_consec)
    genolg1[genolg1$individual_id==i,][,3]=tempgeno$cumulative_count
}

recombination_counts=as.data.frame(matrix(ncol=5,nrow=1))
recombination_counts[is.na(recombination_counts)]=""
colnames(recombination_counts)=c("individual","recomb_event_id","position","event_type","num_recomb")
counter=1
for(i in unique(as.numeric(genolg1$individual_id))){
    tempgeno=filter(genolg1, individual_id==i)
    for(g in 1:(nrow(tempgeno))){
    if(g != nrow(tempgeno)){
        recombination_counts$individual[counter]=i
        recombination_counts$recomb_event_id[counter]=g
        recombination_counts$position[counter]=(as.numeric(tempgeno$cumulative_count[g])+1)
        recombination_counts$event_type[counter]=paste(tempgeno$genotype[g], "to", tempgeno$genotype[(g+1)])
        recombination_counts=add_row(recombination_counts)
        counter=counter+1
    }
    }
}

hist(as.numeric(recombination_counts$position))



#     ____________________________________________________________________________
#     PARSE GENOTYPE INFO. FOR 1 INDIVIDUAL (useful for debugging loop above) ####


#setup dataframe - for counting the number of consecutive genotypes that are the same
genotypelist=as.data.frame(matrix(nrow=1, ncol=5))
genotypelist[is.na(genotypelist)]=""
colnames(genotypelist)=c("genotype","number_consec","lg", "individual", "individual_id")

#setup init. variables
individual.index.number=5
geno_count=1
first=quadrant.to.analyse[individual.index.number,3]
firstlg=quadrant.to.analyse[1,3]
genotypelist$genotype[geno_count]=first
genotypelist$number_consec[geno_count]=1
genotypelist$lg[geno_count]=quadrant.to.analyse[1,3]
genotypelist$individual[geno_count]=quadrant.to.analyse[individual.index.number,2]
genotypelist$individual_id[geno_count]=individual.index.number

listof.lg=unique(as.numeric(quadrant.to.analyse[1, 3:ncol(quadrant.to.analyse)]))

for(i in 4:ncol(quadrant.to.analyse)){
    if(quadrant.to.analyse[individual.index.number,i] == first & quadrant.to.analyse[1,i] == firstlg){
        genotypelist$number_consec[geno_count]=(as.numeric(genotypelist$number_consec[geno_count])+1)
    } else {
        firstlg=quadrant.to.analyse[1,i]
        genotypelist=add_row(genotypelist)
        genotypelist[nrow(genotypelist),]=""
        
        geno_count=geno_count+1
        first=quadrant.to.analyse[individual.index.number,i]
        
        genotypelist$genotype[geno_count]=first
        genotypelist$number_consec[geno_count]=1
        genotypelist$lg[geno_count]=quadrant.to.analyse[1,i]
        genotypelist$individual[geno_count]=quadrant.to.analyse[individual.index.number,2]
        genotypelist$individual_id[geno_count]=individual.index.number
    }
}

listof.lg=unique(as.numeric(genotypelist$lg))
listof.lg
