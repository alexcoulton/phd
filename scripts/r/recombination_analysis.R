# setwd("~/Google Drive/R/scripts/rotation1scripts/") #mac
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts/") #windows
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)

alexq1=read.csv("alexq1_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)
alexq2=read.csv("alexq2_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)
alexq3=read.csv("alexq3_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)
alexq4=read.csv("alexq4_man_cur_sk_map_genogeno.csv", header=T, stringsAsFactors = F)

# testing whether the number of markers in pairwiseLGcomparisons and these maps are the same (they are)
# testcomparison=read.csv("pairwiseLGcomparisons/q1_q2/onebq1_combq1_q2_LG15_3D.csv", header=T, stringsAsFactors = F)
# length(alexq1[1,which(alexq1[1,]==15)])
# length(testcomparison$marker1)

# all(colnames(alexq1)==colnames(alexq2)) #check these genotype files are based on the same genetic map

quadrant.to.analyse=alexq1
quadrant.to.analyse.name="q1"

#setup dataframe - for counting the number of consecutive genotypes that are the same
makegenotypelist=function(quadrant.to.analyse){
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
    #remove last empty row
    genolist2=genolist2[-nrow(genolist2),]
    
    #concatenate consecutive genotypes within LGs within individuals that are the same (i.e. previously separated by NoCalls / "-")
    #NOTE: using only first LG at the moment to test
    counter=1
    while(counter < 20){ #run concatenation loop multiple times to remove multiple identical consecutive genotypes
        for(i in 1:(nrow(genolist2)-1)){
            if(genolist2[i,1] == genolist2[(i+1),1]
                 & genolist2[i,3] == genolist2[(i+1),3] #if lg is the same
                 & genolist2[i,4] == genolist2[(i+1),4]){ #if individual is the same
                genolist2[i,2]=(as.numeric(genolist2[i,2])+as.numeric(genolist2[(i+1),2]))
                genolist2[(i+1),]="del"
            }
        }
        if(length(which(genolist2[,1]=="del")) != 0){
            genolist2=genolist2[-which(genolist2[,1]=="del"),]
        }
        counter=counter+1
    }
    return(genolist2)
}

q1genotypelist=makegenotypelist(alexq1)
q2genotypelist=makegenotypelist(alexq2)
q3genotypelist=makegenotypelist(alexq3)
q4genotypelist=makegenotypelist(alexq4)

makehistlist=function(linkagegroup, genolist2, return_recombinationcounts){ 
    genolg1=filter(genolist2, lg==linkagegroup) #replace this to examine a different linkage group
    
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
    
    #assign number of recombinations to respective event types (a transition from B to A or vice versa is counted as 2 events)
    recombination_counts$num_recomb[which(recombination_counts$event_type=="B to A")]=2
    recombination_counts$num_recomb[which(recombination_counts$event_type=="A to B")]=2
    recombination_counts$num_recomb[which(recombination_counts$event_type=="H to A")]=1
    recombination_counts$num_recomb[which(recombination_counts$event_type=="H to B")]=1
    recombination_counts$num_recomb[which(recombination_counts$event_type=="A to H")]=1
    recombination_counts$num_recomb[which(recombination_counts$event_type=="B to H")]=1
    
    
    #produce weighted histogram based on these new values
    hist_list=list()
    for(i in 1:(nrow(recombination_counts)-1)){
        hist_list=c(hist_list, rep(recombination_counts$position[i], recombination_counts$num_recomb[i]))    
    }
    if(return_recombinationcounts == T){
        return(recombination_counts)
    } else {
        return(hist_list)
    }
}

length(unique(makehistlist(20, q1genotypelist, T)$individual))
length(unique(alexq1$probeset_id[4:nrow(alexq1)]))

getaverageposition=function(lg, genotypelist){
    avgevents=makehistlist(lg,genotypelist,T)
    listofavg.pos=vector()
    for(i in unique(avgevents$individual)){
        fil.avgevents=filter(avgevents, individual==i)
        fil.avgevents$position=as.numeric(as.character(fil.avgevents$position)) #convert position column to numeric
        averageposition=(sum(fil.avgevents$position)/nrow(fil.avgevents))
        listofavg.pos=c(listofavg.pos, averageposition)
    }
    return(listofavg.pos)
}

conv.vec.to.df=function(vector, quad, temp){
    df1=as.data.frame(vector)
    df1$quad=quad
    df1$temp=temp
    return(df1)
}

setupanova=function(lgtotest, return_df){
    avpos1=getaverageposition(lgtotest, q1genotypelist)
    avpos2=getaverageposition(lgtotest, q2genotypelist)
    avpos3=getaverageposition(lgtotest, q3genotypelist)
    avpos4=getaverageposition(lgtotest, q4genotypelist)
    
    df1=conv.vec.to.df(avpos1, 1, 10)
    df2=conv.vec.to.df(avpos2, 2, 26)
    df3=conv.vec.to.df(avpos3, 3, 14)
    df4=conv.vec.to.df(avpos4, 4, 28)
    
    alldf=rbind(df1,df2,df3,df4)
    alldf$quad=as.factor(alldf$quad)
    alldf$temp=as.factor(alldf$temp)
    alldf.lm=lm(alldf$vector ~ alldf$quad)
    
    if(return_df==T){
        return(alldf)
    } else {
    return(alldf.lm)
    }
}

lg1anova=setupanova(15, T)

lg1anova

hist(filter(lg1anova, quad==1)$vector[1:80])
hist(filter(lg1anova, quad==2)$vector[1:(nrow(filter(lg1anova,quad==2))-1)])
hist((filter(lg1anova, quad==3)$vector[1:(nrow(filter(lg1anova,quad==3))-1)]))
hist(filter(lg1anova, quad==4)$vector[1:(nrow(filter(lg1anova,quad==4))-1)])



a1=aov(lg1anova$vector ~ lg1anova$temp)
TukeyHSD(x=a1, 'lg1anova$temp', conf.level=0.95)


hist(as.numeric(makehistlist(linkagegroup = 15, genolist2 = q4genotypelist, return_recombinationcounts = F)))
mean(as.numeric(makehistlist(linkagegroup = 15, genolist2 = q1genotypelist, return_recombinationcounts = F)))





#counts the number of recombination events within a specified bin in the LG. Needs to be updated to take into account
#multiple recombination events, e.g. a transition from A to B = 2 events etc. 
countevents=function(lg, genotypelist, position_start, position_end, quadrant_name){
    lg20_recombinationcounts=makehistlist(lg, genotypelist, T)
    lg20_recombinationcounts$position=as.numeric(as.character(lg20_recombinationcounts$position))
    recomb_bin=filter(lg20_recombinationcounts, position>position_start & position<position_end)
    
    nounique_total=length(unique(quadrant_name[4:nrow(alexq1),2]))
    nounique_within_bin=length(unique(recomb_bin$individual))
    no_zero_events=(nounique_total-nounique_within_bin)
    
    listofevents=vector()
    for(i in unique(recomb_bin$individual)){
        temp_filter=filter(recomb_bin, individual==i)
        listofevents=c(listofevents, nrow(temp_filter))
    }
    for(i in 1:no_zero_events){
        listofevents=c(listofevents, 0)
    }
    return(listofevents)
}

noq1.20=countevents(20, q1genotypelist, 18, 50, alexq1)
noq2.20=countevents(20, q2genotypelist, 18, 50, alexq2)
noq3.20=countevents(20, q3genotypelist, 18, 50, alexq3)
noq4.20=countevents(20, q4genotypelist, 18, 50, alexq4)
?t.test()
t.test(noq3.20, noq2.20)

noq1.10=countevents(10, q1genotypelist, 57, 91, alexq1)
noq2.10=countevents(10, q2genotypelist, 57, 91, alexq2)
noq3.10=countevents(10, q3genotypelist, 57, 91, alexq3)
noq4.10=countevents(10, q4genotypelist, 57, 91, alexq4)
noq2.10
noq3.10

t.test(noq1.10, noq2.10)

#LGs that represent near-full chromosomes:
# LG1
# LG3 (although quite a large chunk is missing at the start of the chromosome)
# LG7 (again, quite a large chunk missing)
# LG9
# LG10 (represents one arm of the chromosome)
# LG15 (one arm)
# LG17
# LG20 (very big gap)
# LG23 (quite a large gap)



#     ____________________________________________________________________________
#     MAKE MULTILINE PLOT                                                                                                         ####


filelist=list.files("pairwiseLGcomparisons/allq_allq/")

lglist=c(1,3,7,9,10,15,17,20,23)
for(lgtomakehistlistfor in lglist){
    blastfile=read.csv(paste("pairwiseLGcomparisons/allq_allq/",filelist[grep(paste("_LG",lgtomakehistlistfor,"_", sep=""), filelist)], sep=""), header=T, stringsAsFactors = F)
    blastfile$blast_percentage[is.na(blastfile$blast_percentage)==T]=""
    blastfile$blast_percentage=as.numeric(blastfile$blast_percentage)
    blastfile$blast_percentage=round(blastfile$blast_percentage, 2)
    
    chromosomeforhere=blastfile[1, 17]
    
    
    q1hist=makehistlist(lgtomakehistlistfor, q1genotypelist)
    q2hist=makehistlist(lgtomakehistlistfor, q2genotypelist)
    q3hist=makehistlist(lgtomakehistlistfor, q3genotypelist)
    q4hist=makehistlist(lgtomakehistlistfor, q4genotypelist)
    
    q1histtable=as.data.frame(table(as.vector(as.numeric(q1hist))))
    q2histtable=as.data.frame(table(as.vector(as.numeric(q2hist))))
    q3histtable=as.data.frame(table(as.vector(as.numeric(q3hist))))
    q4histtable=as.data.frame(table(as.vector(as.numeric(q4hist))))
    
    q1histtable$temp=10
    q2histtable$temp=26
    q3histtable$temp=14
    q4histtable$temp=28
    
    allhisttable=rbind(q1histtable, q2histtable, q3histtable, q4histtable)
    colnames(allhisttable)=c("x","value","temp")
    allhisttable$temp=as.factor(allhisttable$temp)
    allhisttable$x=as.numeric(as.character(allhisttable$x))
    
    pdf(file=paste("multilineplots/recombination_analysis/new_custom/customprinttest_lg", lgtomakehistlistfor,".pdf", sep=""), width=15, height=10)
    graph=ggplot(data=allhisttable, aes(x=x, y=value, group=temp))
    print(graph + geom_line(aes(color=temp, linetype=temp), size=1.2) + 
        theme_bw() + ggtitle(paste("linkage group", lgtomakehistlistfor, "chromosome", chromosomeforhere)) + 
        scale_linetype_manual(values=c("twodash","dotted","solid","dashed")) + 
        scale_color_brewer(palette = "Spectral") +
        scale_x_continuous(breaks=seq(1,length(blastfile$blast_percentage),1), labels=blastfile$blast_percentage) +
        theme(axis.text.x=element_text(angle=90, hjust=1)) +
        xlab("Physical distance (%)") + 
        ylab ("# Recombination events"))
    dev.off()
}





#make 4x4 histogram

list_of_lg=as.numeric(unique(q1genotypelist$lg))[-36]
for(x in list_of_lg){
    q1hist=makehistlist(x, q1genotypelist)
    q2hist=makehistlist(x, q2genotypelist)
    q3hist=makehistlist(x, q3genotypelist)
    q4hist=makehistlist(x, q4genotypelist)
    
    q1histtable=as.data.frame(table(as.vector(as.numeric(q1hist))))
    q2histtable=as.data.frame(table(as.vector(as.numeric(q2hist))))
    q3histtable=as.data.frame(table(as.vector(as.numeric(q3hist))))
    q4histtable=as.data.frame(table(as.vector(as.numeric(q4hist))))
    
    q1histtable$quad=1
    q2histtable$quad=2
    q3histtable$quad=3
    q4histtable$quad=4
    
    allhisttable=rbind(q1histtable, q2histtable, q3histtable, q4histtable)
    allhisttable$quad=as.factor(allhisttable$quad)
    colnames(allhisttable)=c("x","value","quad")
    graph=ggplot(data=allhisttable, aes(x=x, y=value, group=quad))
    # windows()
    png(filename=paste("multilineplots/recombination_analysis/multiline_lg_",x,".png",sep=""), width=2000,height=2000)
    print(graph + geom_line(aes(color=quad, linetype=quad), size=1.2) + theme_bw())
    dev.off()
    
    
    listofbiggest=list()
    biggest1=max(as.numeric(q1hist))
    biggest2=max(as.numeric(q2hist))
    biggest3=max(as.numeric(q3hist))
    biggest4=max(as.numeric(q4hist))
    truebiggest=max(as.numeric(c(biggest1, biggest2, biggest3, biggest4)))

    
    #png(filename=paste("histogram_recombination_counts/",quadrant.to.analyse.name, "/", "hist_",quadrant.to.analyse.name,"_lg_",x,".png",sep=""))
    png(filename=paste("histogram_recombination_counts/all_quadrants/newhist_lg_",x,".png",sep=""), width=1000, height=1000)
    par(mfrow=c(2,2))
    hist(as.numeric(q1hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    hist(as.numeric(q2hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    hist(as.numeric(q3hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    hist(as.numeric(q4hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    dev.off()
}

?hist
?axis



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
