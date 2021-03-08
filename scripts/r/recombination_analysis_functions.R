# setwd("~/Google Drive/R/scripts/rotation1scripts/") #mac
# setwd("C:/Users/ac14037/project.phd.main/") #windows
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)

list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)
list.of.one.arm.lgs.to.test = c(10, 21)

#setup dataframe - for counting the number of consecutive genotypes that are the same
makegenotypelist=function(quadrant.to.analyse){
    #args:
    #quadrant.to.analyse - a dataframe containing genotype information in rQTL format #see alexq1_man_cur_sk_map_genogeno.csv for format
    genotypelist=as.data.frame(matrix(nrow=1, ncol=5))
    genotypelist[is.na(genotypelist)]=""
    colnames(genotypelist)=c("genotype","number_consec","lg", "individual", "individual_id")
    
    geno_count=1
    for(p in 2:nrow(quadrant.to.analyse)){
        
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
    if(length(which(genotypelist$genotype == "-")) != 0){
        genolist2=genotypelist[-which(genotypelist$genotype=="-"),]
    } else {
        genolist2 = genotypelist
    }
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


makehistlist=function(linkagegroup, genolist2, return_recombinationcounts, show_nonevents, homozygous){ 
    #args:
    # linkagegroup - give the linkage group to parse as an integer (e.g. 1)
    # genolist2 - provide genotypelist made with makegenotypelist
    # return_recombinationcounts - return dataframe if true, list if false
    # show_nonevents - show individuals that do not have any recombination events
    # homozygous - boolean, is this population homozygous (e.g. F7, double haploid)? if so, don't double the number of rows for A -> B transitions in the dataframe
    
    if(missing(homozygous)) homozygous = F
    
    if(missing(show_nonevents)){
        show_nonevents = F
    }
    genolg1=filter(genolist2, lg==linkagegroup) #replace this to examine a different linkage group
    # browser()
    
    #count recombination events and note where they occur
    #first, setup a cumulative count column
    genolg1=add_column(genolg1, .after="number_consec", cumulative_count="")
    for(i in unique(as.numeric(genolg1$individual_id))){
        tempgeno=filter(genolg1, individual_id==i)
        tempgeno$cumulative_count=cumsum(tempgeno$number_consec)
        genolg1[genolg1$individual_id==i,][,3]=tempgeno$cumulative_count
        
    }
    
    recombination_counts=as.data.frame(matrix(ncol=6,nrow=1))
    recombination_counts[is.na(recombination_counts)]=""
    colnames(recombination_counts)=c("individual","recomb_event_id","position","event_type","num_recomb", "is.recomb.event")
    no_recomb_events = newdf(c("individual","recomb_event_id","position","event_type","num_recomb", "is.recomb.event"))
    counter=1
    counter2 = 1
    for(i in unique(as.numeric(genolg1$individual_id))){
        tempgeno=filter(genolg1, individual_id==i)
        for(g in 1:(nrow(tempgeno))){
            if(nrow(tempgeno) == 1){
                
                no_recomb_events$individual[counter2]=i
                no_recomb_events$recomb_event_id[counter2]=g
                no_recomb_events$is.recomb.event[counter2] = "N"
                no_recomb_events = add_row(no_recomb_events)
                
                counter2 = counter2 + 1
            } else if(g != nrow(tempgeno)){
                recombination_counts$individual[counter]=i
                recombination_counts$recomb_event_id[counter]=g
                recombination_counts$position[counter]=(as.numeric(tempgeno$cumulative_count[g])+1)
                recombination_counts$event_type[counter]=paste(tempgeno$genotype[g], "to", tempgeno$genotype[(g+1)])
                recombination_counts$is.recomb.event[counter] = "Y"
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
    
    #duplicate rows with double events
    if(homozygous != T){
        double.events = filter(recombination_counts, num_recomb == 2)
        single.events = filter(recombination_counts, num_recomb == 1)
        
        no_recomb_events[is.na(no_recomb_events)] = ""
        no_recomb_events = no_recomb_events[-nrow(no_recomb_events),]
        
        recombination_counts = rbind(no_recomb_events, single.events, double.events, double.events)
    } else {
        double.events = filter(recombination_counts, num_recomb == 2)
        single.events = filter(recombination_counts, num_recomb == 1)
        
        no_recomb_events[is.na(no_recomb_events)] = ""
        no_recomb_events = no_recomb_events[-nrow(no_recomb_events),]
        
        recombination_counts = rbind(no_recomb_events, single.events, double.events)
    }
    
    recombination_counts = recombination_counts[sort(as.numeric(recombination_counts$individual), index.return = T)[[2]],]
    rownames(recombination_counts) = 1:nrow(recombination_counts)
    
    
    
    b = unique(recombination_counts$individual)
    newsorted = newdf(colnames(recombination_counts))
    for(i in b){
        x = filter(recombination_counts, individual == i)
    x = x[sort(x$recomb_event_id, index.return = T)[[2]], ]
    newsorted = rbind(newsorted, x)
    }
    
    recombination_counts = newsorted
    recombination_counts = recombination_counts[-1,]
    rownames(recombination_counts) = 1:nrow(recombination_counts)
    
    if(show_nonevents == F){
        g = which(recombination_counts$is.recomb.event == "N")
        g1 = which(colnames(recombination_counts) == "is.recomb.event")
        if(length(g) > 0) recombination_counts = recombination_counts[-g,]
        if(length(g1) > 0) recombination_counts = recombination_counts[-g1, ]
        
    } else {}
    
    #produce weighted histogram based on these new values
    
    # browser()
    
    if(return_recombinationcounts == T){
        return(recombination_counts)
    } else {
        hist_list=vector()
        if(nrow(recombination_counts) > 0){
            for(i in 1:(nrow(recombination_counts)-1)){
                hist_list=c(hist_list, rep(recombination_counts$position[i], recombination_counts$num_recomb[i]))    
            }
        }
        hist_list = as.numeric(hist_list)
        return(hist_list)
    }
}


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


#takes a numeric vector, with quadrant and temperature specified, and converts it into a dataframe
conv.vec.to.df=function(vector, quad, temp){
    df1=as.data.frame(vector)
    df1$quad=quad
    df1$temp=temp
    return(df1)
}

#this function is perhaps too specific in its aims: it assumes a particular approach which I have now changed (not binning LGs)
setupanova=function(lgtotest, return_df, average_events){
    if(average_events == T){
        avpos1=getaverageposition(lgtotest, q1genotypelist)
        avpos2=getaverageposition(lgtotest, q2genotypelist)
        avpos3=getaverageposition(lgtotest, q3genotypelist)
        avpos4=getaverageposition(lgtotest, q4genotypelist)
    } else {
        avpos1=as.numeric(makehistlist(lgtotest, q1genotypelist, F))
        avpos2=as.numeric(makehistlist(lgtotest, q2genotypelist, F))
        avpos3=as.numeric(makehistlist(lgtotest, q3genotypelist, F))
        avpos4=as.numeric(makehistlist(lgtotest, q4genotypelist, F))
    }
        
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

#counts the number of recombination events per individual within a specified bin in the LG. Needs to be updated to take into account
#multiple recombination events, e.g. a transition from A to B = 2 events etc. 
countevents=function(lg, genotypelist, position_start, position_end, quadrant_name){
    lg20_recombinationcounts=makehistlist(lg, genotypelist, T)
    lg20_recombinationcounts$position=as.numeric(as.character(lg20_recombinationcounts$position))
    recomb_bin=filter(lg20_recombinationcounts, position>position_start & position<position_end)
    
    nounique_total=length(unique(quadrant_name[4:nrow(initial.data[[1]]),2]))
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


setupglmdf=function(genotypelist, dataset, populationsize, temp, list.of.intervals){
    #Creates a new dataframe containing the number of recombination events (count data), for each interval for a particular temperature
    #list.of.intervals is a list of intervals for which to count in
    testdf=newdf(nrow = 1, ncol=3, c("rf","interval","temp"))
    counter=1
    for(i in 1:length(list.of.intervals)){
        eventslist=countevents(LG.TO.ANALYSE, genotypelist, list.of.intervals[[i]][1], list.of.intervals[[i]][2], dataset)
        for(p in 1:length(eventslist)){
            testdf$rf[counter]=eventslist[p]
            testdf$interval[counter]=i
            counter=counter+1
            testdf=add_row(testdf)
        }
        testdf=add_row(testdf)
    }
    testdf$temp=temp
    testdf=testdf[-nrow(testdf),]
    return(testdf)
}


setup.glm.df.rf=function(genotypelist, dataset, populationsize, temp, list.of.intervals){
    #Creates a new dataframe containing the recombination frequencies for each interval for a particular temperature
    #list.of.intervals is a list of intervals for which to calculate recombination frequency
    testdf=newdf(nrow = 1, ncol=3, c("rf","interval","temp"))
    counter=1
    for(i in 1:length(list.of.intervals)){
        testdf$rf[counter]=(sum(countevents(15, genotypelist, list.of.intervals[[i]][1], list.of.intervals[[i]][2], dataset)))#/populationsize)
        testdf$interval[counter]=i
        counter=counter+1
        testdf=add_row(testdf)
    }
    testdf$temp=temp
    testdf=testdf[-nrow(testdf),]
    return(testdf)
}

#https://stackoverflow.com/questions/26754745/replacing-nas-in-r-numeric-vectors-with-values-calculated-from-neighbours
replacena <- function(l) {
    stopifnot(is.numeric(l))
    indx <- is.na(l)
    l[indx] <- vapply(which(indx), function(x) {
        counter = 1
        beforex = l[x - 1]
        afterx = l[x + 1]
        
        #search through neighbours until nearest non-NA value before x is found
        while(is.na(beforex) == T & counter < 20){
            beforex = l[x - counter]
            counter = counter + 1
        }
        
        counter = 1
        #search through neighbours until nearest non-NA value after x is found
        while(is.na(afterx) == T & counter < 20){
            afterx = l[x + counter]
            counter = counter + 1
        }
        
        if(is.na(beforex) == T | is.na(afterx) == T){
            return(NA)
        } else return(mean(c(beforex, afterx)))
    }, FUN.VALUE = double(1))
    
    l
}

return.nonmonotonic.indexes = function(x){
    #pass a numeric vector, returns a logical vector indicating which indices are non-monotonic
    monotonic.vector = vector()
    for(i in 1:(length(x)) - 1){
        monotonic.vector = c(monotonic.vector, (x[i] < x[i + 1]))
    }
    monotonic.vector = !monotonic.vector
    return(monotonic.vector)
}

#https://stackoverflow.com/questions/26754745/replacing-nas-in-r-numeric-vectors-with-values-calculated-from-neighbours
#THIS IS NOT FINISHED
# replace.non.monotonic.values <- function(l) {
#     stopifnot(is.numeric(l))
#     indx <- return.nonmonotonic.indexes(l)
#     l[indx] <- vapply(which(indx), function(x) {
#         counter = 1
#         beforex = l[x - 1]
#         afterx = l[x + 1]
#         
#         #search through neighbours until nearest non-NA value before x is found
#         while(is.na(beforex) == T){ #NEED TO CHANGE THIS FROM is.na TO is.monotonic (which I have not yet written)
#             beforex = l[x - counter]
#             counter = counter + 1
#         }
#         
#         counter = 1
#         #search through neighbours until nearest non-NA value after x is found
#         while(is.na(afterx) == T){ #NEED TO CHANGE THIS FROM is.na TO is.monotonic
#             afterx = l[x + counter]
#             counter = counter + 1
#         }
#         
#         mean(c(beforex, afterx))
#     }, FUN.VALUE = double(1))
#     
#     l
# }
# 
