#script for comparing Allen's genetic maps to the maps I produce with rQTL
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
library(qtl)


#     ____________________________________________________________________________
#     COMPARISON OF LINKAGE GROUPS                                                                                        ####


allenaxp=read.csv("Allen Characterization of Wheat breeders array suppl. docs/Allen genetic map AxP.csv")
allenaxp$marker=gsub('-','.',allenaxp$marker) #fix affycode field so matches with flipped_data

#myaxpq1=t(read.csv("allchrpost_quad1.csv", header=F))
myaxpq1=read.csv("alexq2_q2_0.31_postsplitgeno.csv", header=F)
myaxpq1=add_column(myaxpq1, g="",.before=1)
myaxpq1=t(myaxpq1)


myaxpq1[,1]
chrlist=c("1A","2A","3A","4A","5A","6A","7A","1B","2B","3B","4B","5B","6B","7B","1D","2D","3D","4D","5D","6D","7D")
mylgnum=1:35 #mylgnum does actually reflect LG# here

#do some comparisons on the datasets
nrow(allenaxp)
nrow(myaxpq1)
#my map has more markers than Sascha's...
length(which(myaxpq1[,1] %in% allenaxp[,1]))
#both maps only share 1987 markers in common (2/3s) ... 

myaxpq1=myaxpq1[which(myaxpq1[,1] %in% allenaxp[,1]),] #make a reduced dataset that only contains markers found in both maps (my map has more markers than Saschas, which will alter similarity scores)
nrow(myaxpq1)

go=data.frame("mylg","saschalg","size_mylg","size_saschalg","simcalc_clust (%)", "no. markers shared")
similaritylist=list()
for(i in chrlist){
    for(p in mylgnum){
        subt=allenaxp[allenaxp[,2]==i,]
        subt2=myaxpq1[myaxpq1[,2]==p,]
        mylgsize=nrow(subt2)
        comp=length(which(subt2[,1] %in% subt[,1])) #find markers in my LG that are also in Sascha's LG
        noncomp=subt2[which(!(subt2[,1] %in% subt[,1])),] #find the markers that arn't...
        
        if(comp>10){
            glob=which(as.vector(noncomp[,1]) %in% allenaxp[,1])
            m=allenaxp[match(as.vector(noncomp[,1]), allenaxp[,1]),2] #other LG in Allen data that markers match to
            #glob2=which(as.vector(noncomp[,1]) %in% allenaxp[allenaxp[,2]==g,][,1])
            if(length(glob)>4){print(glob); print(i); print(p); print(m)}
            simcalc=((comp/mylgsize)*100) #perhaps not the best measure of similarity, as doesn't take into account the differing sizes of LGs
                                                                        #instead, lets measure how many markers in my LG are also in Sacha's LG (i.e. the comp variable in raw form)
            
            l=c(p,i,nrow(subt2),nrow(subt),simcalc,comp)
            pp=t(data.frame(l,stringsAsFactors = F))
            colnames(pp)=colnames(go)
            go=rbind(go,pp)
        }
    }
}

unique(allenaxp[,2])

#     ____________________________________________________________________________
#     ORDINALITY TESTING                                                                                                            ####


#test - need to assess similarity of ordinality within clusters of similar markers between mine and sascha's LGs
#first isolate only the markers that are in both datasets
subt=allenaxp[allenaxp[,2]=="1A",]
subt2=myaxpq1[myaxpq1[,2]==8,] #NB.: NEED TO INSPECT go AND FIND MATCHING LG-CHROMOSOME PAIR

which(subt[,1] %in% subt2[,1])

only1=subt[which(subt[,1] %in% subt2[,1]),]
only2=subt2[which(subt2[,1] %in% subt[,1]),]

#make dataframe of the index at which marker from 1 dataset occurs in the other dataset
orderlist=list()
for(i in 1:nrow(only1)){
numsame=as.numeric(which(only1[,1][i] == only2[,1]))
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

#format allen genetic map for analysis with rqtl
allenaxp_rqtl=t(allenaxp)
ncol(allenaxp_rqtl)
uniq=unique(allenaxp_rqtl[2,])
for(i in 1:length(uniq)){
    allenaxp_rqtl[allenaxp_rqtl==uniq[i]]=i
}

allenaxp_rqtl[allenaxp_rqtl=="AA"]="A"
allenaxp_rqtl[allenaxp_rqtl=="BB"]="B"
allenaxp_rqtl[allenaxp_rqtl=="AB"]="H"
rownames(allenaxp_rqtl)[1]=""
rownames(allenaxp_rqtl)[2]=""
rownames(allenaxp_rqtl)[3]=""
write.table(allenaxp_rqtl, "allenaxp_rqtlformat.csv", row.names=T, col.names=F, sep=",")

allendata = read.cross("csv", "./", "allenaxp_rqtlformat.csv", estimate.map=F) #NOTE THE THIRD ARGUMENT — CHANGE BETWEEN flipped_data2.csv AND flipped_data2_pruned.csv ACCORDINGLY

#plot RF graph for allendata... seems to be ordering much better than rQTL.
plotRF(allendata, alternate.chrid = T)
