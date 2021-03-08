#parsing tabular blast search of 35k probes against 
#wheat genome

setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts/")
blastg=read.table("allprobegenome.blast", stringsAsFactors = F)
blastg[,1]=gsub('-','\\.', blastg[,1])
blastg_orig = blastg

pro="AX.94524011"
blastg[blastg==pro,]

class(blastg$V11[1])
blastg$V11[1]

blastg_filtered=blastg[which(blastg$V11 < 1e-5),]
nrow(blastg_orig)
nrow(blastg_filtered)


#add probes with only one blast hit to a new list
listofunique=list()
for(i in 1:nrow(blastg_filtered)){
    if(blastg_filtered$V1[i] %in% blastg_filtered$V1[-i] == F){
        listofunique=c(listofunique,blastg_filtered$V1[i])
    }
}

#have a look at unique markers
length(listofunique)
dataofunique=t(data.frame(listofunique))
head(listofunique)

listofunique
rowenaprobes=read.csv("flipped_data2.csv", header=T)
rowenaprobes2=colnames(rowenaprobes)
head(rowenaprobes2)
rowenaprobes2=rowenaprobes2[-c(1,2)]
rowenaprobes2=gsub('\\.','-',rowenaprobes2) #fix affycode field so matches with flipped_data
?gsub


length(which(rowenaprobes2 %in% listofunique))
