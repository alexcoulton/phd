setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
allprobes=read.csv("rowenaallprobes.csv", stringsAsFactors = F,header=F)
allprobes=t(allprobes)
colnames(allprobes)=allprobes[1,]

?grep
grep("Apogee",allprobes[,1])
apogeeall=allprobes[grep("Apogee",allprobes[,1]),]
paragonall=allprobes[grep("Paragon", allprobes[,1]),]

apogeeall=data.frame(apogeeall,stringsAsFactors = F)
paragonall=data.frame(paragonall,stringsAsFactors = F)

b=ncol(apogeeall)-1
c=length(which(apogeeall[2,]==apogeeall[3,]))
c/b*100

h=length(which(paragonall[4,]==paragonall[6,]))

h/b*100

#convert to character matrix
paragonall2=as.data.frame(as.matrix(paragonall), stringsAsFactors = F)
paragonall2[,1]

apogeeall2=as.data.frame(as.matrix(apogeeall), stringsAsFactors = F)
apogeeall2[,1]

pairwiseapogee=matrix(nrow=4,ncol=4)
pairwiseapogee=data.frame(pairwiseapogee)
pairwiseparagon=matrix(nrow=4,ncol=4)
pairwiseparagon=data.frame(pairwiseparagon)


pairwiseparagon[2:4,1]=paragonall2[,1]
pairwiseparagon[1,2:4]=paragonall2[,1]
pairwiseapogee[2:4,1]=apogeeall2[,1]
pairwiseapogee[1,2:4]=apogeeall2[,1]
for(i in 1:3){
    for(p in 1:3){
        c=length(which(paragonall[i,]==paragonall[p,]))
        pairwiseparagon[i+1,p+1]=(c-1)/b*100
    }
}

for(i in 1:3){
    for(p in 1:3){
        d=length(which(apogeeall[i,]==apogeeall[p,]))
        pairwiseapogee[i+1,p+1]=(d-1)/b*100
    }
}


h/b*100


write.csv(apogeeall2, "rowenaapogeeall.csv",row.names=F,col.names=F)
write.csv(paragonall2, "rowenaparagonall.csv", row.names=F,col.names=F)
