############IMPORT ALL PROBES FROM ALEX DATA, MAKE PAIRWISE COMPARISONS OF PARENTS##############

setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#allprobes=read.table("alexdataallprobes.txt", stringsAsFactors = F)
#allprobes=read.table("alexax_corrected.txt", stringsAsFactors = F) #lets try pairwise comparisons with renamed files
allprobes=read.table("alexaxnew(A8_incl).txt", stringsAsFactors = F) #lets try pairwise comparisons with renamed files
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

pairwiseapogee=matrix(nrow=6,ncol=6)
pairwiseapogee=data.frame(pairwiseapogee)
pairwiseparagon=matrix(nrow=7,ncol=7)
pairwiseparagon=data.frame(pairwiseparagon)


pairwiseparagon[2:7,1]=paragonall2[,1]
pairwiseparagon[1,2:7]=paragonall2[,1]
pairwiseapogee[2:6,1]=apogeeall2[,1]
pairwiseapogee[1,2:6]=apogeeall2[,1]
for(i in 1:6){
    for(p in 1:6){
        c=length(which(paragonall[i,]==paragonall[p,]))
        pairwiseparagon[i+1,p+1]=(c-1)/b*100
    }
}

for(i in 1:5){
    for(p in 1:5){
        d=length(which(apogeeall[i,]==apogeeall[p,]))
        pairwiseapogee[i+1,p+1]=(d-1)/b*100
    }
}


h/b*100

rowenaapogeeall=read.csv("rowenaapogeeall.csv",header=F)
rowenaapogeeall=rowenaapogeeall[-1,]

yo=length(which((apogeeall2[3,]==rowenaapogeeall[1,])))
yo/b*100

rowenaparagonall=read.csv("rowenaparagonall.csv",header=F)
rowenaparagonall=rowenaparagonall[-1,]

bla=length(which((paragonall2[3,]==rowenaparagonall[1,])))
bla/b*100

ncol(allprobes2)
ncol(rowenaparagonall)

# paragon ids in allprobes2: 31 242 244
# apogee ids in allprobes2: 286, 308, 328
length(which(allprobes2[328,]==rowenaparagonall[1,]))/ncol(allprobes2)
length(which(allprobes2[308,]==rowenaapogeeall[1,]))/ncol(allprobes2)
