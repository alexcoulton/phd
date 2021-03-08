setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
icimap=read.csv("flipped_data2.csv", header=T)
??transpose
icimap=t(icimap)
icimap=icimap[,-1]
icimap=icimap[-1,]
icimap=icimap[-1,]
icimap[icimap=='-']='X'

icimap[icimap=='A']=2
icimap[icimap=='B']=0
icimap[icimap=='H']=1
icimap[icimap=='X']=-1

write.table(icimap, file="icimap.csv", col.names=FALSE, sep=",")

ncol(icimap)
print(icimap[icimap !='X' & icimap!='A' & icimap!='H'])
