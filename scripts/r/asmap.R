#install.packages("ASMap")
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
library(ASMap)
library(qtl)



axiomdata = read.cross("csv", "./", "alexq1_man_cur_sk_map_genogeno.csv", genotypes=c("A", "H", "B"),estimate.map=F)
#axiomdata = read.cross("csv", "./", "alexq2_man_cur_sk_map_genogeno.csv", genotypes=c("A", "H", "B"),estimate.map=F)

summary(axiomdata)

axiomdata=convert2bcsft(axiomdata,BC.gen=0,F.gen=2)

axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1, anchor=T)



q1_asmap_map=pull.map(axiomdata_map, as.table=T)
#q2_asmap_map=pull.map(axiomdata_map, as.table=T)

#heatMap(axiomdata_map)

#do some plotting
#chr_num=20
#plot(q1_asmap_map[q1_asmap_map$chr==chr_num,2], q2_asmap_map[q2_asmap_map$chr==chr_num,2], xlim=c(0,100), ylim=c(0,100)); abline(0,1)


r=rownames(q1_asmap_map)
q1_asmap_map=add_column(q1_asmap_map, ske="")    
q1_asmap_map$marker=r
q1_asmap_map[,1]=as.character(q1_asmap_map[,1])

for(i in 1:nrow(q1_asmap_map)){
    q1_asmap_map[i,1]=paste("LG",q1_asmap_map[i,1],sep="")    
}

colnames(q1_asmap_map)=c("lg","pos","ske","marker")


write.csv(q1_asmap_map, "q1_asmap_map.csv", row.names=F)


#     ____________________________________________________________________________
#     REPEAT OF CODE FOR Q2                                                                                                     ####

axiomdata = read.cross("csv", "./", "alexq2_man_cur_sk_map_genogeno.csv", genotypes=c("A", "H", "B"),estimate.map=F)

summary(axiomdata)

axiomdata=convert2bcsft(axiomdata,BC.gen=0,F.gen=2)

axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1, anchor=T)
?mstmap

q2_asmap_map=pull.map(axiomdata_map, as.table=T)

#heatMap(axiomdata_map)

#do some plotting
#chr_num=20
#plot(q1_asmap_map[q1_asmap_map$chr==chr_num,2], q2_asmap_map[q2_asmap_map$chr==chr_num,2], xlim=c(0,100), ylim=c(0,100)); abline(0,1)


r=rownames(q2_asmap_map)
q2_asmap_map=add_column(q2_asmap_map, ske="")    
q2_asmap_map$marker=r
q2_asmap_map[,1]=as.character(q2_asmap_map[,1])

for(i in 1:nrow(q2_asmap_map)){
    q2_asmap_map[i,1]=paste("LG",q2_asmap_map[i,1],sep="")    
}

colnames(q2_asmap_map)=c("lg","pos","ske","marker")


write.csv(q2_asmap_map, "q2_asmap_map.csv", row.names=F)


#     ____________________________________________________________________________
#     QUADRANT 3                                                                                                                            ####

axiomdata = read.cross("csv", "./", "alexq3_man_cur_sk_map_genogeno.csv", genotypes=c("A", "H", "B"),estimate.map=F)

summary(axiomdata)

axiomdata=convert2bcsft(axiomdata,BC.gen=0,F.gen=2)

axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1)
?mstmap

q3_asmap_map=pull.map(axiomdata_map, as.table=T)

#heatMap(axiomdata_map)

#do some plotting
#chr_num=20
#plot(q1_asmap_map[q1_asmap_map$chr==chr_num,2], q2_asmap_map[q2_asmap_map$chr==chr_num,2], xlim=c(0,100), ylim=c(0,100)); abline(0,1)


r=rownames(q3_asmap_map)
q3_asmap_map=add_column(q3_asmap_map, ske="")    
q3_asmap_map$marker=r
q3_asmap_map[,1]=as.character(q3_asmap_map[,1])

for(i in 1:nrow(q3_asmap_map)){
    q3_asmap_map[i,1]=paste("LG",q3_asmap_map[i,1],sep="")    
}

colnames(q3_asmap_map)=c("lg","pos","ske","marker")


write.csv(q3_asmap_map, "q3_asmap_map.csv", row.names=F)


#     ____________________________________________________________________________
#     QUADRANT 4                                                                                                                            ####

axiomdata = read.cross("csv", "./", "alexq4_man_cur_sk_map_genogeno.csv", genotypes=c("A", "H", "B"),estimate.map=F)

summary(axiomdata)

axiomdata=convert2bcsft(axiomdata,BC.gen=0,F.gen=2)

axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1)
?mstmap
?mstmap


q4_asmap_map=pull.map(axiomdata_map, as.table=T)

#heatMap(axiomdata_map)

#do some plotting
#chr_num=20
#plot(q1_asmap_map[q1_asmap_map$chr==chr_num,2], q2_asmap_map[q2_asmap_map$chr==chr_num,2], xlim=c(0,100), ylim=c(0,100)); abline(0,1)


r=rownames(q4_asmap_map)
q4_asmap_map=add_column(q4_asmap_map, ske="")    
q4_asmap_map$marker=r
q4_asmap_map[,1]=as.character(q4_asmap_map[,1])

for(i in 1:nrow(q4_asmap_map)){
    q4_asmap_map[i,1]=paste("LG",q4_asmap_map[i,1],sep="")    
}

colnames(q4_asmap_map)=c("lg","pos","ske","marker")


write.csv(q4_asmap_map, "q4_asmap_map.csv", row.names=F)

#     ____________________________________________________________________________
#     ALL QUADRANTS                                                                                                                     ####

axiomdata = read.cross("csv", "./", "alexq1_man_cur_sk_map_genogeno_all.csv", genotypes=c("A", "H", "B"),estimate.map=F)

summary(axiomdata)

axiomdata=convert2bcsft(axiomdata,BC.gen=0,F.gen=2)

axiomdata_map=mstmap(axiomdata, trace=T, dist.fun="kosambi", id="probeset_id", p.value=1)
?mstmap
?mstmap


allq_asmap_map=pull.map(axiomdata_map, as.table=T)

#heatMap(axiomdata_map)

#do some plotting
#chr_num=20
#plot(q1_asmap_map[q1_asmap_map$chr==chr_num,2], q2_asmap_map[q2_asmap_map$chr==chr_num,2], xlim=c(0,100), ylim=c(0,100)); abline(0,1)


r=rownames(allq_asmap_map)
allq_asmap_map=add_column(allq_asmap_map, ske="")    
allq_asmap_map$marker=r
allq_asmap_map[,1]=as.character(allq_asmap_map[,1])

for(i in 1:nrow(allq_asmap_map)){
    allq_asmap_map[i,1]=paste("LG",allq_asmap_map[i,1],sep="")    
}

colnames(allq_asmap_map)=c("lg","pos","ske","marker")


write.csv(allq_asmap_map, "allq_asmap_map.csv", row.names=F)

