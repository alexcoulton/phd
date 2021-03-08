#going to plot snp density to infer centromere position
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts_v2/")

allprobegenomeblast=read.table("original_data/allprobegenome.blast")
nrow(allprobegenomeblast)
allprobegenomeblast.filtered=filter(allprobegenomeblast, V11<1e-20)

#grab best hits (by smallest e value) for each probe #THIS TAKES A VERY LONG TIME --- read the file instead (see write.csv below)
# allprobe.fil.2 = newdf(paste("V", 1:12, sep = ""))
# for(i in unique(allprobegenomeblast.filtered$V1)){
#     tempdf = filter(allprobegenomeblast.filtered, V1 == i)
#     tempdf$V11 = as.numeric(tempdf$V11)
#     smallest.e.val = min(tempdf$V11)
#     tempdf = filter(tempdf, V11 == smallest.e.val)
#     if(nrow(tempdf) > 1){
#         tempdf = tempdf[1,]
#     }
#     allprobe.fil.2 = rbind(allprobe.fil.2, tempdf)
# }

#write.csv(allprobe.fil.2, "processed_data/allprobegenomeblast_unique_best_hits.csv", row.names = F)
allprobe.fil.2 = read.csv("processed_data/allprobegenomeblast_unique_best_hits.csv", header = T, stringsAsFactors = F)

allprobe.fil.2$V9 = as.numeric(allprobe.fil.2$V9)

chromosomecounts=read.table("original_data/chromosomecounts.txt")

listofchromosomes=unique(chromosomecounts$V2)
listofchromosomes=as.character(listofchromosomes)
listofchromosomes

hitsforchromo=filter(allprobe.fil.2, V2==listofchromosomes[17])

hitsforchromo$percentage_start = (hitsforchromo$V9 / chromosomecounts$V1[1])*100
hitsforchromo$mb = hitsforchromo$V9 / 10^6

hist(hitsforchromo$mb, main="Chromosome 6B", xlim=c(0, 750), breaks = 1000, xaxt = "n", xlab = "Physical Distance (Mb)")    
axis(side=1, at=seq(0, 750, 25), labels=seq(0, 750, 25), xlab = "Physical Distance (Mb)")

windows()
par(mfrow=c(1,2))

plot2 = ggplot(data = hitsforchromo, aes(hitsforchromo$mb)) + geom_histogram(binwidth = 1) + theme_classic() + 
    scale_x_continuous(breaks = c(seq(0, 800, 25)), expand = c(0.01, 0.01)) + xlab("Physical Position (Mb)") + ylab("Frequency")

?hist

count = 1
for(i in listofchromosomes){
    hitsforchromo=filter(allprobe.fil.2, V2 == i)
    hitsforchromo$percentage_start = (hitsforchromo$V9 / chromosomecounts$V1[count])*100
    count = count + 1
    pdf(file = paste("plots/snpdensityplots/chromosome_", i, ".pdf", sep = ""))
    #png(filename=paste("plots/snpdensityplots/chromosome_",i,".png",sep=""))
    hist(hitsforchromo$percentage_start, main=i, xlim=c(0, 100), breaks = 1000, xaxt = "n")    
    axis(side=1, at=seq(0, 100, 5), labels=seq(0, 100, 5))
    dev.off()
}

possible_centromeric_regions=as.data.frame(matrix(ncol=3,nrow=1))
colnames(possible_centromeric_regions)=c("chromosome","region","frequency")
index=1
for(i in listofchromosomes){
    hitsforchromo=filter(allprobegenomeblast.filtered, V2 == i)
    info=filter(as.data.frame(table(cut(hitsforchromo$V9, breaks=20, ordered_result = T))), Freq==min(Freq))
    possible_centromeric_regions$chromosome[index]=i
    possible_centromeric_regions$region[index]=as.character(info[1,1])
    possible_centromeric_regions$frequency[index]=info[1,2]
    possible_centromeric_regions=add_row(possible_centromeric_regions)
    index=index+1
}

moot=as.data.frame(strsplit(possible_centromeric_regions$region, ","))
moot=t(moot)
rownames(moot)=c()
moot[,1]=substr(moot[,1],2,9)
moot[,2]=substr(moot[,2],1,nchar(moot[,2])-1)
moot=moot[-nrow(moot),]

moot[,1]=as.numeric(moot[,1])
moot[,2]=as.numeric(moot[,2])


possible_centromeric_regions=possible_centromeric_regions[-nrow(possible_centromeric_regions),]

possible_centromeric_regions=cbind(possible_centromeric_regions, moot)
possible_centromeric_regions$fulllength=chromosomecounts$V1
possible_centromeric_regions$fulllength=as.numeric(as.character(possible_centromeric_regions$fulllength))
possible_centromeric_regions[,4]=as.numeric(as.character(possible_centromeric_regions[,4]))
possible_centromeric_regions[,5]=as.numeric(as.character(possible_centromeric_regions[,5]))
possible_centromeric_regions$percentage_start=(possible_centromeric_regions[,4]/possible_centromeric_regions$fulllength)*100
possible_centromeric_regions$percentage_end=(possible_centromeric_regions[,5]/possible_centromeric_regions$fulllength)*100


write.csv(possible_centromeric_regions, "possible_centromeric_regions.csv", row.names=F)
