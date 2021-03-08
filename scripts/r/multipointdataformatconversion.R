setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#setwd("~/Google Drive/R/scripts/rotation1scripts/")

quadrant="q4"

#rawdata=read.csv(paste("alex",quadrant,".csv",sep=""), header=T)
rawdata=read.csv("alexq4_mancur_ske.csv", header=T) #updated 07/12/17
# rawdata = c.x.af2 #from c.x.aprocessing2.R
rawdata=rawdata[-1,-c(1:2)]
rawdata=t(rawdata)
#write.csv(rawdata, "alexq2_multipointformat.csv")

b=paste(c(rawdata[1,]), sep="", collapse="")
b
?write

for(i in 1:nrow(rawdata)){
    #write(paste(rownames(rawdata)[i], paste(c(rawdata[i,]), sep="", collapse="")), file=paste("alex",quadrant,"_multipointformat.txt",sep=""), append=T)
    write(paste(rownames(rawdata)[i], paste(c(rawdata[i,]), sep="", collapse="")), file="rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxaf2.multpoint.txt", append=T)
}

rownames(rawdata)[1]

