#getting probe sequences from 35k array in suitable format for blast

setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts/")
kprobes=read.csv("C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/35k_probe_set.csv",header=T)
which(regexpr("\\[...\\]", kprobes$Sequence)==-1) ###check if any of the sequences don't follow the standard format

#conversion of factors to characters
kprobes2=kprobes
kprobes2$Sequence=as.character(kprobes2$Sequence)
class(kprobes2$Sequence[1])
kprobes2$SNPType=as.character(kprobes2$SNPType)
class(kprobes2$SNPType[1])

#replace [G\A] format in sequences with SNP ambiguity codes
for(i in 1:length(kprobes2$Sequence)){
kprobes2$Sequence[i]=gsub("\\[...\\]", kprobes2$SNPType[i], kprobes2$Sequence[i])
}

#write sequences to fasta file
for(i in 1:nrow(kprobes2)){
    write(c(paste(">", kprobes2$ï..35K.SNPId[i], sep=""), kprobes2$Sequence[i]), file="35kprobes.fa", append=T)
}


