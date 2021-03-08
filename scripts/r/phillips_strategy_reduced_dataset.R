#generate reduced dataset al la Phillips et al., 2015
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts/")
library(dplyr)
alldata=read.csv("alexq1_man_cur_sk_map_genogeno_all.csv", header=T, stringsAsFactors = F)

lgtopick=1
total=alldata[,c(1,2)]
for(i in 1:35){
    lgtopick=i
    b=alldata[,c(1,2,which(alldata[1,]==lgtopick))]
    ncol(b)
    divider=ncol(b)/10
    round(divider)
    b_reduced=b[,c(seq(3,ncol(b), round(divider)))]
    total=cbind(total, b_reduced)
}
unique(as.numeric(unique(total[1,])))

#total is the reduced dataset with a limited number of markers
q1reduced=total[c(1,grep("Q1", total[,2])),2:ncol(total)]
q2reduced=total[c(1,grep("Q2", total[,2])),2:ncol(total)]
q3reduced=total[c(1,grep("Q3", total[,2])),2:ncol(total)]
q4reduced=total[c(1,grep("Q4", total[,2])),2:ncol(total)]


write.csv(q1reduced, "reduced_datasets/q1reduced.csv", row.names=F)
write.csv(q2reduced, "reduced_datasets/q2reduced.csv", row.names=F)
write.csv(q3reduced, "reduced_datasets/q3reduced.csv", row.names=F)
write.csv(q4reduced, "reduced_datasets/q4reduced.csv", row.names=F)
