filetoload="itparam5.csv"
param=read.csv(filetoload,header=F, stringsAsFactors = F)
#rfstart
param[1,2]=0.2
#rfend
param[2,2]=0.3
#lodstart
param[3,2]=20
#lodend
param[4,2]=20
write.table(param, "itparam10.csv", col.names = F, row.names = F, sep=",")
