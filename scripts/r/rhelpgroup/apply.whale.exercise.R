source("E:/phd.project.main/rotation1scripts_v4/scripts/r/functions.R")
library(ggplot2)
library(reshape2)

dat1 <- apply(fin.dd[, 3:4], 1, function(x){
    if(x[1] == 2) x[2] <- x[2] * 100
    if(x[1] == 3) x[2] <- x[2] * 30.48
    x
})
fin.dd[, 3:4] <- t(dat1)


plotdata1 <- as.data.frame(t(tapply(fin.dd$Length, list(fin.dd$Species, fin.dd$Year), mean)))
plotdata1$year <- rownames(plotdata1)




plotdata1 <- melt(plotdata1)
ggplot(plotdata1, aes(x = year, y = value, group = variable, color = variable)) + geom_line(size = 1)


