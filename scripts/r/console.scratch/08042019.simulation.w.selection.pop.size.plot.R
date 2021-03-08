setwd("E:/phd.project.main/")
load("rotation1scripts_v4/saved.objects/recomb.sim/all.pop.seg.ratios")

library(ggplot2)
library(gridExtra)

v(all.pop.seg.ratios[[1]])
v(all.pop.seg.ratios[[2]])
all.pop.seg.ratios


all.pop.seg.ratios2 = sapply(all.pop.seg.ratios, function(x){
    x = as.data.frame(x)
    unlist(lapply(x, function(y){
        mean(y)
    }))
})

g = do.call(cbind, all.pop.seg.ratios2)

all.pop.seg.ratios.sd = sapply(all.pop.seg.ratios, function(x){
    x = as.data.frame(x)
    unlist(lapply(x, function(y){
        sd(y)
    }))
})

g1 = do.call(cbind, all.pop.seg.ratios.sd)

g2 = cbind(g, g1)

g2 = as.data.frame(g2)


colnames(g2) = c("96mean", "300mean", "1000mean", "10000mean", "96sd", "300sd", "1000sd", "10000sd")

y2 = c(0.35, 0.65)

p1 = plotter2(g2, g2$`96mean`, g2$`96sd`, title = "Population Size: 96", ylimits = y2, hline = T)
p2 = plotter2(g2, g2$`300mean`, g2$`300sd`, title = "Population Size: 300", ylimits = y2, hline = T)
p3 = plotter2(g2, g2$`1000mean`, g2$`1000sd`, title = "Population Size: 1000", ylimits = y2, hline = T)
p4 = plotter2(g2, g2$`10000mean`, g2$`10000sd`, title = "Population Size: 10000", ylimits = y2, hline = T)




grid.arrange(p1, p2, p3, p4)

par(mfrow = c(2, 2))

ylimits = c(0.4, 0.6)
plot(all.pop.seg.ratios2[[1]], ylim = ylimits)
plot(all.pop.seg.ratios2[[2]], ylim = ylimits)
plot(all.pop.seg.ratios2[[3]], ylim = ylimits)
plot(all.pop.seg.ratios2[[4]], ylim = ylimits)




