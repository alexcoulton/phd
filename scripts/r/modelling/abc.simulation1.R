library(dplyr)
setwd("E:/phd.project.main/")
sim1_cleaned = readLines("rotation1scripts_v4/original_data/abc.simulation/sim1_cleaned.txt")

sim1_cleaned = read_delim("rotation1scripts_v4/original_data/abc.simulation/sim1_cleaned.txt", " ")

sim1.2 = arrange(sim1_cleaned, simstat1)
sim1.3 = sim1.2[1:150, ]

sim5b = read_delim("rotation1scripts_v4/original_data/abc.simulation/sim5b.cleaned.txt", " ")

#arrange by summary statistic column and pick the top 10% of hits
sim5b2 = arrange(sim5b, simstat1)
sim5b3 = sim5b2[1:500, ]

boxplot(sim5b3$simstat1)
boxplot(sim5b$simstat1)

library(abc)

abc(sim5b3$simstat1)


par(mfrow=c(2, 2))
plot(table(sim1_cleaned$selection.strength), xlim = c(0, 50))
plot(table(sim1_cleaned$selection.position), xlim = c(0, 250))

plot(table(sim1.3$selection.strength), xlim = c(0, 50))
plot(table(sim1.3$selection.position), xlim = c(0, 250))

quantile(sim5b$selection.strength)
quantile(sim5b3$selection.strength)

library(readr)
library(dplyr)
library(tibble)

dist.pos = lapply(sort(unique(sim5b$selection.position)), function(x){
    g = filter(sim5b, selection.position == x)
    mean(g$simstat1)
})

dist.pos = unlist(dist.pos)
plot(1 - dist.pos)

dist.strength = lapply(4:28, function(x){
    g = filter(sim5b, selection.strength == x)
    g = sort(as.numeric(g$simstat1))
    mean(g)
})

dist.strength = unlist(dist.strength)
plot(1 - dist.strength)


median(sim5b$simstat1[which(sim5b$selection.strength == 33)])
