setwd("E:/phd.project.main/")

yield = read.csv("rotation1scripts_v4/original_data/yield.data/global.csv")

yield2 = split(yield, yield$Item)

library(dplyr)

yield = bind_rows(yield2[[3]], yield2[[4]], yield2[[2]], yield2[[1]])
yield$Item = factor(yield$Item, levels = unique(yield$Item))
colnames(yield)[[8]] = "Crop"

library(ggplot2)
library(viridis)
plot1 = ggplot(yield, aes(x = Year, y = (Value / 1000000), color = Crop)) + geom_line(size = 1.5) + 
    scale_color_viridis(discrete = T, option = "E") + theme_bw() + ggtitle("Global crop production") +
    ylab("Total produced (Millions of tonnes)")

svg("rotation1scripts_v4/plots/yield/global.crop.production.svg", 6, 4)
plot1
dev.off()

png("rotation1scripts_v4/plots/yield/global.crop.production.png", 1500, 1200,
        res = 300)
plot1
dev.off()


genome.sizes = data.frame(c("Wheat", "Barley", "Human", "Rice", "Arabidopsis"), c(17000, 5300, 3200, 430, 135))
colnames(genome.sizes) = c("Crop", "Size")
genome.sizes$Crop = factor(genome.sizes$Crop, levels = genome.sizes$Crop)
plot2 = ggplot(genome.sizes, aes(x = Crop, y = Size)) + geom_bar(fill = "#000000", stat = "identity") +
    theme_classic() + ylab("Genome Size (Mb)")

png("rotation1scripts_v4/plots/yield/genome.size.comp.png", 1500, 1200,
        res = 300)
plot2
dev.off()



