genome.sizes = data.frame(c("Wheat", "Barley", "Human", "Rice", "Arabidopsis"), c(17000, 5300, 3000, 430, 135))

colnames(genome.sizes) = c("Organism", "Size")
genome.sizes$Organism = factor(genome.sizes$Organism, levels = genome.sizes$Organism)
ggplot(genome.sizes, aes(x = Organism, y = Size)) + geom_bar(stat = "identity", fill = "#000000") +
    theme_classic() + ylab("Genome Size (Mb)") 

svg("E:/ac14037.backup/Google Drive/University/PhD/presentations/round table talk/figures/genome.size.bar.svg", 4, 4)
ggplot(genome.sizes, aes(x = Organism, y = Size)) + geom_bar(stat = "identity", fill = "#000000") +
    theme_classic() + ylab("Genome Size (Mb)") 
dev.off()