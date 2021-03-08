library(readxl)
library(ggplot2)
library(reshape2)
sim.stats1.old = read_excel("E:/ac14037.backup/Google Drive/University/PhD/Seg dist simulation/list of simulations v3.xlsx")

sim.stats1 = read.csv("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/all.analysesv2.csv")
sim.stats1 = cbind(1:nrow(sim.stats1), sim.stats1)
sim.stats1 = convert.to.character.data.frame(sim.stats1)


sim.stats.selec = sim.stats1[c(3, 25:34), ]
sim.stats.selec = convert.to.character.data.frame(sim.stats.selec)
sim.stats.selec[1, 5:7] = c(T, 200, 0)
# sim.stats.selec[1, 3] = T

sim.stats.selec = as.data.frame(sim.stats.selec)
sim.stats.selec[, 7] = as.numeric(sim.stats.selec[, 7])
sim.stats.selec[2:nrow(sim.stats.selec), 7] = 1 / sim.stats.selec[2:nrow(sim.stats.selec), 7]

sim.stats.selec2 = sim.stats.selec[, c(7, 11:15)]
colnames(sim.stats.selec2) = c("Selection strength", "0.05", "0.01", "0.001", "0.05 FDR", "0.05 Bonferroni")

p.data1 = melt(sim.stats.selec2, "Selection strength")
colnames(p.data1) = c("s.strength", "var", "val")
p.data1$val = as.numeric(p.data1$val)
p.data1$val = p.data1$val / 10
plot1 = ggplot(p.data1, aes(x = s.strength, y = val, group = var, color = var, shape = var)) + scale_x_continuous(labels = c("0", "0.05", "0.1", "0.15", "0.2", "0.25"), breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) + geom_line(size = 0.7, alpha = 0.8) + theme_classic() + xlab("Selection strength") +
    ylab("Proportion of sim. w/ sig. distorted markers (%)") + scale_color_manual(labels = c("0.05", "0.01", "0.001", "0.05 FDR", "0.05 Bonf."), values = c("#050517", "#FF8C42", "#FF3C38", "#372549", "#A23E48")) + scale_shape_manual(labels = c("0.05", "0.01", "0.001", "0.05 FDR", "0.05 Bonf."), values = c(15, 16, 17, 18, 3)) + coord_cartesian(xlim = c(0, 0.25)) + guides(color = guide_legend(title = "Threshold"), shape = guide_legend(title = "Threshold")) + geom_point(size = 2) + ggtitle("")
plot1

pdf("rotation1scripts_v4/plots/simulation/pedsim/threshold.comparisonv2.pdf", 8, 8)
plot1
dev.off()


sim.stats2 = read.csv("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/all.analysesv2.csv", stringsAsFactors = F, header = T)
# sim.stats2 = read_excel("E:/ac14037.backup/Google Drive/University/PhD/Seg dist simulation/list of simulations v3.xlsx")
sim.stats2 = sim.stats2[1:44, ]
sim.stats2 = cbind(1:nrow(sim.stats2), sim.stats2)


colnames(sim.stats2)[7] = "selec.str"
sim.stats2 = as.data.frame(sim.stats2)
sim.stats2$selec.str = as.numeric(sim.stats2$selec.str)

sim.stats3 = sim.stats2

sim.stats3$selec.str = (1 / sim.stats3$selec.str)
sim.stats3$selec.str[1:4] = 0
# v(arrange(sim.stats3, Pop, selec.str))

colnames(sim.stats3)[14] = "fdr0.05"
sim.stats3$pop.size = as.factor(sim.stats3$pop.size)

plot2 = ggplot(sim.stats3, aes(x = selec.str, y = fdr0.05 / 10, group = pop.size, color = pop.size, shape = pop.size)) + geom_point(size = 2, alpha = 0.7) + 
    geom_line(size = 0.7) + geom_abline(intercept = 950, size = 1) + xlab("Selection strength") + ylab("Proportion of sim. w/ sig. distorted markers (%)") +
    guides(shape = guide_legend(title = "Pop. Size"), color = guide_legend(title = "Pop. Size")) + theme_classic() + ggtitle("")


pdf("rotation1scripts_v4/plots/simulation/pedsim/pop.size.vs.selec.strengthv2.pdf", 8, 8)
plot2
dev.off()

pdf("rotation1scripts_v4/plots/simulation/pedsim/proportion.sims.combv4.pdf", 12, 6)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()




library(extrafont)
loadfonts(dev = "win")
windowsFonts()

postscript("rotation1scripts_v4/plots/simulation/pedsim/proportion.sims.combv4.eps", horizontal = F, onefile = F, paper = "special")
plot1
dev.off()

grid.arrange(plot1, plot2, ncol = 2)
g = arrangeGrob(plot1, plot2, ncol = 2)
plot1
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig5.eps", g, width = 17.4, height = 9, units = "cm", device = cairo_pdf)

#### NEW FIGURES FOR REVISED MANUSCRIPT ####
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig5a.eps", plot1, width = 8.7, height = 9, units = "cm", device = cairo_pdf)
ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/Fig5b.eps", plot2, width = 8.7, height = 9, units = "cm", device = cairo_pdf)
