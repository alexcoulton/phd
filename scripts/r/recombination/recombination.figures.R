#figure 1
![**`r fig("clusterplot")`** Cluster plot showing a marker that has potentially been assigned an erroneous genotype. The marker is located outside of the posterior cluster of heterozygotes and also has a different genotype to both flanking markers in the genetic mapping data. This marker can therefore be reassigned as a NoCall so as not to incorrectly inflate the recombination frequency and distribution estimates from this genotyping data.    \label{figurelabel}](customfigs/clusterplot.correction.png)

#figure 2
tiff("rotation1scripts_v4/plots/recombination.paper/figure3v2.tiff", units = "cm", width = 10, height = 10, res = 450)
histogram.average.distances.of.recomb.axp
dev.off()

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure2.eps", histogram.average.distances.of.recomb.axp, width = 10, height = 10, units = "cm", device = cairo_pdf)


#figure 3
tiff("rotation1scripts_v4/plots/recombination.paper/figure3v3.tiff", units = "cm", width = 18, height = 18, res = 300, pointsize = 7.5)
# do.call(grid.arrange, plots2)
grid.arrange(plots2[[10]], plots2[[5]], ncol = 2)
dev.off()
gen.plot1 = grid.arrange(plots2[[10]], plots2[[5]], ncol = 2)
ggsave("rotation1scripts_v4/plots/recombination.paper/figure3v3.eps", )

ggsave("rotation1scripts_v4/plots/recombination.paper/figure3v3.eps", gen.plot1, width = 18, height = 18, units = "cm", device = cairo_pdf)


ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure3.eps", grid.arrange(plots2[[10]], plots2[[5]], ncol = 2), width = 19.05, height = 10, units = "cm", device = cairo_pdf)

fig3.plot.eps1 = plots2[[10]]
fig3.plot.eps2 = plots2[[5]]

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure3a.eps", fig3.plot.eps1, width = 15, height = 15, units = "cm", device = cairo_pdf)
ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure3b.eps", fig3.plot.eps2, width = 18, height = 18, units = "cm", device = cairo_pdf)


#figure 4
jpeg("rotation1scripts_v4/plots/recombination.paper/figure4v2.jpg", units = "cm", width = 35, height = 25, res = 450, pointsize = 2)
# plot_grid(axp.genedistplots[[4]][[5]], 
#                     axp.genedistplots[[3]][[5]], 
#                     axp.genedistplots[[2]][[5]], 
#                     axp.genedistplots[[1]][[5]], 
#                     genedisthist3b, ncol = 1,
#                     align = "v", axis = "l")
plot_grid(axp.genedistplots[[4]][[num1]],
                    axp.genedistplots[[3]][[num1]],
                    axp.genedistplots[[2]][[num1]],
                    axp.genedistplots[[1]][[num1]],
                    genedisthist3b,
                    ncol = 1,
                    rel_heights = c(0.4, 0.2, 0.2, 0.2),
                    align = "v", axis = "l")
dev.off()

num1 = 2
jpeg("rotation1scripts_v4/plots/recombination.paper/figure4v3.w.labels.jpg", units = "cm", width = 35, height = 25, res = 450, pointsize = 2)
plot_grid(axp.genedistplotsv2[[4]][[num1]],
                    axp.genedistplotsv2[[3]][[num1]],
                    axp.genedistplotsv2[[2]][[num1]],
                    axp.genedistplotsv2[[1]][[num1]],
                    genedisthist2a, 
                    align = "v", axis = "l", ncol = 1,
                    rel_heights = c(0.325, 0.2, 0.2, 0.2))
dev.off()

num1 = 2
tiff("rotation1scripts_v4/plots/recombination.paper/figure4v3.tiff", units = "cm", width = 18, height = 18, res = 450, pointsize = 2)
grid.arrange(arrangeGrob(recomb.plot1, left = y.grob))
dev.off()

recomb.plot1.save = grid.arrange(arrangeGrob(recomb.plot1, left = y.grob))

ggsave("rotation1scripts_v4/plots/recombination.paper/figure5new.eps", recomb.plot1, width = 17.4, height = 8, units = "cm", device = cairo_pdf)

save_plot("rotation1scripts_v4/plots/recombination.paper/figure5new.eps", recomb.plot1, ncol = 1, nrow = 5, base_height = 1.5, base_width = 8)


png("E:/ac14037.backup/Google Drive/University/PhD/presentations/round table talk/figures/tempshift.png", 4500, 2500,
        res = 400)
plot_grid(axp.genedistplots[[4]][[5]], 
                    axp.genedistplots[[3]][[5]], 
                    axp.genedistplots[[2]][[5]], 
                    axp.genedistplots[[1]][[5]], 
                    genedisthist3b, ncol = 1,
                    align = "v", axis = "l")
dev.off()


ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure4.eps", 
             plot_grid(axp.genedistplots[[4]][[5]], 
                                 axp.genedistplots[[3]][[5]], 
                                 axp.genedistplots[[2]][[5]], 
                                 axp.genedistplots[[1]][[5]], 
                                 genedisthist3b, ncol = 1,
                                 align = "v", axis = "l"),
             width = 16, height = 16, units = "cm", device = cairo_pdf)


#figure 5
tiff("rotation1scripts_v4/plots/recombination.paper/figure6.tiff", units = "cm", width = 12, height = 10, res = 450)
axp.chromo.maplengths.plot
dev.off()

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure5.eps", axp.chromo.maplengths.plot, width = 12, height = 10, units = "cm", device = cairo_pdf)


#figure 6
tiff("rotation1scripts_v4/plots/recombination.paper/figure7.tiff", units = "cm", width = 10, height = 10, res = 450)
axp.recomb.freq.plot1
dev.off()

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure6.eps", axp.recomb.freq.plot1, width = 10, height = 10, units = "cm", device = cairo_pdf)


#figure 7
tiff("rotation1scripts_v4/plots/recombination.paper/figure7.tiff", units = "cm", width = 10, height = 10, res = 450)
grid.arrange(hs.all.hist1, hs.all.hist2, hs.all.hist3, hs.all.hist4)
dev.off()

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure7.eps", grid.arrange(hs.all.hist1, hs.all.hist2, hs.all.hist3, hs.all.hist4), width = 12, height = 12, units = "cm", device = cairo_pdf)

#figure 8
tiff("rotation1scripts_v4/plots/recombination.paper/figure8.tiff", units = "cm", width = 10, height = 10, res = 450)
histogram.sim.means1
dev.off()

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure8.eps", histogram.sim.means1, width = 10, height = 10, units = "cm", device = cairo_pdf)


#figure 9 qtl - not using anymore??
tiff("rotation1scripts_v4/plots/recombination.paper/figure6.tiff", units = "cm", width = 19.05, height = 19.05, res = 350)
do.call(grid.arrange, qtl.plots)
dev.off()

ggsave(file = "rotation1scripts_v4/plots/recombination.paper/figure9.eps", do.call(grid.arrange, qtl.plots), width = 19.05, height = 19.05, units = "cm", device = cairo_pdf)


#figure 10
tiff("rotation1scripts_v4/plots/recombination.paper/figure10.tiff", units = "cm", width = 14, height = 12, res = 350)
par(mfrow = c(1, 2))
boxplot(alab.tree2, random.tree2, names = c("MGS", "RGS"), ylab = "Tree length")
boxplot(alab.tree.nobar2, random.tree.nobar2, names = c("MGS", "RGS"), ylab = "Tree length")
dev.off()




#figure 11
tiff("rotation1scripts_v4/plots/recombination.paper/figure11.tiff", units = "cm", width = 15, height = 15, res = 300)
par(mfrow = c(1, 2))
plot(alab.tree)
add.scale.bar(x = 0.12, y = 2.5, length = 0.05)
title("Meiosis Gene Set")
plot(random.tree)
add.scale.bar(x = 0.12, y = 2.5, length = 0.05)
title("Random Gene Set")
dev.off()

#figure 12a
tiff("rotation1scripts_v4/plots/recombination.paper/figure12a.tiff", units = "cm", width = 10, height = 10, res = 300)
plot(unroot(alab.tree.nobar), type = "unrooted", cex = 0.8)
title("(a) Meiosis Gene Set")
add.scale.bar(length = 0.0005)
dev.off()

#figure 12b
tiff("rotation1scripts_v4/plots/recombination.paper/figure12b.tiff", units = "cm", width = 13, height = 10, res = 300)
plot(unroot(random.tree.nobar), type = "unrooted", cex = 0.8)
title("(b) Random Gene Set")
add.scale.bar(length = 0.0005)
dev.off()

