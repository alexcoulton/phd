# Load of libraries

library(reshape2)
library(ggplot2)

# load("rotation1scripts_v4/saved.objects/heatmap.sim.df2")
# heatmap.sim.df2 = heatmap.sim.df[which(heatmap.sim.df$selec.strengths < 20), ]
# heatmap.sim.df2 = heatmap.sim.df2[which(heatmap.sim.df2$pop.sizes < 1001), ]

load("rotation1scripts_v4/saved.objects/pedsim.heatmap.data")

pedsim.heatmap.data = pedsim.heatmap.data[which(pedsim.heatmap.data$s.str < 21), ]

# Elaboration of heatmap (white - steelblue)


plotheat = ggplot(pedsim.heatmap.data, aes(pop, s.str, fill = dev )) +
    geom_tile() +
    scale_fill_gradientn(colors = c("#1a6aed", "#f2ff00", "#ff0000")) +
    ylab("Selection strength (1/y)") +
    xlab("Population size") +
    theme(legend.title = element_text(size = 10),
                legend.text = element_text(size = 12),
                plot.title = element_text(size=16),
                axis.title=element_text(size=14,face="bold"),
                axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = "Distortion") + scale_y_reverse(breaks = seq(20, 2, -2)) + scale_x_continuous(breaks = seq(100, 1000, 100)) + theme_classic() + 
    theme(axis.text.x = element_text(angle = 45), axis.title.x = element_text(vjust = 5))

pdf("rotation1scripts_v4/plots/simulation/pedsim/heatmap.pdf", 10, 10)
plotheat
dev.off()




ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/plotheat.", g, width = 17.4, height = 9, units = "cm", device = cairo_pdf)

tiff("rotation1scripts_v4/plots/simulation/pedsim/Fig4.tiff", units = "cm", width = 8.4, height = 7, res = 600)
plotheat
dev.off()