library(readr)
library(ASMap)
new1 = read_csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.noaxc.csv")

convert.rqtl.to.asmap.format = function(rqtl.df){
    #Convert a dataframe of genotyping data in the format of alexq1_man_cur_sk_map_genogeno.csv for use with mstmap.data.frame (ASMap function)
    rqtl.df = rqtl.df[-1, ]
    rqtl.df = rqtl.df[, -1]
    rqtl.df = as.data.frame(t(rqtl.df))
    rqtl.df = convert.to.character.data.frame(rqtl.df)
    colnames(rqtl.df) = rqtl.df[1, ]
    rqtl.df = rqtl.df[-1, ]
    rqtl.df[rqtl.df == "H"] = "X"
    as.data.frame(rqtl.df)
} 

g = convert.rqtl.to.asmap.format(new1)


g2 = mstmap.data.frame(g, "RIL2", p.value = 1e-25)

g2 = mstmap.data.frame(g1, "RIL2", p.value = 1e-25)

g1 = read.csv("rotation1scripts_v4/temp/asmap.not.work.csv", stringsAsFactors = F, header = T)

rownames(g1) = g1[, 1]
g1 = g1[, -1]

write.csv(g, "rotation1scripts_v4/temp/asmap.not.work.csv")
write.csv(gts, "rotation1scripts_v4/temp/asmap.work.csv")

row.names(g)
row.names(gts)


axiom.data = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.noaxc.csv", genotypes = c("A", "H", "B"), 
                                                estimate.map = F)

#Specify the Filial generation of the cross with the F.gen parameter
axiomdata = convert2bcsft(axiom.data, BC.gen = 0, F.gen = 2, estimate.map = F)

axiomdata_map = mstmap(axiomdata, trace = T, dist.fun = "kosambi", id = "probeset_id", p.value = 1e-25, anchor = F)

axiomdata_map

map1 = pull.map(axiomdata_map, as.table = T)
max(map1[which(map1$chr == "1.8"), ]$pos)


heatMap(axiomdata_map, "1.8") #linkage group 1.8 is chromosome 1A for c.x.a.f2.noaxc.csv
nrow(q2[[20]])

#### HEATMAP FOR SIMULATED DATA ####

load("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb1a.pop96")

#### HEATMAP FOR SIM. DATA #### 
gen.asmap.for.sim = function(sim.dat1, p.value1){
    if(missing(p.value1)) p.value1 = 4
    sim1.rqtl2 = sim.dat1
    sim1.rqtl2[sim1.rqtl2 == "H"] = "X"
    sim1.rqtl2 = as.data.frame(t(sim1.rqtl2))
    sim1.rqtl2 = convert.to.character.data.frame(sim1.rqtl2)
    sim1.rqtl3 = mstmap.data.frame(sim1.rqtl2, "RIL2", "kosambi", p.value = p.value1)
    sim1.rqtl3
}

all.asmap.sim = lapply(list.of.first.gen.recomb1a.pop96, gen.asmap.for.sim)
all.asmap.sim.lengths = lapply(all.asmap.sim, function(x){
    g = pull.map(x, as.table = T)
    max(g$pos)
})

mean(unlist(all.asmap.sim.lengths))
sd(unlist(all.asmap.sim.lengths))




#make heatmap of simulated data for rf and lod scores 
asmap1 = gen.asmap.for.sim(list.of.first.gen.recomb1a.pop96[[100]])




pdf("rotation1scripts_v4/plots/simulation/rf.heatmap.sim.dat.1a.pdf", 10, 10)
heatMap(asmap1, what = "rf")
dev.off()


#### 96 individual genetic map ####

#first read in data using rQTL function
axiom.data = read.cross(format = "csv", dir = "./", file = "rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/cxa651_redone.csv", genotypes = c("A", "H", "B"), 
                                                estimate.map = F)

#Specify the Filial generation of the cross with the F.gen parameter
axiomdata = convert2bcsft(axiom.data, BC.gen = 0, F.gen = 2, estimate.map = F)


axiomdata_map = mstmap(axiomdata, trace = T, dist.fun = "kosambi", id = "probeset_id", p.value = 1e-13, anchor = F)

realdat.heatmap = heatMap(axiomdata_map, "1.8")

pdf("rotation1scripts_v4/plots/simulation/rf.heatmap.real.dat.1a.pdf", 10, 10)
heatMap(axiomdata_map, "1.8", what = "rf")
dev.off()

ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/tester.eps", g, width = 17.4, height = 9, units = "cm", device = cairo_pdf)

g = rownames(pull.map(axiomdata_map, "1.8", as.table = T))

g == colnames(rf.dat1)

#trying to make a recombination frequency heatmap with ggplot... doesn't seem to work
rf.dat1 = pull.rf(axiomdata_map, what = "rf", "1.8")
rf.dat1 = as.data.frame(rf.dat1)

rf.dat1[is.na(rf.dat1)] = 0

ggplot(rf.dat1) + geom_tile()
library(reshape2)
rf.dat2 = melt(rf.dat1)

rf.dat2$variable = unlist(lapply(1:224, function(x) rep(x, 224)))
rf.dat2$var2 = rep(1:224, 224)
colnames(rf.dat2)[2] = "RF"

plot3 = ggplot(rf.dat2, aes(x = variable, y = var2)) + geom_tile(aes(fill = RF)) + 
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() + ggtitle("(a)") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.line = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust = 0.05))
plot3


pdf("rotation1scripts_v4/plots/test1.pdf", 10, 10)
plot3
dev.off()



sim.heatmap

grid.arrange(realdat.heatmap, sim.heatmap)

nrow(q2[[20]])

#length of real data in centimorgans
max(pull.map(axiomdata_map, as.table = T, chr = "1.8")$pos)

#load pedigreeSim data
load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2.noselec.pop300")

sim1.1 = all.sim.files[[1]]

sim1.1 = convert.sim.data.to.rqtl.format(sim1.1)
write.csv(sim1.1, "rotation1scripts_v4/temp/sim1.1.csv", row.names = F)
sim1.2 = rqtl.read("rotation1scripts_v4/temp/sim1.1.csv")

sim1.2rf = pull.rf(sim1.2, what = "rf")
sim1.2rf2 = melt(sim1.2rf)
sim1.2rf2[is.na(sim1.2rf2)] = 0
sim1.2rf2$Var1 = unlist(lapply(1:ncol(sim1.2rf), function(x) rep(x, ncol(sim1.2rf))))
sim1.2rf2$Var2 = rep(1:ncol(sim1.2rf), ncol(sim1.2rf))
colnames(sim1.2rf2)[3] = "RF"

plot4 = ggplot(sim1.2rf2, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = RF)) + 
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() + ggtitle("(b)") +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.line = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust = 0.05)) 
    
plot4

head(rf.dat2)
rf.dat2.1 = rf.dat2[, c(1, 3, 2)]
colnames(rf.dat2.1) = colnames(sim1.2rf2)

rf.dat2.1$type = "(a)"
sim1.2rf2$type = "(b)"

all.rf.dat = rbind(rf.dat2.1, sim1.2rf2)


all.rf.plot = ggplot(all.rf.dat, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = RF)) + facet_grid(. ~ type) +
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.line = element_blank(), axis.title = element_blank(), plot.title = element_text(hjust = 0.05)) +
    theme(strip.text = element_text(size = 13, hjust = 0.04), strip.background = element_rect(colour = "white", fill = "#ffffff"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

all.rf.plot

all.rf.plot2 = ggplot(all.rf.dat, aes(x = Var1, y = Var2)) + geom_tile(aes(fill = RF)) + facet_grid(. ~ type) +
    scale_fill_gradientn(colours = c("#AF0E0E", "#FFFC84", "#045ECC")) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.05)) +
    theme(strip.text = element_text(size = 13, hjust = 0.04), strip.background = element_rect(colour = "white", fill = "#ffffff"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_continuous(expand = c(0.01, 0.01)) + scale_y_continuous(expand = c(0.01, 0.01)) +
    xlab("Marker number") + ylab("Marker number")

all.rf.plot2

tiff("rotation1scripts_v4/plots/simulation/pedsim/Fig1v3_plos.tiff", units = "cm", width = 17.4, height = 7.3, res = 410)
all.rf.plot2
dev.off()






g2 = arrangeGrob(plot3, plot4, ncol = 2)

ggsave(file = "rotation1scripts_v4/plots/simulation/pedsim/heatmap.real.vs.sim.pdf", g2, width = 17.4, height = 9, units = "cm", device = cairo_pdf)

heatmap.rf.plots1 = list(plot3, plot4)
s(heatmap.rf.plots1, "rotation1scripts_v4/saved.objects/heatmap.rf.plots1", "26032019.asmap.heatmap.a.x.c.f2.R")

tiff("rotation1scripts_v4/plots/simulation/pedsim/heatmap.rf.real.vs.sim.tiff", units = "cm", width = 17.4, height = 7.3, res = 1000)
grid.arrange(plot3, plot4, ncol = 2)
dev.off()

tiff("heatmap.real.vs.simv2.tiff", units = "cm", width = 17.4, height = 9, res = 100)
plot3
dev.off()



length(colnames(sim1.2rf))




pdf("rotation1scripts_v4/plots/simulation/pedsim/heatmap.f2.300pop.no.selec.pdf", 10, 10)
heatMap(sim1.2, what = "rf")
dev.off()
