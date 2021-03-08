load("rotation1scripts_v4/saved.objects/all.selec.plots")
for(i in 1:1000){
    png(p("rotation1scripts_v4/plots/simulation/selec/plot", i, ".png"), 1000, 1000)
    print(all.selec.plots[[i]])
    dev.off()
}

load("rotation1scripts_v4/saved.objects/peak.dist")


for(i in 1:length(peak.dist)){
    pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.peak.dist/hist", i, ".pdf"), 8, 8)
    hist(peak.dist[[i]], breaks = 100, xlim = c(0, 225), main = "Histogram of peak segregation distortion in 1000 simulations", xlab = "Locus")
    dev.off()
}


load("rotation1scripts_v4/saved.objects/peak.dist2")


for(i in 1:length(peak.dist2)){
    pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.peak.dist/hist_noselec_", i, ".pdf"), 8, 8)
    hist(peak.dist2[[i]], breaks = 100, xlim = c(0, 225), main = "Histogram of peak segregation distortion in 1000 simulations", xlab = "Locus")
    dev.off()
}

pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.peak.dist/hist_noselec_", 6, ".pdf"), 8, 8)
hist(peak.dist2[[6]], breaks = 100, xlim = c(0, 399), main = "Histogram of peak segregation distortion in 1000 simulations", xlab = "Locus")
dev.off()


load("rotation1scripts_v4/saved.objects/magnitude.dist")

for(i in 1:length(magnitude.dist)){
    pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.mag.dist/hist_mag_", i, ".pdf"), 8, 8)
    hist(magnitude.dist[[i]], breaks = 100, xlim = c(0, 0.4), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")
    dev.off()
}

hist(magnitude.dist[[1]], breaks = 100, xlim = c(0, 0.3), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")

load("rotation1scripts_v4/saved.objects/magnitude.dist2")

for(i in 1:length(magnitude.dist2)){
    pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.mag.dist/hist_mag_noselec_", i, ".pdf"), 8, 8)
    hist(magnitude.dist2[[i]], breaks = 100, xlim = c(0, 0.4), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")
    dev.off()
}

hist(magnitude.dist2[[7]], breaks = 100, xlim = c(0, 0.3), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")



load("rotation1scripts_v4/saved.objects/list.of.selec.10k.stat")
load("rotation1scripts_v4/saved.objects/list.of.selec.1k.stat")

#grab magnitude of distortion for 1k selection populations
mag.dist3 = lapply(list.of.selec.1k.stat, function(x) x[[5]])

#grab magnitude of distortion for 10k selection populations
mag.dist4 = lapply(list.of.selec.10k.stat, function(x) x[[5]])

hist(mag.dist3[[6]], breaks = 100, xlim = c(0, 0.3), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")
hist(mag.dist4[[1]], breaks = 100, xlim = c(0, 0.3), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")

for(i in 1:length(mag.dist3)){
    pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.mag.dist/hist_mag_selec_pop1k_", i, ".pdf"), 8, 8)
    hist(mag.dist3[[i]], breaks = 100, xlim = c(0, 0.4), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")
    dev.off()
}

for(i in 1:length(mag.dist4)){
    pdf(p("rotation1scripts_v4/plots/simulation/selec/hist.mag.dist/hist_mag_selec_pop10k_", i, ".pdf"), 8, 8)
    hist(mag.dist4[[i]], breaks = 100, xlim = c(0, 0.4), main = "Histogram of magnitude of peak segregation distortion in 1000 simulations", xlab = "Magnitude")
    dev.off()
}



magnitude.dist.no.selec = magnitude.dist2
mag.vs.popsize = data.frame(c(mean(magnitude.dist.no.selec[[1]]), mean(magnitude.dist.no.selec[[2]]), mean(magnitude.dist.no.selec[[3]])), c(300, 1000, 10000))
colnames(mag.vs.popsize) = c("mean.mag", "pop.size")


load("rotation1scripts_v4/saved.objects/recomb.sim/simulations.recomb.any.pop.50.to.10000")

simulations.recomb.any.pop.50.to.10000

check.magnitude.of.distortion = function(x){
    #what is the magnitude of distortion in the simulations at the peak
    #args:
    #x - a list of genotyping dataframes 
    g = unlist(lapply(x, function(y){
        g = convert.recomb.to.seg(y)
        g1 = abs(g - 0.5)
        max(g1)
    }))
    
    g
}

simulations.recomb.any.pop.50.to.10000.mag = mclapply(simulations.recomb.any.pop.50.to.10000, function(x){
    check.magnitude.of.distortion(x)
}, mc.cores = 30)

load("rotation1scripts_v4/saved.objects/simulations.recomb.any.pop.50.to.10000.mag")

lapply(simulations.recomb.any.pop.50.to.10000.mag, function(x){
    mean(na.omit(x))
})


load("rotation1scripts_v4/saved.objects/simulations.recomb.any.pop.50.to.1000.mag")
load("rotation1scripts_v4/saved.objects/simulations.recomb.any.pop.1050.to.4000.mag")

all.sim.mags1 = c(simulations.recomb.any.pop.50.to.1000.mag, simulations.recomb.any.pop.1050.to.4000.mag)

all.sim.mags.means = unlist(lapply(all.sim.mags1, mean))

plot(seq(50, 4000, 50), all.sim.mags.means)


