#execution of recombination / segregation distortion simulation
setwd("~/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

load("rotation1scripts_v4/saved.objects/recombination.profiles1a")

#### ANALYSIS ####

library(gridExtra)

any.1000 = gen.data(185, 300)

g = convert.recomb.to.seg2(any.1000)

#make 8 different F2 generations
list.of.first.gen = mclapply(1:8, function(x){
    g = gen.data(185, 300)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob = mclapply(1:8, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.recomb.prof = mclapply(1:8, function(x){
    g = gen.data(120, 300, recombination.profiles = recombination.profiles, selection.pos = 60)
    g
}, mc.cores = corestouse)


list.of.first.gen.w.selec = mclapply(1:8, function(x){
    g = gen.data(185, 300, selection = 10, selection.pos = 170)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 10, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec2 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 8, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec3 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 7, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec4 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 6, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec5 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 5, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec6 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 4, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec7 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 3, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec8 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recomb.diffs, selection = 2, selection.pos = 110)
    g
}, mc.cores = corestouse)

list.of.first.gen.w.prob.w.selec2 = mclapply(1:1000, function(x){
    g = gen.data(120, 300, recombination.profiles = recombination.profiles, selection = 10, selection.pos = 60)
    g
}, mc.cores = corestouse)


seg.ratios1 = lapply(list.of.first.gen.w.prob.w.selec, convert.recomb.to.seg)
markers.w.max.seg = lapply(seg.ratios1, function(x){
    distx = abs(x - 1)
    round(mean(which(distx == max(distx))))
})

tab1 = sort(table(unlist(markers.w.max.seg)))

sum(as.numeric(tab1[which(as.numeric(names(tab1)) < 120 & as.numeric(names(tab1)) > 100)]))

plots1 = make.plots(list.of.first.gen.w.prob.w.selec3[1:8])
pdf("rotation1scripts_v4/plots/seg.dist.simulations/marker120ind300/f2.recomb.prof.1-7selecpos100v2.pdf", 16, 11)
do.call(grid.arrange, plots1)
dev.off()

#iterate new generations from an F2, simulates SSD
list.of.gen = lapply(1:8, function(x){
    a2 = gen.next.generation(list.of.first.gen[[1]], 185)
    # print(which(a2 == "H")) #heterozygosity decreases each generation
    any.1000 <<- a2
}, mc.cores = corestouse)


list.of.gen.no.initial.seg = mclapply(1:8, function(x){
    a2 = gen.next.generation(list.of.first.gen[[2]], 185)
    # print(which(a2 == "H")) #heterozygosity decreases each generation
    any.1000 <<- a2
}, mc.cores = corestouse)


list.of.gen.w.prob = mclapply(1:8, function(x){
    a2 = gen.next.generation(list.of.first.gen.w.prob[[1]], 120, recombination.profiles = recomb.diffs, selection = 10, selection.pos = 110)
    # print(which(a2 == "H")) #heterozygosity decreases each generation
    any.1000 <<- a2
}, mc.cores = corestouse)


test1a.1 = gen.data(225, 94, recombination.profiles = recombination.profiles1a)

testers1 = mclapply(1:8, function(x){
    gen.data(225, 94, recombination.profiles = recombination.profiles1a, selection = 3, selection.pos = 100)
}, mc.cores = corestouse)

testers1.plots = make.plots(testers1)
do.call(grid.arrange, testers1.plots)

any.1000 = testers1[[2]]
count1 = make.counter()
list.of.gen.no.initial.seg.w.prob = lapply(1:8, function(x){
    print(count1())
    a2 = gen.next.generation(any.1000, 225, recombination.profiles = recombination.profiles1a, selection = 3, selection.pos = 100)
    # print(which(a2 == "H")) #heterozygosity decreases each generation
    any.1000 <<- a2
})

g123 = make.plots(list.of.gen.no.initial.seg.w.prob)
do.call(grid.arrange, g123)

list.first.gen.recomb.anywhere = mclapply(1:1000, function(x){
    g = gen.data(225, 300)
    g
}, mc.cores = corestouse)

which.seg = unlist(lapply(list.first.gen.recomb.anywhere, function(x){
    if(length(check.for.distorted.markers(x)) > 0) return(T)
    return(F)
}))

length(which(which.seg))

lapply(list.of.gen.no.initial.seg.w.prob, function(x){
    sum(unlist(lapply(x, function(y){
        length(which(y == "H"))
    })))
})

write.csv(list.of.gen.no.initial.seg.w.prob[[1]], "rotation1scripts_v4/original_data/simulated.genotype.data/1gen.csv", row.names = F)
write.csv(list.of.gen.no.initial.seg.w.prob[[2]], "rotation1scripts_v4/original_data/simulated.genotype.data/2gen.csv", row.names = F)

plots2 = make.plots(list.of.gen.w.prob)
pdf("rotation1scripts_v4/plots/seg.dist.simulations/marker120ind300/fseries.recomb.5a.recomb.profile.selec.110v2.pdf", 16, 11)
do.call(grid.arrange, plots2)
dev.off()

# s(recomb.anywhere.df, "rotation1scripts_v4/saved.objects/recomb.anywhere.df", "recombination.modelling.R")


#### LARGER DATASETS ####

list.of.first.gen.recomb.any = mclapply(1:1000, function(x){
    g = gen.data(185, 300)
    g
}, mc.cores = corestouse)

list.of.first.gen.recomb.prob = mclapply(1:1000, function(x){
    g = gen.data(225, 300, recombination.profiles = recombination.profiles1a)
    g
}, mc.cores = corestouse)

g = find.recomb.w.seg.dist(list.of.first.gen.recomb.any)
length(which(g))

g = find.recomb.w.seg.dist(list.of.first.gen.recomb.prob)
length(which(g))



load("rotation1scripts_v4/saved.objects/list.of.first.gen.recomb.5b.225.300")
plot123 = make.plots(list.of.first.gen.recomb.5b.225.300)



#### ANIMATION ####

testdata1 = list.of.first.gen[[1]]
for(i in 1:nrow(testdata1)){
    g = testdata1[1:i, ]
    seg.ratios = unlist(lapply(g, function(x){
        a = length(which(x == "A"))
        b = length(which(x == "B"))
        
        a / b
    }))
    
    
    newdata1 = data.frame(seg.ratios, 1:length(seg.ratios), stringsAsFactors = F)
    colnames(newdata1) = c("seg.ratio", "marker")
    jim = ggplot(newdata1, aes(x = marker, y = seg.ratio)) + geom_point() + coord_cartesian(ylim = c(0.2, 1.8)) + geom_hline(yintercept = 1)
    
    png(p("rotation1scripts_v4/plots/seg.dist.simulations/animation/", i, ".png"), 1000, 900)
    print(jim)
    dev.off()
    
}


#### NEW DATASETS ####

corestouse = 20

list.of.first.gen.recomb.anywherepop300 = mclapply(1:1000, function(x){
    g = gen.data(225, 300)
    g
}, mc.cores = corestouse)

save(list.of.first.gen.recomb.anywherepop300, "rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop300")

list.of.first.gen.recomb.anywherepop1000 = mclapply(1:1000, function(x){
    g = gen.data(225, 1000)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.recomb.anywherepop1000, "rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop1000", "recombination.sim.execution.R")

list.of.first.gen.recomb.anywherepop10000 = mclapply(1:1000, function(x){
    g = gen.data(225, 10000)
    g
}, mc.cores = corestouse)

s(list.of.first.gen.recomb.anywherepop10000, "rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.anywherepop10000", "recombination.sim.execution.R")

p1223 = make.plots(list.of.first.gen.recomb.anywherepop10000[1:8])
do.call(grid.arrange, p1223)

list.of.first.gen.recomb.prob.10000pop = mclapply(1:1000, function(x){
    g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a)
    g
}, mc.cores = corestouse)

load("rotation1scripts_v4/saved.objects/list.of.first.gen.recomb.1a.225.300")
load("rotation1scripts_v4/saved.objects/list.of.first.gen.recomb.5b.225.300")

g2 = make.plots(list.of.first.gen.recomb.5b.225.300[1:8])
do.call(grid.arrange, g2)


#### LOAD GNU PARALLEL FILES ####

all.recomb.prob2 = list()
for(i in 1:50){
    
    load(p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob", i))
    all.recomb.prob2 = c(all.recomb.prob2, list.of.first.gen.recomb.prob)
    
}

all.recomb.prob.pop10000 = list()
for(i in 1:50){
    
    load(p("rotation1scripts_v4/saved.objects/recomb.sim/list.of.first.gen.recomb.prob10000pop", i))
    all.recomb.prob.pop10000 = c(all.recomb.prob.pop10000, list.of.first.gen.recomb.prob.10000pop)
    
}


length(all.recomb.prob.pop10000)
v(all.recomb.prob.pop10000[[1]])

all10000 = make.plots(all.recomb.prob.pop10000[1:9])
do.call(grid.arrange, all10000)


#### 



#check how many of the recombination profiles are unique






