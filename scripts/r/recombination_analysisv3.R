# setwd("~/Google Drive/R/scripts/rotation1scripts/") #mac

setwd("E:/phd.project.main/rotation1scripts_v4/")

library(gridExtra)
library(ggsignif)
library(dunn.test)
library(foreign)
library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(reshape2)

source("rotation1scripts_v4/scripts/r/recombination_analysis_functions.R")
source("scripts/r/functions.R")

#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####

#grab mean distance from centre for each quad
grab.mean.dist = function(list.of.df, variable){
    #args:
    # list.of.df = a list of dataframes
    # variable = the variable in each dataframe to calculate the mean for
    as.character(lapply(list.of.df, function(i) mean(as.numeric(na.omit(i[[variable]])))))
}

grab.sd.dist = function(list.of.df, variable){
    #args:
    # list.of.df = a list of dataframes
    # variable = the variable in each dataframe to calculate the mean for
    as.character(lapply(list.of.df, function(i) sd(as.numeric(na.omit(i[[variable]])))))
}

grab.se.dist = function(list.of.df, variable){
    #args:
    # list.of.df = a list of dataframes
    # variable = the variable in each dataframe to calculate the mean for
    as.character(lapply(list.of.df, function(i) sd(as.numeric(na.omit(i[[variable]])))/sqrt((nrow(i)-1))))
}

get.central.lg.position = function(lg){
    #get the position of the central marker in a linage group (e.g. for an lg of 100 markers, returns 50)
    filtered.genotype.df=initial.data[[1]][1:3, which(initial.data[[1]][1,]==lg)]
    
    central.position=ncol(filtered.genotype.df)/2
    # browser()
    return(central.position)
}



add.dist.col = function(histlist, lg){
    #adds distance, measured as number of markers, of recombination event from central marker
    histlist$dist.from.centre = abs(as.numeric(list.of.central[lg])-as.numeric(histlist$position))
    return(histlist)
}

grab.info = function(histlist){ #this could probably just be written as an anonymous function
    info.needed = histlist[, c(5, 7)]
    return(info.needed)
}


perform.stats = function(lg.to.test){
    # browser()
    histlists=lapply(genotypelists, makehistlist, linkagegroup=lg.to.test, return_recombinationcounts=T)
    histlists = lapply(histlists, add.dist.col, lg=lg.to.test)
    # browser()
    
    #filter histlist dataframe for number of recombination events and distance of recombination event from central marker.
    anova.df.sep=lapply(histlists, grab.info)
    
    all.anova.df = newdf(colnames = c("num_recomb","dist.from.centre","quad"))
    for(i in 1:length(anova.df.sep)){
        anova.df.sep[[i]]$quad = i
        all.anova.df = rbind(all.anova.df, anova.df.sep[[i]])
    }
    
    
    all.anova.df$quad=as.factor(all.anova.df$quad)
    
    all.anova.df=all.anova.df[-nrow(all.anova.df),]
    all.anova.df$dist.from.centre = as.numeric(all.anova.df$dist.from.centre)
    
    ano = lm(all.anova.df$dist.from.centre ~ all.anova.df$quad)
    summary(ano)
    
    md = grab.mean.dist(anova.df.sep, "dist.from.centre")
    sd = grab.sd.dist(anova.df.sep, "dist.from.centre")
    
    
    a1=aov(all.anova.df$dist.from.centre ~ all.anova.df$quad)
    summary(a1)
    TukeyHSD(x=a1, 'all.anova.df$quad', conf.level=0.95)
    
    #do non-parametric test
    kt = kruskal.test(all.anova.df$dist.from.centre ~ all.anova.df$quad)
    # dt = dunn.test(all.anova.df$dist.from.centre, all.anova.df$quad, method = "bonferroni")
    # results = list(md, kt, sd, dt)
    results = list(md, kt, sd)
    return(results)
}

#### SETUP ####

list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)

load.ind.to.remove = function(x){
    read.csv(paste("rotation1scripts_v4/processed_data/listofindividualsover.recomb.thres.q", x, 
                                 ".csv", sep = ""), header = T,
                     stringsAsFactors = F)
}

ind.to.remove = lapply(1:4, load.ind.to.remove)

load.data = function(x) read.csv(paste("rotation1scripts_v4/processed_data/genotypes/alexq", x, "_man_cur_sk_map_genogeno.csv", sep = ""), 
                                                                 header = T, stringsAsFactors = F)

initial.data = lapply(1:4, load.data)

genotypelists = lapply(initial.data, makegenotypelist)

#remove individuals over threshold from genotypelists
for(quad in 1:4){
    for(i in 1:nrow(ind.to.remove[[quad]])){
        rem = which(genotypelists[[quad]]$individual_id == ind.to.remove[[quad]]$individual[i] & genotypelists[[quad]]$lg == ind.to.remove[[quad]]$lg[i])
        genotypelists[[quad]] = genotypelists[[quad]][-rem, ]
    }
}

list.of.central=1:35
#get the middle marker position for each linkage group
list.of.central=as.numeric(as.character(lapply(list.of.central, get.central.lg.position)))



#     ____________________________________________________________________________
#     TEST FOR DIFFERENCES IN MARKER POSITION                                                                ####


#setting up the ANOVA to test distance from centre of LG between quads -- change lg.to.test as appropriate
list.of.one.arm.lgs.to.test = c(10, 21)
lg.to.test = list.of.full.lgs.to.test[[1]]

# histlists.old = lapply(list.of.full.lgs.to.test, function(x){
#     lapply(genotypelists, function(y){
#         makehistlist(x, y, T)
#     })
# })


all.stats.results = lapply(list.of.full.lgs.to.test, perform.stats)

resultsdf = newdf(colnames = c("LG", "mean_dist_q1", "mean_dist_q2", "mean_dist_q3", "mean_dist_q4", "sd_dist_q1", "sd_dist_q2", "sd_dist_q3", "sd_dist_q4", "kw_df", "kw_csq", "kw_p"))

for(i in seq_along(all.stats.results)){
    resultsdf$LG[i] = list.of.full.lgs.to.test[i]
    resultsdf$mean_dist_q1[i] = all.stats.results[[i]][[1]][1]
    resultsdf$mean_dist_q2[i] = all.stats.results[[i]][[1]][2]
    resultsdf$mean_dist_q3[i] = all.stats.results[[i]][[1]][3]
    resultsdf$mean_dist_q4[i] = all.stats.results[[i]][[1]][4]
    resultsdf$sd_dist_q1[i] = all.stats.results[[i]][[3]][1]
    resultsdf$sd_dist_q2[i] = all.stats.results[[i]][[3]][2]
    resultsdf$sd_dist_q3[i] = all.stats.results[[i]][[3]][3]
    resultsdf$sd_dist_q4[i] = all.stats.results[[i]][[3]][4]
    resultsdf$kw_df[i] = all.stats.results[[i]][[2]][[2]]
    resultsdf$kw_csq[i] = all.stats.results[[i]][[2]][[1]]
    resultsdf$kw_p[i] = all.stats.results[[i]][[2]][[3]]
    resultsdf = add_row(resultsdf)
}

resultsdf$kw_p = as.numeric(resultsdf$kw_p)
resultsdf$kw_p_rounded = round(resultsdf$kw_p, 5)
resultsdf = resultsdf[-nrow(resultsdf),]
write.csv(resultsdf, "processed_data/statsresults.csv", row.names = F)

#setup dataframe for bar plot
results.df.melted = melt(resultsdf, id = "LG")
results.df.melted = cbind(results.df.melted[1:40,], results.df.melted[41:80,])
results.df.melted = results.df.melted[,-c(4,5)]
colnames(results.df.melted) = c("LG", "variable", "mean", "sd")
results.df.melted = filter(results.df.melted, LG == 17)
results.df.melted$mean = as.numeric(results.df.melted$mean)
results.df.melted$sd = as.numeric(results.df.melted$sd)

popsizes = c(83, 73, 78, 81)

results.df.melted$se = results.df.melted$sd / popsizes

results.df.melted$temp = c(10, 26, 14, 28)
results.df.melted = results.df.melted[sort(results.df.melted$temp, index.return = T)$ix,]
results.df.melted$temp = as.factor(results.df.melted$temp)

#do barplot
pdf(file = "plots/mean_distance.pdf")
plot = ggplot(data = results.df.melted, aes(x = temp, y = mean)) + geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .2, size = 1.2, position = position_dodge(.9)) + 
    geom_signif(comparisons = list(c("10", "26"), c("26", "14"),
                                                                 c("10", "28")),
                            map_signif_level = T, annotations = c("***", "*", "**"),
                            y_position = c(24, 26, 28)) +
    theme_classic() + scale_y_continuous(name = "Mean distance of recombination event from centre of LG (# markers)", 
                                                                             limits = c(0, 30)) +
    scale_x_discrete(name = "Temperature treatment ((deg)C)", labels = c("10", "14", "26", "28"))
print(plot)
dev.off()



#     ____________________________________________________________________________
#     RECOMBINATION FREQUENCY ANALYSIS                                                                                ####



test.recombination.frequency = function(lg.to.test2) {
    test.geno.list = lapply(genotypelists, makehistlist, linkagegroup = lg.to.test2, return_recombinationcounts = T)    
    test.geno.list = lapply(test.geno.list, function(x){
        recomb.freq = newdf(c("individual_id", "num_recomb_events", "chromosome"))
        for(i in 1:length(unique(x$individual))){
            individual.events = filter(x, individual == unique(x$individual)[i])
            recomb.freq$individual_id[i] = individual.events$individual[1]
            recomb.freq$num_recomb_events[i] = nrow(individual.events)
            recomb.freq$chromosome[i] = lg.to.test2
            recomb.freq = add_row(recomb.freq)
        }
        return(recomb.freq)
    })
    
    for(i in 1:length(test.geno.list)){
        test.geno.list[[i]]$quad = i
        test.geno.list[[i]] = test.geno.list[[i]][-nrow(test.geno.list[[i]]),]
    }
    
    all.recomb.freq = rbind(test.geno.list[[1]], test.geno.list[[2]], test.geno.list[[3]], test.geno.list[[4]])
    all.recomb.freq$quad = as.factor(all.recomb.freq$quad)
    all.recomb.freq$num_recomb_events = as.numeric(all.recomb.freq$num_recomb_events)
    return(all.recomb.freq)
}

all.freq.tests = lapply(list.of.full.lgs.to.test, test.recombination.frequency)

all.freq.tests2 = lapply(all.freq.tests, function(x){
    kt = kruskal.test(x$num_recomb_events ~ x$quad)
    dt = dunn.test(x$num_recomb_events, x$quad, method = "bonferroni")    
    results = list(kt, dt)
    return(results)
})

frequency.results = newdf(c("lg", "kw_csq", "kw_df", "kw_p", "dunn_comp1", "dunn_comp1_adj_p", "dunn_comp2", "dunn_comp2_adj_p", 
                                                        "dunn_comp3", "dunn_comp3_adj_p", "dunn_comp4", "dunn_comp4_adj_p",
                                                        "dunn_comp5", "dunn_comp5_adj_p", "dunn_comp6", "dunn_comp6_adj_p"))

for(i in 1:length(all.freq.tests2)){
    frequency.results$lg[i] = list.of.full.lgs.to.test[i]
    frequency.results$kw_csq[i] = all.freq.tests2[[i]][[1]][[1]]
    frequency.results$kw_df[i] = all.freq.tests2[[i]][[1]][[2]]
    frequency.results$kw_p[i] = all.freq.tests2[[i]][[1]][[3]]
    frequency.results$dunn_comp1[i] = all.freq.tests2[[i]][[2]][[5]][[1]]
    frequency.results$dunn_comp1_adj_p[i] = all.freq.tests2[[i]][[2]][[4]][[1]]
    frequency.results$dunn_comp2[i] = all.freq.tests2[[i]][[2]][[5]][[2]]
    frequency.results$dunn_comp2_adj_p[i] = all.freq.tests2[[i]][[2]][[4]][[2]]
    frequency.results$dunn_comp3[i] = all.freq.tests2[[i]][[2]][[5]][[3]]
    frequency.results$dunn_comp3_adj_p[i] = all.freq.tests2[[i]][[2]][[4]][[3]]
    frequency.results$dunn_comp4[i] = all.freq.tests2[[i]][[2]][[5]][[4]]
    frequency.results$dunn_comp4_adj_p[i] = all.freq.tests2[[i]][[2]][[4]][[4]]
    frequency.results$dunn_comp5[i] = all.freq.tests2[[i]][[2]][[5]][[5]]
    frequency.results$dunn_comp5_adj_p[i] = all.freq.tests2[[i]][[2]][[4]][[5]]
    frequency.results$dunn_comp6[i] = all.freq.tests2[[i]][[2]][[5]][[6]]
    frequency.results$dunn_comp6_adj_p[i] = all.freq.tests2[[i]][[2]][[4]][[6]]
    frequency.results = add_row(frequency.results)
}

write.csv(frequency.results, "processed_data/frequencyresults.csv", row.names = F)


#     ____________________________________________________________________________
#     cM vs Mb plots                                                                                                                    ####

#NB. REQUIRES estimated.maps from asmapv2.R and list.of.physical.maps from below

windows()
par(mfrow=c(2,2))
grabq = function(x) c(x, x+10, x+20, x+30)
lgindex = 5


for(i in grabq(lgindex)){
    plot(list.of.physical.maps[[lgindex]], estimated.maps[[i]]$pos, ylim=c(0,200))
}

dfphys.cm = newdf(c("phys.pos","cM", "quad"))

for(i in 1:length(list.of.physical.maps[[1]])){
    dfphys.cm$phys.pos[i] = list.of.physical.maps[[1]][i]
    dfphys.cm = add_row(dfphys.cm)
}


lg = 5
g = data.frame(c(list.of.physical.maps[[lg]], list.of.physical.maps[[lg]], list.of.physical.maps[[lg]], list.of.physical.maps[[lg]]), 
                     c(estimated.maps[[lg]]$pos, estimated.maps[[lg+10]]$pos, estimated.maps[[lg+20]]$pos, estimated.maps[[lg+30]]$pos),
                     c(rev((estimated.maps[[lg]]$pos/max(estimated.maps[[lg]]$pos))*100), 
                         rev((estimated.maps[[lg+10]]$pos/max(estimated.maps[[lg+10]]$pos)*100)), 
                         rev((estimated.maps[[lg+20]]$pos/max(estimated.maps[[lg+20]]$pos)*100)), 
                         rev((estimated.maps[[lg+30]]$pos/max(estimated.maps[[lg+30]]$pos))*100)),
                     c(rep(1,length(list.of.physical.maps[[lg]])), rep(2, length(list.of.physical.maps[[lg]])), 
                         rep(3, length(list.of.physical.maps[[lg]])), rep(4, length(list.of.physical.maps[[lg]]))))

colnames(g) = c("phys", "cm", "cm.norm", "quad")
g$quad = as.factor(g$quad)

conv.q.to.temp = function(x){
    if(x == 1) y = 10
    if(x == 2) y = 26
    if(x == 3) y = 14
    if(x == 4) y = 28
    return(y)
}

g$temp = lapply(g$quad, conv.q.to.temp)
g$temp = as.factor(as.character(g$temp))
colnames(g) = c("phys", "cm", "cm.norm", "quad", "Temperature")

#THIS PLOT IS TO BE USED IN CONJUNCTION WITH plot2 IN snpdensityplot.R
plot1 = ggplot(data = g, aes(x = phys, y = cm.norm, group = Temperature, color = Temperature, linetype = Temperature)) + geom_line(size = 1.2) + theme_classic() +
    xlab("Physical Position (Mb)") + ylab("Normalized Genetic Distance (cM)") + scale_x_continuous(breaks = c(seq(0, 800, 25)), expand = c(0.045, 0.045)) +
    scale_colour_discrete(name = "Temperature ((deg)C)", breaks = c(10, 14, 26, 28), labels = c(10, 14, 26, 28)) + 
    scale_shape_discrete(name = "Temperature ((deg)C)", breaks = c(10, 14, 26, 28), labels = c(10, 14, 26, 28)) +
    scale_linetype_discrete(name = "Temperature ((deg)C)", breaks = c(10, 14, 26, 28), labels = c(10, 14, 26, 28)) +
    theme(legend.position = c(0.85, 0.3), legend.background = element_rect(linetype = "solid", size = 0.5, color = "black"))

pdf(file = "plots/cm.vs.physical.distance/chromosome6b.w.snp.density.v2.pdf", width = 10, height = 10)
grid.arrange(plot1, plot2, nrow = 2) #SEE snpdensityplot.R for plot2
dev.off()

ggplot(data = g, aes(x = phys, y = cm, group = quad, color = quad, linetype = quad)) + geom_line(size = 1.2) + theme_classic() #change y from cm to cm.norm and back for normalized cM distances --- good for examining slope



#     ____________________________________________________________________________
#     ANOVA PHYSICAL DISTANCES                                                                                                ####

#ANOVA PHASE 2
#setting up anova with physical distances
lg.to.test = 15
list.of.corresponding.chromosomes.to.lgs = c("2A", "3A", "5A", "7B", "6B", "1A", "3B", "7A", "2D", "4A")
list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)

library(zoo)

list.of.physical.maps = Map(function(chromosome, LG){ #ALTER return.bp FLAG ACCORDINGLY
    # browser()
    phys1 = genetictophysical(chromosome, colnames(initial.data[[1]][1,which(initial.data[[1]][1,]==LG)]), ".", 90, return.bp = F)
    # replacena(phys1)
    # length(na.approx(phys1))
    # length(phys1)
    na.approx(phys1, na.rm = F)
}, list.of.corresponding.chromosomes.to.lgs, list.of.full.lgs.to.test)



do.phys.stats = function(lg.to.test){
    
    #list.of.physical.maps = lapply(list.of.physical.maps, function(x) x/10^6)
    
    #lg.to.test = list.of.full.lgs.to.test[[6]]
    
    phys.positions = list.of.physical.maps[[which(list.of.full.lgs.to.test == lg.to.test)]]
    
    histlists=lapply(genotypelists, makehistlist, linkagegroup=lg.to.test, return_recombinationcounts=T)
    histlists = lapply(histlists, add.dist.col, lg=lg.to.test)
    
    histlists2 = lapply(histlists, function(x){
        x$phy.position = phys.positions[as.numeric(x$position)]
        return(x)
    })
    
    histlists2 = lapply(histlists2, function(x){
        x$dist.from.fifty = abs(50 - x$phy.position)
        return(x)
    })
    
    for(i in 1:length(histlists2)){
        histlists2[[i]]$quad = i
    }
    
    all.anova.phys.df = rbind(histlists2[[1]], histlists2[[2]], histlists2[[3]], histlists2[[4]])
    all.anova.phys.df$quad = as.factor(all.anova.phys.df$quad)
    all.anova.phys.df$dist.from.fifty = as.numeric(as.character(unlist(all.anova.phys.df$dist.from.fifty)))
    
    md = grab.mean.dist(histlists2, "dist.from.fifty")
    sd = grab.se.dist(histlists2, "dist.from.fifty")
    #se = grab.se.dist(histlists2, "dist.from.fifty")
    
    
    #try transformation
    # all.anova.phys.df$log.dist = log(all.anova.phys.df$dist.from.fifty)
    kt = kruskal.test(all.anova.phys.df$dist.from.fifty ~ all.anova.phys.df$quad)
    # dt = dunn.test(all.anova.phys.df$dist.from.fifty, all.anova.phys.df$quad, method = "bonferroni")
    # results = list(md, kt, sd, dt)
    results = list(md, kt, sd)
    return(results)
}

all.phys.stats.results = lapply(list.of.full.lgs.to.test, do.phys.stats)

resultsdf.phys = newdf(colnames = c("LG", "mean_dist_q1", "mean_dist_q2", "mean_dist_q3", "mean_dist_q4", "sd_dist_q1", "sd_dist_q2", "sd_dist_q3", "sd_dist_q4", "kw_df", "kw_csq", "kw_p"))

for(i in seq_along(all.phys.stats.results)){
    resultsdf.phys$LG[i] = list.of.full.lgs.to.test[i]
    resultsdf.phys$mean_dist_q1[i] = all.phys.stats.results[[i]][[1]][1]
    resultsdf.phys$mean_dist_q2[i] = all.phys.stats.results[[i]][[1]][2]
    resultsdf.phys$mean_dist_q3[i] = all.phys.stats.results[[i]][[1]][3]
    resultsdf.phys$mean_dist_q4[i] = all.phys.stats.results[[i]][[1]][4]
    resultsdf.phys$sd_dist_q1[i] = all.phys.stats.results[[i]][[3]][1]
    resultsdf.phys$sd_dist_q2[i] = all.phys.stats.results[[i]][[3]][2]
    resultsdf.phys$sd_dist_q3[i] = all.phys.stats.results[[i]][[3]][3]
    resultsdf.phys$sd_dist_q4[i] = all.phys.stats.results[[i]][[3]][4]
    resultsdf.phys$kw_df[i] = all.phys.stats.results[[i]][[2]][[2]]
    resultsdf.phys$kw_csq[i] = all.phys.stats.results[[i]][[2]][[1]]
    resultsdf.phys$kw_p[i] = all.phys.stats.results[[i]][[2]][[3]]
    resultsdf.phys = add_row(resultsdf.phys)
}

resultsdf.phys$kw_p = as.numeric(resultsdf.phys$kw_p)
resultsdf.phys$kw_p_rounded = round(resultsdf.phys$kw_p, 5)
resultsdf.phys = resultsdf.phys[-nrow(resultsdf.phys),]
write.csv(resultsdf.phys, "processed_data/statsresults.phys.csv", row.names = F)


resultsdf.phys.melted = melt(resultsdf.phys, "LG")
resultsdf.phys.melted = resultsdf.phys.melted[1:40, ]

resultsdf.phys.melted2 = melt(resultsdf.phys, "LG")
resultsdf.phys.melted2 = resultsdf.phys.melted2[41:80, ]



conv.q.to.temp2 = function(x){
    if(x == "mean_dist_q1") y = 10
    if(x == "mean_dist_q2") y = 26
    if(x == "mean_dist_q3") y = 14
    if(x == "mean_dist_q4") y = 28
    return(y)
}

resultsdf.phys.melted$temp = lapply(resultsdf.phys.melted$variable, conv.q.to.temp2)
resultsdf.phys.melted$temp = as.factor(as.character(resultsdf.phys.melted$temp))
resultsdf.phys.melted$value = as.numeric(resultsdf.phys.melted$value)
resultsdf.phys.melted$chr = lapply(resultsdf.phys.melted$LG, function(x) list.of.corresponding.chromosomes.to.lgs[which(list.of.full.lgs.to.test == x)])
resultsdf.phys.melted$chr = as.factor(as.character(resultsdf.phys.melted$chr))

resultsdf.phys.melted = cbind(resultsdf.phys.melted, resultsdf.phys.melted2$value)
colnames(resultsdf.phys.melted) = c(colnames(resultsdf.phys.melted)[-6], "sd.value")
resultsdf.phys.melted$sd.value2 = as.numeric(as.character(resultsdf.phys.melted$sd.value))





#plot clustered histogram for all LGs showing mean distance from 50 % of the chromosome
# pdf(file = "plots/clustered.histogram.phys.distance/clusteredphys.pdf", width = 10, height = 10)
ggplot(resultsdf.phys.melted, aes(chr, y = value, ymin = value, ymax = value+sd.value2, fill = temp)) + 
    geom_bar(color = "black", stat = "identity", position = "dodge") +
    theme_classic() + coord_cartesian(ylim=c(25,41)) + scale_fill_brewer(palette = 1) +
    xlab("Chromosome") + ylab("Mean physical distance from the centre of chromosome (%)") +
    geom_errorbar(position = "dodge")
# dev.off()
    


a1=aov(all.anova.phys.df$dist.from.fifty ~ all.anova.phys.df$quad)
summary(a1)
TukeyHSD(x=a1, 'all.anova.phys.df$quad', conf.level=0.95)

####calculating standardised residuals
stres<- (a1$residuals - mean(a1$residuals))/ sd(a1$residuals)

shapiro.test(stres)
ks.test(stres, rnorm)

####Checking model assumptions:
#1. Residuals are normally distributed
windows (6,6); par (mfrow = c (2,2))    # for windows users
hist (stres)
qqnorm(stres, cex = 1.8, pch = 20); qqline(stres, lty=2, lwd = 2) 

#2. Variability of residuals should be uniform across the range of fitted values 
plot(stres ~ a1$fitted.values, pch = 20, cex = 2, cex.lab = 1.5) # a graphical check
fligner.test(stres ~ a1$fitted.values)                                                     # a quantitative test for the question - is there a pattern? 

#3. No strongly influential observations
windows (6,6); par (mfrow = c (2,2))    # for windows users
#quartz (6,6); par (mfrow = c (2,2))        # for Mac users
plot(a1)    #    4th panel shows the influence statistics


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### SETUP RF GENERAL LINEAR MODEL ANALYSIS                                                                    ####
popq1=83
popq2=73
popq3=78
popq4=81

LG.TO.ANALYSE=1

initial.data[[1]][,which(initial.data[[1]][1,]==LG.TO.ANALYSE)]
howmany=length(which(initial.data[[1]][1,]==LG.TO.ANALYSE)) #how many markers in the LG?

#make list of intervals for this particular LG
s1=seq(1,howmany,round(howmany/4))
s2=Map(c, s1[-length(s1)], s1[-1]) #https://stackoverflow.com/questions/41761999/create-a-sequence-of-intervals-in-r
s3=seq(1,4,2)
for(i in 1:4){
    s2[[i]][2]=(s2[[i]][2]-1)
}
s2 #list of intervals

(sum(countevents(15, genotypelists[[1]], s2[[1]][1], s2[[1]][2], initial.data[[1]]))/popq1)


g1=setup.glm.df.rf(genotypelists[[1]], initial.data[[1]], popq1, 10, s2)
g2=setup.glm.df.rf(genotypelists[[2]], initial.data[[2]], popq2, 26, s2)
g3=setup.glm.df.rf(genotypelists[[3]], initial.data[[3]], popq3, 14, s2)
g4=setup.glm.df.rf(genotypelists[[4]], initial.data[[4]], popq4, 28, s2)
allg=rbind(g1,g2,g3,g4)

allg=filter(allg, is.na(allg$interval)==F)
allg$interval=as.factor(allg$interval)
allg$temp=as.factor(allg$temp)

#note: this model doesn't yield any p-values. This is because there is not enough data / degrees of freedom
# see https://stats.stackexchange.com/questions/94078/why-do-i-not-get-a-p-value-from-this-anova-in-r for explanantion
model1=lm(allg$rf ~ allg$interval + allg$temp + allg$interval:allg$temp)
model1=lm(allg$rf ~ allg$interval + allg$temp) #works without the interaction term    
anova(model1)

#negative binomial regression model
model2=glm.nb(allg$rf ~ allg$interval + allg$temp + allg$interval:allg$temp)
summary(model2)

                                                                    ####calculating standardised residuals
stres<- (model1$residuals - mean(model1$residuals))/ sd(model1$residuals)



####Checking model assumptions:
#1. Residuals are normally distributed
windows (6,6); par (mfrow = c (2,2))    # for windows users
hist (stres)
qqnorm(stres, cex = 1.8, pch = 20); qqline(stres, lty=2, lwd = 2) 

#2. Variability of residuals should be uniform across the range of fitted values 
plot(stres ~ model1$fitted.values, pch = 20, cex = 2, cex.lab = 1.5) # a graphical check
fligner.test(stres ~ model1$fitted.values)                                                     # a quantitative test for the question - is there a pattern? 

#3. No strongly influential observations
windows (6,6); par (mfrow = c (2,2))    # for windows users
#quartz (6,6); par (mfrow = c (2,2))        # for Mac users
plot(model1)    #    4th panel shows the influence statistics

g1=setupglmdf(genotypelists[[1]], initial.data[[1]], popq1, 10)
g2=setupglmdf(genotypelists[[2]], initial.data[[2]], popq2, 26)
g3=setupglmdf(genotypelists[[3]], initial.data[[3]], popq3, 14)
g4=setupglmdf(genotypelists[[4]], initial.data[[4]], popq4, 28)

allg=rbind(g1,g2,g3,g4)




res.aov=aov(allg$rf ~ allg$interval * allg$temp)
summary(res.aov)






sum(countevents(20, genotypelists[[1]], 0, 25, initial.data[[1]]))/popq1
sum(countevents(20, genotypelists[[2]], 0, 25, initial.data[[2]]))/popq2
sum(countevents(20, genotypelists[[3]], 0, 25, initial.data[[3]]))/popq3
sum(countevents(20, genotypelists[[4]], 0, 25, initial.data[[4]]))/popq4

noq1.20=countevents(20, genotypelists[[1]], 0, 25, initial.data[[1]])
noq2.20=countevents(20, genotypelists[[2]], 18, 50, initial.data[[2]])
noq3.20=countevents(20, genotypelists[[3]], 18, 50, initial.data[[3]])
noq4.20=countevents(20, genotypelists[[4]], 18, 50, initial.data[[4]])
?t.test()
t.test(noq3.20, noq2.20)

noq1.10=countevents(10, genotypelists[[1]], 57, 91, initial.data[[1]])
noq2.10=countevents(10, genotypelists[[2]], 57, 91, initial.data[[2]])
noq3.10=countevents(10, genotypelists[[3]], 57, 91, initial.data[[3]])
noq4.10=countevents(10, genotypelists[[4]], 57, 91, initial.data[[4]])
noq2.10
noq3.10

t.test(noq1.10, noq2.10)

lg=9
startbin=1
endbin=15
noq1.15=countevents(lg, genotypelists[[1]], startbin,endbin, initial.data[[1]])
noq2.15=countevents(lg, genotypelists[[2]], startbin,endbin, initial.data[[2]])
noq3.15=countevents(lg, genotypelists[[3]], startbin,endbin, initial.data[[3]])
noq4.15=countevents(lg, genotypelists[[4]], startbin,endbin, initial.data[[4]])

df1=conv.vec.to.df(noq1.15, 1, 10)
df2=conv.vec.to.df(noq2.15, 2, 26)
df3=conv.vec.to.df(noq3.15, 3, 14)
df4=conv.vec.to.df(noq4.15, 4, 28)

alldf=rbind(df1,df2,df3,df4)
alldf$quad=as.factor(alldf$quad)
alldf$temp=as.factor(alldf$temp)
alldf.lm=lm(alldf$vector ~ alldf$quad)

a1=aov(alldf$vector ~ alldf$temp)
TukeyHSD(x=a1, 'alldf$temp', conf.level=0.95)



noq1.15
noq2.15

t.test(noq1.15, noq2.15)






#LGs that represent near-full chromosomes:
# LG1
# LG3 (although quite a large chunk is missing at the start of the chromosome)
# LG7 (again, quite a large chunk missing)
# LG9
# LG10 (represents one arm of the chromosome)
# LG15 (one arm)
# LG17
# LG20 (very big gap)
# LG23 (quite a large gap)



#     ____________________________________________________________________________
#     MAKE MULTILINE PLOT                                                                                                         ####


library(zoo) #loads na.approx function --- for approximating missing blast percentage values
source("scripts/r/functions.R")
source("scripts/r/recombination_analysis_functions.R")

filelist=list.files("processed_data/pairwiseLGcomparisons/allq_allq/")

lglist=c(1,3,5,7,9,10,15,17,20,23)
list.of.full.lgs.to.test = c(1, 3, 5, 7, 9, 10, 15, 17, 20, 23)
list.of.chromosome.assignments.for.lgs = c("2A", "3A", "5A", "7B", "6B", "1A", "3B", "7A", "2D", "4A")
for(lgtomakehistlistfor in lglist){
lgtomakehistlistfor=1
    blastfile=read.csv(paste("processed_data/pairwiseLGcomparisons/allq_allq/",filelist[grep(paste("_LG",lgtomakehistlistfor,"_", sep=""), filelist)], sep=""), header=T, stringsAsFactors = F)
    blastfile$blast_percentage[is.na(blastfile$blast_percentage)==T]=""
    blastfile$blast_percentage=as.numeric(blastfile$blast_percentage)
    blastfile$blast_percentage=round(blastfile$blast_percentage, 2)
    
    blastfile$blast_percentage = replacena(genetictophysical(chromosome = list.of.chromosome.assignments.for.lgs[which(lglist == lgtomakehistlistfor)], 
                                                                                                 listofmarkers = grab.list.of.markers(lgtomakehistlistfor), #see recombination_diagnostics.R for this function 
                                                                                                 markerformat = ".", threshold = 70))
    blastfile$blast_percentage = round(blastfile$blast_percentage, 2)
    blastfile$blast_percentage = as.numeric(blastfile$blast_percentage)
    blastfile$blast_per_filled=blastfile$blast_percentage
    #blastfile$blast_per_filled=replacena(blastfile$blast_per_filled)
    #blastfile$blast_per_filled=replacena(replacena(round(genetictophysical("3B", blastfile$marker1, "."),2))) #for LG15 (this is wrongly assigned as 3D, in reality it is 3B)
    chromosomeforhere=blastfile[1, 17]
    
    print("making histlists")
    
    q1hist=makehistlist(lgtomakehistlistfor, genotypelists[[1]], F)
    q2hist=makehistlist(lgtomakehistlistfor, genotypelists[[2]], F)
    q3hist=makehistlist(lgtomakehistlistfor, genotypelists[[3]], F)
    q4hist=makehistlist(lgtomakehistlistfor, genotypelists[[4]], F)
    print("made histlists")
    
    q1histtable=as.data.frame(table(as.vector(as.numeric(q1hist))))
    q2histtable=as.data.frame(table(as.vector(as.numeric(q2hist))))
    q3histtable=as.data.frame(table(as.vector(as.numeric(q3hist))))
    q4histtable=as.data.frame(table(as.vector(as.numeric(q4hist))))
    
    q1histtable$temp=10
    q2histtable$temp=26
    q3histtable$temp=14
    q4histtable$temp=28
    
    allhisttable=rbind(q1histtable, q2histtable, q3histtable, q4histtable)
    colnames(allhisttable)=c("x","value","temp")
    allhisttable$temp=as.factor(allhisttable$temp)
    allhisttable$x=as.numeric(as.character(allhisttable$x))
    allhisttable$phys=blastfile$blast_per_filled[allhisttable$x]
    allhisttable$phys=as.numeric(as.character(allhisttable$phys))
    
    temperatures = c(10, 14, 26, 28)
    phys.values = lapply(temperatures, function(x){ 
        as.numeric(allhisttable[allhisttable$temp == x, ]$phys)
        }
    )
    
    grab.unique = function(x){
        b = c(1, 2, 3, 4)
        b = b[-which(b == x)]
        phys.values[[x]][which(!phys.values[[x]] %in% union(union(phys.values[[b[1]]], phys.values[[b[2]]]), phys.values[[b[3]]]))]
    }
    
    grab.shared.unique = function(x, y){
        b = 1:4
        b = b[-which(b == x | b == y)]
        g = as.numeric(phys.values[[x]]); g2 = as.numeric(phys.values[[y]])
        intersect.x.y = intersect(g, g2)
        intersect.x.y[which(!intersect.x.y %in% as.numeric(union(phys.values[[b[1]]], phys.values[[b[[2]]]])))]
        
    }
    
    list.of.shared.unique = grab.shared.unique(3, 4)
    
    
    
    list.of.unique = lapply(1:4, grab.unique)
    phys.unique.df = newdf(c("unique.val", "temp"))
    
    conv.i.to.temp = function(x){
        if(x == 1) y = 10
        if(x == 2) y = 14
        if(x == 3) y = 26
        if(x == 4) y = 28
        return(y)
    }
    
    for(i in 1:4){
        for(p in seq_along(list.of.unique[[i]])){
            phys.unique.df[nrow(phys.unique.df), 1] = list.of.unique[[i]][[p]]
            phys.unique.df[nrow(phys.unique.df), 2] = conv.i.to.temp(i)
            phys.unique.df = add_row(phys.unique.df)
        }
    }
    phys.unique.df = phys.unique.df[-nrow(phys.unique.df), ]
    phys.unique.df$count = ""
    for(i in 1:nrow(phys.unique.df)){
        phys.unique.df$count[i] = allhisttable$value[which(allhisttable$phys == phys.unique.df$unique.val[i] & allhisttable$temp == phys.unique.df$temp[i])]
    }
    phys.unique.df$unique.val = as.numeric(phys.unique.df$unique.val)
    phys.unique.df = phys.unique.df[sort(phys.unique.df$unique.val, index.return = T)[[2]],]
    phys.unique.df[] = lapply(phys.unique.df, as.numeric)
    
    phys.unique.df$unique.val = as.factor(phys.unique.df$unique.val)
    
    pdf(file=paste("plots/unique.values/uniquevalues.lg", lgtomakehistlistfor, ".pdf", sep = ""),
            width = 10, height = 10)
    print(ggplot(data = phys.unique.df) + geom_bar(aes(x = unique.val, y = count, group = temp), stat = 'identity') +
        facet_wrap(~temp, nc = 1))
    dev.off()
    
    
    print("writing pdf")
    pdf(file=paste("plots/multilineplots/recombination_analysis/stacked_bar/phys/customprinttest_lg", lgtomakehistlistfor,".pdf", sep=""), width=15, height=10)
    # windows()
    # graph=ggplot(data=allhisttable, aes(x=phys, y=value, group=temp)) #CHANGE FROM x=x to x=phys ACCORDINGLY - ONE PLOTS THE PHYSICAL DISTANCE ACCROSS THE X-AXIS
    # print(graph + geom_line(aes(color=temp, linetype=temp), size=1.2) + 
    #     theme_bw() + ggtitle(paste("Linkage group", lgtomakehistlistfor, "Chromosome", chromosomeforhere)) + 
    #     scale_linetype_manual(name = "Temperature ((deg)C)", values=c("twodash","dotted","solid","dashed")) + 
    #     scale_color_brewer(name = "Temperature ((deg)C)", palette = "Spectral") +
    #     scale_x_continuous(breaks=seq(0, 100, 5), labels=seq(0, 100, 5)) +
    #     theme(axis.text.x=element_text(angle=90, hjust=1)) +
    #     theme_classic() +
    #     guides(fill=guide_legend(title="Temperature ((deg)C)")) +
    #     xlab("Physical distance (%)") + 
    #     ylab ("# Recombination events"))
    
    print(ggplot(data = allhisttable) + geom_bar(aes(x = x, y = value, group = temp), stat = 'identity') +
                    facet_wrap(~temp, nc = 1) +
            xlab("Marker position in linkage group") +
            ylab ("# Recombination events"))
    
    dev.off()
    print("done writing pdf")
}




#make 4x4 histogram

list_of_lg=as.numeric(unique(genotypelists[[1]]$lg))[-36]
for(x in list_of_lg){
    q1hist=makehistlist(x, genotypelists[[1]])
    q2hist=makehistlist(x, genotypelists[[2]])
    q3hist=makehistlist(x, genotypelists[[3]])
    q4hist=makehistlist(x, genotypelists[[4]])
    
    q1histtable=as.data.frame(table(as.vector(as.numeric(q1hist))))
    q2histtable=as.data.frame(table(as.vector(as.numeric(q2hist))))
    q3histtable=as.data.frame(table(as.vector(as.numeric(q3hist))))
    q4histtable=as.data.frame(table(as.vector(as.numeric(q4hist))))
    
    q1histtable$quad=1
    q2histtable$quad=2
    q3histtable$quad=3
    q4histtable$quad=4
    
    allhisttable=rbind(q1histtable, q2histtable, q3histtable, q4histtable)
    allhisttable$quad=as.factor(allhisttable$quad)
    colnames(allhisttable)=c("x","value","quad")
    graph=ggplot(data=allhisttable, aes(x=x, y=value, group=quad))
    # windows()
    png(filename=paste("multilineplots/recombination_analysis/multiline_lg_",x,".png",sep=""), width=2000,height=2000)
    print(graph + geom_line(aes(color=quad, linetype=quad), size=1.2) + theme_bw())
    dev.off()
    
    
    listofbiggest=list()
    biggest1=max(as.numeric(q1hist))
    biggest2=max(as.numeric(q2hist))
    biggest3=max(as.numeric(q3hist))
    biggest4=max(as.numeric(q4hist))
    truebiggest=max(as.numeric(c(biggest1, biggest2, biggest3, biggest4)))

    
    #png(filename=paste("histogram_recombination_counts/",quadrant.to.analyse.name, "/", "hist_",quadrant.to.analyse.name,"_lg_",x,".png",sep=""))
    png(filename=paste("histogram_recombination_counts/all_quadrants/newhist_lg_",x,".png",sep=""), width=1000, height=1000)
    par(mfrow=c(2,2))
    hist(as.numeric(q1hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    hist(as.numeric(q2hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    hist(as.numeric(q3hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    hist(as.numeric(q4hist), ylim=c(0,(truebiggest+10)), xlim=c(0,120), axes=F)
    axis(side=1, at=seq(0,120,10), labels=seq(0,120,10))
    axis(side=2, at=seq(0,truebiggest+10, 5), labels=seq(0,truebiggest+10, 5))
    dev.off()
}

?hist
?axis



#     ____________________________________________________________________________
#     PARSE GENOTYPE INFO. FOR 1 INDIVIDUAL (useful for debugging loop above) ####


#setup dataframe - for counting the number of consecutive genotypes that are the same
genotypelist=as.data.frame(matrix(nrow=1, ncol=5))
genotypelist[is.na(genotypelist)]=""
colnames(genotypelist)=c("genotype","number_consec","lg", "individual", "individual_id")

#setup init. variables
individual.index.number=5
geno_count=1
first=quadrant.to.analyse[individual.index.number,3]
firstlg=quadrant.to.analyse[1,3]
genotypelist$genotype[geno_count]=first
genotypelist$number_consec[geno_count]=1
genotypelist$lg[geno_count]=quadrant.to.analyse[1,3]
genotypelist$individual[geno_count]=quadrant.to.analyse[individual.index.number,2]
genotypelist$individual_id[geno_count]=individual.index.number

listof.lg=unique(as.numeric(quadrant.to.analyse[1, 3:ncol(quadrant.to.analyse)]))

for(i in 4:ncol(quadrant.to.analyse)){
    if(quadrant.to.analyse[individual.index.number,i] == first & quadrant.to.analyse[1,i] == firstlg){
        genotypelist$number_consec[geno_count]=(as.numeric(genotypelist$number_consec[geno_count])+1)
    } else {
        firstlg=quadrant.to.analyse[1,i]
        genotypelist=add_row(genotypelist)
        genotypelist[nrow(genotypelist),]=""
        
        geno_count=geno_count+1
        first=quadrant.to.analyse[individual.index.number,i]
        
        genotypelist$genotype[geno_count]=first
        genotypelist$number_consec[geno_count]=1
        genotypelist$lg[geno_count]=quadrant.to.analyse[1,i]
        genotypelist$individual[geno_count]=quadrant.to.analyse[individual.index.number,2]
        genotypelist$individual_id[geno_count]=individual.index.number
    }
}

listof.lg=unique(as.numeric(genotypelist$lg))
listof.lg
