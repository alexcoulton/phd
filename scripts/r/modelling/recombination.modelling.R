#### INITIAL DATA ####

corestouse = 4

source("rotation1scripts_v4/scripts/r/functions.R")

#modelling recombination / segregation distortion

#take avalon cadenza 5A recombination map to use for probabilities in model
load("rotation1scripts_v4/saved.objects/c651.hist")

# c651.hist = as.numeric(c651.hist)

# c651.table = table(factor(c651.hist, levels = 1:185)) #there are 185 markers in c.x.a 5A

# tot.events = sum(c651.table)
# probability.vec = as.numeric(c651.table) / tot.events
# probability.vec = probability.vec[-length(probability.vec)]

# probability.vec = pxwmaps2$cM
# probability.vec = diff(probability.vec)

# probability.vec = probability.vec / sum(probability.vec)
# plot(probability.vec)

# prob.cent = (length(probability.vec) / 2)
#make central markers very low prob
# probability.vec[(prob.cent - 25):(prob.cent + 25)] = rnorm(51, 0.006, sd = 0.001)

#scale back up to one
# probability.vec = probability.vec * (1 / sum(probability.vec))
# plot(probability.vec)


#--------------- make low prob center prob vec

# RANGE1 = 4
# POWER1 = 10
# probability.vec2 = probability.vec

# prob.cent = (length(probability.vec2) / 2)
#make central markers very low prob
# probability.vec2[(prob.cent - RANGE1):(prob.cent + RANGE1)] = (probability.vec2[(prob.cent - RANGE1):(prob.cent + RANGE1)] * (c(seq(1, RANGE1), 10, seq(RANGE1, 1)) * POWER1))

#scale back up to one
# probability.vec2 = probability.vec2 * (1 / sum(probability.vec2))



#--------------- make medium prob center prob vec
# probability.vec3 = probability.vec

# prob.cent = (length(probability.vec3) / 2)
#make central markers very low prob
# probability.vec3[(prob.cent - RANGE1):(prob.cent + RANGE1)] = (probability.vec3[(prob.cent - RANGE1):(prob.cent + RANGE1)] * (c(seq(1, RANGE1), 10, seq(RANGE1, 1)) * (POWER1 * 2)))

#scale back up to one
# probability.vec3 = probability.vec3 * (1 / sum(probability.vec3))



# plot(probability.vec)
# plot(probability.vec2)
# plot(probability.vec3)

# prob.vec.equal = rep((1 / 76), 76)

# recomb.qtl.dat.sim = gen.data(77, 1000, prob.vec = prob.vec.equal)
# r2 = gen.next.generation(recomb.qtl.dat.sim, 77, prob.vec = prob.vec.equal)
# r3 = gen.next.generation(r2, 77, prob.vec = prob.vec.equal)
# r4 = gen.next.generation(r3, 77, prob.vec = prob.vec.equal)
# r5 = gen.next.generation(r4, 77, prob.vec = prob.vec.equal)
# r6 = gen.next.generation(r5, 77, prob.vec = prob.vec.equal)

# evaluate.hets(r6)

# conv.sim2 = function(x){
#     ro3 = convert.sim.to.rqtl.ske.format(x)
    
#     colnames(ro3)[3:ncol(ro3)] = pxwmaps2$marker
#     ro3[1, 3:ncol(ro3)] = "6B"
#     ro3[1, 1] = NA
#     colnames(ro3)[2] = "probeset_id"    
#     ro3
# }

# simf2 = conv.sim2(recomb.qtl.dat.sim)
# simf3 = conv.sim2(r2)

# evaluate.hets(simf3)

# simf2.test = qtlanalysis.all(simf2, pxwmaps2, "RIL7", geno.prob.step = 0, fgen = 2, marker.format = "-", qtl.algorithm = "em")
# simf3.test = qtlanalysis.all(simf3, pxwmaps2, "RIL7", geno.prob.step = 0, fgen = 3, marker.format = "-", qtl.algorithm = "em")

# simf2.test[[1]][[1]][[2]] #plot
# simf3.test[[1]][[1]][[2]] #plot


# simf2.test[[1]][[1]][[5]]
# simf3.test[[1]][[1]][[5]]

# lapply(simf3.test[[1]][[1]][[6]], function(x){
#     x$effect.means$Means
# })



# ro3.qtl2.wo.count1 = qtlanalysis.all(ro3, pxwmaps2, "RIL7", geno.prob.step = 0, fgen = 7, marker.format = "-", qtl.algorithm = "em")
# ro3.qtl3.wo.count1 = qtlanalysis.all(ro3, pxwmaps2, "RIL7", geno.prob.step = 0, fgen = "riself", marker.format = "-", qtl.algorithm = "em")

# ro3.qtl2.wo.test = qtlanalysis.all(ro3, pxwmaps2, "RIL7", geno.prob.step = 0, fgen = 7, marker.format = "-", qtl.algorithm = "em")

#crosstype bcfst, w/ dist qtl
# ro3.qtl2[[1]][[1]][[2]]
# ro3.qtl2[[1]][[2]][[2]]
#crosstype riself, w/ dist qtl
# ro3.qtl3[[1]][[1]][[2]]
# ro3.qtl3[[1]][[2]][[2]]

#crosstype bcfst, no dist qtl
# ro3.qtl2.wo[[1]][[1]][[2]]
# ro3.qtl2.wo[[1]][[2]][[2]]
#crosstype riself, no dist qtl
# ro3.qtl3.wo[[1]][[1]][[2]]
# ro3.qtl3.wo[[1]][[2]][[2]]

# v(detect.recombination(ro3))

# ro3.qtl2.wo.count1[[1]][[1]]
# ro3.qtl3.wo.count1[[1]][[1]]



# load("rotation1scripts_v4/saved.objects/recombination.profiles")

#can't just use recombination profiles from individuals as these are made up of two gametes
#hence are doubled in the amount of events they contain.
#here i will trim down the recombination profiles by removing events with the smallest difference in position
#from adjacent events
# recomb.diffs = lapply(recombination.profiles, function(x){
#     diffs = diff(x)
#     diffs.sort = sort(diffs, index.return = T)$ix
#     # browser()
#     g = x
#     if(length(x) > 1){
#         g = x[-diffs.sort[1:round((length(x) / 2))]]
#     }
#     g
# })

# mean(unlist(lapply(recomb.diffs, function(x) length(x))))


#### FUNCTIONS ####

gen.hap = function(no.markers, haplotype.one, haplotype.two, selection, selection.pos, prob.vec, recombination.profiles){
    #args:
    #recomb.param: integer from 1:9 indicating the upper range of possible recombination position
    #no.markers: integer, number of markers to generate haplotype for
    #haplotype.one: character vector, series of "A" and "B" representing first haplotype
    #haplotype.two: character vector, series of "A" and "B" representing second haplotype
    #selection: integer specifying the denominator of the selection factor (1/selection gametes are killed)
    #selection.pos: integer specifying the marker at which to apply selection, if missing, defaults to central marker
    #prob.vec: numeric vector of probabilities for recombination distribution
    #recombination.profiles: list of numeric vectors, each vector containing positions of recombination events in an individuals of a mapping population
    
    if(missing(haplotype.one)) haplotype.one = rep("A", no.markers)
    if(missing(haplotype.two)) haplotype.two = rep("B", no.markers)
    if(missing(selection)) selection = F
    if(missing(selection.pos)) selection.pos = round(no.markers / 2)
    if(missing(prob.vec)) prob.vec = F
    if(missing(recombination.profiles)) recombination.profiles = F
    #recombination can occur anywhere
    recomb.param = no.markers - 1
    
    
    
    make.haplotype = function(recomb.param, no.markers, haplotype.one, haplotype.two, prob.vec, recombination.profiles){
        
        if(missing(haplotype.one)) haplotype.one = rep("A", no.markers)
        if(missing(haplotype.two)) haplotype.two = rep("B", no.markers)
        if(missing(prob.vec)) prob.vec = F
        if(missing(recombination.profiles)) recombination.profiles = F
 
        if(length(prob.vec) == 1){ #if we are not using a probability vector
            recomb.pos = sample(recomb.param, 1)    
            all.recomb.pos = recomb.pos
        } else {
            #get number of recombination events from poisson distribution with mean 1.3
            #number.events1 = rpois(1, 1.3) #- this includes zeros, not good as can't have zero recombination events
            
            number.events1 = sample((sapply(1:1000, function(x) rpois(1, 0.5)) + 1), 1) #mean of this distribution is ~ 1.4
            #get positions of recombination events
            # qtl.position = 30
            # if(haplotype.one[qtl.position] == "A" & haplotype.two[qtl.position] == "A"){
            #     print("jaja")
            #     all.recomb.pos = sample(recomb.param, number.events1, prob = probability.vec)
            # }
            # if(haplotype.one[qtl.position] == "A" & haplotype.two[qtl.position] == "B"){
            #     print("jaja2")
            #     all.recomb.pos = sample(recomb.param, number.events1, prob = probability.vec2)
            # }
            # if(haplotype.one[qtl.position] == "B" & haplotype.two[qtl.position] == "A"){
            #     print("jaja3")
            #     all.recomb.pos = sample(recomb.param, number.events1, prob = probability.vec2)
            # }
            # if(haplotype.one[qtl.position] == "B" & haplotype.two[qtl.position] == "B"){
            #     print("jaja4")
            #     all.recomb.pos = sample(recomb.param, number.events1, prob = probability.vec3)
            # }

            all.recomb.pos = sample(recomb.param, number.events1, prob = prob.vec)
            
            #old method of selecting 2nd event
            # recomb.pos = sample(recomb.param, 1, prob = prob.vec) #recomb w. probability vector
            # all.recomb.pos = recomb.pos
            # second.event.dice = sample(100, 1)
            
            # if(second.event.dice >= 50){
            #     prob.vec2 = prob.vec
            #     #set the probability of the second event being in the same position as the first event to zero
            #     prob.vec2[recomb.pos] = 0
            #     #lower the probabilities of flanking markers from the first event
            # 
            #     if((recomb.pos + 10) > recomb.param){
            #         forward.flank = (recomb.param - recomb.pos)
            #     } else {
            #         forward.flank = 10
            #     }
            # 
            #     if((recomb.pos - 10) < 1){
            #         backward.flank = (recomb.pos - 1)
            #     } else {
            #         backward.flank = 10
            #     }
            # 
            #     #reduce probabilities of 10 flanking markers (or all flanking markers if close to the end) of previous recombination postition
            #     if(forward.flank != 0){
            #         prob.vec2[(recomb.pos + 1):(recomb.pos + forward.flank)] = prob.vec2[(recomb.pos + 1):(recomb.pos + forward.flank)] / length(prob.vec2[(recomb.pos + 1):(recomb.pos + forward.flank)]):1
            #     }
            # 
            #     if(backward.flank != 0){
            #         prob.vec2[(recomb.pos - 1):(recomb.pos - backward.flank)] = prob.vec2[(recomb.pos - 1):(recomb.pos - backward.flank)] / length(prob.vec2[(recomb.pos + 1):(recomb.pos + backward.flank)]):1
            #     }
            # 
            # 
            #     #scale probabilities so they sum to 1 again
            #     prob.vec2 = prob.vec2 * (1 / sum(prob.vec2))
            # 
            #     #get second recombination position
            #     recomb.pos2 = sample(recomb.param, 1, prob = prob.vec2)
            # 
            # 
            #     #combine both recombination positions and sort them so smallest happens first
            #     all.recomb.pos = sort(c(recomb.pos, recomb.pos2))
            # }
        }
        
        perform.recombination = function(hap1, hap2, recomb.pos){
            recombined.hap1 = c(hap1[1:recomb.pos], hap2[(recomb.pos + 1): no.markers])
            recombined.hap2 = c(hap2[1:recomb.pos], hap1[(recomb.pos + 1): no.markers])
            list(recombined.hap1, recombined.hap2)
        }
        
        
        #if we are using recombination profiles rather than a probability vector
        if(length(recombination.profiles) > 1){
            #randomly pick a recombination profile to use
            profile.to.use = recombination.profiles[[sample(length(recombination.profiles), 1)]]
            
            for(i in 1:length(profile.to.use)){
                recomb1 = perform.recombination(haplotype.one, haplotype.two, profile.to.use[i])
                #update haplotypes (input to perform.recombination) to new recombined haplotypes
                haplotype.one = recomb1[[1]]
                haplotype.two = recomb1[[2]]
            }
            
            recombined.hap1 = haplotype.one
            recombined.hap2 = haplotype.two
            
        } else { #if we are using a probability vector
            
            # recombined.hap1 = c(haplotype.one[1:recomb.pos], haplotype.two[(recomb.pos + 1): no.markers])
            # recombined.hap2 = c(haplotype.two[1:recomb.pos], haplotype.one[(recomb.pos + 1): no.markers])
            
            profile.to.use = all.recomb.pos
            
            if(length(all.recomb.pos) > 0){
                for(i in 1:length(profile.to.use)){
                    recomb1 = perform.recombination(haplotype.one, haplotype.two, profile.to.use[i])
                    #update haplotypes (input to perform.recombination) to new recombined haplotypes
                    haplotype.one = recomb1[[1]]
                    haplotype.two = recomb1[[2]]
                }
            }
            
            
            recombined.hap1 = haplotype.one
            recombined.hap2 = haplotype.two
        }
        
 
        #generate gametes
        
        
        # browser()
        
        both.haps = list(recombined.hap1, recombined.hap2)
        #pick one of the gametes at random
        selected.hap = both.haps[[sample(2, 1)]]
        
        
        return(selected.hap)
    }
    
    # if(length(prob.vec) != 224) browser()
    # if(recomb.param != 224) browser()
    
    haplotype = make.haplotype(recomb.param, no.markers, haplotype.one, haplotype.two, prob.vec, recombination.profiles)
    
    
    if(selection != F){
        #implement selection at a particular marker
        dice = sample(selection, 1)
        
        #make sure that selection is possible in this line - e.g. the line is not homozygous
        if(haplotype.one[selection.pos] != haplotype.two[selection.pos]){
            #1/2 chance of deleting this gamete. make new haplotype until there is not an "A" at marker 2
            while(dice == 1 & haplotype[selection.pos] == "A"){
                haplotype = make.haplotype(recomb.param, no.markers, haplotype.one, haplotype.two, prob.vec)
                # browser()
            }    
        }
        
    }
    
    
    
    
    return(haplotype)
}

gen.zygote = function(no.markers, hap1, hap2, selection, selection.pos, prob.vec, recombination.profiles){
    #if haplotypes not provided, generate some
    if(missing(selection)) selection = F
    if(missing(selection.pos)) selection.pos = round(no.markers / 2)
    if(missing(prob.vec)) prob.vec = F
    if(missing(recombination.profiles)) recombination.profiles = F
    if(missing(hap1)) hap1 = gen.hap(no.markers, selection = selection, selection.pos = selection.pos, prob.vec = prob.vec, recombination.profiles = recombination.profiles)
    if(missing(hap2)) hap2 = gen.hap(no.markers, selection = selection, selection.pos = selection.pos, prob.vec = prob.vec, recombination.profiles = recombination.profiles)
    
    
    both.hap = data.frame(hap1, hap2)
    
    hap3 = ""
    for(i in 1:nrow(both.hap)){
        if(both.hap[i, 1] == "A" & both.hap[i, 2] == "A") hap3 = c(hap3, "A")
        if(both.hap[i, 1] == "B" & both.hap[i, 2] == "B") hap3 = c(hap3, "B")
        if(both.hap[i, 1] == "B" & both.hap[i, 2] == "A") hap3 = c(hap3, "H")
        if(both.hap[i, 1] == "A" & both.hap[i, 2] == "B") hap3 = c(hap3, "H")
    }
    
    hap3 = hap3[2:(no.markers+1)]
    # attr(hap3, "hap1") = hap1
    # attr(hap3, "hap2") = hap2
    
    
    return(hap3)
}

library(parallel)
gen.data = function(re, no.individuals, selection, selection.pos, prob.vec, recombination.profiles){
    #generate genotyping data
    #args:
    # re - an integer, the number of markers in the dataset
    # no.individuals - integer, the number of individuals in the dataset
    # selection - integer, strength of selection
    # selection.pos - integer, position of selection strength
    # prob.vec - numeric vector of recombination probabilities
    # recombination.profiles - list of numeric vectors, possible recombination positions in gametes
    
    if(missing(selection)) selection = F
    if(missing(selection.pos)) selection.pos = round(re / 2)
    if(missing(prob.vec)) prob.vec = F
    if(missing(recombination.profiles)) recombination.profiles = F
    #simulate recombination data for 300 individuals
    if(missing(no.individuals)) no.individuals = 300
    
    # print(prob.vec)
    
    recomb.data = lapply(1:no.individuals, function(x) gen.zygote(re, selection = selection, selection.pos = selection.pos, prob.vec = prob.vec, recombination.profiles = recombination.profiles))
    recomb.data = as.data.frame(t(as.data.frame(recomb.data)))
    recomb.data = reset.rownames(recomb.data)
    convert.to.character.data.frame(recomb.data)
}

generate.haplotypes.from.zygote = function(zygote){
    g = zygote
    g1 = g
    g2 = g
    
    for(i in 1:length(g)){
        if(g[i] == "A"){
            g1[i] = "A"
            g2[i] = "A"
        } 
        if(g[i] == "B"){
            g1[i] = "B"
            g2[i] = "B"
        }
        if(g[i] == "H"){
            g1[i] = "A"
            g2[i] = "B"
        }
    }
    
    return(list(as.character(g1), as.character(g2)))
}

gen.next.generation = function(recombination.data, no.markers, selection, selection.pos, prob.vec, recombination.profiles){
    if(missing(selection)) selection = F
    if(missing(selection.pos)) selection.pos = round(no.markers / 2)
    if(missing(prob.vec)) prob.vec = F
    if(missing(recombination.profiles)) recombination.profiles = F
    # browser()
    
    recombination.data = as.data.frame(t(recombination.data), stringsAsFactors = F)
    # browser()
    # print(ncol(recombination.data))
    c12 = make.counter()
    new.recomb = lapply(recombination.data, function(x){
        new.haps = generate.haplotypes.from.zygote(x)
        
        hapf3_1 = gen.hap(no.markers, new.haps[[1]], new.haps[[2]], selection, selection.pos, prob.vec, recombination.profiles)
        hapf3_2 = gen.hap(no.markers, new.haps[[1]], new.haps[[2]], selection, selection.pos, prob.vec, recombination.profiles)
        
        gen.zygote(no.markers, hapf3_1, hapf3_2, selection, selection.pos)
        # print(p("whooo", c12()))    
    })#, mc.cores = corestouse)
    
    # browser()
    
    
    
    new.recomb = as.data.frame(t(as.data.frame(new.recomb)))
    
    
    return(new.recomb)
}



perform.1000.simulations.fgen.keep = function(pop.size, selection.strength, selection.pos, number.of.markers, recombination.profile, corestouse, fgen){
    #keeps every filial generation
    load("rotation1scripts_v4/saved.objects/recombination.profiles1a")
    
    if(missing(corestouse)) corestouse = 35
    if(missing(selection.pos)) selection.pos = 200
    if(missing(number.of.markers)) number.of.markers = 225
    if(missing(recombination.profile)) recombination.profile = F
    if(missing(selection.strength)) selection.strength = F
    
    if(fgen < 3){
        return("fgen must be at least 3")
    }
    
    if(selection.strength == F){
        if(class(recombination.profile) != "list"){
            simulations1 = mclapply(1:1000, function(x){
                g = list(gen.data(number.of.markers, pop.size))
                
                all.gens = lapply(2:fgen, function(x){
                    g <<- c(g, list(gen.next.generation(g[[length(g)]], number.of.markers)))
                    return()
                })
                
                g
                
            }, mc.cores = corestouse)     
        } else {
            simulations1 = mclapply(1:1000, function(x){
                g = list(gen.data(number.of.markers, pop.size, recombination.profiles = recombination.profile))
                
                all.gens = lapply(2:fgen, function(x){
                    g <<- c(g, list(gen.next.generation(g[[length(g)]], number.of.markers, recombination.profiles = recombination.profile)))
                    return()
                })
                
                g
            }, mc.cores = corestouse) 
        }
        
    } else {
        if(class(recombination.profile) != "list"){
            simulations1 = mclapply(1:1000, function(x){
                g = list(gen.data(number.of.markers, pop.size, selection = selection.strength, selection.pos = selection.pos))
                
                all.gens = lapply(2:fgen, function(x){
                    g <<- c(g, list(gen.next.generation(g[[length(g)]], number.of.markers, selection = selection.strength, 
                                                                        selection.pos = selection.pos)))
                    return()
                })
                
                g
            }, mc.cores = corestouse)
        } else {
            simulations1 = mclapply(1:1000, function(x){
                g = list(gen.data(number.of.markers, pop.size, recombination.profiles = recombination.profile, selection = selection.strength, selection.pos = selection.pos))
                
                all.gens = lapply(2:fgen, function(x){
                    g <<- c(g, list(gen.next.generation(g[[length(g)]], number.of.markers, selection = selection.strength, 
                                                                        selection.pos = selection.pos, recombination.profiles = recombination.profile)))
                    return()
                })
                
                
                g
            }, mc.cores = corestouse)
        }
        
    }
    
    
    
    simulations1
    
}

perform.1000.simulations.fgen = function(pop.size, selection.strength, selection.pos, number.of.markers, recombination.profile, corestouse, fgen, number.simulations){
    load("rotation1scripts_v4/saved.objects/recombination.profiles1a")
    
    if(missing(corestouse)) corestouse = 35
    if(missing(selection.pos)) selection.pos = 200
    if(missing(number.of.markers)) number.of.markers = 225
    if(missing(recombination.profile)) recombination.profile = F
    if(missing(selection.strength)) selection.strength = F
    if(missing(number.simulations)) number.simulations = 1000
    
    if(fgen < 3){
        return("fgen must be at least 3")
    }
    
    if(selection.strength == F){
        if(class(recombination.profile) != "list"){
            simulations1 = mclapply(1:number.simulations, function(x){
                g = gen.data(number.of.markers, pop.size)
                
                all.gens = lapply(3:fgen, function(x){
                    g <<- gen.next.generation(g, number.of.markers)
                    return()
                })
                
                g
                
            }, mc.cores = corestouse)     
        } else {
            simulations1 = mclapply(1:number.simulations, function(x){
                g = gen.data(number.of.markers, pop.size, recombination.profiles = recombination.profile)
                
                all.gens = lapply(3:fgen, function(x){
                    g <<- gen.next.generation(g, number.of.markers, recombination.profiles = recombination.profile)
                    return()
                })
                
                g
            }, mc.cores = corestouse) 
        }
        
    } else {
        if(class(recombination.profile) != "list"){
            simulations1 = mclapply(1:number.simulations, function(x){
                g = gen.data(number.of.markers, pop.size, selection = selection.strength, selection.pos = selection.pos)
                
                all.gens = lapply(3:fgen, function(x){
                    g <<- gen.next.generation(g, number.of.markers, selection = selection.strength, 
                                                                        selection.pos = selection.pos)
                    return()
                })
                
                g
            }, mc.cores = corestouse)
        } else {
            simulations1 = mclapply(1:number.simulations, function(x){
                g = gen.data(number.of.markers, pop.size, recombination.profiles = recombination.profile, selection = selection.strength, selection.pos = selection.pos)
                
                all.gens = lapply(3:fgen, function(x){
                    g <<- gen.next.generation(g, number.of.markers, selection = selection.strength, 
                                                                        selection.pos = selection.pos, recombination.profiles = recombination.profile)
                    return()
                })
                
                
                g
            }, mc.cores = corestouse)
        }
        
    }
    
    
    
    simulations1
    
}





perform.1000.simulations = function(pop.size, selection.strength, selection.pos, number.of.markers, recombination.profile, corestouse){
    load("rotation1scripts_v4/saved.objects/recombination.profiles1a")
    
    if(missing(corestouse)) corestouse = 35
    if(missing(selection.pos)) selection.pos = 200
    if(missing(number.of.markers)) number.of.markers = 225
    if(missing(recombination.profile)) recombination.profile = F
    if(missing(selection.strength)) selection.strength = F
    
    
    
    if(selection.strength == F){
        if(class(recombination.profile) != "list"){
            simulations1 = mclapply(1:1000, function(x){
                g = gen.data(number.of.markers, pop.size)
                g
            }, mc.cores = corestouse)     
        } else {
            simulations1 = mclapply(1:1000, function(x){
                g = gen.data(number.of.markers, pop.size, recombination.profiles = recombination.profile)
                g
            }, mc.cores = corestouse) 
        }
        
    } else {
        if(class(recombination.profile) != "list"){
            simulations1 = mclapply(1:1000, function(x){
                g = gen.data(number.of.markers, pop.size, selection = selection.strength, selection.pos = selection.pos)
                g
            }, mc.cores = corestouse)
        } else {
            simulations1 = mclapply(1:1000, function(x){
                g = gen.data(number.of.markers, pop.size, recombination.profiles = recombination.profile, selection = selection.strength, selection.pos = selection.pos)
                g
            }, mc.cores = corestouse)
        }
        
    }
    
    
    
    simulations1
    
}



test.data = gen.data(5, 5)
g = gen.next.generation(test.data, 5)
g = gen.next.generation(g, 5)


#### ANALYSIS FUNCTIONS 1 ####


convert.recomb.to.seg = function(recomb.data){
    #args: 
    # recomb.data - dataframe containing genotyping information
    seg.data = unlist(lapply(recomb.data, function(x){
        # browser()
        a = length(which(x == "A"))
        b = length(which(x == "B"))
        a / (a + b)
    })) 
    seg.data
}

convert.recomb.to.seg2 = function(recomb.data, p.value1){
    #args:
    # recomb.data - a dataframe containing genotype information
    # p.value1 - boolean flag indicating whether to only return p-values from chi-square tests
    
    if(missing(p.value1)) p.value1 = F
    
    seg.data = lapply(recomb.data, function(x){
        # browser()
        a = length(which(x == "A"))
        b = length(which(x == "B"))
        chisq.test(c(a, b))
    })
 
    if(p.value1 == T){
        return(unlist(lapply(seg.data, function(x) x[[3]])))
    } else {
        return(seg.data)
    }
    
}

convert.recomb.to.seg.f2 = function(recomb.data, p.value1){
    #args:
    # recomb.data - a dataframe containing genotype information
    # p.value1 - boolean flag indicating whether to only return p-values from chi-square tests
    
    if(missing(p.value1)) p.value1 = F
    
    seg.data = lapply(recomb.data, function(x){
        # browser()
        a = length(which(x == "A"))
        h = length(which(x == "H"))
        b = length(which(x == "B"))
        chisq.test(c(a, h, b), p = c(1/4, 1/2, 1/4))
    })
    
    if(p.value1 == T){
        return(unlist(lapply(seg.data, function(x) x[[3]])))
    } else {
        return(seg.data)
    }
}


# lapply(selected1, function(x){
#     # browser()
#     a = length(which(x == "A"))
#     b = length(which(x == "B"))
#     chisq.test(c(a, b))
# })



check.for.distorted.markers = function(recomb.data, thresholdlevel, F2){
    #returns the coordinates of markers that show significant seg. dist. according to chi-square test
    #args:
    # recomb.data - a dataframe containing genotyping data
    # thresholdlevel - an integer specifying the level of significance, defaults to 0.05
    if(missing(thresholdlevel)) thresholdlevel = 0.05
    if(missing(F2)) F2 = F
    
    if(F2 == T){
        g = convert.recomb.to.seg.f2(recomb.data)
    } else {
        g = convert.recomb.to.seg2(recomb.data)
    }
    
    g2 = which(unname(unlist(lapply(g, function(x) x[3]))) < thresholdlevel)
    return(g2)
}

check.for.distorted.markers.fdr = function(recomb.data, thresholdlevel, method1, F2){
    #returns the coordinates of markers that show significant seg. dist. according to chi-square test
    #args:
    # recomb.data - a dataframe containing genotyping data
    # thresholdlevel - an integer specifying the level of significance, defaults to 0.05
    if(missing(thresholdlevel)) thresholdlevel = 0.05
    if(missing(method1)) method1 = "BH"
    if(missing(F2)) F2 = F
    
    if(F2 == T){
        g = convert.recomb.to.seg.f2(recomb.data)
    } else {
        g = convert.recomb.to.seg2(recomb.data)
    }
    
    g2 = which(p.adjust(unname(unlist(lapply(g, function(x) x[3]))), method1) < thresholdlevel)
    return(g2)
}

prepare.gg.data = function(recomb.data){
    new.data = data.frame(convert.recomb.to.seg(recomb.data), "", stringsAsFactors = F)
    new.data[, 2][check.for.distorted.markers(recomb.data)] = "T"
    new.data[, 2][check.for.distorted.markers(recomb.data, 0.01)] = "T1"
    new.data[, 2][check.for.distorted.markers(recomb.data, 0.001)] = "T2"
    colnames(new.data) = c("seg.ratio", "sig")
    # new.data$chromo = unique(recomb.data[1, 3:ncol(recomb.data)])
    new.data
}

make.plots = function(list.of.recomb.data, plot.titles){
    #args:
    # list.of.recomb.data - a list of genotyping dataframes
    count1 = make.counter()
    if(missing(plot.titles)) plot.titles = 1:length(list.of.recomb.data)

    lapply(list.of.recomb.data, function(x){
        g = prepare.gg.data(x)
        # browser()
        ggplot(g, aes(x = 1:nrow(g), y = seg.ratio, group = sig, color = sig)) + geom_point() + 
            geom_hline(yintercept = 0.5) + ggtitle(plot.titles[count1()])# + coord_cartesian(ylim = c(0.5, 1.6))
    })
}

find.recomb.w.seg.dist = function(list.of.recomb.data){
    unlist(lapply(list.of.recomb.data, function(x){
        g = which(unlist(lapply(convert.recomb.to.seg2(x), function(x2) x2[3])) < 0.05)
        if(length(g) > 0) return(T)
        if(length(g) == 0) return(F)
    }))
}

#### ANALYSIS FUNCTIONS 2 ####

#functions to evaluate lists of genotype data simulations

number.within.10 = function(x, selection.position){
    #grab number of simulations in which the peak of segregation distortion is within 10 markers of the point of selection
    #args:
    #x - a list of genotyping dataframes 

    if(missing(selection.position)) selection.position = 200

    g = unlist(lapply(x, function(y){
        g = convert.recomb.to.seg(y)
        g1 = abs(g - 0.5)
        mean(which(g1 == max(g1)))
    }))

    length(which(g > (selection.position - 10) & g < (selection.position + 10)))
}

number.dist.markers = function(y, threshold, F2){
    #check number of simulations with distorted markers
    #args:
    #y - a list of genotyping dataframes 
    #threshold - integer, p-value threshold, defaults to 0.05
    if(missing(threshold)) threshold = 0.05
    if(missing(F2)) F2 = F
    
    length(which(unlist(lapply(y, function(x){
    g = check.for.distorted.markers(x, threshold, F2 = F2)
    if(length(g) > 0) return(T)
    return(F)
    }))))
}


number.dist.markers.fdr = function(y, threshold, method1, F2){
    #check number of simulations with distorted markers
    #args:
    #y - a list of genotyping dataframes 
    #threshold - integer, p-value threshold, defaults to 0.05
    if(missing(threshold)) threshold = 0.05
    if(missing(method1)) method1 = "BH"
    if(missing(F2)) F2 = F
    
    length(which(unlist(lapply(y, function(x){
    g = check.for.distorted.markers.fdr(x, threshold, method1, F2 = F2)
    if(length(g) > 0) return(T)
    return(F)
    }))))
}

mean.sd.num.distorted = function(y, threshold, F2){
    #grab the mean and standard deviations of number of distorted markers for each group of 1000 simulations
    #args:
    #y - a list of genotyping dataframes 

    if(missing(threshold)) threshold = 0.05
    if(missing(F2)) F2 = F
    q2 = unlist(lapply(y, function(x){
        g = check.for.distorted.markers(x, threshold, F2 = F2)
        length(g)
    }))
    
    g1 = list(mean(q2[which(q2 != 0)]), sd(q2[which(q2 != 0)]))
    names(g1) = c("mean", "sd")
    g1
}

mean.sd.num.distorted.fdr = function(y, threshold, method1, F2){
    #grab the mean and standard deviations of number of distorted markers for each group of 1000 simulations
    #args:
    #y - a list of genotyping dataframes 
    if(missing(threshold)) threshold = 0.05
    if(missing(method1)) method1 = "BH"
    if(missing(F2)) F2 = F

    q2 = unlist(lapply(y, function(x){
        g = check.for.distorted.markers.fdr(x, threshold, method1, F2 = F2)
        length(g)
    }))
    
    g1 = list(mean(q2[which(q2 != 0)]), sd(q2[which(q2 != 0)]))
    names(g1) = c("mean", "sd")
    g1
}


peak.dist = function(x){
    #what is the mean peak of distortion and what is the sd of the peak of distortion between simulations
    #args:
    #x - a list of genotyping dataframes 
    g = unlist(lapply(x, function(y){
        g = convert.recomb.to.seg(y)
        g1 = abs(g - 0.5)
        mean(which(g1 == max(g1)))
    }))

    g
}

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



#ns = number.sims.w.dist.markers
#mn = mean.number.distorted.markers.in.sims.w.distortion
#sdn = sd.number.distorted.markers.in.sims.w.distortion
#nsp = number.of.simulations.where.peak.of.distortion.is.within.10.markers.of.selection.locus
#mmp = mean.magnitude.of.peak.distortion
#mmpsd = sd.magnitude.of.peak.distortion
#maxmag = max.magnitude.of.peak.distortion
#minmag = min.magnitude.of.peak.distortion


perform.analysis = function(geno.list, selection, selection.position, selection.strength, recombination.position, F2, fgen){
    #args:
    # geno.list - a list of genotyping dataframes
    # selection - boolean flag - was selection used in this simulation?
    # selection.position - integer, the position of selection
    # selection.strength - integer, the strength of selection
    # recombination.position - string, the recombination profile (if any) used
    # F2 - boolean flag indicating whether to use f2 type test of segregation distortion, otherwise assumes 1:1 ratio of homozygotes as distribution (ignores hets)
    # fgen - integer specifying the filial generation of the simulation
    

    if(missing(selection)) selection = "-"
    if(missing(selection.position)) selection.position = "-"
    if(missing(selection.strength)) selection.strength = "-"    
    
    if(missing(F2)) F2 = F
    if(missing(recombination.position)){
        return("need to enter recombination.position")
    } else {
    
    number.simulations = length(geno.list)
    pop.size = nrow(geno.list[[1]])
    num.markers = ncol(geno.list[[1]])    
    ns0.05 = number.dist.markers(geno.list, 0.05, F2 = F2)    
    ns0.01 = number.dist.markers(geno.list, 0.01, F2 = F2)
    ns0.001 = number.dist.markers(geno.list, 0.001, F2 = F2)
    ns.fdr.0.05 = number.dist.markers.fdr(geno.list, 0.05, F2 = F2)
    ns.bon.0.05 = number.dist.markers.fdr(geno.list, 0.05, "bonferroni", F2 = F2)
    mnsd = mean.sd.num.distorted(geno.list, 0.05)
    mn.0.05 = mnsd[[1]]
    sdn.0.05 = mnsd[[2]]
    mnsd = mean.sd.num.distorted(geno.list, 0.01)
    mn.0.01 = mnsd[[1]]
    sdn.0.01 = mnsd[[2]]
    mnsd = mean.sd.num.distorted(geno.list, 0.001)
    mn.0.001 = mnsd[[1]]
    sdn.0.001 = mnsd[[2]]
    mnsd = mean.sd.num.distorted.fdr(geno.list, 0.05)
    mn.fdr.0.05 = mnsd[[1]]
    sdn.fdr.0.05 = mnsd[[2]]
    mnsd = mean.sd.num.distorted.fdr(geno.list, 0.05, "bonferroni")

    mn.bon.0.05 = mnsd[[1]]
    sdn.bon.0.05 = mnsd[[2]]
    if(selection == "-" | selection.position == "-"){
        nsp = "-"
    } else {
        nsp = number.within.10(geno.list, selection.position)    
    }
    
    mag.distortion = check.magnitude.of.distortion(geno.list)
    mmp = mean(mag.distortion)
    mmpsd = sd(mag.distortion)
    maxmag = max(mag.distortion)
    minmag = min(mag.distortion)

    summarydf1 = newdf(c("R.object.name", "pop.size", "fgen", "selection", "selection.position", "selection.strength", "recombination.position", "number.simulations", "number.markers", "ns.0.05", "ns.0.01", "ns.0.001", "ns.fdr.0.05", "ns.bon.0.05", "mn.0.05", "mn.0.01", "mn.0.001", "mn.fdr.0.05", "mn.bon.0.05", "sdn.0.05", "sdn.0.01", "sdn.0.001", "sdn.fdr.0.05", "sdn.bon.0.05", "nsp", "mmp", "mmpsd", "maxmag", "minmag"))    


    summarydf1[1, ] = c(deparse(substitute(geno.list)), pop.size, fgen, selection, selection.position, selection.strength, recombination.position, number.simulations, num.markers, ns0.05, ns0.01, ns0.001, ns.fdr.0.05, ns.bon.0.05, mn.0.05, mn.0.01, mn.0.001, mn.fdr.0.05, mn.bon.0.05, sdn.0.05, sdn.0.01, sdn.0.001, sdn.fdr.0.05, sdn.bon.0.05, nsp, mmp, mmpsd, maxmag, minmag)

    summarydf1

    # summarydf1$R.object.name = 
    # summarydf1$pop.size = pop.size
    # summarydf1$selection = selection
    # summarydf1$selection.position = selection.position
    # summarydf1$selection.strength = selection.strength
    # summarydf1$recombination.position = recombination.position
    # summarydf1$

    }

    
}

