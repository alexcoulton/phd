#array hybridization simulation


#NB array here refers to the programmatical sense, not the biological sense. I should of called this a "matrix"
generate.array = function(population.size, number.of.loci, static){
    #args:
    #static: if true, distribute genotypes in exactly 1:2:1 ratio, if false, use a multinomial distribution
    if(missing(static)) static = F
    array = as.data.frame(matrix(nrow = population.size, ncol = number.of.loci))
    
    if(static == T){
        make.genotypes = function(){
            first.a = sample(population.size, (population.size / 4))
            first.b = sample(as.vector(1:population.size)[-first.a], (population.size / 4))
            first.h = sample(as.vector(1:population.size)[-c(first.a, first.b)], (population.size / 2))
            return(list(first.a, first.b, first.h))
        }
        
        all.geno = lapply(1:number.of.loci, function(x) make.genotypes())
            
        for(i in 1:length(all.geno)){
            array[all.geno[[i]][[1]], i] = "a"
            array[all.geno[[i]][[2]], i] = "b"
            array[all.geno[[i]][[3]], i] = "h"
        }
    } else {
        make.genotypes = function(pop){
            geno = unlist(lapply(1:pop, function(x) which(rmultinom(1, 1, c(0.25, 0.25, 0.5)) == 1)))
            geno[geno == 1] = "a"
            geno[geno == 2] = "b"
            geno[geno == 3] = "h"
            return(geno)
        }
        
        all.geno = lapply(1:number.of.loci, function(x) make.genotypes(population.size))
        
        array = as.data.frame(all.geno)
        colnames(array) = c("loc.1", "loc.2", "loc.3")
    }
    
    return(array)
}

gen.seg.dist = function(array, prob){
    #args:
    #prob: a vector of probabilities (for binomial, do c(0.5, 0.5))
    array$pick = ""
    for(i in 1:nrow(array)){
        geno.coord = which(rmultinom(1, 1, prob) == 1)
        # binary = rbinom(1, 1, prob) # code for binomial distribution, not needed anymore as this can be acheived w/ multinomial dist.
        array$pick[i] = array[i, geno.coord]
        
    }
    
    table(array$pick)
}


gen.chi.sq.distribution = function(pop.size, is.static){
    if(missing(is.static)) is.static = F
    if(is.static == F){
        array1 = generate.array(pop.size, 3)
    } else {
        array1 = generate.array(pop.size, 3, T)
    }
    
    return.p.value = function(prob.vector){
        d = as.data.frame(gen.seg.dist(array1, prob.vector))
        d = rbind(d, d[3, ])
        d[3:4, 2] = (d[3, 2] / 2)
        # browser()
        chisq.test(d$Freq)[3]
    }
    
    gen.p.dist = function(prob.vector2) lapply(1:100, function(x) return.p.value(prob.vector2))
    
    p.equal.prob = gen.p.dist(c(1/3, 1/3, 1/3))
    p.no.int = gen.p.dist(c(1, 0, 0))
    p.small.int = gen.p.dist(c(0.8, 0.1, 0.1))
    p.medium.int = gen.p.dist(c(0.7, 0.2, 0.1))
    
    
    get.stats = function(p.dist){
        q1 = mean(unname(unlist(p.dist)))
        q2 = sd(unname(unlist(p.dist)))
        q3 = data.frame(q1, q2)
        return(q3)
    }
    
    d.list = lapply(list(p.no.int, p.small.int, p.medium.int, p.equal.prob), get.stats)
    d = rbind(d.list[[1]], d.list[[2]], d.list[[3]], d.list[[4]])
    library(dplyr)
    d = add_column(d, cat = "")
    d$cat = factor(c("1, 0, 0", "0.8, 0.1, 0.1", "0.7, 0.2, 0.1", "1/3, 1/3, 1/3"))
    d$cat = factor(d$cat, levels = c("1, 0, 0", "0.8, 0.1, 0.1", "0.7, 0.2, 0.1", "1/3, 1/3, 1/3"))
    
    attr(d, "array") = array1
    
    return(d)
}

pop1000 = lapply(1:4, function(x) gen.chi.sq.distribution(1000))
pop100 = lapply(1:10, function(x) gen.chi.sq.distribution(100))
pop60 = lapply(1:4, function(x) gen.chi.sq.distribution(60))
pop30 = lapply(1:4, function(x) gen.chi.sq.distribution(30))

c10 = make.counter()
pop.strings = c(1000, 100, 60, 30)

make.plot = function(chi.sq.df, pop.size){
    ggplot(chi.sq.df, aes(x = cat, y = q1)) + geom_bar(stat = "identity", width = 0.5) + 
        geom_errorbar(aes(ymin = q1 - q2, ymax = q1 + q2, width = 0.1)) +
        xlab("Probability distribution") + ylab("Chi-square p-value") + 
        geom_hline(yintercept = 0.05) +
        labs(title = p("Pop. size: ", pop.size)) 
}

plots = lapply(list(pop1000, pop100, pop60, pop30), make.plot)    

library(gridExtra)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]])

