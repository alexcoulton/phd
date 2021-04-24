options(show.error.locations = T)

options(error = function() { traceback(2); if(!interactive()) quit('no', status = 1, runLast = FALSE) })


#### FUNCTIONS ####

perform.dqc.analysis = function(ind.avg.df, var.name1, dqc.info1){
    #ind.avg.df - MRD df produced by ind.avg
    #var.name1 - variable name to test (phys.dist.long, phys.dist.short, or phys.dist)
    #dqc.info1 dqc info data frame from read_delim
    
    axp.avg.chromo = arrange(ind.avg.df, chromo)
    axp.avg.chromo2 = split.df.into.list.of.dfs.by.column(axp.avg.chromo, "chromo")
    
    #match up phys.dist and dqc/qc information
    axp.avg.chromo3 = lapply(axp.avg.chromo2, function(x){
        x$individual = multi.str.split(as.character(x$individual), "_call_code", 1)
        lowqc.temp = dqc.info1[-which(!dqc.info1$`Sample Filename` %in% x$individual), ]
        lowqc.temp2 = bind_cols(lowqc.temp, x[match(lowqc.temp$`Sample Filename`, x$individual), ])
        lowqc.temp2[c("Sample Filename", "DQC", "QC call_rate", "phys.dist", "phys.dist.long", "phys.dist.short", "marker.dist", "chromo")]
    })
    
    #(6)
    #linear models of MRD vs DQC
    lm1 = lapply(axp.avg.chromo3, function(x){
        # browser()
        summary(lm(x[[var.name1]] ~ x$DQC))
    })
    
    lm2 = lapply(axp.avg.chromo3, function(x){
        # browser()
        summary(lm(x[[var.name1]] ~ x$`QC call_rate`))
    })
    
    
    dqcp1 = p.adjust(sapply(lm1, function(x){
        x[[4]][2, 4]
    }), method = "fdr")
    
    
    dqcr1 = round(unlist(sapply(lm1, function(x){
        x[8]
    })), digits = 4)
    
    qcp2 = p.adjust(sapply(lm2, function(x){
        x[[4]][2, 4]
    }), method = "fdr")
    
    qcr2 = round(unlist(sapply(lm2, function(x){
        x[8]
    })), digits = 4)
    
    axp.dqc.df = data.frame(dqcp1, dqcr1, qcp2, qcr2)
    axp.dqc.df
}

detect.recombination = function(genotype.data, conservative, liss.df){
    #detects recombination events in genotyping data by looking for transitions between genotypes. the transition.coord column represents the coordinate of the probe before the recombination event in the genotyping dataframe. the num.events column is either 1 or 2 depending on the transition type (a transition from homozygous genotype of parent 1 to homozygous parent 2 counts as two events for an F2 population).
    #args:
    # genotype.data - a genotyping dataframe in rQTL format, the same format as alexq1_man_cur_sk_map_genogeno.csv
    # conservative - Boolean, if TRUE, removes recombination events in which the middle marker in a sequence of three markers is different from the flanking markers, i.e. A B A
        #n.b. the conservative argument has now been repurposed to function with the variable booboo1 to only remove double events that are within 30 (or x) Mb of each other.
    if(missing(conservative)) conservative = F
    genotype.data[1, 1] = NA
    orig.data = genotype.data
    genotype.data = as.data.frame(t(genotype.data))
    genotype.data = convert.to.character.data.frame(genotype.data)
    
    
    colnames(genotype.data) = genotype.data[2, ]
    colnames(genotype.data)[1] = "chromo"
    
    genotype.data2 = split.df.into.list.of.dfs.by.column(genotype.data, "chromo", F)
    
    count.chromo = 3
    all.recomb = lapply(genotype.data2[3:length(genotype.data2)], function(x){
        
        # browser()
        
        count1 = 2
        
        # counter.debug = 1
        recomb.individuals = lapply(x[, 2:ncol(x)], function(y){
            # counter.debug <<- counter.debug + 1
            
            #if all of the markers are nocalls, change them all to the same genotype
            if(all(y == "-")) y = "B"
            
            # special case, deal with genotypes that have nocalls at the first marker:
            if(y[1] == "-"){
                
                y2 = which(y == "-")
                
                #find start and end nocalls
                y3 = paste(sapply(split(y2,cumsum(c(1,diff(y2)-1))),deparse),collapse=", ")
                
                start.nocalls = regmatches(y3, regexpr("^1\\:[0-9]*", y3))
                
                if(length(start.nocalls) == 0){ #if only the first marker is a nocall (not the second or third)
                    y[1] = y[2]
                } else {
                    start.nocalls = strsplit(start.nocalls, ":")
                    start.nocalls = as.numeric(start.nocalls[[1]])
                    start.nocalls = seq(start.nocalls[1], start.nocalls[2], 1)
                    y[start.nocalls] = y[(max(start.nocalls) + 1)]
                }
            }
            
            # tryCatch(if((y[length(y)] == "-") != T & (y[length(y)] == "-") != F) browser(), error = function(e) browser())
            
            #special case, deal with genotypes that have nocalls at the last marker:
            if(y[length(y)] == "-"){
                
                y2 = which(y == "-")
                
                #find start and end nocalls
                y3 = paste(sapply(split(y2,cumsum(c(1,diff(y2)-1))),deparse),collapse=", ")
                
                end.nocalls = regmatches(y3, regexpr(paste0("[0-9]*\\:", length(y), "$"), y3))
                
                if(length(end.nocalls) == 0){
                    y[length(y)] = y[(length(y) - 1)]
                } else {
                    end.nocalls = strsplit(end.nocalls, ":")
                    end.nocalls = as.numeric(end.nocalls[[1]])
                    end.nocalls = seq(end.nocalls[1], end.nocalls[2], 1)
                    y[end.nocalls] = y[(min(end.nocalls) - 1)]    
                }
            }
            
            
            #begin main analysis
            g1 = y
            g2 = y
            g1[g1 == "A"] = 1
            g1[g1 == "H"] = 10
            g1[g1 == "B"] = 100
            g1[g1 == "-"] = 100000
            g1 = as.numeric(g1)
            #get all transition points
            transition.coord = which(diff(g1) != 0 & abs(diff(g1)) < 10000)
            
            a1 = g2[transition.coord]
            a2 = g2[transition.coord + 1]
            a3 = g2[transition.coord + 2]
            
            geno.before = as.character(a1)
            geno.after = as.character(a2)
            geno.after2 = as.character(a3)
            
            recombination.events = data.frame(geno.before, geno.after, geno.after2, transition.coord, stringsAsFactors = F)
            
            if(nrow(recombination.events) != 0){
                recombination.events$individual = colnames(x)[count1]
                recombination.events$individual.id = count1
                recombination.events$chromo = names(genotype.data2)[count.chromo]    
                recombination.events$probe.before.transition = rownames(x)[recombination.events$transition.coord]
                recombination.events$probe.after.transition = rownames(x)[(recombination.events$transition.coord + 1)]
                recombination.events$probe.after.transition2 = rownames(x)[(recombination.events$transition.coord + 2)]
                recombination.events$nocall.trans = F
                
            }
            
            
            #dealing with special case A - B or A - - B or A - - - B etc (regardless of how long the run of nocalls in between genotypes is).
             #i.e. where a legitimate recombination event occurs (change in parental genotype phase) that is seperated by nocalls
            transition.coord.nocall = which(diff(g1) > 10000)
            transition.coord.nocall2 = which(diff(g1) < -10000)
            
            before.nocall = g2[transition.coord.nocall]
            after.nocall = g2[transition.coord.nocall2 + 1]
            after.nocall2 = g2[transition.coord.nocall2 + 2]
            
            # if(length(before.nocall) != length(after.nocall)) browser()
            
            nocall.transitions = which(before.nocall != after.nocall)
            
            if(length(nocall.transitions) > 0){
                recombination.events2 = data.frame(before.nocall[nocall.transitions], after.nocall[nocall.transitions], after.nocall2[nocall.transitions], transition.coord.nocall[nocall.transitions], stringsAsFactors = F)
                colnames(recombination.events2) = c("geno.before", "geno.after", "geno.after2", "transition.coord")
                recombination.events2$individual = colnames(x)[count1]
                recombination.events2$individual.id = count1
                recombination.events2$chromo = names(genotype.data2)[count.chromo]
                recombination.events2$probe.before.transition = rownames(x)[transition.coord.nocall[nocall.transitions]]
                recombination.events2$probe.after.transition = rownames(x)[transition.coord.nocall2[nocall.transitions] + 1]
                recombination.events2$probe.after.transition2 = rownames(x)[transition.coord.nocall2[nocall.transitions] + 2]
                recombination.events2$nocall.trans = T
                
                if(nrow(recombination.events) > 0 & nrow(recombination.events2) > 0){
                    recombination.events = as.data.frame(rbind(recombination.events, recombination.events2))
                    
                }
                
            }
            
            
            
            
            # if(length(transition.coord.nocall2) > 0) browser()
            
            
            
            
            count1 <<- count1 + 1
            # browser()
            recombination.events
        })
        
        
        count.chromo <<- count.chromo + 1
        
        all.recomb.ind = combine.list.of.data.frames(recomb.individuals)
        all.recomb.ind
        
    })
    
    all.recomb2 = combine.list.of.data.frames(all.recomb)
    all.recomb2$num.events = 1
    w1 = which(all.recomb2$geno.before %in% c("A", "B") & all.recomb2$geno.after %in% c("A", "B"))

    if(length(w1) > 0){
        all.recomb2[w1, ]$num.events = 2
    }

    all.recomb2$individual.trunc = multi.str.split(all.recomb2$individual, "_call_code", 1)
    
    
    if(conservative == T){
        #remove recombination events of A B A or B A B etc.
        remove.double.events = function(x){
            #this now splits the recombination dataframe in two, with non-double events retained in all.recomb2, and double events sent to booboo1
            g = split(x, x$chromo)
            
            over.chromosomes = lapply(g, function(y){
                g2 = split(y, y$individual)
                
                over.individuals = lapply(g2, function(z){
                    #get double event coords
                    coord1 = which(diff(z[which(z$nocall.trans == F), ]$transition.coord) == 1)
                    # if(length(coord1) > 0) browser()
                    
                    z.test = z[coord1, ]
                    if(nrow(z.test) > 0){
                        g = lapply(1:nrow(z.test), function(qqq){
                            g1 = liss.df[[2]][which(liss.df[[2]]$marker == z.test$probe.before.transition[qqq]), ]$phys.bp
                            g2 = liss.df[[2]][which(liss.df[[2]]$marker == z.test$probe.after.transition[qqq]), ]$phys.bp
                            g3 = liss.df[[2]][which(liss.df[[2]]$marker == z.test$probe.after.transition2[qqq]), ]$phys.bp
                            
                            (g3 - g1) / 1000000
                        })
                        
                        z.test$dist.double = unlist(g)
                        booboo1 <<- c(booboo1, list(z.test))
                    }
                    
                    #remove double event coords from main (these are now in booboo1)
                    if(length(coord1) > 0) z = z[-c(coord1, (coord1 + 1)), ]
                    
                    z
                })
                
                bind_rows(over.individuals)
            })
            bind_rows(over.chromosomes)
        }
        
        all.recomb2 = remove.double.events(all.recomb2)
    }
    
    all.recomb2
}

makehistlist.new = function(detect.recomb.df, physical){
    #counts recombination events 
    #args:
    # detect.recomb.df - recombination dataframe produced by detect.recombination()
    # physical - boolean, if true return physical positions of markers rather than marker positions in LG
    
    if(missing(physical)) physical = F
    
    recombination2 = split.df.into.list.of.dfs.by.column(detect.recomb.df, "chromo", F)
    
    if(physical == F){
        lapply(recombination2, function(x){
            if(length(which(x$num.events == 2)) > 0){
                two.events = rep(x[which(x$num.events == 2), ]$transition.coord, 2)
                
                one.events = x[which(x$num.events == 1), ]$transition.coord
                
                return(c(two.events, one.events))
                
                
            } else {
                return(x$transition.coord)
            }
        })    
    } else {
        lapply(recombination2, function(x){
            if(length(which(x$num.events == 2)) > 0){
                two.events = rep(x[which(x$num.events == 2), ]$phys.pos, 2)
                
                # browser()
                
                one.events = x[which(x$num.events == 1), ]$phys.pos
                
                return(c(two.events, one.events))
                
                
            } else {
                return(x$phys.pos)
            }
        })
    }
}

grab.recombinations = function(recomb.df, geno.df, chromo1, position){
    #shows the recombination events as they occur in the genotyping dataframe
    #args:
    # recomb.df - recombination dataframe produced by detect.recombination
    # geno.df - genotyping dataframe in rQTL format
    # chromo1 - string, name of the chromosome or linkage group
    # position - integer, coordinate position of marker on the chromosome of choice
    
    recomb.df2 = recomb.df[which(recomb.df$transition.coord == position & recomb.df$chromo == chromo1), ]
    geno.df = geno.df[, c(1, 2, which(geno.df[1, ] == chromo1))]
    
    geno.df[which(geno.df$probeset_id %in% recomb.df2$individual), c(2, (position + 2):(position + 3))]
}

examine.posteriors = function(axiom.folder){
    #grabs the areas of posterior cluster assignments
    #args:
    # axiom.folder - the folder containing the posterior file (usually C:/Users/Public/Documents/AxiomAnalysisSuite/Output/PROJECT.NAME/)
    
    sub_plot_cov2 <- function(vx, vy, cov, x=0, y=0, col.ellipse, np=100, ...) {
        # browser()
        theta <- 0.5 * atan2(cov*2, vx-vy)
        alpha <- 2 * pi * (0:np)/np
        sint <- sin(theta)
        cost <- cos(theta)
        sina <- sin(alpha)
        cosa <- cos(alpha)
        a <- 2*sqrt(vx*cost*cost + vy*sint*sint + cov*2*sint*cost)
        b <- 2*sqrt(vx*sint*sint + vy*cost*cost - cov*2*sint*cost)
        x <- x + a * cosa * cost - b * sina * sint
        y <- y + a * cosa * sint + b * sina * cost
        
        output1 = list(x, y)
        names(output1) = c("x", "y")
        output1
        # lines(x, y, type="l", col=col.ellipse, ...)
        # browser()
    }
    
    #parse cluster posteriors file
    posteriors.raw = read.delim(paste0(axiom.folder, "AxiomGT1.snp-posteriors.txt"), comment.char = "#", sep = "\t", stringsAsFactors = F)
    posteriors.raw
    
    bb = lapply(posteriors.raw$BB, function(x){
        g = as.data.frame(t(as.data.frame(strsplit(x, ","))))
        convert.to.character.data.frame(g)
        colnames(g) = c("x", "vx", "unk1", "unk2", "y", "vy", "cov")
        g
    })
    
    bb2 = bind_rows(bb)
    
    ab = lapply(posteriors.raw$AB, function(x){
        g = as.data.frame(t(as.data.frame(strsplit(x, ","))))
        convert.to.character.data.frame(g)
        colnames(g) = c("x", "vx", "unk1", "unk2", "y", "vy", "cov")
        g
    })
    
    ab2 = bind_rows(ab)
    
    aa = lapply(posteriors.raw$AA, function(x){
        g = as.data.frame(t(as.data.frame(strsplit(x, ","))))
        convert.to.character.data.frame(g)
        colnames(g) = c("x", "vx", "unk1", "unk2", "y", "vy", "cov")
        g
    })
    
    aa2 = bind_rows(aa)
    
    bb2 = cbind(posteriors.raw$id, bb2)
    ab2 = cbind(posteriors.raw$id, ab2)
    aa2 = cbind(posteriors.raw$id, aa2)
    
    posteriors2 = list(bb2, ab2, aa2)
    posteriors2
    
    post.assign = c("bb", "ab", "aa")
    
    for(i in 1:3){
        posteriors2[[i]]$geno = post.assign[i]
    }
    
    post3 = bind_rows(posteriors2)
    
    post3$x = as.numeric(post3$x)
    post3$vx = as.numeric(post3$vx)
    post3$unk1 = as.numeric(post3$unk1)
    post3$unk2 = as.numeric(post3$unk2)
    post3$y = as.numeric(post3$y)
    post3$vy = as.numeric(post3$vy)
    post3$cov = as.numeric(post3$cov)
    
    post4 = transpose(post3)
    rownames(post4) = colnames(post3)
    
    post5 = lapply(post4, function(x){
        as.matrix(as.data.frame(sub_plot_cov2(as.numeric(x[3]), as.numeric(x[7]), as.numeric(x[8]), as.numeric(x[2]), as.numeric(x[6]))))
    })
    
    post6 = sapply(post5, function(x){
        slot(Polygon(x), "area")
    })
    
    post3$area = post6
    post3    
    
}


evaluate.cluster.plots = function(recomb.df, axiom.folder){
    #get recombination events which immediately change back
    #args:
    # recomb.df - recombination dataframe produced with detect.recombination
    # axiom.folder - path to axiom folder containing calls, posteriors and summary txt files (usually C:/Users/Public/Documents/AxiomAnalysisSuite/Output/PROJECT.NAME/)
    
    sub_plot_cov2 <- function(vx, vy, cov, x=0, y=0, col.ellipse, np=100, ...) {
        # browser()
        theta <- 0.5 * atan2(cov*2, vx-vy)
        alpha <- 2 * pi * (0:np)/np
        sint <- sin(theta)
        cost <- cos(theta)
        sina <- sin(alpha)
        cosa <- cos(alpha)
        a <- 2*sqrt(vx*cost*cost + vy*sint*sint + cov*2*sint*cost)
        b <- 2*sqrt(vx*sint*sint + vy*cost*cost - cov*2*sint*cost)
        x <- x + a * cosa * cost - b * sina * sint
        y <- y + a * cosa * sint + b * sina * cost
        
        output1 = list(x, y)
        names(output1) = c("x", "y")
        output1
        # lines(x, y, type="l", col=col.ellipse, ...)
        # browser()
    }
    
    recomb.df2 = recomb.df[which(recomb.df$geno.before == recomb.df$geno.after2), ]
    recomb.df2 = arrange(recomb.df2, probe.before.transition)
    
    change.mid.to.nocall = as.numeric()
    
    calls.raw = read.delim(paste0(axiom.folder, "AxiomGT1.calls.txt"), comment.char = "#", sep = "\t", stringsAsFactors = F)
    
    #parse cluster posteriors file
    posteriors.raw = read.delim(paste0(axiom.folder, "AxiomGT1.snp-posteriors.txt"), comment.char = "#", sep = "\t", stringsAsFactors = F)
    posteriors.raw
    
    bb = lapply(posteriors.raw$BB, function(x){
        g = as.data.frame(t(as.data.frame(strsplit(x, ","))))
        convert.to.character.data.frame(g)
        colnames(g) = c("x", "vx", "unk1", "unk2", "y", "vy", "cov")
        g
    })
    
    bb2 = bind_rows(bb)
    
    ab = lapply(posteriors.raw$AB, function(x){
        g = as.data.frame(t(as.data.frame(strsplit(x, ","))))
        convert.to.character.data.frame(g)
        colnames(g) = c("x", "vx", "unk1", "unk2", "y", "vy", "cov")
        g
    })
    
    ab2 = bind_rows(ab)
    
    aa = lapply(posteriors.raw$AA, function(x){
        g = as.data.frame(t(as.data.frame(strsplit(x, ","))))
        convert.to.character.data.frame(g)
        colnames(g) = c("x", "vx", "unk1", "unk2", "y", "vy", "cov")
        g
    })
    
    aa2 = bind_rows(aa)
    
    bb2 = cbind(posteriors.raw$id, bb2)
    ab2 = cbind(posteriors.raw$id, ab2)
    aa2 = cbind(posteriors.raw$id, aa2)
    
    posteriors2 = list(bb2, ab2, aa2)
    
    #parse signals file
    signals.raw = read.delim(paste0(axiom.folder, "AxiomGT1.summary.txt"), comment.char = "#", sep = "\t", stringsAsFactors = F)
    
    for(i in 1:nrow(recomb.df2)){
        #write sample file for Ps_Visualization()
        sample.name = recomb.df2$individual.trunc[[i]]
        sample.name2 = gsub(" ", "\\.", sample.name)
        sample.name2 = gsub("-", "\\.", sample.name2)
        
        probes = switch.affy.format(c(recomb.df2$probe.before.transition[[i]], 
                                                                    recomb.df2$probe.after.transition[[i]], recomb.df2$probe.after.transition2[[i]]))
        
        calls2 = calls.raw[na.omit(match(probes, calls.raw$probeset_id)), which(colnames(calls.raw) == sample.name2)]
        
        genotype.positions = as.numeric()
        
        for(i2 in 1:length(probes)){
            
            signals2 = signals.raw[grep(probes[i2], signals.raw$probeset_id), c(1, which(colnames(signals.raw) == sample.name2))]
            
            sample.x.coord = log2(signals2[1, 2]) - log2(signals2[2, 2])
            sample.y.coord = ((log2(signals2[1, 2]) + log2(signals2[2, 2])) / 2)
            
            posteriors3 = lapply(posteriors2, function(z){
                z = z[which(z[, 1] == probes[i2]), ]
                z
            })
            
            if(calls2[i2] == 2) num1 = 1
            if(calls2[i2] == 1) num1 = 2
            if(calls2[i2] == 0) num1 = 3
            
            oval.data = sub_plot_cov2(as.numeric(posteriors3[[num1]]$vx), 
                                                                as.numeric(posteriors3[[num1]]$vy), 
                                                                as.numeric(posteriors3[[num1]]$cov), 
                                                                as.numeric(posteriors3[[num1]]$x), 
                                                                as.numeric(posteriors3[[num1]]$y), "#000000", np = 100)
            
            genotype.positions = c(genotype.positions, point.in.polygon(sample.x.coord, sample.y.coord, oval.data$x, oval.data$y))
            
        }
        
        if(genotype.positions[2] == 0){
            change.mid.to.nocall = c(change.mid.to.nocall, 1)
        } else {
            change.mid.to.nocall = c(change.mid.to.nocall, 0)
        }
        
        print(paste0("done row ", i, " of ", nrow(recomb.df2)))
        
        
    }
    
    change.mid.to.nocall
    
    recomb.df3 = cbind(recomb.df2, change.mid.to.nocall)
    
    recomb.df4 = recomb.df3[which(recomb.df3$change.mid.to.nocall == 1), ]
    change.to.nocall = recomb.df4[c("individual", "probe.after.transition")]
    change.to.nocall
    
}

examine.cluster.plots = function(recomb.df){
    #performs plotting to pdf
    
    
    
    for(i in 1:nrow(recomb.df)){
        
        # browser()
        sample.name = recomb.df$individual.trunc[[i]]
        samp.file = data.frame(sample.name, "black")
        colnames(samp.file) = c("sample", "color")
        write.table(samp.file, "samp.file.txt", sep = "\t", quote = F)
        
        pid.file = c("probeset_id", switch.affy.format(c(recomb.df$probe.before.transition[[i]], 
                                                                                                         recomb.df$probe.after.transition[[i]], recomb.df$probe.after.transition2[[i]])))
        writeLines(pid.file, "pid.file.txt")
        
        cluster.data = Ps_Visualization(pidFile = "pid.file.txt",
                                                                        sampleFile = "samp.file.txt", output.File = paste0("plots/plot", i, ".pdf"),
                                                                        summaryFile = "AxiomGT1.summary.txt", callFile = "AxiomGT1.calls.txt",
                                                                        posteriorFile = "AxiomGT1.snp-posteriors.txt",
                                                                        col.AA = "#d6d6d6",
                                                                        col.AB = "#d6d6d6",
                                                                        col.BB = "#d6d6d6",
                                                                        confidenceFile = "AxiomGT1.confidences.txt", plot.prior = T,
                                                                        plot.ref = F, plot.intensity = F, plot.type = "pdf", plot.width = 15,
                                                                        plot.height = 10, num.cols = 3, num.rows = 1)
        
    }
}

dist.from.cent = function(histlist, physical){
    #calculate distance of recombination positions from central marker
    #args:
    # histlist - output from makehistlist.new
    
    if(missing(physical)) physical = F
    
    if(physical == F){
        manntest1 = lapply(histlist, function(x){
            
            center = round(max(x) / 2)
            g = sort(x)
            abs(g - center)
            
        })
        
    } else {
        manntest1 = lapply(histlist, function(x){
            
            center = sort(x - 50)
            abs(center)
            
        })
    }
    
    # browser()
    manntest1
}

perform.test = function(list.of.histlists, treatments.comp, physical){
    #UPDATE: INVALID STATISTICAL ANALYSIS. RECOMBINATION EVENTS FROM THE SAME INDIVIDUAL ARE NOT INDEPENDENT - see notebook
    #performs mann-whitney comparison of distance of recombination events from centre
    #of chromosome between chromosomes of different treatments
    #args:
    # list.of.histlists - list of outputs from makehistlist.new for each treatment / population
    # treatments.comp - numeric vector with two elements, each corresponding the population to be compared (e.g. c(1, 2))
    
    if(missing(physical)) physical = F
    
    
    comp1 = treatments.comp
    
    axc.test = lapply(list.of.histlists, function(x){
        dist.from.cent(x, physical)
    })
    
    #tranpose list of lists
    z2 = lapply(1:length(list.of.histlists[[1]]), function(a){
        lapply(axc.test, function(a1){
            a1[[a]]
        })
    })
    
    axc.test2 = lapply(z2, function(x){
        
        wilcox.test(x[[comp1[1]]], x[[comp1[2]]])
        
    })
    
    axc.test3 = unlist(lapply(axc.test2, function(x){
        x[3]
    }))
    
    axc.test3
}

#version of ASMap::mstmap.data.frame() function without column name spaces bug
m2 = function (object, pop.type, dist.fun = "kosambi", objective.fun = "COUNT", 
                             p.value = 1e-06, noMap.dist = 15, noMap.size = 0, miss.thresh = 1, 
                             mvest.bc = FALSE, detectBadData = FALSE, as.cross = TRUE, 
                             return.imputed = FALSE, trace = FALSE, ...) 
{
    if (trace) 
        trace <- "MSToutput.txt"
    if (is.character(trace)) {
        ftrace <- file(trace, "w")
        sink(trace, type = "output", append = FALSE)
        on.exit(sink(type = "output"))
        on.exit(close(ftrace), add = TRUE)
        trace <- TRUE
    }
    if (!all(sapply(object, is.character)) & !all(sapply(object, 
                                                                                                             is.numeric))) 
        stop("Columns of the marker data.frame must be all of type \"character\" or all of \"numeric\".")
    RILN <- paste("RIL", 1:20, sep = "")
    allow.pop <- c("BC", "DH", "ARIL", RILN)
    if (!(pop.type %in% allow.pop)) 
        stop("Population type needs to be \"BC\",\"DH\",\"ARIL\" or \"RILn\" (see ?mstmap.data.frame).")
    if (pop.type %in% RILN) 
        rnum <- substring(pop.type, 4, nchar(pop.type))
    if (!(dist.fun %in% c("haldane", "kosambi"))) 
        stop("Distance function needs to be \"haldane\" or \"kosambi\" (see ?mstmap.data.frame).")
    if (!(objective.fun %in% c("COUNT", "ML"))) 
        stop("Objective function needs to be \"COUNT\" or \"ML\" (see ?mstmap.data.frame).")
    alleles <- unique(unlist(lapply(object, unique)))
    allow.list <- c("A", "a", "B", "b", "-", "U")
    ptype <- pop.type
    
    
    if (pop.type %in% c("BC", "DH", "ARIL")) {
        if (all(sapply(object, is.numeric))) {
            if (any(apply(object, 2, function(el) length(el[is.na(el)])))) 
                stop("Numeric input cannot contain missing values")
            if (any(apply(object, 2, function(el) el < 0 | el > 
                                        1))) 
                stop("Ony values between 0 and 1 are allowed if input is numeric.")
            if (as.cross) 
                stop("Cross object cannot be returned for numeric input.")
        }
        else {
            if (!all(alleles %in% allow.list)) 
                stop("Non-allowable allele encodings in DH marker set")
            cons <- c("1", "1", "2", "2", NA, NA, NA)
        }
        pop.type <- "DH"
    }    else {
        
        if (is.numeric(object)) 
            stop("Numeric input is not available for RILn populations.")
        allow.list <- c(allow.list, "X")
        if (!all(alleles %in% allow.list)) 
            stop("Non-allowable allele encodings in RIL marker set")
        cons <- c("1", "1", "3", "3", NA, NA, "2")
    }
    
    # browser()
    
    if (length(grep(" ", rownames(object)))) {
        
        
        
        rownames(object) <- gsub(" ", "-", rownames(object))
        warning("Replacing spaces in genotype names with a " - 
                            " separator\n")
    }
    
    
    
    if (length(grep(" ", names(object)))) {
        # browser()
        names(object) <- gsub(" ", "&", names(object))
        # warning("Replacing spaces in marker nxames with a - separator\n")
    }
    
    dist.fun <- match.arg(dist.fun)
    objective.fun <- match.arg(objective.fun)
    param <- list(population_type = pop.type, distance_function = dist.fun, 
                                cut_off_p_value = as.double(p.value), no_map_dist = as.double(noMap.dist), 
                                no_map_size = as.integer(noMap.size), missing_threshold = as.double(miss.thresh), 
                                estimation_before_clustering = as.integer(mvest.bc), 
                                detect_bad_data = as.integer(detectBadData), objective_function = objective.fun, 
                                trace = as.integer(trace))
    mst <- .Call("mst", param, object)
    map <- lapply(mst, function(el) el$map)
    marks <- unlist(lapply(map, names))
    ordm <- pmatch(marks, rownames(object))
    if (any(is.na(ordm))) 
        omit.list <- t(object[is.na(ordm), drop = FALSE])
    ordm <- ordm[!is.na(ordm)]
    object <- object[ordm, ]
    chn <- paste("L", 1:length(map), sep = "")
    if (ptype == "ARIL") {
        imf <- switch(dist.fun, haldane = imf.h, kosambi = imf.k)
        mf <- switch(dist.fun, haldane = mf.h, kosambi = mf.k)
        map <- lapply(map, function(el, imf, mf) {
            nams <- names(el)
            rf <- mf(diff(el))
            rf <- (rf/2)/(1 - rf)
            el <- cumsum(c(0, imf(rf)))
            names(el) <- nams
            el
        }, imf, mf)
    }
    if (as.cross) {
        co <- list()
        spl <- rep(1:length(map), times = sapply(map, length))
        object <- as.matrix(object)
        al <- allow.list %in% alleles
        cons <- cons[al]
        al <- allow.list[al]
        for (i in 1:length(al)) object[object == al[i]] <- cons[i]
        dn <- dimnames(object)[[1]]
        object <- apply(object, 2, function(el) as.numeric(el))
        dimnames(object)[[1]] <- dn
        datm <- lapply(split.data.frame(object, spl), t)
        mo <- lapply(1:length(map), function(el, datm, map) {
            temp <- list()
            temp$data <- datm[[el]]
            temp$map <- map[[el]]
            class(temp) <- "A"
            temp
        }, datm, map)
        co$geno <- mo
        names(co$geno) <- chn
        if (return.imputed) {
            co$imputed.geno <- lapply(mst, function(el) {
                names(el)[2] <- "data"
                el$data <- t(el$data)
                el$map <- el$map[names(el$map) %in% dimnames(el$data)[[2]]]
                el$data <- el$data[, pmatch(names(el$map), dimnames(el$data)[[2]]), 
                                                     drop = FALSE]
                el
            })
            names(co$imputed.geno) <- chn
        }
        co$pheno <- data.frame(Genotype = factor(dimnames(object)[[2]]))
        wp <- (1:23)[allow.pop %in% ptype]
        class(co) <- c(c("bc", "dh", "riself", rep("f2", 20))[wp], 
                                     "cross")
        if (wp %in% 4:23) 
            co <- convert2bcsft(co, F.gen = wp - 3, estimate.map = FALSE)
        object <- co
    }
    else {
        do <- list()
        chrv <- rep(chn, times = sapply(map, length))
        dist <- unlist(map)
        do$geno <- cbind.data.frame(markers = marks, chr = chrv, 
                                                                dist = dist, object)
        rownames(do$geno) <- NULL
        if (return.imputed) {
            imp <- lapply(mst, function(el) {
                nm <- names(el$map)[names(el$map) %in% dimnames(el$imputed_values)[[1]]]
                pm <- pmatch(nm, dimnames(el$imputed_values)[[1]])
                el$imputed_values <- el$imputed_values[pm, , 
                                                                                             drop = FALSE]
                el$imputed_values
            })
            chri <- rep(chn, times = sapply(imp, function(el) dim(el)[1]))
            marki <- unlist(lapply(imp, rownames))
            disti <- dist[pmatch(marki, names(dist))]
            do$imputed.geno <- cbind.data.frame(markers = marki, 
                                                                                    chr = chri, disti = disti, do.call("rbind", imp))
            rownames(do$imputed.geno) <- NULL
        }
        object <- do
    }
    if (exists("omit.list")) 
        object$omit <- omit.list
    object
}



extract.centimorgan.from.genotypes = function(geno.df, pop.type){
    #assumes markers in the genotyping dataframe are already in the correct order and assigned to chromosomes
    #uses ASMap to convert data frame to rQTL cross object, then pulls the pariwise RF to retain marker order and 
    #calculates map distance directly from RFs using BatchMap::kosambi()
    #args:
    # geno.df - genotyping dataframe in rQTL (man_cur_sk) format
    if(missing(pop.type)) pop.type = "RIL2"
    if(!"ASMap" %in% (.packages())) library(ASMap)
    if(!"BatchMap" %in% (.packages())) library(BatchMap)
    t1 = transpose(geno.df)
    rownames(t1) = colnames(geno.df)
    colnames(t1) = t1[2, ]
    t1 = t1[-(1:2), ]
    colnames(t1)[1] = "chromo"
    t2 = split.df.into.list.of.dfs.by.column(t1, "chromo", F)
    
    t3 = lapply(t2, function(x){
        x = x[, -1]
        x[x == "H"] = "X"
        x
    })
    
    t4 = lapply(t3, function(x){
        maps1 = m2(x, pop.type = pop.type, dist.fun = "kosambi", p.value = 2)
        # browser()
        
        
        rf1 = as.data.frame(pull.rf(maps1))
        rf2 = rf1[match(rownames(x), rownames(rf1)) , match(rownames(x), colnames(rf1))]
        
        #remove first row to make consecutive marker pair RFs the diagonal of the matrix for extraction with diag()
        rf2 = rf2[, -1]
        rf2 = as.matrix(rf2)
        
        # as1 = pull.map(maps1, as.table = T)
        # as2 = abs(as1$pos - max(as1$pos))
        # as3 = as2[length(as2):1]
        
        
	rf2.values = diag(rf2)
	rf2.coords = which(rf2.values >= 0.5)
	if(length(rf2.coords) > 0) rf2.values[rf2.coords] = 0.49
        c(0, cumsum(kosambi(rf2.values)))
        
        
        
    })
    
    # browser()
    
    count = 1
    t5 = lapply(t4, function(x){
        x2 = as.data.frame(x)
        x2$chromo = names(t4)[count]
        count <<- count + 1
        colnames(x2) = c("cm", "chromo")
        x2
    })
    
    t5 = bind_rows(t5)
    
    t5$marker = colnames(geno.df)[3:ncol(geno.df)]
    
    t5
}





prepare.genetic.map.comparison.plots = function(list.of.geno.dfs, list.of.centimorgans, list.of.populations, labels.on){
    if(missing(list.of.populations)) list.of.populations = 1:length(list.of.geno.dfs)
    if(missing(labels.on)) labels.on = F
    geno.centi2 = bind_cols(list.of.centimorgans)
    
    geno.centi3 = melt(geno.centi2, "chromo")
    
    geno.centi4 = geno.centi3[which(geno.centi3$variable %in% c("cm", "cm1", "cm2", "cm3")), ]
    geno.centi4$marker = rep(geno.centi2$marker, 4)
    
    geno.centi5 = split.df.into.list.of.dfs.by.column(geno.centi4, "chromo", do.sort = F)
    
    geno.centi6 = lapply(geno.centi5, function(x){
        x$xaxis = rep(1:(nrow(x) / 4), 4)
        x$yaxis = c(rep(1, (nrow(x) / 4)), rep(2, (nrow(x) / 4)), rep(3, (nrow(x) / 4)), rep(4, (nrow(x) / 4)))
        x$value = as.numeric(x$value)
        x
    })
    
    
    geno.centi6 = lapply(geno.centi6, function(x){
        genovar1 = unique(x$variable)
        x$variable = as.character(x$variable)
        for(i in 1:length(genovar1)){
            x$variable[which(x$variable == genovar1[i])] = list.of.populations[i]
        }
        x$variable = factor(x$variable, levels = unique(x$variable))    
        x
    })
    
    geno.centi6 = lapply(geno.centi6, function(x){
        genovar1 = unique(x$yaxis)
        x$yaxis = as.character(x$yaxis)
        for(i in 1:length(genovar1)){
            x$yaxis[which(x$yaxis == genovar1[i])] = list.of.populations[i]
        }
        x$yaxis = factor(x$yaxis, levels = unique(x$yaxis))    
        x
    })
    
    # ggplot(geno.centi6[[1]], aes(x = yaxis, y = value, group = variable, color = variable)) + geom_point() + geom_line(aes(group = xaxis)) +
    #     theme_bw()
    
    
    ##some debugging
    # browser()
    # list.of.centimorgans[[1]][which(list.of.centimorgans[[1]]$chromo == "4B"), ]
    
    chromos1 = unique(geno.centi2$chromo)
    

    
    count = 1
    
    # browser()
    
    plots1 = lapply(geno.centi6, function(x){
        # ggplot(x, aes(x = yaxis, y = value, group = variable, color = variable)) + geom_line()    
        
        x$centromere = ""    
        
        if(x$chromo[1] == "7B"){
            
            x[which(x$xaxis == 19), ]$centromere = T
        }
        
        if(x$chromo[1] == "3B"){
            
            x[which(x$xaxis == 15), ]$centromere = T
        }
        
        
        blankx1 = x[1:(nrow(x) / 4), ]
        blankx1$yaxis = ""
        blankx1$value = NA
        
        x2 = bind_rows(blankx1, x)
        
        # x$yaxis = factor(x$yaxis, levels = unique(x$yaxis)) 
        buffer.val1 = 300
        # browser()
        
        # x2.2 = split(x2, x2$variable)
        # x2.3 = lapply(x2.2, function(x){
        #     x$value = ((x$value / max(na.omit(x$value))) * 100)
        #     x
        # })
        # 
        # x2.4 = bind_rows(x2.3)
        # if(x$chromo[1] == "3B") browser()
        
        if(labels.on == T){
            
            x.labeldf = x2[which(x2$yaxis == "10°C"), ]
            max.val1 = max(na.omit(x2$value))
            x.labeldf2 = data.frame(x.labeldf$marker, seq(0, max.val1, (max.val1 / nrow(x.labeldf)))[1:nrow(x.labeldf)])
            colnames(x.labeldf2) = c("marker", "y")
            x.labeldf2$x = 1
            x.labeldf2$variable = ""
            x.labeldf2$centromere = ""
            
            x.labeldf3 = x.labeldf2
            
            x.labeldf3.2 = x2[which(x2$yaxis == "10°C"), ]
            x.labeldf3$y = x.labeldf3.2$value
            
            x.labeldf2$group = 1:nrow(x.labeldf2)
            x.labeldf3$group = 1:nrow(x.labeldf3)
            x.labeldf3$x = 1.9 # line position data end
            
            x.linedf = bind_rows(x.labeldf3, x.labeldf2)
            x.linedf$x[x.linedf$x == 1] = 1.5 #line position label end
            
            p1 = ggplot(x2, aes(x = yaxis, y = value, group = variable, color = centromere, shape = centromere)) + geom_point() + geom_line(aes(group = xaxis)) +
                theme_classic(base_size = 22) + ggtitle(chromos1[count]) + labs(color = "Population") + 
                ylab("Genetic Distance (cM)") + xlab("Population") + theme(legend.position = "none") +
                scale_y_continuous(breaks = seq(0, (max(na.omit(x2$value)) + 40), 50)) + coord_cartesian(ylim = c(0, max(na.omit(x2$value)))) +
                scale_x_discrete(breaks = c("10°C", "14°C", "26°C", "28°C")) + scale_y_reverse(breaks = seq(0, 200, 50), labels = seq(0, 200, 50)) +
                geom_text(data = x.labeldf2, aes(x = x, y = y, label = marker), nudge_x = -0.05, color = "#000000", size = 5) + 
                geom_line(data = x.linedf, aes(x = x, y = y, group = group), color = "#ababab", size = 0.2)
            
            p1
        } else {
            p1 = ggplot(x2, aes(x = yaxis, y = value, group = variable, color = centromere, shape = centromere)) + geom_point() + geom_line(aes(group = xaxis)) +
                theme_classic(base_size = 22) + ggtitle(chromos1[count]) + labs(color = "Population") + 
                ylab("Genetic Distance (cM)") + xlab("Population") + theme(legend.position = "none", axis.title = element_text(size = 10)) +
                scale_y_continuous(breaks = seq(0, (max(na.omit(x2$value)) + 40), 50)) + #coord_cartesian(ylim = c(-buffer.val1, (max(na.omit(x2$value)) + buffer.val1))) +
                scale_x_discrete(breaks = c("10°C", "14°C", "26°C", "28°C"))
        }
        
        
        
        
        count <<- count + 1
        p1
    })
    
    plots1
}


ind.avg = function(recomb.df, manual.center, use.centromere, scale.by.arm){
    #compute the average distance of recombination events from the center of the chromosome for each individual
    #args:
    # recomb.df - a recombination data.frame produced using detect.recombination() that also contains physical distances of markers as a percentage in a phys.pos column
    # manual.center - integer, set the position of the center of the LG manually
    
    if(missing(manual.center)) manual.center = 0
    if(missing(use.centromere)) use.centromere = F
    if(missing(scale.by.arm)) scale.by.arm = F
    
    for.ind1 = lapply(unique(recomb.df$individual), function(x){
        recomb.df2 = recomb.df[which(recomb.df$individual == x), ]
        
        for.chromo1 = lapply(unique(recomb.df2$chromo), function(y){
            recomb.df3 = recomb.df2[which(recomb.df2$chromo == y), ]
            
            if(manual.center == 0){
                lg.center = round(axc.lg.lengths[which(names(axc.lg.lengths) == y)] / 2)    
            } else {
                lg.center = manual.center
            }
            
            double.events = which(recomb.df3$num.event == 2)
            single.events = which(recomb.df3$num.event == 1)
            if(length(double.events) > 0){
                recomb.df3 = recomb.df3[c(single.events, double.events, double.events), ]
                
            } 
            
            phys.midpoint = 50
            
            if(use.centromere == T){
                chromosomecounts.v2 = chromosomecounts[-22, ]
                chromosomecounts.v2$centro = centromeres$pos.bp
                chromosomecounts.v2$centro.per = (chromosomecounts.v2$centro / chromosomecounts.v2$V1) * 100
                phys.midpoint = chromosomecounts.v2[which(chromosomecounts.v2$V2 == paste0("chr", y)), ]$centro.per
                phys.midpoint.bp = chromosomecounts.v2[which(chromosomecounts.v2$V2 == paste0("chr", y)), ]$centro
                chromo.length = chromosomecounts.v2[which(chromosomecounts.v2$V2 == paste0("chr", y)), ]$V1
                
                arm1.length = chromo.length - phys.midpoint.bp
                arm2.length = chromo.length - arm1.length
                arms = c(arm1.length, arm2.length)
                arms = sort(arms, decreasing = T)
                names(arms) = c("long", "short")
            }
            
            if(y == "7D"){
                #the centromeric position for 7D is located at 53%, making the positions smaller than the centromere the long arm; whereas all the others are before 50%
                recomb.df3.short.arm = recomb.df3[which(recomb.df3$phys.pos > phys.midpoint), ]
                recomb.df3.long.arm = recomb.df3[which(recomb.df3$phys.pos <= phys.midpoint), ]    
            } else {
                recomb.df3.short.arm = recomb.df3[which(recomb.df3$phys.pos <= phys.midpoint), ]
                recomb.df3.long.arm = recomb.df3[which(recomb.df3$phys.pos > phys.midpoint), ]    
            }
            
            if(scale.by.arm == T){
                
                if(nrow(recomb.df3.long.arm) > 0){
                    recomb.df3.long.arm$phys.pos = recomb.df3.long.arm$phys.pos - phys.midpoint
                    l.arm.len = 100 - phys.midpoint
                    l.arm.scaler = 100 / l.arm.len
                    recomb.df3.long.arm$phys.pos = recomb.df3.long.arm$phys.pos * l.arm.scaler    
                }
                
                if(nrow(recomb.df3.short.arm) > 0){
                    s.arm.scaler = 100 / phys.midpoint
                    recomb.df3.short.arm$phys.pos = recomb.df3.short.arm$phys.pos * s.arm.scaler 
                    #need to reverse the positions so a decreasing percentage approaches the centromere, as in the long arm
                    recomb.df3.short.arm$phys.pos = abs(recomb.df3.short.arm$phys.pos - 100)
                }
                
            }
            
            # ((recomb.df3$phys.bp - phys.midpoint.bp) / arms[which.max(arms)]) * 100
            
            
            marker.dist = mean(abs(na.omit(recomb.df3$transition.coord) - lg.center))
            phys.dist = mean(abs(na.omit(recomb.df3$phys.pos) - phys.midpoint))
            # phys.dist.short = mean(abs(na.omit(recomb.df3.short.arm$phys.pos) - phys.midpoint))
            # phys.dist.long = mean(abs(na.omit(recomb.df3.long.arm$phys.pos) - phys.midpoint))
            
            phys.dist.short = mean(abs(na.omit(recomb.df3.short.arm$phys.pos)))
            phys.dist.long = mean(abs(na.omit(recomb.df3.long.arm$phys.pos)))
            
            
            # if(is.na(phys.dist)) browser()
            distances1 = c(marker.dist, phys.dist, y, phys.dist.short, phys.dist.long)
            
            names(distances1) = c("marker.dist", "phys.dist", "chromo", "phys.dist.short", "phys.dist.long")
            distances1
            
        })
        
        for.chromo2 = transpose(as.data.frame(for.chromo1))
        colnames(for.chromo2) = c("marker.dist", "phys.dist", "chromo", "phys.dist.short", "phys.dist.long")
        for.chromo2$individual = x
        for.chromo2
        
    })
    
    final.df = bind_rows(for.ind1)
    final.df$phys.dist = as.numeric(final.df$phys.dist)
    final.df$phys.dist.short = as.numeric(final.df$phys.dist.short)
    final.df$phys.dist.long = as.numeric(final.df$phys.dist.long)
    final.df$marker.dist = as.numeric(final.df$marker.dist)
    final.df
}

printplot = function(plotname, plotobject){
    pdf(paste0("rotation1scripts_v4/plots/recombination.distribution11072019/", plotname, ".pdf"), 6, 6)
    print(plotobject)
    dev.off()
}

exportresults = function(resultsdf, resultsname){
    write.csv(resultsdf, paste0("rotation1scripts_v4/original_data/recombination.distribution.chap.results/", resultsname, ".csv"), row.names = T)
}

comb.treatments = function(x, ntreatments){
    #adds treatment column to each dataframe in list of dataframes and combines all dataframes together. this enables linear model / anova 
    #(use treatment as independent variable)
    if(missing(ntreatments)) ntreatments = 4
    for(i in 1:ntreatments){
        x[[i]]$treatment = i
    }
    
    bind_rows(x)
}

report_kruskal = function(kruskal.sum){
    chi1 = round(kruskal.sum$statistic, digits = 2)
    df1 = round(kruskal.sum$parameter, digits = 2)
    p.val = round(kruskal.sum$p.value, digits = 5)
    paste0("chi-square = ", chi1, ", d.f. = ", df1, ", p = ", p.val)
}

report_t = function(t.test1){
    t1 = round(t.test1$statistic, digits = 2)
    df1 = round(t.test1$parameter, digits = 2)
    p.val = round(t.test1$p.value, digits = 5)
    paste0("t = ", t1, ", d.f. = ", df1, ", p = ", p.val)
}

report_mann = function(t.test1){
    t1 = round(t.test1$statistic, digits = 2)
    p.val = round(t.test1$p.value, digits = 5)
    paste0("W = ", t1, ", p = ", p.val)
}


get.phys.for.chromos = function(geno.df, interpolate.na, include.chromo, marker.format){
    if(missing(interpolate.na)) interpolate.na = T
    if(missing(include.chromo)) include.chromo = F
    if(missing(marker.format)) marker.format = "."
    
    chromosomes1 = unique(as.character(geno.df[1, ]))[-(1:2)]
    ini.axc = convert.to.character.data.frame(as.data.frame(t(geno.df)))
    colnames(ini.axc)[1] = "marker"
    
    ini.axc2 = ini.axc[-(1:2), ]
    ini.axc3 = split.df.into.list.of.dfs.by.column(ini.axc2, "marker", F)
    # browser()
    ini.axc.phys = lapply(ini.axc3, function(x){
        
        
        # browser()
        phys1 = genetictophysical(unique(x$marker), rownames(x), marker.format)
        phys.bp = attr(phys1, "bp")
        # browser()
        
        if(interpolate.na == T){
            phys1 = na.approx(phys1, na.rm = F)
        }
        
        # rownames(x)
        data.frame(rownames(x), phys1, phys.bp)
    })
    
    if(include.chromo == T){
        for(i in 1:length(ini.axc.phys)){
            ini.axc.phys[[i]]$chromo = chromosomes1[i]
        }    
    }
    
    ini.axc.phys2 = bind_rows(ini.axc.phys)
    colnames(ini.axc.phys2)[1:3] = c("marker", "phys.pos", "phys.bp")
    
    ini.axc.phys2$marker = ini.axc.phys2$marker
    ini.axc.phys2
}

#### SIMULATION ####

process.sim.recomb = function(geno, sample.size){
    rd1 = geno
    rd1 = as.data.frame(rbind(rd1[1, ], rd1))
    rd1 = convert.to.character.data.frame(rd1)
    rd1[1, ] = 1
    
    tobind1 = rd1[, 1:2]
    colnames(tobind1) = c("tobind1", "tobind2")
    
    rd1 = cbind(tobind1, rd1)
    rd1 = as.data.frame(rd1)
    rd1 = convert.to.character.data.frame(rd1)
    rd1[, 2] = rownames(rd1)
    rd1[, 1] = c("", 1:(nrow(rd1) - 1))
    
    rd1[1, 2] = ""
    rd1[2, 2] = "F2_1_1"
    
    rd1 = rd1[c(1, sample(2:nrow(rd1), sample.size)), ]
    x = detect.recombination(rd1)
    phys.pos.before = axp.3b.map[match(x$probe.before.transition, axp.3b.map$marker), ]$phys
    phys.pos.after = axp.3b.map[match(x$probe.after.transition, axp.3b.map$marker), ]$phys
    phys.pos.midpoint = (phys.pos.before + phys.pos.after) / 2
    x$phys.pos = phys.pos.midpoint
    x$chromo = "3B"
    x
}

convert.sim.to.rqtl.ske.format = function(sim.geno.df){
    rd1 = sim.geno.df
    rd1 = as.data.frame(rbind(rd1[1, ], rd1))
    rd1 = convert.to.character.data.frame(rd1)
    rd1[1, ] = 1
    
    rd1 = cbind(rd1[, 1:2], rd1)
    rd1 = as.data.frame(rd1)
    rd1 = convert.to.character.data.frame(rd1)
    rd1[, 2] = rownames(rd1)
    rd1[, 1] = c("", 1:(nrow(rd1) - 1))
    colnames(rd1)[1:2] = c("t1", "t2")
    
    rd1[1, 2] = ""
    rd1[2, 2] = "F2_1_1"
    rd1
}

#### MORE FUNCTIONS #### 

Ps_Visualization2 <- function(pidFile, output.File="Clusters.pdf", output.dir=NULL, summaryFile="AxiomGT1.summary.txt", callFile=NULL, posteriorFile=NULL, multiallele.posteriorFile=NULL, confidenceFile=NULL, sampleFile=NULL, labelsFile=NULL, specialSNPsFile=NULL, reportFile=NULL, temp.dir="Temp", keep.temp.dir=TRUE, use.temp.dir=FALSE, plot.intensity=TRUE, plot.genotype=TRUE, plot.calls=TRUE, plotNoCalls=TRUE, plot.confs=FALSE, plot.ref=TRUE, refFile=NULL, refNoCalls=TRUE, plot.prior=FALSE, plot.prior.ref=FALSE, priorFile=NULL, multiallele.priorFile=NULL, plot.posterior=TRUE, plot.posterior.ref=FALSE, multiallele.order="decreasing", multiallele.all=TRUE, transform="MvA", K=2, autotetraploid=FALSE, confidences.cutoff=NULL, confs.thresh=c(0,0.01,0.05,0.10,0.15), confs.col=c("green","gold","darkorange1","firebrick1"), col.AA=NULL, col.AB=NULL, col.BB=NULL, col.A=NULL, col.B=NULL, col.CN0=NULL, col.OTV=NULL, col.NC=NULL, col.OTV1=NULL, col.NC1=NULL, xlim.intensity=NULL, ylim.intensity=NULL, xlim.genotype=NULL, ylim.genotype=NULL, x.axis.mirror=TRUE, vertical.line=FALSE, max.num.SNP.draw=5000, cex.main=1.2, cex.lab=1.2, num.rows=NULL, num.cols=NULL, plot.width=15.5, plot.height=12, plot.type="pdf", plot.all=TRUE, batch=F){
    
    ### checking if perl is installed on the system
    hold <- system("perl -e1")
    if(hold>=1){
        stop("perl is not installed. Please install perl and retry.\n\n")
    }
    
    ### set the global option to keep characters strings instead of changing to factors
    options(stringsAsFactors=FALSE) 
    
    ### check that the required arguments have been supplied
    ### do all checks that involve a "stop" statement first before reading anything in
    if(missing(pidFile)){
        stop("\n\npidFile argument is missing. Please supply a SNP list.")    
    }
    if(length(pidFile)>1){
        stop("\n\npidFile argument should be only one file. Please check and retry.")
    }
    pid <- read.delim(pidFile, header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")[,1]
    if(length(pid)==0){
        stop("\n\nThe pidFile is empty.Please check and retry.")
    }
    
    
    ### check that the output file type is one of the three options
    if(!plot.type%in%c("pdf","png","svg")){
        stop("\n\nplot.type must be either pdf, png, or svg. Please check and retry.")
    }
    if(sum(strsplit(output.File,".p")[[1]][length(strsplit(output.File,".p")[[1]])]%in%c("df","ng"),strsplit(output.File,".s")[[1]][length(strsplit(output.File,".s")[[1]])]%in%"vg")!=1){
        stop("\n\noutput.File name must end in pdf, png, or svg. Please check and retry.")
    }
    if(plot.type=="pdf" && strsplit(output.File,".p")[[1]][length(strsplit(output.File,".p")[[1]])]!="df"){
        stop("\n\nplot.type is pdf but output.File does not end in pdf. Please check and retry.")
    }
    if(plot.type=="png" && strsplit(output.File,".p")[[1]][length(strsplit(output.File,".p")[[1]])]!="ng"){
        stop("\n\nplot.type is png but output.File does not end in png. Please check and retry.")
    }
    if(plot.type=="svg" && strsplit(output.File,".s")[[1]][length(strsplit(output.File,".s")[[1]])]!="vg"){
        stop("\n\nplot.type is svg but output.File does not end in svg. Please check and retry.")
    }
    
    
    ### check that the plotting options are TRUE/FALSE only
    if(!(plot.all%in%c(TRUE,FALSE))){
        stop("\n\nplot.all must be either TRUE or FALSE.    Please check and retry.")
    } 
    if(!(plot.intensity%in%c(TRUE,FALSE))){
        stop("\n\nplot.intensity must be either TRUE or FALSE.    Please check and retry.")
    }    
    if(!(plot.genotype%in%c(TRUE,FALSE))){
        stop("\n\nplot.genotype must be either TRUE or FALSE.    Please check and retry.")
    }
    if(!(plot.calls%in%c(TRUE,FALSE))){
        stop("\n\nplot.calls must be either TRUE or FALSE.    Please check and retry.")
    }
    if(!(plot.confs%in%c(TRUE,FALSE))){
        stop("\n\nplot.confs must be either TRUE or FALSE.    Please check and retry.")
    }
    if(!(plot.ref%in%c(TRUE,FALSE))){
        stop("\n\nplot.ref must be either TRUE or FALSE.    Please check and retry.")
    }
    if(!(plot.prior%in%c(TRUE,FALSE))){
        stop("\n\nplot.prior must be either TRUE or FALSE.    Please check and retry.")
    }
    if(!(plot.posterior%in%c(TRUE,FALSE))){
        stop("\n\nplot.posterior must be either TRUE or FALSE.    Please check and retry.")
    }
    if(plot.confs==TRUE){
        if(!(length(confs.thresh)==5)){
            stop("\n\nconfs.thresh must be 5 numbers.    Please check and retry.")
        }
        if(!(is.numeric(confs.thresh))){
            stop("\n\nconfs.thresh must be numeric.    Please check and retry.")
        }
        if(!(length(confs.col)==4)){
            stop("\n\nconfs.col must be 4 colors.    Please check and retry.")
        }
    }
    
    ### check that some plotting options have been selected
    if(plot.intensity==FALSE && plot.genotype==FALSE){
        stop("\n\nPlease select intensities or genotypes or both to be plotted.\n\n")
    }
    if(plot.confs==FALSE && plot.calls==FALSE){
        stop("\n\nPlease select confidences or calls or both to be plotted.\n\n")
    }
    
    ### set up plotting order using plotting type options
    plot.order <- NA
    if(plot.intensity==TRUE){
        if(plot.calls==TRUE){
            plot.order <- c(plot.order,"intensity.calls")
        }
        if(plot.ref==TRUE){
            plot.order <- c(plot.order,"intensity.ref")
        }
        if(plot.confs==TRUE){
            plot.order <- c(plot.order,"intensity.confs")
        }
    }
    if(plot.genotype==TRUE){
        if(plot.calls==TRUE){
            plot.order <- c(plot.order,"genotype.calls")
        }
        if(plot.ref==TRUE){
            plot.order <- c(plot.order,"genotype.ref")
        }
        if(plot.confs==TRUE){
            plot.order <- c(plot.order,"genotype.confs")
        }
    }
    plot.order <- plot.order[-1]
    if(length(plot.order)==0){
        stop("\n\nNo plotting options selected.    Please check that one of intensities, genotypes, calls, confidences, or references has been selected.")
    }
    
    ### check that the multiallele plotting order is either alphabetical or in order of which allele is most common
    if(!(multiallele.order%in%c("decreasing","alphabetical"))){
        stop("\n\nmultiallele.order must be decreasing or alphabetical.    Please check and retry.\n\n")
    }
    ### check that option for only plotting all samples with multiallele probesets is either TRUE or FALSE
    if(!multiallele.all%in%c(TRUE,FALSE)){
        stop("\n\nmultiallele.all must be either TRUE or FALSE. Please check and retry.\n\n")
    }
    
    ### check that the transformation is one of the 7 options
    if(!(transform%in%c("log2","MvA","CES","CSS","RvT","MvA-sqrt","polar"))){
        stop("\n\ntransform must be MvA, MvA-sqrt, log2, CES, CSS, RvT, or polar.    Please check and retry.\n\n")
    }
    
    ### check that the K value supplied for the non-MvA transforms is in the correct range
    if((transform%in%c("CES","CSS"))){
        if(!(K%in%c(1,2,3,4,5))){
            stop("\n\nK must be an integer from 1 to 5.    Please check and retry.\n\n")
        }
        K.val <- K
    }
    
    ### check that the autotetraploid option is TRUE or FALSE
    if(!(autotetraploid%in%c(TRUE,FALSE))){
        stop("\n\nautotetraploid must be either TRUE or FALSE.    Please check and retry.")
    }
    
    ### check that the user-supplied x and y limits make sense
    if(!is.null(xlim.intensity)){
        if(length(xlim.intensity)!=2){
            stop("\n\nxlim.intensity must be two numbers.    Please check and retry.\n\n")
        }
        if(xlim.intensity[1]>xlim.intensity[2]){
            stop("\n\nxlim.intensity[1] must be smaller than xlim.intensity[2].    Please check and retry.\n\n")
        }
    }
    if(!is.null(ylim.intensity)){
        if(length(ylim.intensity)!=2){
            stop("\n\nylim.intensity must be two numbers.    Please check and retry.\n\n")
        }
        if(ylim.intensity[1]>ylim.intensity[2]){
            stop("\n\nylim.intensity[1] must be smaller than ylim.intensity[2].    Please check and retry.\n\n")
        }
    }
    if(!is.null(xlim.genotype)){
        if(length(xlim.genotype)!=2){
            stop("\n\nxlim.genotype must be two numbers.    Please check and retry.\n\n")
        }
        if(xlim.genotype[1]>xlim.genotype[2]){
            stop("\n\nxlim.genotype[1] must be smaller than xlim.genotype[2].    Please check and retry.\n\n")
        }
    }
    if(!is.null(ylim.genotype)){
        if(length(ylim.genotype)!=2){
            stop("\n\nylim.genotype must be two numbers.    Please check and retry.\n\n")
        }
        if(ylim.genotype[1]>ylim.genotype[2]){
            stop("\n\nylim.genotype[1] must be smaller than ylim.genotype[2].    Please check and retry.\n\n")
        }
    }
    
    ### check that the X axis mirror option is either TRUE or FALSE
    if(!(x.axis.mirror%in%c(FALSE,TRUE))){
        stop("\n\nx.axis.mirror must be TRUE or FALSE.    Please check and retry.\n\n")
    }
    
    ### check that the vertical line option is either TRUE or FALSE
    if(!(vertical.line%in%c(FALSE,TRUE))){
        stop("\n\nvertical.line must be TRUE or FALSE.    Please check and retry.\n\n")
    }
    
    
    ### if the special SNPs file is not supplied, then only autosomal probesets should be supplied
    ### specialSNPsFile should only be one file, not one per batch
    if(is.null(specialSNPsFile) && !is.null(reportFile)){
        stop("\n\nSpecial SNPs file is missing.    Please supply a special SNPs file to plot gender separately for non-PAR X and Y SNPs.\n\n")
    }
    if(!is.null(specialSNPsFile) && length(specialSNPsFile)>1){
        stop("\n\nThere should only be one special SNPs file provided. Please check and retry.\n\n")
    }
    if(!(is.null(specialSNPsFile)) && is.null(reportFile)){
        stop("\n\nReport file is missing.    Please supply a report file to plot gender separately for non-PAR X and Y SNPs.\n\n")
    }
    
    ### if plotting confidences, must have a confidences file
    if((plot.confs==TRUE) && (is.null(confidenceFile))){
        stop("\n\nConfidence plotting is turned on.    Please supply a confidences file and retry.")
    }
    ### can't use confidences cutoff without a confidences file
    if((is.null(confidenceFile)) && (!(is.null(confidences.cutoff)))){
        stop("\n\nThe confidences cutoff value can only be used with a confidences file.    Please supply the name of a confidences file and retry.")
    }
    ### if a confidences file is supplied but a confidences.cutoff is not, 
    ### then a default of 1 (i.e. no changes) is used 
    if((!(is.null(confidenceFile))) && ((is.null(confidences.cutoff))) && plot.confs==TRUE){
        confidences.cutoff <- 1
        cat("\nconfidences.cutoff has not been supplied. No changes will be made to the calls based on confidence values.\n\n")
    }
    ### done with all basic "stop" checks
    
    ### fix temp.dir name so that it does not have a "/" at the end so perl and R 
    ### can make and write to the temp.dir regardless of OS
    for(i in 1:length(temp.dir)){
        if(!is.na(temp.dir[i])){
            hold <- strsplit(temp.dir[i],"")[[1]]
            while((hold[length(hold)]=="/")==TRUE){
                hold <- hold[-length(hold)]
            }
            temp.dir[i] <- paste(hold,collapse="")
        } else{
            temp.dir[i] <- "."
        }
    }
    ### fix output.dir name so that it does not have a "/" at the end so the complete output file name will be correct
    if(!is.null(output.dir)){
        for(i in 1:length(output.dir)){
            if(!is.na(output.dir[i])){
                hold <- strsplit(output.dir[i],"")[[1]]
                while((hold[length(hold)]=="/")==TRUE){
                    hold <- hold[-length(hold)]
                }
                output.dir[i] <- paste(hold,collapse="")
            } else{
                output.dir[i] <- "."
            }
        }
    }
    
    ### Ps_Vis can be run with a list of summary files, call files, and what temporary directories to use
    ### the input lists must be of the same length
    ### non-batch potting usage is a list of length 1 (i.e. only one file each)
    ### if (MA) posterior, (MA) priors, and/or report files are provided, they must be of the same length too
    l <- length(summaryFile)
    l1 <- length(callFile)
    l2 <- length(temp.dir)
    if(l!=l1){
        stop("\n\nList of summary files and calls files must have the same number of elements.
                 Please check the length of the lists and retry.")
    }
    if(l!=l2){
        stop("\n\nList of call files and temporary directories must have the same number of elements.
                 Please check the length of the lists and retry.")
    }
    if(!is.null(output.dir)){
        if(l!=length(output.dir) && plot.all==FALSE){
            stop("\n\nList of output directories and summary files must have the same number of elements. 
                     Please check the length of the lists and retry.")
        }
        }
    if(plot.ref==TRUE){
        if(l!=length(refFile) && !is.null(refFile)){
            stop("\n\nList of summary files and reference files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        }
    if(plot.posterior==TRUE){
        if(l!=length(posteriorFile) && !is.null(posteriorFile)){
            stop("\n\nList of summary files and posteriors files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        if(l!=length(multiallele.posteriorFile) && !is.null(multiallele.posteriorFile)){
            stop("\n\nList of summary files and multiallelic posteriors files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        }
    if(plot.prior==TRUE){
        if(l!=length(priorFile) && !is.null(priorFile)){
            stop("\n\nList of summary files and priors files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        if(l!=length(multiallele.priorFile) && !is.null(multiallele.priorFile)){
            stop("\n\nList of summary files and multiallelic priors files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        }
    if(!is.null(reportFile)){
        if(l!=length(reportFile)){
            stop("\n\nList of summary files and report files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        }
    if(!is.null(labelsFile)){
        if(l!=length(labelsFile)){
            stop("\n\nList of summary files and labels files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        }
    if(!is.null(sampleFile)){
        if(l!=length(sampleFile)){
            stop("\n\nList of summary files and sample files must have the same number of elements.
                     Please check the length of the lists and retry.")
        }
        }
    
    ### if more than 1 batch has been supplied, then batch==TRUE
    if(l>1){
        batch <- TRUE
    }
    
    ### pull the calls from the calls file APT header
    ### if a line starts with "#%call-code-", check if it's an older file without copy number aware genotyping
    ### check if the "use copy num calls" flag has been set to FALSE and if so, update the posteriors file when necessary
    call.codes <- rep(NA,3)
    use_copy_num_calls <- TRUE
    for(j in 1:length(callFile)){
        i <- 1
        if(tools::file_ext(callFile[j])%in%c("gz","gzip")){
            calls.comments <- read.table(gzfile(callFile[j]),stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
        } else{
            calls.comments <- read.table(callFile[j],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
        }
        while(length(calls.comments)>0){
            ### stop when you get to the end of the APT header and hit the calls header which starts with probeset_id
            if(length(grep("probeset_id",calls.comments[[1]]))>0){
                break
            }
            ### pull all of the call codes out of the headers
            if(length(which(unlist(lapply(strsplit(calls.comments[[1]],"call-code-"),length))==2))>0){
                call.codes <- rbind(call.codes,strsplit(strsplit(strsplit(calls.comments[[1]],"call-code-")[[1]][2],"=")[[1]][2],":")[[1]])
            }
            i <- i+1
            if(tools::file_ext(callFile[j])%in%c("gz","gzip")){
                calls.comments <- read.table(gzfile(callFile[j]),stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
            } else{
                calls.comments <- read.table(callFile[j],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
            }
        }
        if(!is.null(dim(call.codes))){
            call.codes <- call.codes[!duplicated(call.codes[,1]),]
        }
    }
    ### check if the call headers had the call codes
    ### if not, set them to the original diploid values
    if(is.null(dim(call.codes))){
        call.codes <- c("NoCall","-1","2")
        call.codes <- rbind(call.codes,c("AA","0","2"))
        call.codes <- rbind(call.codes,c("AB","1","2"))
        call.codes <- rbind(call.codes,c("BB","2","2"))
        call.codes <- rbind(call.codes,c("OTV","-2","2"))
        call.codes <- rbind(call.codes,c("ZeroCN","3","0"))
        call.codes <- rbind(call.codes,c("A","4","1"))
        call.codes <- rbind(call.codes,c("B","5","1"))
        rownames(call.codes)[1] <- ""
    } else{
        call.codes <- call.codes[-1,]
    }
    ### check if any batches have used the same call codes for different genotype call assignments
    if(length(unique(call.codes[,2]))<length(call.codes[,2])){
        stop(c("\n\nThe same call code has been used for different genotype calls. The following calls have been assigned multiple times:\n",paste(unique(call.codes[duplicated(call.codes[,2]),2])," "),"\nIf multiple batches have been provided, call codes may be different across calls files.    Please check and retry."))
    }
    
    ### if there's reference calls, check that they have the same call codes as the calls file
    ### if not, map the two sets of call codes and store them
    if(!is.null(refFile)){
        calls.ref.match <- rep(NA,3)
        ref.codes <- rep(NA,3)
        for(j in 1:length(refFile)){
            if(!is.null(refFile[j])){
                i <- 1
                if(tools::file_ext(refFile[j])%in%c("gz","gzip")){
                    ref.comments <- read.table(gzfile(refFile[j]),stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                } else{
                    ref.comments <- read.table(refFile[j],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                }
                while(length(ref.comments)>0){
                    ### stop when you get to the end of the APT header and hit the calls header which starts with probeset_id
                    if(length(grep("probeset_id",ref.comments[[1]]))>0){
                        break
                    }
                    ### pull all of the call codes out of the headers
                    if(length(which(unlist(lapply(strsplit(ref.comments[[1]],"call-code-"),length))==2))>0){
                        ref.codes <- rbind(ref.codes,strsplit(strsplit(strsplit(ref.comments[[1]],"call-code-")[[1]][2],"=")[[1]][2],":")[[1]])
                    }
                    i <- i+1
                    if(tools::file_ext(refFile[j])%in%c("gz","gzip")){
                        ref.comments <- read.table(gzfile(refFile[j]),stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                    } else{
                        ref.comments <- read.table(refFile[j],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                    }
                }
                if(!is.null(dim(ref.codes))){
                    ref.codes <- ref.codes[!duplicated(ref.codes[,1]),]
                }
            }
        }
        if(!is.null(dim(ref.codes))){
            ### check if there are either codes in the reference file that don't appear in the call codes or 
            ### if there are mismatches between the calls and references
            ref.codes <- ref.codes[-1,]
            for(j in 1:dim(ref.codes)[1]){
                if(sum(ref.codes[j,1]%in%call.codes[,1])==0){
                    call.codes <- rbind(call.codes,ref.codes[j,])
                } else if(ref.codes[j,1]%in%call.codes[,1] && ref.codes[j,2]!=call.codes[call.codes[,1]==ref.codes[j,1],2]){
                    calls.ref.match <- rbind(calls.ref.match,c(ref.codes[j,1],call.codes[call.codes[,1]==ref.codes[j,1],2],ref.codes[j,2]))
                }
            }
            if(!is.null(dim(calls.ref.match))){
                calls.ref.match <- calls.ref.match[-1,]
            }
            if(is.null(dim(calls.ref.match))){
                calls.ref.match <- as.data.frame(calls.ref.match)
            }
            call.codes <- call.codes[order(as.numeric(call.codes[,2])),]
        }
    }
    
    
    ### set the colors
    ### update user-supplied values in each probeset so that multiallelic probesets keep the correct colors
    geno.col <- cbind((-4:5),c("#ABFFE3","#D9D9D9","#66FFCC","#BEBEBE","#FF0000","#FFD700","#0000FF","#FFFFFF","#FFCACA","#9EEFE3"),c(23,22,23,22,24,21,25,5,24,25))
    ### call codes 3, 4, and 5 may be used in dynamic call assignment
    ### remove calls 3, 4, and/or 5 from geno.col if they appear in the call codes assignments with calls that are different from A, B, CN0
    ### if 4 and 5 don't appear at all then add in the default colors in case haploid posteriors are provided but haploid calls aren't
    if("3"%in%call.codes[,2]){
        if(call.codes[call.codes[,2]=="3",1]!="ZeroCN"){
            geno.col <- geno.col[-which(geno.col[,1]==3),]
        }
    }
    if("4"%in%call.codes[,2]){
        if(call.codes[call.codes[,2]=="4",1]!="A"){
            geno.col <- geno.col[-which(geno.col[,1]==4),]
        }
    } else{
        if(call.codes[call.codes[,2]=="4",1]=="A" && length(which(geno.col[,1]==4))==0){
            geno.col <- rbind(geno.col,c("4","#FFCACA","24"))
        }
    }
    if("5"%in%call.codes[,2]){
        if(call.codes[call.codes[,2]=="5",1]!="B"){
            geno.col <- geno.col[-which(geno.col[,1]==5),]
        }
    } else{
        if(call.codes[call.codes[,2]=="5",1]=="B" && length(which(geno.col[,1]==5))==0){
            geno.col <- rbind(geno.col,c("5","#9EEFE3","25"))
        }
    }
    
    ### use call codes to see if there's a possibility that multiallelic calls will appear
    multiallele <- FALSE
    if(dim(call.codes[-which(call.codes[,1]%in%c("OTV","NoCall","OTV_1","NoCall_1","ZeroCN","AA","AB","BB","A","B")),])[1]>0){
        multiallele <- TRUE
    }
    
    ### set the multi-allele colors and shapes
    ### users cannot change the multi-allele colors and shapes
    multi.col <- c("#00FF00","#FF00FF","#000080","#006400","#BA55D3","#B8860B","#FF1493","#556B2F","#FF8C00","#DC143C","#000000","#8B4513")
    multi.pch <- c(24,21,25,23,22,24,21,25,23,22,24,21,25,23,22)
    
    ### read in the multi posteriors headers
    mean.codes <- rep(NA,2)
    var.codes <- rep(NA,2)
    if(!is.null(multiallele.posteriorFile)){
        for(j in 1:length(multiallele.posteriorFile)){
            i <- 1
            if(!is.na(multiallele.posteriorFile[j])){
                if(tools::file_ext(multiallele.posteriorFile[j])%in%c("gz","gzip")){
                    ma.comments <- read.table(gzfile(multiallele.posteriorFile[j]),stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                } else{
                    ma.comments <- read.table(multiallele.posteriorFile[j],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                }
                while(length(ma.comments)>0){
                    ### stop when you get to the end of the APT header and hit the MA posteriors header which starts with probeset_id
                    if(sum(ma.comments[1]=="probeset_id")>0){
                        break
                    }
                    ### pull all of the allele and mean/variance codes out of the headers
                    if(length(unlist(strsplit(ma.comments[[1]],"data-order-mean-nalleles-")))==2){
                        mean.codes <- rbind(mean.codes,c(strsplit(strsplit(ma.comments[[1]],"data-order-mean-nalleles-")[[1]][2],"=")[[1]][1],strsplit(strsplit(ma.comments[[1]],"data-order-mean-nalleles-")[[1]][2],"=")[[1]][2]))
                    }
                    if(length(unlist(strsplit(ma.comments[[1]],"data-order-covariance-nalleles-")))==2){
                        var.codes <- rbind(var.codes,c(strsplit(strsplit(ma.comments[[1]],"data-order-covariance-nalleles-")[[1]][2],"=")[[1]][1],strsplit(strsplit(ma.comments[[1]],"data-order-covariance-nalleles-")[[1]][2],"=")[[1]][2]))
                    }
                    i <- i+1
                    if(tools::file_ext(multiallele.posteriorFile[j])%in%c("gz","gzip")){
                        ma.comments <- read.table(gzfile(multiallele.posteriorFile[j]),stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                    } else{
                        ma.comments <- read.table(multiallele.posteriorFile[j],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=(i-1),check.names=F,comment.char = "")
                    }
                }
            }
        }
        if(is.null(dim(mean.codes))){
            cat("\nDefault values for the multiallelic mean order will be used.\n\n")
            mean.codes <- cbind(2:12,rbind(c("A,B"),c("A,B,C"),c("A,B,C,D"),c("A,B,C,D,E"),c("A,B,C,D,E,F"),c("A,B,C,D,E,F,G"),c("A,B,C,D,E,F,G,H"),c("A,B,C,D,E,F,G,H,I"),c("A,B,C,D,E,F,G,H,I,J"),c("A,B,C,D,E,F,G,H,I,J,K"),c("A,B,C,D,E,F,G,H,I,J,K,L")))
        }
        if(is.null(dim(var.codes))){
            cat("\nDefault values for the multiallelic covariance order will be used.\n\n")
            var.codes <- cbind(2:12,rbind(c("varA,covAB,varB"),c("varA,covAB,covAC,varB,covBC,varC"),c("varA,covAB,covAC,covAD,varB,covBC,covBD,varC,covCD,varD"),c("varA,covAB,covAC,covAD,covAE,varB,covBC,covBD,covBE,varC,covCD,covCE,varD,covDE,varE"),c("varA,covAB,covAC,covAD,covAE,covAF,varB,covBC,covBD,covBE,covBF,varC,covCD,covCE,covCF,varD,covDE,covDF,varE,covEF,varF"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,varB,covBC,covBD,covBE,covBF,covBG,varC,covCD,covCE,covCF,covCG,varD,covDE,covDF,covDG,varE,covEF,covEG,varF,covFG,varG"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,varB,covBC,covBD,covBE,covBF,covBG,covBH,varC,covCD,covCE,covCF,covCG,covCH,varD,covDE,covDF,covDG,covDH,varE,covEF,covEG,covEH,varF,covFG,covFH,varG,covGH,varH"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,varC,covCD,covCE,covCF,covCG,covCH,covCI,varD,covDE,covDF,covDG,covDH,covDI,varE,covEF,covEG,covEH,covEI,varF,covFG,covFH,covFI,varG,covGH,covGI,varH,covHI,varI"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,varD,covDE,covDF,covDG,covDH,covDI,covDJ,varE,covEF,covEG,covEH,covEI,covEJ,varF,covFG,covFH,covFI,covFJ,varG,covGH,covGI,covGJ,varH,covHI,covHJ,varI,covIJ,varJ"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,varE,covEF,covEG,covEH,covEI,covEJ,covEK,varF,covFG,covFH,covFI,covFJ,covFK,varG,covGH,covGI,covGJ,covGK,varH,covHI,covHJ,covHK,varI,covIJ,covIK,varJ,covJK,varK"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,covAL,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,covBL,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,covCL,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,covDL,varE,covEF,covEG,covEH,covEI,covEJ,covEK,covEL,varF,covFG,covFH,covFI,covFJ,covFK,covFL,varG,covGH,covGI,covGJ,covGK,covGL,varH,covHI,covHJ,covHK,covHL,varI,covIJ,covIK,covIL,varJ,covJK,covJL,varK,covKL,varL")))
        }
        mean.codes <- mean.codes[!duplicated(mean.codes[,1]),]
        if(dim(mean.codes)[1]>1){
            mean.codes <- mean.codes[-1,]
            ### add in the defaults for any allele numbers not in the headers
            hold <-    cbind(2:12,rbind(c("A,B"),c("A,B,C"),c("A,B,C,D"),c("A,B,C,D,E"),c("A,B,C,D,E,F"),c("A,B,C,D,E,F,G"),c("A,B,C,D,E,F,G,H"),c("A,B,C,D,E,F,G,H,I"),c("A,B,C,D,E,F,G,H,I,J"),c("A,B,C,D,E,F,G,H,I,J,K"),c("A,B,C,D,E,F,G,H,I,J,K,L")))
            mean.codes <- rbind(mean.codes,hold[-which(hold[,1]%in%mean.codes[,1]),])
        } else{
            cat("\nDefault values for the multiallelic mean order will be used.\n\n")
            mean.codes <- cbind(2:12,rbind(c("A,B"),c("A,B,C"),c("A,B,C,D"),c("A,B,C,D,E"),c("A,B,C,D,E,F"),c("A,B,C,D,E,F,G"),c("A,B,C,D,E,F,G,H"),c("A,B,C,D,E,F,G,H,I"),c("A,B,C,D,E,F,G,H,I,J"),c("A,B,C,D,E,F,G,H,I,J,K"),c("A,B,C,D,E,F,G,H,I,J,K,L")))
        }
        ### dimension handling if mean.codes is only one line
        if(is.null(dim(mean.codes))){
            mean.codes <- as.matrix(t(mean.codes))
        }
        var.codes <- var.codes[!duplicated(var.codes[,1]),]
        if(dim(var.codes)[1]>1){
            var.codes <- var.codes[-1,]
            ### add in the defaults for any allele numbers not in the headers
            hold <- cbind(2:12,rbind(c("varA,covAB,varB"),c("varA,covAB,covAC,varB,covBC,varC"),c("varA,covAB,covAC,covAD,varB,covBC,covBD,varC,covCD,varD"),c("varA,covAB,covAC,covAD,covAE,varB,covBC,covBD,covBE,varC,covCD,covCE,varD,covDE,varE"),c("varA,covAB,covAC,covAD,covAE,covAF,varB,covBC,covBD,covBE,covBF,varC,covCD,covCE,covCF,varD,covDE,covDF,varE,covEF,varF"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,varB,covBC,covBD,covBE,covBF,covBG,varC,covCD,covCE,covCF,covCG,varD,covDE,covDF,covDG,varE,covEF,covEG,varF,covFG,varG"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,varB,covBC,covBD,covBE,covBF,covBG,covBH,varC,covCD,covCE,covCF,covCG,covCH,varD,covDE,covDF,covDG,covDH,varE,covEF,covEG,covEH,varF,covFG,covFH,varG,covGH,varH"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,varC,covCD,covCE,covCF,covCG,covCH,covCI,varD,covDE,covDF,covDG,covDH,covDI,varE,covEF,covEG,covEH,covEI,varF,covFG,covFH,covFI,varG,covGH,covGI,varH,covHI,varI"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,varD,covDE,covDF,covDG,covDH,covDI,covDJ,varE,covEF,covEG,covEH,covEI,covEJ,varF,covFG,covFH,covFI,covFJ,varG,covGH,covGI,covGJ,varH,covHI,covHJ,varI,covIJ,varJ"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,varE,covEF,covEG,covEH,covEI,covEJ,covEK,varF,covFG,covFH,covFI,covFJ,covFK,varG,covGH,covGI,covGJ,covGK,varH,covHI,covHJ,covHK,varI,covIJ,covIK,varJ,covJK,varK"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,covAL,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,covBL,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,covCL,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,covDL,varE,covEF,covEG,covEH,covEI,covEJ,covEK,covEL,varF,covFG,covFH,covFI,covFJ,covFK,covFL,varG,covGH,covGI,covGJ,covGK,covGL,varH,covHI,covHJ,covHK,covHL,varI,covIJ,covIK,covIL,varJ,covJK,covJL,varK,covKL,varL")))
            var.codes <- rbind(var.codes,hold[-which(hold[,1]%in%var.codes[,1]),])
        } else{
            cat("\nDefault values for the multiallelic covariance order will be used.\n\n")
            var.codes <- cbind(2:12,rbind(c("varA,covAB,varB"),c("varA,covAB,covAC,varB,covBC,varC"),c("varA,covAB,covAC,covAD,varB,covBC,covBD,varC,covCD,varD"),c("varA,covAB,covAC,covAD,covAE,varB,covBC,covBD,covBE,varC,covCD,covCE,varD,covDE,varE"),c("varA,covAB,covAC,covAD,covAE,covAF,varB,covBC,covBD,covBE,covBF,varC,covCD,covCE,covCF,varD,covDE,covDF,varE,covEF,varF"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,varB,covBC,covBD,covBE,covBF,covBG,varC,covCD,covCE,covCF,covCG,varD,covDE,covDF,covDG,varE,covEF,covEG,varF,covFG,varG"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,varB,covBC,covBD,covBE,covBF,covBG,covBH,varC,covCD,covCE,covCF,covCG,covCH,varD,covDE,covDF,covDG,covDH,varE,covEF,covEG,covEH,varF,covFG,covFH,varG,covGH,varH"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,varC,covCD,covCE,covCF,covCG,covCH,covCI,varD,covDE,covDF,covDG,covDH,covDI,varE,covEF,covEG,covEH,covEI,varF,covFG,covFH,covFI,varG,covGH,covGI,varH,covHI,varI"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,varD,covDE,covDF,covDG,covDH,covDI,covDJ,varE,covEF,covEG,covEH,covEI,covEJ,varF,covFG,covFH,covFI,covFJ,varG,covGH,covGI,covGJ,varH,covHI,covHJ,varI,covIJ,varJ"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,varE,covEF,covEG,covEH,covEI,covEJ,covEK,varF,covFG,covFH,covFI,covFJ,covFK,varG,covGH,covGI,covGJ,covGK,varH,covHI,covHJ,covHK,varI,covIJ,covIK,varJ,covJK,varK"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,covAL,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,covBL,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,covCL,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,covDL,varE,covEF,covEG,covEH,covEI,covEJ,covEK,covEL,varF,covFG,covFH,covFI,covFJ,covFK,covFL,varG,covGH,covGI,covGJ,covGK,covGL,varH,covHI,covHJ,covHK,covHL,varI,covIJ,covIK,covIL,varJ,covJK,covJL,varK,covKL,varL")))
        }
        ### dimension handling if var.codes is only one line
        if(is.null(dim(var.codes))){
            var.codes <- as.matrix(t(var.codes))
        }
    } else if(multiallele==TRUE){
        ### if multiallelic call codes appear but a multiallelic posteriors file was not provided, use the default values
        mean.codes <- cbind(2:12,rbind(c("A,B"),c("A,B,C"),c("A,B,C,D"),c("A,B,C,D,E"),c("A,B,C,D,E,F"),c("A,B,C,D,E,F,G"),c("A,B,C,D,E,F,G,H"),c("A,B,C,D,E,F,G,H,I"),c("A,B,C,D,E,F,G,H,I,J"),c("A,B,C,D,E,F,G,H,I,J,K"),c("A,B,C,D,E,F,G,H,I,J,K,L")))
        var.codes <- cbind(2:12,rbind(c("varA,covAB,varB"),c("varA,covAB,covAC,varB,covBC,varC"),c("varA,covAB,covAC,covAD,varB,covBC,covBD,varC,covCD,varD"),c("varA,covAB,covAC,covAD,covAE,varB,covBC,covBD,covBE,varC,covCD,covCE,varD,covDE,varE"),c("varA,covAB,covAC,covAD,covAE,covAF,varB,covBC,covBD,covBE,covBF,varC,covCD,covCE,covCF,varD,covDE,covDF,varE,covEF,varF"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,varB,covBC,covBD,covBE,covBF,covBG,varC,covCD,covCE,covCF,covCG,varD,covDE,covDF,covDG,varE,covEF,covEG,varF,covFG,varG"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,varB,covBC,covBD,covBE,covBF,covBG,covBH,varC,covCD,covCE,covCF,covCG,covCH,varD,covDE,covDF,covDG,covDH,varE,covEF,covEG,covEH,varF,covFG,covFH,varG,covGH,varH"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,varC,covCD,covCE,covCF,covCG,covCH,covCI,varD,covDE,covDF,covDG,covDH,covDI,varE,covEF,covEG,covEH,covEI,varF,covFG,covFH,covFI,varG,covGH,covGI,varH,covHI,varI"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,varD,covDE,covDF,covDG,covDH,covDI,covDJ,varE,covEF,covEG,covEH,covEI,covEJ,varF,covFG,covFH,covFI,covFJ,varG,covGH,covGI,covGJ,varH,covHI,covHJ,varI,covIJ,varJ"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,varE,covEF,covEG,covEH,covEI,covEJ,covEK,varF,covFG,covFH,covFI,covFJ,covFK,varG,covGH,covGI,covGJ,covGK,varH,covHI,covHJ,covHK,varI,covIJ,covIK,varJ,covJK,varK"),c("varA,covAB,covAC,covAD,covAE,covAF,covAG,covAH,covAI,covAJ,covAK,covAL,varB,covBC,covBD,covBE,covBF,covBG,covBH,covBI,covBJ,covBK,covBL,varC,covCD,covCE,covCF,covCG,covCH,covCI,covCJ,covCK,covCL,varD,covDE,covDF,covDG,covDH,covDI,covDJ,covDK,covDL,varE,covEF,covEG,covEH,covEI,covEJ,covEK,covEL,varF,covFG,covFH,covFI,covFJ,covFK,covFL,varG,covGH,covGI,covGJ,covGK,covGL,varH,covHI,covHJ,covHK,covHL,varI,covIJ,covIK,covIL,varJ,covJK,covJL,varK,covKL,varL")))
    }
    
    ### check issues with pidFile: correct formatting, correct headers
    if (!missing(pidFile)){
        hold <- read.delim(pidFile, header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
        pid.read.in <- pidFile
        if(dim(hold)[2]!=1){
            stop("\n\npidFile must only have only column")
        }
        if(!("probeset_id"%in%names(hold))){
            cat(paste("\nHeader for pidFile is ",names(hold),". This should be probeset_id.\n",sep=""))
            cat("Header will be changed to probeset_id and file written to working directory as pidFile.ps.\n")
            write.table(hold,file="pidFile.ps",quote=FALSE,row.names=FALSE,col.names="probeset_id")
            names(hold) <- "probeset_id"
            pid.read.in <- "pidFile.ps"
        }
    }
    
    ### read in the headers from the calls and summary and reference files 
    ### to match up the sample order
    ### this is used to highlight the correct sample when plotting
    ### track any summary or reference samples that aren't in the calls file for removal later on
    calls.header <- as.list(rep(NA,l))
    calls.order <- as.list(rep(NA,l))
    summary.header <- as.list(rep(NA,l))
    summary.order <- as.list(rep(NA,l))
    hold.data.object <- as.list(rep(NA,l))
    if(!is.null(refFile)){
        ref.header <- as.list(rep(NA,l))
        ref.order <- as.list(rep(NA,l))
    }
    if(!is.null(confidenceFile)){
        conf.header <- as.list(rep(NA,l))
        conf.order <- as.list(rep(NA,l))
    }
    for(i in 1:l){
        line_counter <- 0
        while(sum(is.na(summary.header[[i]]))>0){
            summary.header[[i]] <- as.character(read.table(summaryFile[i],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=line_counter,check.names=F)[1,])
            if(summary.header[[i]][1]!="probeset_id"){
                summary.header[[i]] <- NA
                line_counter <- line_counter+1
            }
        }
        summary.header[[i]] <- summary.header[[i]][-1]
        hold.data.object[[i]] <- summary.header[[i]]
        ### if the summary file doesn't have a header at all, stop with an error
        if(sum(is.na(summary.header[[i]]))>0){
            if(batch==FALSE){
                stop("\n\nNo sample names found in the summary file. Please check and retry.\n\n")
            } else{
                stop(paste("\n\nNo sample names found in the summary file for batch ",i,". Please check and retry.\n\n",sep=""))
            }
        }
        if(!is.null(callFile) && plot.calls==TRUE){
            line_counter <- 0
            while(sum(is.na(calls.header[[i]]))>0){
                calls.header[[i]] <- as.character(read.table(callFile[i],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=line_counter,check.names=F)[1,])
                if(calls.header[[i]][1]!="probeset_id"){
                    calls.header[[i]] <- NA
                    line_counter <- line_counter+1
                }
            }
            calls.header[[i]] <- calls.header[[i]][-1]
            ### if the summary and calls file don't have any samples in common, stop with an error
            if(sum(summary.header[[i]]%in%calls.header[[i]])==0){
                stop(paste("\n\nNo samples from the calls file found in the summary file for batch ",i,".    Please check and retry.",sep=""))
            }
            ### if the samples don't all match, print a warning
            if(sum(summary.header[[i]]%in%calls.header[[i]])<length(calls.header[[i]]) || sum(calls.header[[i]]%in%summary.header[[i]])<length(summary.header[[i]])){
                if(batch==FALSE){
                    cat("\nNot all samples match between calls and summary files. Plots will only include samples that appear in both the calls and summary files.\n\n")
                } else{
                    cat(paste("\nNot all samples match between calls and summary files for batch ",i,". Plots will only include samples that appear in both the calls and summary files.\n\n",sep=""))
                }
            }
            hold.data.object[[i]] <- hold.data.object[[i]][hold.data.object[[i]]%in%calls.header[[i]]]
        }
        ### read in the reference file if supplied
        if(!is.null(refFile) && plot.ref==TRUE){
            if(!is.na(refFile[[i]])){
                line_counter <- 0
                while(sum(is.na(ref.header[[i]]))>0){
                    ref.header[[i]] <- as.character(read.table(refFile[i],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=line_counter,check.names=F)[1,])
                    if(ref.header[[i]][1]!="probeset_id"){
                        ref.header[[i]] <- NA
                        line_counter <- line_counter+1
                    }
                }
                ref.header[[i]] <- ref.header[[i]][-1]
                ### if the reference and calls file don't have any samples in common, put out a warning
                if(!is.null(callFile)){
                    if(sum(ref.header[[i]]%in%calls.header[[i]])==0){
                        if(batch==FALSE){
                            cat("\nNo samples from the calls file found in the reference file. Reference plots will be blank.\n\n")
                        } else{
                            cat(paste("\nNo samples from the calls file found in the reference file for batch ",i,". Reference plots will be blank.\n\n",sep=""))
                        }
                    }
                }
            }
        }
        ### read in the confidences file if supplied
        if(!is.null(confidenceFile) && plot.confs==TRUE){
            line_counter <- 0
            while(sum(is.na(conf.header[[i]]))>0){
                conf.header[[i]] <- as.character(read.table(confidenceFile[i],stringsAsFactors=F,quote="",sep="\t",nrows=1,skip=line_counter,check.names=F)[1,])
                if(conf.header[[i]][1]!="probeset_id"){
                    conf.header[[i]] <- NA
                    line_counter <- line_counter+1
                }
            }
            conf.header[[i]] <- conf.header[[i]][-1]
            ### if the summary and confidences file don't have any samples in common, stop with an error
            if(sum(summary.header[[i]]%in%conf.header[[i]])==0){
                stop(paste("\n\nNo samples from the confidences file found in the summary file for batch ",i,".    Please check and retry.",sep=""))
            }
            if(sum(summary.header[[i]]%in%conf.header[[i]])<length(conf.header[[i]]) || sum(conf.header[[i]]%in%summary.header[[i]])<length(summary.header[[i]])){
                if(batch==FALSE){
                    cat("\nNot all samples match between confidences and summary or calls files. Confidence plots will only include samples that appear in both the confidences and summary or calls files.\n\n")
                } else{
                    cat(paste("\nNot all samples match between confidences and summary or calls files for batch ",i,". Confidence plots will only include samples that appear in both the confidences and summary or calls files.\n\n",sep=""))
                }
            }
            if(is.null(callFile)){
                hold.data.object[[i]] <- hold.data.object[[i]][hold.data.object[[i]]%in%conf.header[[i]]]
            }
        }
        ### match up all of the samples
        summary.order[[i]] <- rep(0,length(summary.header[[i]]))
        for(j in 1:length(summary.header[[i]])){
            if(summary.header[[i]][j]%in%hold.data.object[[i]]){
                summary.order[[i]][j] <- which(hold.data.object[[i]]==summary.header[[i]][j])
            }
        }
        if(!is.null(callFile) && plot.calls==TRUE){
            calls.order[[i]] <- rep(0,length(calls.header[[i]]))
            for(j in 1:length(calls.header[[i]])){
                if(calls.header[[i]][j]%in%hold.data.object[[i]]){
                    calls.order[[i]][j] <- which(hold.data.object[[i]]==calls.header[[i]][j])
                }
            }
        }
        if(!is.null(refFile) && plot.ref==TRUE){
            if(!is.na(refFile[i])){
                ref.order[[i]] <- rep(0,length(ref.header[[i]]))
                for(j in 1:length(ref.header[[i]])){
                    if(ref.header[[i]][j]%in%hold.data.object[[i]]){
                        ref.order[[i]][j] <- which(hold.data.object[[i]]==ref.header[[i]][j])
                    }
                }
            }
        }
        if(!is.null(confidenceFile) && plot.confs==TRUE){
            conf.order[[i]] <- rep(0,length(conf.header[[i]]))
            for(j in 1:length(conf.header[[i]])){
                if(conf.header[[i]][j]%in%hold.data.object[[i]]){
                    conf.order[[i]][j] <- which(hold.data.object[[i]]==conf.header[[i]][j])
                }
            }
        }
    }
    if(!is.null(refFile) && plot.ref==TRUE){
        ### update the refFile name holder to put NULL for any batches without reference files
        ### this keeps the perl script from misfiring
        if(length(refFile)<l){
            refFile <- c(refFile,rep("NULL",length((length(refFile)+1):l)))
        }
        refFile[is.na(refFile)] <- "NULL"
    }
    
    ### read in the special SNPs file
    if(!is.null(specialSNPsFile)){
        specialSNPs <- read.delim(specialSNPsFile, header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
        if(dim(specialSNPs)[2]!=4){
            stop("\n\nSpecial SNPs file should have 4 columns.    Please check and retry.\n")
        }
        if((names(specialSNPs)[1]!="probeset_id")||(names(specialSNPs)[2]!="chr")||(names(specialSNPs)[3]!="copy_male")||(names(specialSNPs)[4]!="copy_female")){
            names(specialSNPs) <- c("probeset_id","chr","copy_male","copy_female")
        }
    } else{
        specialSNPs <- NULL
    }
    
    ### read in the report file
    if(!is.null(reportFile)){
        report <- as.list(rep(NA,l))
        gender <- as.list(rep(NA,l))
        for(i in 1:l){
            if(!is.na(reportFile[i])){
                report[[i]] <- read.delim(reportFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                if(sum(names(report[[i]])%in%c("cel_files","Sample","computed_gender","ComputedGender"))<2){
                    stop("\n\nColumns 'cel_files' and 'computed_gender' or 'Sample' and 'ComputedGender' must be in the report file.    Please check and retry.")
                }
                gender[[i]] <- as.data.frame(cbind(report[[i]]$cel_files,report[[i]]$computed_gender))
                if(sum(dim(gender[[i]]))==0){
                    gender[[i]] <- as.data.frame(cbind(report[[i]]$Sample,report[[i]]$ComputedGender))
                }
                if(sum(dim(gender[[i]]))==0){
                    stop("\n\nReport file must contain sample names and gender.    Please check the column names and retry.\n")
                }
                names(gender[[i]])[1] <- "Sample"
                names(gender[[i]])[2] <- "Gender"
                ### check if the gender data frame is missing any samples and give them a gender of unknown if so
                if(sum(is.na(gender[[i]]$Gender))>0){
                    if(batch==FALSE){
                        cat(paste("\nWarning: ",sum(is.na(gender[[i]]$Gender))," samples of ",dim(gender[[i]])[1]," samples are missing gender data in the report file. Missing genders have been set to unknown.\n",sep=""))
                    } else{
                        cat(paste("\nWarning: ",sum(is.na(gender[[i]]$Gender))," samples    of ",dim(gender[[i]])[1]," samples are missing gender data in the report file for batch ",i,". Missing genders have been set to unknown.\n",sep=""))
                    }
                    gender[[i]]$Gender[which(is.na(gender[[i]]$Gender))] <- "unknown"
                }
                if(sum(gender[[i]]$Sample%in%summary.header[[i]])==0){
                    gender[[i]] <- NA
                    if(batch==FALSE){
                        cat("\nWarning: no samples from the report file found in the summary file. Gender-separated plots will not be made.\n")
                    } else{
                        cat(paste("\nWarning: no samples from the report file found in the summary file for batch ",i,". Gender-separated plots will not be made for X and Y probesets.\n",sep=""))
                    }
                } else if(sum(gender[[i]]$Sample%in%summary.header[[i]]) != length(summary.header[[i]])){
                    gender[[i]]$Gender[-which(gender[[i]]$Sample%in%summary.header[[i]])] <- "unknown"
                    gender[[i]] <- cbind(summary.header[[i]][-which(summary.header[[i]]%in%gender[[i]]$Sample)],rep("unknown",length(summary.header[[i]][-which(summary.header[[i]]%in%gender[[i]]$Sample)])))
                }
                if(sum(is.na(gender[[i]]))==0){
                    if((sum(gender[[i]]$Gender=="unknown")/length(gender[[i]]$Gender))>=0.90){
                        if(batch==TRUE){
                            cat(paste("\nMore than 90% of the samples in batch ",i," have unknown gender. Gender-separated plots will not be made for X and Y probesets in batch ",i,".\n\n"),sep="")
                        } else{
                            cat("\nMore than 90% of the samples have unknown gender. Gender-separated plots will not be made for X and Y probesets.\n\n")
                        }
                    }
                }
            }
        }
    } else{
        report <- NULL
        gender <- NULL
    }
    
    ### read in the labels file to use as main label instead of probeset names if requested
    if(!(is.null(labelsFile))){
        labels <- as.list(rep(NA,l))
        for(i in 1:l){
            if(!is.na(labelsFile[i])){
                labels[[i]] <- read.delim(labelsFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                if(dim(labels[[i]])[2]!=2){
                    stop("\n\nLabels file must have two columns only.    Please check labels file and retry.")
                }
                if(sum(names(labels[[i]])==c("probeset_id","label"))!=2){
                    cat(paste("\nLabels file header line is: ", names(labels[[i]])[1], names(labels[[i]])[2],"\n",sep=" "))
                    cat("Labels file header line will be changed to 'probeset_id label'.\n")
                    names(labels[[i]]) <- c("probeset_id","label")
                }
            }
        }
    } else{
        labels <- NULL
    }
    
    ### read in the sample highlight file if supplied
    if(!is.null(sampleFile)){
        samples <- as.list(rep(NA,l))
        for(i in 1:l){
            if(!is.na(sampleFile[i])){
                samples[[i]] <- read.delim(sampleFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                if(dim(samples[[i]])[2]!=2){
                    stop("\n\nSample highlight file must have two columns only.    Please check samples file and retry.")
                }
                if(sum(names(samples[[i]])==c("sample","color"))!=2){
                    cat(paste("\nSample highlight file header line is: ", names(samples[[i]])[1], names(samples[[i]])[2],"\n",sep=" "))
                    cat("Sample highlight file header line will be changed to 'sample color'.\n")
                    names(samples[[i]]) <- c("sample","color")
                }
            }
        }
    } else{
        samples <- NULL
    }
    
    ### if the user doesn't already have a temp dir with files made, create it and run the perl script now
    if(use.temp.dir==FALSE){
        ### create temp dir
        for(i in 1:l){
            dir.create(temp.dir[i], showWarnings=F, recursive=T)
        }
        
        ### run perl script using summary, calls, confidences, and ref files
        ### output files are named using the probeset name "AX-probeset.txt", one text file per probeset
        ### and consolidates the input summary, calls, confidences, and ref files and are written to the temporary directory
        ### the cat/gsub line is to handle file and folder names with spaces so that they are passed to perl correctly
        for(i in 1:l){
            cmd_pidFile <- paste('"',pid.read.in,'"',sep='')
            cmd_summaryFile <- paste('"',summaryFile[i],'"',sep='')
            cmd_callFile <- ifelse(is.null(callFile[i]), 'NULL', paste('"',callFile[i],'"',sep=''))
            cmd_confidenceFile <- ifelse(is.null(confidenceFile[i]), 'NULL', paste('"',confidenceFile[i],'"',sep=''))
            cmd_refFile <- ifelse(is.null(refFile[i]), 'NULL',paste('"',refFile[i],'"',sep=''))
            cmd_temp.dir <- paste('"',temp.dir[i],'"',sep='')
            cmd <- paste("perl", paste('"', path.package("SNPolisher"),"/Perl/visualization.pl",'"', sep=''), cmd_pidFile, cmd_summaryFile, cmd_callFile, cmd_confidenceFile, cmd_refFile, cmd_temp.dir, sep=" ")
            sapply(cmd, system)
        }
    }
    
    ### let the user know how many probesets were requested and how many were found
    pid <- read.delim(pidFile, header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")[,1]
    if(length(pid)>length(unique(pid))){
        cat(length(pid), " probesets in the pid file. Only ",length(unique(pid)), " unique probesets are in the pid file and will be plotted. ")
        pid <- unique(pid)
    }
    cat(length(pid), " probesets in the pid file.    First ", ifelse(length(pid)<max.num.SNP.draw,length(pid),max.num.SNP.draw), " requested for plotting. ",sep="")
    cat(paste("Generating the cluster plots for ",min(length(pid),max.num.SNP.draw)," SNPs/probesets.\n\n", sep=""))
    if(length(pid)>max.num.SNP.draw){
        pid <- pid[1:max.num.SNP.draw]
    }
    if(length(pid)>=500){
        cat("Beginning to plot clusters. For larger pid lists (500 or more SNPs), this may produce a very large file.\n")
    }
    
    ### track which probesets are in the special SNPs list (if supplied)
    ### there will never be a chromosome named "A", so all non-X, non-Y, or non-MT chr get that value
    gender_keep <- cbind(pid,rep("A",length(pid)))
    if(!is.null(specialSNPs)){
        if((sum(pid%in%specialSNPs$probeset_id[specialSNPs$chr%in%c("X","Y")])>0) && (sum(unlist(lapply(gender, is.na)))<length(gender))){
            gender_keep <- cbind(pid,rep("A",length(pid)))
            hold <- specialSNPs[specialSNPs$probeset_id%in%pid,1:2]
            gender_keep[unlist(apply(hold,1,function(x) which(gender_keep[,1]==x[1]))),2] <- hold[,2]
        }
    }
    gender_keep <- as.data.frame(gender_keep)
    names(gender_keep) <- c("probeset_id","chr")
    
    ### read in the prior and posterior data if supplied
    postdata <- as.list(NA)
    post_names <- as.list(NA)
    priordata <- as.list(NA)
    prior_names <- as.list(NA)
    multipostdata <- as.list(NA)
    multi_names <- as.list(NA)
    multipriordata <- as.list(NA)
    multiprior_names <- as.list(NA)
    for(i in 1:l){
        ### read in posterior data, if user has supplied it
        ### get list of probeset names without ":1" for matching up data
        if(plot.posterior==TRUE){
            if(!is.null(posteriorFile[i]) && !is.na(posteriorFile[i])){
                if(tools::file_ext(posteriorFile[i])%in%c("gz","gzip")){
                    postdata[[i]] <- read.delim(gzfile(posteriorFile[i]), header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                } else{
                    postdata[[i]] <- read.delim(posteriorFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                }
                ### get rid of any extra blank lines in the file
                postdata[[i]] <- postdata[[i]][which(unlist(lapply((strsplit(postdata[[i]]$id,split="")),length))>0),]
                ### only keep the probesets in the pid file
                keep <- c(which(postdata[[i]]$id%in%pid),which(postdata[[i]]$id%in%paste(pid,":1",sep="")),which(postdata[[i]]$id%in%"GENERIC"))
                keep <- keep[order(keep)]
                postdata[[i]] <- postdata[[i]][keep,]
                ### hardcode in the MvA generic if it's not present in an MvA file
                ### all other transforms and Eureka data must supply their own generic labelled as GENERIC
                ### the posteriors/priors sub-function will look for a line named GENERIC if a probeset isn't present in the posteriors
                if(transform=="MvA" && sum(postdata[[i]]$id=="GENERIC")==0){
                    postdata[[i]] <- rbind(postdata[[i]],postdata[[i]][1,])
                    postdata[[i]][dim(postdata[[i]])[1],1] <- "GENERIC"
                    postdata[[i]][dim(postdata[[i]])[1],2] <- paste(-2,0.06,0,0,10.5,0.08,0,sep=",")
                    postdata[[i]][dim(postdata[[i]])[1],3] <- paste(0,0.06,0,0,11.0,0.08,0,sep=",")
                    postdata[[i]][dim(postdata[[i]])[1],4] <- paste(2,0.06,0,0,10.5,0.08,0,sep=",")
                    postdata[[i]][dim(postdata[[i]])[1],5] <- paste(-0.1,-0.1,-0.1,-0.05,-0.1,-0.05,0,0,0,0,0,0,sep=",")
                }
                post_names[[i]] <- unlist(strsplit(postdata[[i]]$id,split=":1"))
            } else{
                postdata[[i]] <- NA
                post_names[[i]] <- NA
                if(transform=="MvA"){
                    postdata[[i]] <- as.data.frame(t(c(id="GENERIC",BB=paste(-2,0.06,0,0,10.5,0.08,0,sep=","),AB=paste(0,0.06,0,0,11.0,0.08,0,sep=","),AA=paste(2,0.06,0,0,10.5,0.08,0,sep=","),CV=paste(-0.1,-0.1,-0.1,-0.05,-0.1,-0.05,0,0,0,0,0,0,sep=","))))
                    post_names[[i]] <- "GENERIC"
                }
            }
            ### read in multiallele posterior data, if user has supplied it
            if(!is.null(multiallele.posteriorFile[i]) && !is.na(multiallele.posteriorFile[i])){
                if(tools::file_ext(multiallele.posteriorFile[i])%in%c("gz","gzip")){
                    multipostdata[[i]] <- read.delim(gzfile(multiallele.posteriorFile[i]), header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                } else{
                    multipostdata[[i]] <- read.delim(multiallele.posteriorFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                }
                ### only keep the probesets in the pid file
                keep <- c(which(multipostdata[[i]]$probeset_id%in%pid),which(multipostdata[[i]]$probeset_id=="GENERIC"))
                keep <- keep[order(keep)]
                multipostdata[[i]] <- multipostdata[[i]][keep,]
                multi_names[[i]] <- unique(multipostdata[[i]]$probeset_id)
            } else{
                multipostdata[[i]] <- NA
                multi_names[[i]] <- NA
            }
        } else{
            postdata[[i]] <- NA
            multipostdata[[i]] <- NA
            multi_names[[i]] <- NA
        }
        ### read in prior data, if user has supplied it
        ### get list of probeset names without ":1" for matching up data
        if(plot.prior==TRUE){
            if(!is.null(priorFile[i]) && !is.na(priorFile[i])){
                if(tools::file_ext(priorFile[i])%in%c("gz","gzip")){
                    priordata[[i]] <- read.delim(gzfile(priorFile[i]), header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                } else{
                    priordata[[i]] <- read.delim(priorFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                }
                ### get rid of any extra blank lines in the file
                priordata[[i]] <- priordata[[i]][which(unlist(lapply((strsplit(priordata[[i]]$id,split="")),length))>0),]
                ### only keep the probesets in the pid file
                keep <- c(which(priordata[[i]]$id%in%pid),which(priordata[[i]]$id%in%paste(pid,":1",sep="")),which(priordata[[i]]$id%in%"GENERIC"))
                keep <- keep[order(keep)]
                ### hardcode in the MvA generic if it's not present in an MvA file
                ### all other transforms and Eureka data must supply their own generic labelled as GENERIC
                ### the posteriors/priors sub-function will look for a line named GENERIC if a probeset isn't present in the priors
                if(transform=="MvA" && sum(priordata[[i]]$id=="GENERIC")==0){
                    priordata[[i]] <- rbind(priordata[[i]],priordata[[i]][1,])
                    priordata[[i]][dim(priordata[[i]])[1],1] <- "GENERIC"
                    priordata[[i]][dim(priordata[[i]])[1],2] <- paste(-2,0.06,0,0,10.5,0.08,0,sep=",")
                    priordata[[i]][dim(priordata[[i]])[1],3] <- paste(0,0.06,0,0,11.0,0.08,0,sep=",")
                    priordata[[i]][dim(priordata[[i]])[1],4] <- paste(2,0.06,0,0,10.5,0.08,0,sep=",")
                    priordata[[i]][dim(priordata[[i]])[1],5] <- paste(-0.1,-0.1,-0.1,-0.05,-0.1,-0.05,0,0,0,0,0,0,sep=",")
                }
                prior_names[[i]] <- unlist(strsplit(priordata[[i]]$id,split=":1"))
            } else{
                priordata[[i]] <- NA
                prior_names[[i]] <- NA
                if(transform=="MvA"){
                    priordata[[i]] <- as.data.frame(t(c(id="GENERIC",BB=paste(-2,0.06,0,0,10.5,0.08,0,sep=","),AB=paste(0,0.06,0,0,11.0,0.08,0,sep=","),AA=paste(2,0.06,0,0,10.5,0.08,0,sep=","),CV=paste(-0.1,-0.1,-0.1,-0.05,-0.1,-0.05,0,0,0,0,0,0,sep=","))))
                    prior_names[[i]] <- "GENERIC"
                }
            }
            ### read in multiallele prior data, if user has supplied it
            if(!is.null(multiallele.priorFile[i]) && !is.na(multiallele.priorFile[i])){
                if(tools::file_ext(multiallele.priorFile[i])%in%c("gz","gzip")){
                    multipriordata[[i]] <- read.delim(gzfile(multiallele.priorFile[i]), header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                } else{
                    multipriordata[[i]] <- read.delim(multiallele.priorFile[i], header=T, check.names=F, stringsAsFactors=F, comment.char="#", sep="\t")
                }
                ### only keep the probesets in the pid file
                keep <- c(which(multipriordata[[i]]$probeset_id%in%pid),which(multipriordata[[i]]$probeset_id=="GENERIC"))
                keep <- keep[order(keep)]
                multipriordata[[i]] <- multipriordata[[i]][keep,]
                multiprior_names[[i]] <- unique(multipriordata[[i]]$probeset_id)
            } else{
                multipriordata[[i]] <- NA
                multiprior_names[[i]] <- NA
            }
        } else{
            priordata[[i]] <- NA
            multipriordata[[i]] <- NA
            multiprior_names[[i]] <- NA
        }
    }
    
    ### check if the posteriors file is for autotetraploid data or if autotetraploid is TRUE
    ### if so, check and update the call codes and colors and shape assignments
    if(sum(is.na(postdata))==0){
        post.names <- names(postdata[[1]])
        if(length(postdata)>1){
            for(j in 2:length(postdata)){
                post.names <- unique(c(post.names,names(postdata[[j]])))
            }
        }
        post.names <- post.names[-which(post.names%in%c("id","CV","OTV"))]
        if(length(post.names)>0 && sum(post.names%in%call.codes[,1])!=length(post.names)){
            ### see if the posteriors names are for autotetraploids
            if(sum(post.names%in%c("AA","AB","BB"))==0 && sum(post.names%in%c("BBBB","ABBB","AABB","AAAB","AAAA"))==5){
                call.codes <- rep(NA,3)
                call.codes <- c("OTV","-2","4")
                call.codes <- rbind(call.codes,c("NoCall","-1","4"))
                call.codes <- rbind(call.codes,c("AAAA","0","4"))
                call.codes <- rbind(call.codes,c("AAAB","1","4"))
                call.codes <- rbind(call.codes,c("AABB","2","4"))
                call.codes <- rbind(call.codes,c("ABBB","3","4"))
                call.codes <- rbind(call.codes,c("BBBB","4","4"))
                rownames(call.codes)[1] <- ""
                geno.col <- cbind(cbind((-2:4),c("cyan","gray","red","gold","blue","darkgreen","purple"),c(23,22,24,21,25,3,4)))
            }
        }
    }
    ### if there's no generic value, use a hard coded generic
    if(autotetraploid==TRUE && sum(is.na(postdata))>0){
        call.codes <- rep(NA,3)
        call.codes <- c("OTV","-2","4")
        call.codes <- rbind(call.codes,c("NoCall","-1","4"))
        call.codes <- rbind(call.codes,c("AAAA","0","4"))
        call.codes <- rbind(call.codes,c("AAAB","1","4"))
        call.codes <- rbind(call.codes,c("AABB","2","4"))
        call.codes <- rbind(call.codes,c("ABBB","3","4"))
        call.codes <- rbind(call.codes,c("BBBB","4","4"))
        rownames(call.codes)[1] <- ""
        geno.col <- cbind(cbind((-2:4),c("cyan","gray","red","gold","blue","darkgreen","purple"),c(23,22,24,21,25,3,4)))
        if(sum(is.na(postdata))>0){
            for(i in 1:l){
                postdata[[i]] <- NA
                post_names[[i]] <- NA
                postdata[[i]] <- as.data.frame(t(c(id="GENERIC",BBBB=paste(1,0.001,20,20,10.0,0.025,0,sep=","),ABBB=paste(0.5,0.001,20,20,10.25,0.025,0,sep=","),AABB=paste(0,0.001,20,20,10.5,0.025,0,sep=","),AAAB=paste(-0.5,0.001,20,20,10.25,0.025,0,sep=","),AAAA=paste(-1,0.001,20,20,10.0,0.025,0,sep=","),CV=paste(-0.1,-0.1,-0.1,-0.05,-0.1,-0.05,0,0,0,0,0,0,sep=","))))
                post_names[[i]] <- "GENERIC"
            }
        }
    }
    
    
    ### beginning to make plots
    ### if all plots are going in the same file, make that file now
    if(plot.all==TRUE){
        if(!is.null(output.dir)){
            output.File <- paste(output.dir,"/",output.File,sep="")
        }
        dir.create(dirname(output.File), showWarnings = F, recursive = T)
        if(plot.type=="png"){
            png(output.File, width=plot.width, height=plot.height, units="in", res=72)
        } else if(plot.type=="svg"){
            svg(output.File, width=plot.width, height=plot.height)
        } else{
            # pdf(output.File, width=plot.width, height=plot.height)
        }
    }
    mar <- c(5,4,3,1)
    ### need to change area of ylabel for RvT b/c label takes up more space
    if(transform%in%c("RvT","MvA-sqrt") || batch==T){
        mar <- c(5,5,3,1)
    }
    ### setting the number of rows and columns based on plotting types unless already specified by the user
    ### for 1, 2, or 4 options, there are 4 rows without blanks
    ### for 3 or 6 options, there are 4 rows with 3 plots each and a blank plot
    ### for 5 options, there are 5 rows without blanks
    if(is.null(num.rows)){
        num.rows <- 4
        if(length(plot.order)==5){
            num.rows <- 5
        }
    }
    if(is.null(num.cols)){
        num.cols <- 6
    }
    ### run the par command
    par(mfcol=c(num.rows,num.cols), mar=mar)
    
    ### holder variable to track if a new page needs to be made before plotting a probeset
    new_page <- FALSE
    
    ### track the number of plots made
    ### if the plots don't fill an entire column at the end of the batch/probeset, fill out the plots with extra empty plots
    num.plots <- 0
    
    ### set up the pid tracking variable to stop plotting when all probesets in the pid list have been seen
    for(q in 1:length(pid)){
        ### holder variable to track if all loops need to be broken up to the "q" loop to start a new probeset
        ### this should only be done if a probeset is missing calls or summary data i.e. it can't be plotted 
        next_probeset <- FALSE
        ### set up id variable for tracking what probeset is being plotted
        current_id <- pid[q]
        ### set up variable for holding data per probeset for plotting
        current_data <- as.list(rep(NA,l))
        ### holder variable to track if all batches have hit a new probeset name
        num_next <- as.list(rep(NA,l))
        line_counter <- as.list(rep(2,l))
        for(i in 1:l){
            ### current data object takes the smallest number of samples per batch
            current_data[[i]] <- as.data.frame(hold.data.object[[i]])
            names(current_data[[i]])[1] <- "samples"
        }
        ### loop through the lines in the data files reading each line into the obj object
        obj <- as.list(rep(NA,l))
        ### set the first probeset id read in
        if(is.na(current_id)){
            stop(paste("\n\npid number", q, "is missing.    Please check and retry.\n\n",sep=" "))
        }
        ### for a multiallele probeset, we need to store all of the calls for plotting the legend at the end
        legend.calls <- NA
        ### track if reading data in and plotting should continue
        while(sum(unlist(lapply(obj,length)))>0){
            ### loop through the batches
            for(j in 1:l){
                ### read in data for batches that haven't hit the next probeset yet
                if(is.na(num_next[[j]])){
                    if(!file.exists(paste(temp.dir[j], '/',current_id,'.txt', sep=''))){
                        next_probeset <- TRUE
                        cat("\nNo file found for",current_id, "in the temporary directory", temp.dir[j],"\n\n")
                        break
                    }
                    ### read the data file in line by line
                    obj[[j]] <- try(data<-scan(paste(temp.dir[j], '/',current_id,'.txt', sep=''),sep="\t",skip=line_counter[[j]]-1,what="list",nlines=1,quiet=T),silent=T)
                    ### check if the output in the temp directory is empty
                    ### i.e. no probesets in the pid list were found in the data files
                    if(line_counter[[j]]==2 && length(obj[[j]])==0){
                        if(batch==TRUE){
                            cat(paste("\nProbeset ",current_id," is missing calls and summary data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("\nProbeset ",current_id," is missing calls and summary data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        num_next[[j]] <- 0
                    }
                    ### check if you've hit the end of the data input file
                    if(length(obj[[j]])==0){
                        num_next[[j]] <- 0
                    } else if(obj[[j]][1]%in%pid){
                        ### check if the probeset is in the pid file
                        ### for the current probeset, store the data
                        if(obj[[j]][1]==current_id){
                            if(obj[[j]][2]=="calls" && plot.calls==TRUE){
                                hold <- rep(-1,dim(current_data[[j]])[1])
                                hold[calls.order[[j]]] <- as.numeric(unlist(strsplit(obj[[j]][3],","))[-1])[which(calls.order[[j]]>0)]
                                current_data[[j]] <- cbind(current_data[[j]],calls=hold)
                            }
                            if(obj[[j]][2]=="summary"){
                                hold <- unlist(strsplit(obj[[j]][3],","))
                                hold2 <- rep(-1,dim(current_data[[j]])[1])
                                hold2[summary.order[[j]]] <- as.numeric(unlist(strsplit(obj[[j]][3],","))[-1])[which(summary.order[[j]]>0)]
                                current_data[[j]] <- cbind(current_data[[j]],hold2)
                                ### pull out the allele name
                                names(current_data[[j]])[length(names(current_data[[j]]))] <- paste("summary",unlist(strsplit(hold[1],"-"))[length(unlist(strsplit(hold[1],"-")))],sep="_")
                            }
                            if(obj[[j]][2]=="confidences" && plot.confs==TRUE){
                                ### set up the confidences to trigger colors to be transparent for any sample that doesn't have a confidence
                                hold <- rep(-1,dim(current_data[[j]])[1])
                                hold[conf.order[[j]]] <- as.numeric(unlist(strsplit(obj[[j]][3],","))[-1])[which(conf.order[[j]]>0)]
                                current_data[[j]] <- cbind(current_data[[j]],confidences=hold)
                            }
                            if(obj[[j]][2]=="reference" && plot.ref==TRUE && !is.null(refFile)){
                                ### add in no calls for any samples that don't have a reference
                                hold <- rep(-1,dim(current_data[[j]])[1])
                                hold[ref.order[[j]]] <- as.numeric(unlist(strsplit(obj[[j]][3],","))[-1])[which(ref.order[[j]]>0)]
                                current_data[[j]] <- cbind(current_data[[j]],reference=hold)
                            }
                        } else{
                            ### mark that this batch has gotten to the next probeset
                            num_next[[j]] <- 0
                        }
                    } else{
                        ### mark that this probeset is not in the pid list
                        num_next[[j]] <- 0
                    }
                }
            }
            ### if you've hit the end of the data for the current id, plot the current id
            if(sum(unlist(lapply(num_next,is.na)))==0){
                ### check that there's calls/confs, and summary data for the probeset for all batches
                for(j in 1:l){
                    if(sum(names(current_data[[j]])=="calls")==0 && plot.calls==TRUE && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," is missing calls data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," is missing calls data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    } else{
                        legend.calls <- c(legend.calls,unique(current_data[[j]]$calls))
                    }
                    if(sum(names(current_data[[j]])=="calls")>1 && plot.calls==TRUE && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," has multiple lines of calls data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," has multiple lines of calls data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(sum(names(current_data[[j]])=="confidences")==0 && plot.confs==TRUE && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," is missing confidences data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," is missing confidences data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(sum(names(current_data[[j]])=="confidences")>1 && plot.confs==TRUE && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," has multiple lines of confidences data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," has multiple lines of confidences data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(sum(unlist(lapply(strsplit(names(current_data[[j]]),"summary"),length))>=2)<=1 && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," is missing summary data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," is missing summary data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(length(unique(names(current_data[[j]])))!=length(names(current_data[[j]])) && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Multiple data lines with the same name found for probeset ",current_id," in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Multiple data lines with the same name found for probeset ",current_id," and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(sum(names(current_data[[j]])=="reference")>1 && plot.ref==TRUE && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," has multiple lines of reference data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," has multiple lines of reference data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(dim(current_data[[j]])[1]==0 && next_probeset==FALSE){
                        if(batch==TRUE){
                            cat(paste("Probeset ",current_id," does not have any data in batch ",j," and will not be plotted.\n\n",sep=""))
                        } else{
                            cat(paste("Probeset ",current_id," does not have any data and will not be plotted.\n\n",sep=""))
                        }
                        next_probeset <- TRUE
                        break
                    }
                    if(sum(names(current_data[[j]])=="reference")==1){
                        legend.calls <- c(legend.calls,unique(current_data[[j]]$reference))
                    }
                }
                ### if there's a problem, break up to the next level
                if(next_probeset==TRUE){
                    break
                }
                ### check if it's a multiallele probeset before plotting
                multiallele <- FALSE
                multi.callcols <- NA
                ### pull the calls from all 4 sources: calls, references, posteriors, priors (posteriors and priors are lumped together)
                ### sources must be kept separate until the end so that the calls are in order of where they come from
                ### this is so that multiallele dynamic color and shape assignments match AxAS
                ### calls.hold is passed to the sub-function sub_multiallele_pairs for multiallele probesets
                calls <- NA
                calls.ref <- NA
                calls.post <- NA
                calls.hold <- NA
                batch_names <- NA
                summary_names <- NA
                for(j in 1:l){
                    batch_names <- unlist(c(batch_names,ifelse(plot.posterior, multi_names[[j]], multiprior_names[[j]])))
                    batch_names <- batch_names[!is.na(batch_names)]
                    summary_names <- c(summary_names,names(current_data[[j]]))
                    summary_names <- summary_names[!is.na(summary_names)]
                    if(!is.null(callFile)){
                        ### if a confidences cut-off value has been supplied, change any calls with a bad confidence value to NoCall
                        if(!is.null(confidences.cutoff)){
                            if("NoCall"%in%call.codes[,1]){
                                current_data[[j]]$calls[current_data[[j]]$confidences>confidences.cutoff] <- call.codes[call.codes[,1]=="NoCall",2]
                            }
                        }
                        calls <- c(calls,unique(current_data[[j]]$calls))
                        calls <- calls[!is.na(calls)]
                        calls.hold <- c(calls.hold,current_data[[j]]$calls)
                        calls.hold <- calls.hold[!is.na(calls.hold)]
                        calls <- unique(calls)
                        if(plot.ref==TRUE && "reference"%in%names(current_data[[j]])){
                            if(sum(is.na(calls.ref.match))==0){
                                for(n in 1:dim(calls.ref.match)[1]){
                                    current_data[[j]]$reference[current_data[[j]]$reference==calls.ref.match[n,3]] <- calls.ref.match[n,2]
                                    legend.calls[legend.calls==calls.ref.match[n,3]] <- calls.ref.match[n,2]
                                }
                            }
                            calls.ref <- c(calls.ref,unique(current_data[[j]]$reference))
                            calls.ref <- calls.ref[!is.na(calls.ref)]
                            calls.ref <- unique(calls.ref)
                            legend.calls <- unique(legend.calls)
                        }
                        ### pull out the call codes for any posterior ovals that must be plotted if a call appears
                        ### ie if BD appears, then the posterior ovals for BB, BD, and BB must also be plotted
                        ### the same set of ovals must appear if priors are plotted so there isn't a separate loop for those
                        if((plot.posterior==TRUE || plot.prior==TRUE) && sum(!is.na(unique(current_data[[j]]$calls)))>0){
                            if(sum(unique(current_data[[j]]$calls)>=0)>0){
                                alleles <- unique(unlist(strsplit(call.codes[call.codes[,2]%in%unique(current_data[[j]]$calls)[unique(current_data[[j]]$calls)>=0],1],"")))
                                if(length(multipostdata)>=j && sum(is.na(multipostdata[[j]]))==0){
                                    alleles <- unique(c(alleles,unlist(strsplit(multipostdata[[j]]$cluster[multipostdata[[j]]$probeset_id==current_id],split=""))))
                                }
                                alleles <- alleles[order(alleles)]
                                if(length(alleles)>=1){
                                    for(n in 1:length(alleles)){
                                        calls.post <- c(calls.post,alleles[n],paste(alleles[n],alleles[n],sep=""))
                                        for(m in (n+1):length(alleles)){
                                            if(n<m && m<=length(alleles)){
                                                calls.post <- c(calls.post,paste(alleles[n],alleles[m],sep=""))
                                            }
                                        }
                                    }
                                    if(length(alleles)==1){
                                        calls.post <- unique(c("AA",paste("A",alleles,sep=""),paste(alleles,alleles,sep=""),calls.post))
                                    }
                                }
                                calls.post <- unique(calls.post)
                                calls.post <- calls.post[!is.na(calls.post)]
                            }
                        }
                        if(plot.ref==TRUE){
                            calls.ref <- calls.ref[-which(calls.ref%in%calls)]
                            if(length(calls.ref)>0){
                                calls <- c(calls,calls.ref)
                            }
                        }
                        if(plot.posterior==TRUE || plot.prior==TRUE){
                            calls.post <- call.codes[call.codes[,1]%in%calls.post,2]
                            calls.post <- calls.post[-which(calls.post%in%calls)]
                            if(length(calls.post)>0){
                                calls <- c(calls,calls.post)
                            }
                        }
                        calls <- calls[!is.na(calls)]
                        calls <- unique(calls)
                    }
                }
                ### check if it has more than 2 summary lines per batch: it's multiallele
                if(sum(lapply(strsplit(summary_names,split="summary_"),length)==2)>(2*l)){
                    multiallele <- TRUE
                    num_multi_alleles <- length(unique(unlist(strsplit(summary_names,split="summary_")[lapply(strsplit(summary_names,split="summary_"),length)==2])[which(unlist(strsplit(summary_names,split="summary_")[lapply(strsplit(summary_names,split="summary_"),length)==2])!="")]))
                    multi_alleles <- unique(unlist(strsplit(summary_names,split="summary_")[lapply(strsplit(summary_names,split="summary_"),length)==2])[which(unlist(strsplit(summary_names,split="summary_")[lapply(strsplit(summary_names,split="summary_"),length)==2])!="")])
                    ### set the plot to a new page
                    new_page <- TRUE
                    par(mfcol=c(num.rows,num.cols), mar=mar)
                    num.plots <- 0
                    ### pull out the multiallelic call codes b/c they will need to have
                    ### colors and shapes dynamically assigned
                    ### MA call codes are all non-diploid, non-haploid, non-NC calls
                    multi.callcols <- call.codes[-which(call.codes[,1]%in%c("AA","AB","BB","A","B","ZeroCN","NoCall","OTV","NoCall_1","OTV_1")),]
                    ### if there's only one line, turn it into a data frame
                    if(is.null(dim(multi.callcols))&&length(multi.callcols)>0){
                        multi.callcols <- as.data.frame(t(multi.callcols))
                    }
                    ### if something has gone wrong (e.g. only no calls) set the multi.callcols object to be the geno.col object
                    if((is.na(multi.callcols)||(dim(multi.callcols)[1]==0))){
                        if(length(batch_names)>0 && pid[q]%in%batch_names){
                            multi.callcols <- geno.col
                        }
                        if(length(batch_names)==0){
                            multi.callcols <- geno.col
                        }
                    }
                }
                ### set the shapes and colors for all calls that appear for this probeset
                probeset.cols <- geno.col
                ### assign colors and shapes to any calls that aren't already in probeset.cols
                if(sum(calls%in%geno.col[,1])<length(calls)){
                    missing <- calls[-which(calls%in%geno.col[,1])]
                    if(sum(calls%in%geno.col[,1])==0){
                        missing <- calls
                    }
                    if(length(missing)>0){
                        ### remove any assigned shapes that have already been used before assigning multiallele shapes
                        hold.pch <- multi.pch[-which(multi.pch%in%unique(probeset.cols[probeset.cols[,1]>=0,3]))]
                        ### check that there are more multiallele shapes and colors than multiallele calls
                        if(length(hold.pch)>0){
                            while(length(unique(missing))>length(hold.pch)){
                                hold.pch <- c(hold.pch,hold.pch)
                            }
                        }
                        hold.col <- multi.col
                        while(length(unique(missing))>length(hold.col)){
                            hold.col <- c(hold.col,hold.col)
                        }
                        ### assign the shapes and colors to the calls without them
                        for(j in 1:length(missing)){
                            probeset.cols <- rbind(probeset.cols,c(missing[j],hold.col[j],hold.pch[j]))
                        }
                    }
                }
                ### if this probeset is biallelic and the previous probeset was multiallelic, set the plot to a new page
                ### new_page back to false
                if(multiallele==FALSE && new_page==TRUE){
                    par(mfcol=c(num.rows,num.cols), mar=mar)
                    new_page <- FALSE
                }
                ### if it's a multiallele probeset, set up the pairwise biallelic combinations for all batches
                biallele.order <- as.list(NA)
                if(multiallele==TRUE){
                    ### get the biallelic pair combinations to be plotted
                    data <- NA
                    for(j in 1:l){
                        data <- c(data,current_data[[j]]$calls)
                    }
                    if(sum(data%in%call.codes[call.codes[,1]%in%c("OTV","NoCall","OTV_1","NoCall_1","ZeroCN"),2])>0){
                        data <- data[-which(data%in%call.codes[call.codes[,1]%in%c("OTV","NoCall","OTV_1","NoCall_1","ZeroCN"),2])]
                    }
                    data <- data[!is.na(data)]
                    if(length(data)>0){
                        biallele.order <- sub_multiallele_pairs(data, call.codes, multiallele.order, multiallele.all, as.list(NA))
                    } else{
                        biallele.order[[1]][1] <- "A"
                        biallele.order[[1]][2] <- "B"
                    }
                }
                ### plot the batches
                for(j in 1:l){
                    ### store the x and y limits from calls for reference plots
                    xlim.calls <- NA
                    ylim.calls <- NA
                    if(multiallele==TRUE){
                        ### store the multiallele data so that the current_data[[j]] object can be reduced for biallelic pairwise plotting in the loop below
                        multi_data <- current_data[[j]]
                    }
                    ### if a probeset is multiallele, plot all samples and then loop through the different biallelic combinations
                    ### n starts at 1 to plot all samples on the A/B axes and then loops through all of the biallelic combination with those samples only
                    ### when n=1, the behavior is slightly different than for all other n values
                    for(n in 1:(length(biallele.order)+1)){
                        ### for non-multiallele probesets, only plot the full dataset (i.e. one set of biallelic combination only)
                        if(multiallele==FALSE && n>1){
                            next
                        }
                        ### if it's a multiallele probeset with only one biallele.order entry (i.e. all alleles are A and/or B only), skip any n>1
                        if(multiallele==TRUE && n>1 && length(biallele.order)==1){
                            next
                        }
                        ### for a multiallele probeset that is not the full dataset plot (i.e. n>1) set up the current_data[[j]] object to have only the calls that are in the biallelic combination
                        if(multiallele==TRUE && n>1){
                            ### reset back to the full data set
                            current_data[[j]] <- multi_data
                            ### find the heterozygous call code - if the het call made by pasting together the two alleles is in the wrong order, paste it in the other order
                            check <- call.codes[which(call.codes[,1]==paste(biallele.order[[n-1]][1],biallele.order[[n-1]][2],sep="")),2]
                            if(length(check)==0){
                                check <- call.codes[which(call.codes[,1]==paste(biallele.order[[n-1]][2],biallele.order[[n-1]][1],sep="")),2]
                            }
                            calls <- c(-9,-3,-2,-1,call.codes[which(call.codes[,1]==paste(biallele.order[[n-1]][1],biallele.order[[n-1]][1],sep="")),2],call.codes[which(call.codes[,1]==paste(biallele.order[[n-1]][2],biallele.order[[n-1]][2],sep="")),2])
                            if(length(check)==1){
                                calls <- c(calls,check)
                            }
                            ### set the current plotting data to only the calls for the biallelic combination plus no calls and OTV
                            if(dim(current_data[[j]][current_data[[j]]$calls%in%calls,])[1]>0){
                                current_data[[j]] <- current_data[[j]][current_data[[j]]$calls%in%calls,]
                            } else{
                                current_data[[j]]$calls <- as.numeric(call.codes[call.codes[,1]=="NoCall",2])
                            }
                        }
                        ### split the data for X and Y plots
                        gender_length <- 1
                        if(gender_keep$chr[gender_keep$probeset_id==current_id]%in%c("X","Y")){
                            if(!is.null(reportFile)&&!is.na(gender[[j]])){
                                if((sum(gender[[j]]$Gender=="unknown")/length(gender[[j]]$Gender))<0.90){
                                    if(sum(current_data[[j]]$samples%in%gender[[j]]$Sample)>0){
                                        male_data <- current_data[[j]][current_data[[j]]$samples%in%gender[[j]]$Sample[gender[[j]]$Gender=="male"],]
                                        female_data <- current_data[[j]][current_data[[j]]$samples%in%gender[[j]]$Sample[gender[[j]]$Gender=="female"],]
                                        unknown_data <- current_data[[j]][current_data[[j]]$samples%in%gender[[j]]$Sample[gender[[j]]$Gender=="unknown"],]
                                        gender_length <- 2
                                        ### fill out the end of a column if the gender plots can't all be plotted in the remaining column space
                                        extra.plots <- ifelse(num.plots%%num.rows==0, 0, num.rows - num.plots%%num.rows)
                                        if(extra.plots>0 && extra.plots<(m*gender_length)){
                                            for(m in 1:extra.plots){
                                                plot(0,0,type="n",axes=F,xlab="",ylab="")
                                                num.plots <- num.plots+1
                                            }
                                        }
                                    }
                                } else{
                                    ### set the internal chromosome to be A instead of X or Y so that gender-separated plots aren't made
                                    gender_keep$chr[gender_keep$probeset_id==current_id] <- "A"
                                }
                            }
                        }
                        ### run through the gender options for X and Y probesets
                        for(p in 1:gender_length){
                            gender.title <- NA
                            if(p==1){
                                if(gender_keep$chr[gender_keep$probeset_id==current_id]=="Y"){
                                    gender.title    <- "(all samples)"
                                } else if(gender_keep$chr[gender_keep$probeset_id==current_id]=="X"){
                                    current_data[[j]] <- rbind(male_data, unknown_data)
                                    gender.title <- "(male samples)"
                                }
                            } else{
                                if(gender_keep$chr[gender_keep$probeset_id==current_id]=="Y"){
                                    current_data[[j]] <- rbind(male_data, unknown_data)
                                    gender.title    <- "(male samples)"
                                } else if(gender_keep$chr[gender_keep$probeset_id==current_id]=="X"){
                                    current_data[[j]] <- rbind(female_data, unknown_data)
                                    gender.title <- "(female samples)"
                                }
                            }
                            ### plot the different options per probeset for each batch
                            for(m in 1:length(plot.order)){
                                ### a holder variable for if no calls are removed
                                no_plot <- NA
                                ### set up the type of plot being done
                                type <- "genotype"
                                if(length(unlist(strsplit(plot.order[m],"intensity.")))==2){
                                    type <- "intensity" 
                                }
                                col.by <- "reference"
                                if(length(unlist(strsplit(plot.order[m],".c")))==2){
                                    col.by <- "confidences"
                                    if(unlist(strsplit(plot.order[m],".c"))[2]=="alls"){
                                        col.by <- "calls" 
                                    }
                                }
                                ### set up the data for the plotting options
                                ### if it's intensity, then the data is the original summary values
                                if(type=="intensity"){
                                    current_data[[j]]$x_data <- current_data[[j]]$summary_A
                                    current_data[[j]]$y_data <- current_data[[j]]$summary_B
                                    xlab <- "A"
                                    ylab <- "B"
                                    multi.title <- NA
                                    if(multiallele==TRUE){
                                        ### multiallele probesets have different samples plotted depending on which biallelic combination n is on
                                        if(n>1){
                                            current_data[[j]]$x_data <- current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[n-1]][1],sep="_"))]
                                            current_data[[j]]$y_data <- current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[n-1]][2],sep="_"))]
                                            xlab <- biallele.order[[n-1]][1]
                                            ylab <- biallele.order[[n-1]][2]
                                            multi.title <- paste("(",biallele.order[[n-1]][1]," vs ",biallele.order[[n-1]][2],")",sep="")
                                        }
                                        if(n==1){
                                            ### for the first multialele plot, plot all of the samples in the same space
                                            if(length(biallele.order)>1){
                                                current_data[[j]]$x_data <- current_data[[j]][,which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[1]]
                                                current_data[[j]]$y_data <- current_data[[j]][,which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[2]]
                                                xlab <- unlist(strsplit(names(current_data[[j]])[which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[1]],"summary_"))[2]
                                                ylab <- unlist(strsplit(names(current_data[[j]])[which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[2]],"summary_"))[2]
                                                multi.title <- "(all samples)"
                                            } else{
                                                ### if there's a multiallele probeset but it only has one biallelic combination for all samples (i.e. all samples are BB, BD, or DD)
                                                ### then there's only one plot made and it's using the only biallelic combination
                                                current_data[[j]]$x_data <- current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[1]][1],sep="_"))]
                                                current_data[[j]]$y_data <- current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[1]][2],sep="_"))]
                                                xlab <- biallele.order[[1]][1]
                                                ylab <- biallele.order[[1]][2]
                                                multi.title <- paste("(",biallele.order[[1]][1]," vs ",biallele.order[[1]][2],")",sep="")
                                            }
                                        }
                                    }
                                }
                                ### if it's genotype (i.e. not intensity), then the data is the transformed summary values
                                if(type=="genotype"){
                                    current_data[[j]]$x_data <- log2(current_data[[j]]$summary_A)-log2(current_data[[j]]$summary_B)
                                    current_data[[j]]$y_data <- (log2(current_data[[j]]$summary_A)+log2(current_data[[j]]$summary_B))/2
                                    xlab <- "log2(A)-log2(B)"
                                    ylab <- "(log2(A)+log2(B))/2"
                                    multi.title <- NA
                                    if(transform=='CES'){
                                        current_data[[j]]$x_data <- sinh(K*(current_data[[j]]$summary_A-current_data[[j]]$summary_B)/(current_data[[j]]$summary_A+current_data[[j]]$summary_B))/sinh(K.val)
                                        current_data[[j]]$y_data <- log2(current_data[[j]]$summary_A+current_data[[j]]$summary_B)
                                        xlab <- 'sinh(K(A-B)/(A+B))/sinh(K)'
                                        ylab <- 'log2(A+B)'
                                    } else if(transform=='CSS'){
                                        current_data[[j]]$x_data <- asinh(K*(current_data[[j]]$summary_A-current_data[[j]]$summary_B)/(current_data[[j]]$summary_A+current_data[[j]]$summary_B))/asinh(K.val)
                                        current_data[[j]]$y_data <- log2(current_data[[j]]$summary_A+current_data[[j]]$summary_B)
                                        xlab <- 'asinh(K(A-B)/(A+B))/asinh(K)'
                                        ylab <- 'log2(A+B)'
                                    } else if(transform=='log2'){
                                        current_data[[j]]$x_data <- log2(current_data[[j]]$summary_A)
                                        current_data[[j]]$y_data <- log2(current_data[[j]]$summary_B)
                                        xlab <- 'log2(A)'
                                        ylab <- 'log2(B)'
                                    } else if(transform=='RvT'){
                                        current_data[[j]]$x_data <- atan(current_data[[j]]$summary_A/current_data[[j]]$summary_B)
                                        current_data[[j]]$y_data <- log(sqrt(current_data[[j]]$summary_A^2+current_data[[j]]$summary_B^2))
                                        xlab <- 'atan(A/B)'
                                        ylab <- expression(paste("ln(",sqrt("A"^2+"B"^2),")"))
                                    } else if(transform=="MvA-sqrt"){
                                        current_data[[j]]$x_data <- sqrt(current_data[[j]]$summary_A)-sqrt(current_data[[j]]$summary_B)
                                        current_data[[j]]$y_data <- (sqrt(current_data[[j]]$summary_A)+sqrt(current_data[[j]]$summary_B))/2
                                        xlab <- expression(paste(sqrt("A"),"-",sqrt("B")))
                                        ylab <- expression(paste("(",sqrt("A"),"+",sqrt("B"),")/2"))
                                    } else if(transform=="polar"){
                                        current_data[[j]]$x_data <- (2/pi)*atan(current_data[[j]]$summary_A/current_data[[j]]$summary_B)
                                        current_data[[j]]$y_data <- current_data[[j]]$summary_A+current_data[[j]]$summary_B
                                        xlab <- expression(paste("(2/",pi,")*atan(A/B)"))
                                        ylab <- "A+B"
                                    } else if(multiallele==TRUE){
                                        ### multiallele probesets have different samples plotted depending on which biallelic combination n is on
                                        if(n>1){
                                            current_data[[j]]$x_data <- log2(current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[n-1]][1],sep="_"))])
                                            current_data[[j]]$y_data <- log2(current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[n-1]][2],sep="_"))])
                                            xlab <- paste("log2(",biallele.order[[n-1]][1],")",sep="")
                                            ylab <- paste("log2(",biallele.order[[n-1]][2],")",sep="")
                                            multi.title <- paste("(",biallele.order[[n-1]][1]," vs ",biallele.order[[n-1]][2],")",sep="")
                                        }
                                        if(n==1){
                                            current_data[[j]]$x_data <- log2(current_data[[j]][,which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[1]])
                                            current_data[[j]]$y_data <- log2(current_data[[j]][,which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[2]])
                                            xlab <- paste("log2(",unlist(strsplit(names(current_data[[j]])[which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[1]],"summary_"))[2],")",sep="")
                                            ylab <- paste("log2(",unlist(strsplit(names(current_data[[j]])[which(unlist(lapply(strsplit(names(current_data[[j]]),"summar"),length))==2)[2]],"summary_"))[2],")",sep="")
                                            multi.title <- "(all samples)"
                                        }
                                        ### if there's a multiallele probeset but it only has one biallelic combination for all samples (i.e. all samples are BB, BD, or DD)
                                        ### then there's only one plot made and it's using the only biallelic combination
                                        if(n==1 && length(biallele.order)==1){
                                            current_data[[j]]$x_data <- log2(current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[1]][1],sep="_"))])
                                            current_data[[j]]$y_data <- log2(current_data[[j]][,which(names(current_data[[j]])==paste("summary",biallele.order[[1]][2],sep="_"))])
                                            xlab <- paste("log2(",biallele.order[[1]][1],")",sep="")
                                            ylab <- paste("log2(",biallele.order[[1]][2],")",sep="")
                                            multi.title <- paste("(",biallele.order[[1]][1]," vs ",biallele.order[[1]][2],")",sep="")
                                        }
                                    } else if(multiallele==FALSE && sum(lapply(strsplit(names(current_data[[j]]),split="summary_"),length)==2)>2){
                                        ### if there's a multiallele probeset that only has standard calls so multiallele is set to FALSE,
                                        ### then use the log2 transform and labels but that's it
                                        current_data[[j]]$x_data <- log2(current_data[[j]]$summary_A)
                                        current_data[[j]]$y_data <- log2(current_data[[j]]$summary_B)
                                        xlab <- "log2(A)"
                                        ylab <- "log2(B)"
                                    }
                                }
                                ### set up the bg and pch columns in the data
                                ### plots may be missing all data if it's a gender plot and that gender is not in the data
                                if(dim(current_data[[j]])[1]>0){
                                    current_data[[j]]$bg <- "black"
                                    current_data[[j]]$col <- "black"
                                    current_data[[j]]$pch <- 22
                                    ### fill in the plotting shapes and colors by col.by
                                    current_data[[j]] <- sub_shape_col(current_data[[j]], col.by, probeset.cols, call.codes, refNoCalls, no_plot, plotNoCalls, confs.col, confs.thresh)
                                    ### fill in any user-supplied colors for biallelic plots
                                    if(multiallele==FALSE){
                                        if(sum(c(!is.null(col.AA),!is.null(col.AB),!is.null(col.BB),!is.null(col.A),!is.null(col.B),!is.null(col.CN0),!is.null(col.OTV),!is.null(col.NC),!is.null(col.OTV1),!is.null(col.NC1)))>0){
                                            if(!is.null(col.AA) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="AA",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="AA",2]] <- col.AA
                                            }
                                            if(!is.null(col.AB) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="AB",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="AB",2]] <- col.AB
                                            }
                                            if(!is.null(col.BB) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="BB",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="BB",2]] <- col.BB
                                            }
                                            if(!is.null(col.A) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="A",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="A",2]] <- col.A
                                            }
                                            if(!is.null(col.B) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="B",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="B",2]] <- col.B
                                            }
                                            if(!is.null(col.CN0) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="ZeroCN",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="ZeroCN",2]] <- col.CN0
                                            }
                                            if(!is.null(col.OTV) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="OTV",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="OTV",2]] <- col.OTV
                                            }
                                            if(!is.null(col.NC) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="NoCall",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="NoCall",2]] <- col.NC
                                            }
                                            if(!is.null(col.OTV1) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="OTV_1",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="OTV_1",2]] <- col.OTV1
                                            }
                                            if(!is.null(col.NC1) && length(current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="NoCall_1",2]])>0){
                                                current_data[[j]]$bg[current_data[[j]]$calls==call.codes[call.codes[,1]=="NoCall_1",2]] <- col.NC1
                                            }
                                        }
                                    }
                                    ### fill in the color instead of the bg for any pch that isn't a fill option (pch<21)
                                    current_data[[j]]$col[current_data[[j]]$pch<21] <- current_data[[j]]$bg[current_data[[j]]$pch<21]
                                    ### except CN0 has pch 5 so set that back to black
                                    current_data[[j]]$col[current_data[[j]]$pch==5] <- "black"
                                    ### if it's a gender plot, set the unknown gender samples to be a black asterisk
                                    if(gender_length>1 && gender_keep$chr[gender_keep$probeset_id==current_id]%in%c("X","Y") && dim(unknown_data)[1]>0){
                                        if(sum(current_data[[j]]$samples%in%unknown_data$samples)>0){
                                            current_data[[j]]$bg[which(current_data[[j]]$samples%in%unknown_data$samples)] <- "black"
                                            current_data[[j]]$col[which(current_data[[j]]$samples%in%unknown_data$samples)] <- "black"
                                            current_data[[j]]$pch[which(current_data[[j]]$samples%in%unknown_data$samples)] <- 8
                                        }
                                    }
                                    ### rearrange order so highlighted samples are plotted last
                                    ### change plotting color to highlight color
                                    if(!is.null(samples)){
                                        if(sum(!is.na(samples[[j]]))>0){
                                            if(sum(samples[[j]]$sample%in%current_data[[j]]$samples)>0){
                                                for(k in 1:dim(samples[[j]])[1]){
                                                    current_data[[j]]$bg[which(current_data[[j]]$samples%in%samples[[j]]$sample[k])] <- samples[[j]]$color[k]
                                                }
                                                current_data[[j]] <- current_data[[j]][c((1:length(current_data[[j]]$samples))[!current_data[[j]]$samples%in%samples[[j]]$sample],(1:length(current_data[[j]]$samples))[current_data[[j]]$samples%in%samples[[j]]$sample]),]
                                            }
                                        }
                                    }
                                }
                                ### store the calls used in the legend.calls object if it's a multiallele probeset
                                if(multiallele==TRUE){
                                    legend.calls <- unique(c(legend.calls,unique(current_data[[j]]$calls)))
                                }
                                ### set up the initial x and y data limits to find the limits after the priors and posteriors are handled
                                xlim <- current_data[[j]]$x_data
                                ylim <- current_data[[j]]$y_data
                                ### prepare the posteriors for plotting if requested
                                ### posteriors only plotted for genotype/reference calls, never for intensities or confidences
                                d.post.diploid <- NA
                                d.post.haploid <- NA
                                if(type=="genotype" && dim(current_data[[j]])[1]>0 && ((col.by=="calls"&&plot.posterior==TRUE)||(col.by=="reference"&&plot.posterior.ref==TRUE))){
                                    ### function for producing posteriors/priors for the probeset
                                    ### output is a list whose parts must be assigned to the correct objects
                                    if(multiallele==FALSE && sum(lapply(strsplit(names(current_data[[j]]),split="summary_"),length)==2)<=2){
                                        hold <- sub_prepare_posteriors_priors(postdata[[j]], post_names[[j]], current_data[[j]]$calls, type, col.by, current_id, call.codes, probeset.cols, gender_keep, xlim, ylim, prepare.priors=FALSE, prepare.posteriors=TRUE, specialSNPs)
                                        ### update the color with any user-supplied ones
                                        if(sum(c(!is.null(col.AA),!is.null(col.AB),!is.null(col.BB),!is.null(col.A),!is.null(col.B),!is.null(col.CN0),!is.null(col.OTV),!is.null(col.NC),!is.null(col.OTV1),!is.null(col.NC1)))>0){
                                            if(!is.null(col.AA) && sum(is.na(hold$diploid))==0){
                                                hold$diploid[3,6] <- col.AA
                                            }
                                            if(!is.null(col.AB) && sum(is.na(hold$diploid))==0){
                                                hold$diploid[2,6] <- col.AB
                                            }
                                            if(!is.null(col.BB) && sum(is.na(hold$diploid))==0){
                                                hold$diploid[1,6] <- col.BB
                                            }
                                            if(!is.null(col.A) && sum(is.na(hold$haploid))==0){
                                                hold$haploid[2,6] <- col.A
                                            }
                                            if(!is.null(col.B) && sum(is.na(hold$haploid))==0){
                                                hold$haploid[1,6] <- col.B
                                            }
                                        }
                                    } else{
                                        if(n==1 && length(biallele.order)>1){
                                            biallele <- c("A","B")
                                        } else if(n==1 && length(biallele.order)==1){
                                            biallele <- biallele.order[[1]]
                                        } else{
                                            biallele <- biallele.order[[n-1]]
                                        }
                                        if(sum(is.na(multipostdata[[j]]))==0){
                                            data <- multipostdata[[j]][multipostdata[[j]]$probeset_id%in%c(current_id,"GENERIC"),]
                                        } else{
                                            data <- NA
                                        }
                                        hold <- sub_prepare_posteriors_priors_multi(data, biallele, current_data[[j]]$calls, type, col.by, probeset.cols, multi.col, gender_keep, xlim, ylim, call.codes, mean.codes, var.codes, num_multi_alleles, prepare.priors=FALSE, prepare.posteriors=TRUE)
                                    }
                                    d.post.diploid <- hold$d.diploid
                                    d.post.haploid <- hold$d.haploid
                                    post.diploid <- hold$diploid
                                    post.haploid <- hold$haploid
                                    xlim <- c(xlim,hold$xlim)
                                    ylim <- c(ylim,hold$ylim)
                                }
                                ### prepare the priors for plotting if requested
                                d.prior.diploid <- NA
                                d.prior.haploid <- NA
                                if(type=="genotype" && dim(current_data[[j]])[1]>0 && ((col.by=="calls"&&plot.prior==TRUE)||(col.by=="reference"&&plot.prior.ref==TRUE))){
                                    ### function for producing posteriors/priors for the probeset
                                    ### output is a list whose parts must be assigned to the correct objects
                                    if(multiallele==FALSE && sum(lapply(strsplit(names(current_data[[j]]),split="summary_"),length)==2)<=2){
                                        hold <- sub_prepare_posteriors_priors(priordata[[j]], prior_names[[j]], current_data[[j]]$calls, type, col.by, current_id, call.codes, probeset.cols, gender_keep, xlim, ylim, prepare.priors=TRUE, prepare.posteriors=FALSE, specialSNPs)
                                        ### update the color with any user-supplied ones
                                        if(sum(c(!is.null(col.AA),!is.null(col.AB),!is.null(col.BB),!is.null(col.A),!is.null(col.B),!is.null(col.CN0),!is.null(col.OTV),!is.null(col.NC),!is.null(col.OTV1),!is.null(col.NC1)))>0){
                                            if(!is.null(col.AA) && sum(is.na(hold$diploid))==0){
                                                hold$diploid[3,6] <- col.AA
                                            }
                                            if(!is.null(col.AB) && sum(is.na(hold$diploid))==0){
                                                hold$diploid[2,6] <- col.AB
                                            }
                                            if(!is.null(col.BB) && sum(is.na(hold$diploid))==0){
                                                hold$diploid[1,6] <- col.BB
                                            }
                                            if(!is.null(col.A) && sum(is.na(hold$haploid))==0){
                                                hold$haploid[2,6] <- col.A
                                            }
                                            if(!is.null(col.B) && sum(is.na(hold$haploid))==0){
                                                hold$haploid[1,6] <- col.B
                                            }
                                        }
                                    } else{
                                        if(n==1 && length(biallele.order)>1){
                                            biallele <- c("A","B")
                                        } else if(n==1 && length(biallele.order)==1){
                                            biallele <- biallele.order[[1]]
                                        } else{
                                            biallele <- biallele.order[[n-1]]
                                        }
                                        if(sum(is.na(multipriordata[[j]]))==0){
                                            data <- multipriordata[[j]][multipriordata[[j]]$probeset_id%in%c(current_id,"GENERIC"),]
                                        } else{
                                            data <- NA
                                        }
                                        hold <- sub_prepare_posteriors_priors_multi(data, biallele, current_data[[j]]$calls, type, col.by, probeset.cols, multi.col, gender_keep, xlim, ylim, call.codes, mean.codes, var.codes, num_multi_alleles, prepare.priors=TRUE, prepare.posteriors=FALSE)
                                    }
                                    d.prior.diploid <- hold$d.diploid
                                    d.prior.haploid <- hold$d.haploid
                                    prior.diploid <- hold$diploid
                                    prior.haploid <- hold$haploid
                                    xlim <- c(xlim,hold$xlim)
                                    ylim <- c(ylim,hold$ylim)
                                }
                                ### set up the plot ranges: default range values in case all data is missing, replaced with actual range values for data
                                if(dim(current_data[[j]])[1]>0){
                                    xlim <- range(xlim, finite=TRUE)
                                    ylim <- range(ylim, finite=TRUE)
                                } else{
                                    xlim <- c(ifelse(type=="intensity",0,-3),ifelse(type=="intensity",2500,3))
                                    ylim <- c(ifelse(type=="intensity",0,9),ifelse(type=="intensity",3500,12))
                                }
                                ### mirror the range of the X axis if requested
                                if(x.axis.mirror==TRUE && multiallele==FALSE){
                                    if(type=="genotype" ){
                                        xlim <- c(-abs(xlim)[order(abs(xlim),decreasing=TRUE)][1],abs(xlim)[order(abs(xlim),decreasing=TRUE)][1])
                                    }
                                    if(type=="intensity"){
                                        xlim <- c(min(xlim[1],ylim[1]),max(xlim[2],ylim[2]))
                                        ylim <- xlim
                                    }
                                }
                                ### store calls x and y limits to use with any matching reference plots
                                if(col.by=="calls"){
                                    xlim.calls <- xlim
                                    ylim.calls <- ylim
                                }
                                ### if there are matching x and y limits from a calls plot, set the limits for the reference or confidences plot
                                if(col.by%in%c("reference","confidences")){
                                    if(!is.na(xlim.calls[1])){
                                        xlim <- xlim.calls
                                    }
                                    if(!is.na(ylim.calls[1])){
                                        ylim <- ylim.calls
                                    }
                                }
                                ### update the x and y limits if the user has supplied them
                                if(!is.null(xlim.intensity) && type=="intensity"){
                                    xlim <- xlim.intensity
                                }
                                if(!is.null(ylim.intensity) && type=="intensity"){
                                    ylim <- ylim.intensity
                                }
                                if(!is.null(xlim.genotype) && type=="genotype"){
                                    xlim <- xlim.genotype
                                }
                                if(!is.null(ylim.genotype) && type=="genotype"){
                                    ylim <- ylim.genotype
                                }
                                ### make the plot
                                ### if plot.all is FALSE, set up the name for the new file to be plotted to
                                if(plot.all==FALSE){
                                    output.File <- current_id
                                    if(!is.null(output.dir)){
                                        output.File <- paste(output.dir[[j]],"/",output.File,sep="")
                                    }
                                    if(batch==TRUE){
                                        output.File <- paste(output.File,"_batch_",j,sep="")
                                    }
                                    if(multiallele==TRUE){
                                        output.File <- paste(output.File,"_ma_",n,sep="")
                                    }
                                    if(gender_length>1){
                                        output.File <- paste(output.File,"_gender_",p,sep="")
                                    }
                                    output.File <- paste(output.File,".",plot.order[m],".",plot.type,sep="")
                                    ### create the new file and reset the par to have 1 row and 1 column
                                    dir.create(dirname(output.File), showWarnings = F, recursive = T)
                                    if(plot.type=="png"){
                                        png(output.File, width=plot.width, height=plot.height, units="in", res=72)
                                    } else if(plot.type=="svg"){
                                        svg(output.File, width=plot.width, height=plot.height)
                                    } else{
                                        # pdf(output.File, width=plot.width, height=plot.height)
                                    }
                                    par(mfrow=c(1,1))
                                }
                                ### plotting sub-function takes all of the options that have been set up and puts them onto the plot
                                ### plotting is done regardless of plot.all status
                                
                                # browser()
                                
                                output1 = list(current_data, post.diploid)
                                names(output1) = c("current_data", "post.diploid")
                                return(output1)
                                
                                sub_plot_clusters(current_data[[j]], type, col.by, current_id, xlim, ylim, xlab, ylab, d.prior.diploid, prior.diploid, d.prior.haploid, prior.haploid, d.post.diploid, post.diploid, d.post.haploid, post.haploid, no_plot, labels[[j]], gender.title, multi.title, batch, j, multiallele, vertical.line, cex.main, cex.lab)
                                ### if plot.all is FALSE, close the file that was plotted to
                                if(plot.all==FALSE){
                                    dev.off()
                                }
                                num.plots <- num.plots+1
                            }
                            ### end of plotting options per probeset
                            ### fill out the end of a column if there are fewer plotting options than rows
                            extra.plots <- ifelse(num.plots%%num.rows==0, 0, num.rows - num.plots%%num.rows)
                            if(extra.plots>0 && extra.plots<length(plot.order) && plot.all==TRUE){
                                for(m in 1:extra.plots){
                                    plot(0,0,type="n",axes=F,xlab="",ylab="")
                                    num.plots <- num.plots+1
                                }
                            }
                        }
                        ### end of gender options per probeset
                        ### fill out the end of a column if it's a special SNP
                        extra.plots <- ifelse(num.plots%%num.rows==0, 0, num.rows - num.plots%%num.rows)
                        if(extra.plots>0 && (extra.plots<gender_length && gender_length>1) && plot.all==TRUE){
                            for(m in 1:extra.plots){
                                plot(0,0,type="n",axes=F,xlab="",ylab="")
                                num.plots <- num.plots+1
                            }
                        }
                    }
                    ### end of biallelic pairwise combination plotting
                    ### fill out the end of a column if there isn't enough space left to 
                    ### completely plot the next probeset without going to a new column
                    extra.plots <- ifelse(num.plots%%num.rows==0, 0, num.rows - num.plots%%num.rows)
                    if(extra.plots>0 && extra.plots<length(plot.order) && plot.all==TRUE){
                        for(m in 1:extra.plots){
                            plot(0,0,type="n",axes=F,xlab="",ylab="")
                            num.plots <- num.plots+1
                        }
                    }
                    ### when batch plotting multialleles, automatically reset to a new column for the next batch
                    if(multiallele==TRUE && batch==TRUE && length(biallele.order)>1 && plot.all==TRUE && num.plots%%(num.cols*num.rows)!=0){
                        extra.plots <- ifelse(num.plots%%num.rows==0, 0, num.rows - num.plots%%num.rows)
                        if(extra.plots>0){
                            for(m in 1:extra.plots){
                                plot(0,0,type="n",axes=F,xlab="",ylab="")
                                num.plots <- num.plots+1
                            }
                        }
                    }
                }
                ### finished plotting the batches
                ### if it's a multiallele probeset, plot the legend
                if(multiallele==TRUE && (plot.calls==TRUE || plot.ref==TRUE)){
                    ### fill out the end of a column before plotting the legend if plot.all if TRUE
                    extra.plots <- ifelse(num.plots%%num.rows==0, 0, num.rows - num.plots%%num.rows)
                    if(extra.plots>0 &&    plot.all==TRUE){
                        for(m in 1:extra.plots){
                            plot(0,0,type="n",axes=F,xlab="",ylab="")
                            num.plots <- num.plots+1
                        }
                    }
                    legend.calls <- unique(legend.calls[!is.na(legend.calls)])
                    ### get the matching shapes and colors per call
                    if(length(legend.calls)>0){
                        legend.cols <- probeset.cols[probeset.cols[,1]==legend.calls[1],]
                        if(length(legend.calls)>1){
                            for(m in 2:length(legend.calls)){
                                legend.cols <- rbind(legend.cols,probeset.cols[probeset.cols[,1]==legend.calls[m],])
                            }
                        }
                        if(is.null(dim(legend.cols))){
                            legend.cols <- as.data.frame(t(legend.cols))
                        }
                        ### rearrange the calls to be in alphabetical order
                        legend.order <- rep(0,dim(legend.cols)[1])
                        for(m in 1:dim(legend.cols)[1]){
                            legend.order[m] <- which(call.codes[,2]==legend.cols[m,1])
                        }
                        legend.cols <- legend.cols[order(legend.order),]
                        ### set the axes up to have several more units than the number of calls
                        axes.limit <- ifelse(is.null(dim(legend.cols)),7,(dim(legend.cols)[1]+1))
                        if(axes.limit<7){
                            axes.limit <- 7
                        }
                        ### if plot.all is FALSE, make the new file to hold the legend plot
                        if(plot.all==FALSE){
                            output.File <- paste(current_id,"_legend.",plot.type,sep="")
                            if(!is.null(output.dir)){
                                output.File <- paste(output.dir[[j]],"/",current_id,"_legend.",plot.type,sep="")
                            }
                            ### create the new file and reset the par to have 1 row and 1 column
                            dir.create(dirname(output.File), showWarnings = F, recursive = T)
                            if(plot.type=="png"){
                                png(output.File, width=plot.width, height=plot.height, units="in", res=72)
                            } else if(plot.type=="svg"){
                                svg(output.File, width=plot.width, height=plot.height)
                            } else{
                                pdf(output.File, width=plot.width, height=plot.height)
                            }
                            par(mfrow=c(1,1))
                        }
                        ### make an empty plot
                        plot(1:(axes.limit+1),0:axes.limit,type="n",main="Call Assignments",xlab="",ylab="",xaxt="n",yaxt="n")
                        ### if there's only one call, plot it in the legend
                        if(is.null(dim(legend.cols))){
                            points(2,axes.limit-2,pch=as.numeric(legend.cols[3]),bg=legend.cols[2],col="black")
                            text(4,axes.limit-2,legend.cols[1])
                            text(6,axes.limit-2,call.codes[call.codes[,2]==legend.cols[1],1])
                        } else{
                            ### run through all calls and plot those in the legend
                            ### the vertical offset ensures than the legend starts one unit down from the top and
                            ### end at least one unit up from the bottom
                            for(j in 1:dim(legend.cols)[1]){
                                points(2,axes.limit-j,pch=as.numeric(legend.cols[j,3]),bg=legend.cols[j,2],col="black")
                                text(4,axes.limit-j,legend.cols[j,1])
                                text(6,axes.limit-j,call.codes[call.codes[,2]==legend.cols[j,1],1])
                            }
                        }
                        ### end of plotting the legend if there's at least one call
                        ### close the file if plot.all is FALSE
                        if(plot.all==FALSE){
                            dev.off()
                        }
                    }
                }
                ### if this is a batch plot, set it to the next page in the plot
                if(batch==T) {
                    if(multiallele==FALSE){
                        par(mfcol=c(num.rows,num.cols), mar=mar)
                    }
                    ### reset the x and y limits at the end of the plots
                    xlim.calls <- NA
                    ylim.calls <- NA
                }
            }
            ### if there's a problem, break up one level
            if(next_probeset==TRUE){
                break
            }
            ### increment the line counter
            for(j in 1:l){
                if(length(obj[[j]])>0){
                    if(obj[[j]][1]==current_id){
                        line_counter[[j]] <- line_counter[[j]]+1
                    }
                }
            }
        }
        ### if there's a problem, go on to the next probeset
        if(next_probeset==TRUE){
            next
        }
        ### end of the probeset
    }
    ### close the connection to the output file
    if(plot.all==TRUE){
        dev.off()
    }
    
    ### remove the temporary directory if requested
    if(keep.temp.dir == FALSE){
        ### to handle the permissions error in Windows, delete any files from the
        ### temp directory and then delete it
        ### detect the OS and then run the removal loop
        sysinf <- Sys.info()
        if (!is.null(sysinf)){
            os <- sysinf['sysname']
            if (os == 'Darwin')
                os <- "osx"
        } else {
            os <- .Platform$OS.type
            if (grepl("^darwin", R.version$os))
                os <- "osx"
            if (grepl("linux-gnu", R.version$os))
                os <- "linux"
        }
        if(os=="Linux" || os=="linux"){
            for(i in 1:l){
                system(paste("rm -rf ", temp.dir[i])) 
            }
        } else{
            for(i in 1:l){
                system(paste("rm ", temp.dir[i], "/*",sep=""))
                system(paste("rm -r ", temp.dir[i])) 
            }
        }
    }
    ### end of Ps_Visualization function
    }


sub_plot_cov <- function(vx, vy, cov, x=0, y=0, col.ellipse, np=100, ...) {
    theta <- 0.5 * atan2(cov*2, vx-vy)
    alpha <- 2 * pi * (0:np)/np
    sint <- sin(theta)
    cost <- cos(theta)
    sina <- sin(alpha)
    cosa <- cos(alpha)
    a <- 2*sqrt(vx*cost*cost + vy*sint*sint + cov*2*sint*cost)
    b <- 2*sqrt(vx*sint*sint + vy*cost*cost - cov*2*sint*cost)
    x <- x + a * cosa * cost - b * sina * sint
    y <- y + a * cosa * sint + b * sina * cost
    lines(x, y, type="l", col=col.ellipse, ...)
}


### fill in the shapes and colors for plotting in the current_data[[j]] object
### data is the name of the current_data[[j]] object passed to the subfunction
sub_shape_col <- function(data, col.by, probeset.cols, call.codes, refNoCalls, no_plot, plotNoCalls, confs.col, confs.thresh, ...){
    if(col.by=="reference" && "reference"%in%names(data)){
        data$reference[is.na(data$reference)] <- -1
        data$reference[data$reference==-9] <- -1
        if(sum(call.codes[,1]=="NoCall")==1){
            if(as.numeric(call.codes[call.codes[,1]=="NoCall",2])!=-1){
                data$reference[data$reference==call.codes[call.codes[,1]=="NoCall",2]] <- -1
            }
        }
        if(sum(call.codes[,1]=="OTV")==1){
            if(as.numeric(call.codes[call.codes[,1]=="OTV",2])!=-2){
                data$reference[data$reference==call.codes[call.codes[,1]=="OTV",2]] <- -2
            }
        }
        for(k in 1:dim(probeset.cols)[1]){
            data$bg[data$reference==probeset.cols[k,1]] <- probeset.cols[k,2]
            data$pch[data$reference==probeset.cols[k,1]] <- as.numeric(probeset.cols[k,3])
        }
        ### rearrange order of current_data based on reference calls so haploids are plotted last
        data <- data[order(data$reference,decreasing=F),]
        ### set all colors to white for any no calls in reference data if refNoCalls is false
        if(refNoCalls==FALSE && sum(data$reference%in%c(-1,-9,-3))>0){
            data$bg[data$reference%in%c(-1,-9,-3)] <- "white"
            data$col[data$reference%in%c(-1,-9,-3)] <- "white"
            no_plot <- which(data$reference%in%c(-1,-9,-3))
        }
    } else if(col.by=="calls"){
        data$calls[is.na(data$calls)] <- -1
        data$calls[data$calls==-9] <- -1
        if(sum(call.codes[,1]=="OTV")==1){
            if(as.numeric(call.codes[call.codes[,1]=="NoCall",2])!=-1){
                data$calls[data$calls==call.codes[call.codes[,1]=="NoCall",2]] <- -1
            }
        }
        if(sum(call.codes[,1]=="OTV")==1){
            if(as.numeric(call.codes[call.codes[,1]=="OTV",2])!=-2){
                data$calls[data$calls==call.codes[call.codes[,1]=="OTV",2]] <- -2
            }
        }
        for(k in 1:dim(probeset.cols)[1]){
            data$bg[data$calls==probeset.cols[k,1]] <- probeset.cols[k,2]
            data$pch[data$calls==probeset.cols[k,1]] <- as.numeric(probeset.cols[k,3])
        }
        # rearrange order of current_data based on calls so haploids are plotted last
        data <- data[order(data$calls,decreasing=F),]
        ### set all colors to white for any no calls in calls data if plotNoCalls is false
        if(plotNoCalls==FALSE && sum(data$calls%in%c(-1,-9,-3))>0){
            data$bg[data$calls%in%c(-1,-9,-3)] <- "white"
            data$col[data$calls%in%c(-1,-9,-3)] <- "white"
            no_plot <- which(data$calls%in%c(-1,-9,-3))
        }
    } else if(col.by=="confidences" && "confidences"%in%names(data)){
        data$bg <- confs.col[1]
        data$bg[intersect(which(data$confidences>=confs.thresh[2]),which(data$confidences<confs.thresh[3]))] <- confs.col[2]
        data$bg[intersect(which(data$confidences>=confs.thresh[3]),which(data$confidences<confs.thresh[4]))] <- confs.col[3]
        data$bg[intersect(which(data$confidences>=confs.thresh[4]),which(data$confidences<confs.thresh[5]))] <- confs.col[4]
        data$bg[data$confidences>=confs.thresh[5]] <- "gray"
        ### any sample without a confidence value is set to transparent
        data$bg[data$confidences<0] <- "#FFFFFF00"
        data$col[data$confidences<0] <- "#FFFFFF00"
        data$pch <- 23
        data$pch[intersect(which(data$confidences>=confs.thresh[2]),which(data$confidences<confs.thresh[3]))] <- 24
        data$pch[intersect(which(data$confidences>=confs.thresh[3]),which(data$confidences<confs.thresh[4]))] <- 21
        data$pch[intersect(which(data$confidences>=confs.thresh[4]),which(data$confidences<confs.thresh[5]))] <- 25
        data$pch[data$confidences>confs.thresh[5]] <- 22
    } else{
        data$bg <- "white"
        data$col <- "white"
    }
    return(data)
}


### prepare the posterior and/or prior data for plotting
### split ellipses up into diploid and haploid
### remove any ellipses that shouldn't be plotted
sub_prepare_posteriors_priors <- function(data, names, calls, type, col.by, current_id, call.codes, probeset.cols, gender_keep, xlim, ylim, prepare.priors=FALSE, prepare.posteriors=TRUE, specialSNPs, ...){
    d.diploid <- NA
    d.haploid <- NA
    diploid <- NA
    haploid <- NA    
    if(type=="genotype" && (col.by=="calls"||col.by=="reference") && !is.na(data)){
        posterior <- data[which(names==current_id),]
        d.diploid <- which(posterior$id==current_id)
        d.haploid <- which(posterior$id==paste(current_id,":1",sep=""))
        ### if you haven't found anything, look for a GENERIC line
        if(length(d.diploid)==0 && length(d.haploid)==0){
            posterior <- data[which(names=="GENERIC"),]
            if(gender_keep$chr[which(gender_keep==current_id)]%in%c("Y","MT")){
                d.haploid <- which(posterior$id=="GENERIC")
            } else{
                d.diploid <- which(posterior$id=="GENERIC")
            }
        }
        if(length(d.diploid)==0){
            d.diploid <- NA
        }
        if(length(d.haploid)==0){
            d.haploid <- NA
        }
        if(!is.na(d.diploid)){
            diploid <- do.call('rbind', strsplit(as.character(posterior[d.diploid,c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))]), ','))[,c(1,2,5,6,7)]
            diploid[diploid=="NA"] <- NA
            diploid <- data.frame(apply(diploid, 2, as.numeric))
            ### put colors into the diploid object
            diploid <- cbind(diploid,rep(0,dim(diploid)[1]))
            for(i in 1:dim(diploid)[1]){
                diploid[i,6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]==names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))][i],2],2]
            }
            ### get rid of het ovals for MT and Y
            if(gender_keep$chr[which(gender_keep$probeset_id==current_id)]%in%c("Y","MT")){
                diploid <- diploid[-which(names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))]=="AB"),]
            }
            names(diploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.posterior')
            if(prepare.priors==TRUE){
                names(diploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.prior')
            }
            if(sum(apply(diploid,1,function(x) sum(is.na(x))))>0){
                diploid <- diploid[-which(apply(diploid,1,function(x) sum(is.na(x)))>0),]
            }
            xlim <- c(xlim,diploid$x[!is.na(diploid$x)]+2*sqrt(diploid$vx[!is.na(diploid$vx)]),diploid$x[!is.na(diploid$x)]-2*sqrt(diploid$vx[!is.na(diploid$vx)]))
            ylim <- c(ylim,diploid$y[!is.na(diploid$y)]+2*sqrt(diploid$vy[!is.na(diploid$vy)]),diploid$y[!is.na(diploid$y)]-2*sqrt(diploid$vy[!is.na(diploid$vy)]))
        }
        if(!is.na(d.haploid)){
            haploid <- do.call('rbind', strsplit(as.character(posterior[d.haploid,c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))]), ','))[,c(1,2,5,6,7)]
            haploid[haploid=="NA"] <- NA
            haploid <- data.frame(apply(haploid, 2, as.numeric))
            ### put colors into the diploid object
            haploid <- cbind(haploid,rep(0,dim(haploid)[1]))
            for(i in 1:dim(haploid)[1]){
                if(names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="B"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="A"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))][i]=="BB"){
                    if(length(probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="B",2],2])==1){
                        haploid[i,6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="B",2],2]
                    } else{
                        haploid[i,6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]==names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))][i],2],2]
                    }
                }
                if(names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="B"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="A"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))][i]=="AA"){
                    if(length(probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="A",2],2])==1){
                        haploid[i,6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="A",2],2]
                    } else{
                        haploid[i,6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]==names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))][i],2],2]
                    }
                }
            }
            ### if it's a CN genotyped probeset, there will be colors for calls A and B
            if(length(c(probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="B",2],2],probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="A",2],2]))>0){
                haploid[which(names(posterior)=="BB"),6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="B",2],2]
                haploid[which(names(posterior)=="AA"),6] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]=="A",2],2]
                ### get rid of ovals for any non-populated haploid clusters that aren't special SNPs
                if(!is.null(specialSNPs) && current_id%in%specialSNPs$probeset_id){
                    if(!specialSNPs$chr[specialSNPs$probeset_id==current_id]%in%c("X","Y","MT")){
                        if(sum(calls==as.numeric(call.codes[call.codes[,1]=="A",2]))==0){
                            haploid <- haploid[-2,]
                        }
                        if(sum(calls==as.numeric(call.codes[call.codes[,1]=="B",2]))==0){
                            haploid <- haploid[-1,]
                        }
                    }
                }
            } else{
                ### if it's a Y/X/MT probeset there will not be colors for A and B but AA and BB
                haploid[which(names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))]=="BB"),6] <- "#9EEFE3"
                haploid[which(names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))]=="AA"),6] <- "#FFCACA"
                ### get rid of ovals for any non-populated haploid clusters that aren't special SNPs
                if(!is.null(specialSNPs) && current_id%in%specialSNPs$probeset_id){
                    if(!specialSNPs$chr[specialSNPs$probeset_id==current_id]%in%c("X","Y","MT")){
                        if(sum(calls==as.numeric(call.codes[call.codes[,1]=="AA",2]))==0){
                            haploid <- haploid[-2,]
                        }
                        if(sum(calls==as.numeric(call.codes[call.codes[,1]=="BB",2]))==0){
                            haploid <- haploid[-1,]
                        }
                    }
                }
            }
            ### get rid of het values because haploid doesn't have het
            haploid <- haploid[-which(names(posterior)[c(which(names(posterior)=="BB"),which(names(posterior)=="AB"),which(names(posterior)=="AA"),which(names(posterior)=="BBBB"),which(names(posterior)=="ABBB"),which(names(posterior)=="AABB"),which(names(posterior)=="AAAB"),which(names(posterior)=="AAAA"),which(names(posterior)=="OTV"))]=="AB"),]
            ### finish if there's anything left in haploid, otherwise set d.haploid to NA
            if(dim(haploid)[1]>0){
                names(haploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.posterior')
                if(prepare.priors==TRUE){
                    names(haploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.prior')
                }
                if(sum(apply(haploid,1,function(x) sum(is.na(x))))>0){
                    haploid <- haploid[-which(apply(haploid,1,function(x) sum(is.na(x)))>0),]
                }
                xlim <- c(xlim,haploid$x[!is.na(haploid$x)]+2*sqrt(haploid$vx[!is.na(haploid$vx)]),haploid$x[!is.na(haploid$x)]-2*sqrt(haploid$vx[!is.na(haploid$vx)]))
                ylim <- c(ylim,haploid$y[!is.na(haploid$y)]+2*sqrt(haploid$vy[!is.na(haploid$vy)]),haploid$y[!is.na(haploid$y)]-2*sqrt(haploid$vy[!is.na(haploid$vy)]))
            } else{
                d.haploid <- NA
            }
        }
    }
    ### put the data into a list to return to the plotting function
    d <- list()
    d[["d.diploid"]] <- d.diploid
    d[["d.haploid"]] <- d.haploid
    d[["diploid"]] <- diploid
    d[["haploid"]] <- haploid
    d[["xlim"]] <- xlim
    d[["ylim"]] <- ylim
    return(d)
}


### prepare the posterior and/or prior data for plotting
### split ellipses up into diploid and haploid
### remove any ellipses that shouldn't be plotted
sub_prepare_posteriors_priors_multi <- function(data, biallele, calls, type, col.by, probeset.cols, multi.col, gender_keep, xlim, ylim, call.codes, mean.codes, var.codes, num_multi_alleles, prepare.priors=FALSE, prepare.posteriors=TRUE, ...){
    d.diploid <- NA
    d.haploid <- NA
    diploid <- NA
    haploid <- NA
    calls <- calls[calls>=0]
    if(sum(is.na(data))==0){
        if(sum(data$cluster%in%call.codes[,1])<dim(data)[1]){
            data <- data[which(data$cluster%in%call.codes[,1]),]
        }
    }
    if(type=="genotype" && (col.by=="calls"||col.by=="reference")){
        ### pull out probeset data if it's supplied
        if(sum(is.na(data))==0){
            if(sum(data$probeset_id!="GENERIC")>0){
                data <- data[which(data$probeset_id!="GENERIC"),]
            } else if(sum(data$probeset_id!="GENERIC")==0 && sum(data$probeset_id=="GENERIC")>0){
                data <- data[data$probeset_id=="GENERIC",]
                if(is.na(num_multi_alleles)||is.null(num_multi_alleles)){
                    num_multi_alleles <- 3
                }
                data <- data[data$nAlleles==num_multi_alleles,]
            }
            d.diploid <- data[data$copynumber==2,]
            if(length(d.diploid$probeset_id)>0){
                d.diploid <- d.diploid[d.diploid$cluster%in%c(paste(biallele[1],biallele[1],sep=""),paste(biallele[1],biallele[2],sep=""),paste(biallele[2],biallele[1],sep=""),paste(biallele[2],biallele[2],sep="")),]
            }
            if(length(d.diploid$probeset_id)==0){
                d.diploid <- NA 
            }
            d.haploid <- data[data$copynumber==1,]
            if(length(d.haploid$probeset_id)>0){
                d.haploid <- d.haploid[d.haploid$cluster%in%c(paste(biallele[1],biallele[1],sep=""),paste(biallele[2],biallele[2],sep="")),]
                if(length(d.haploid$probeset_id)==0){
                    d.haploid <- data[data$copynumber==1,]
                    d.haploid <- d.haploid[d.haploid$cluster%in%c(biallele[1],biallele[2]),]
                }
            }
            if(length(d.haploid$probeset_id)==0){
                d.haploid <- NA 
            }
        }
        ### if it's a multiallele probeset with only standard calls, then biallele will be NA
        if(is.na(biallele[1])){
            biallele <- c("A","B")
        }
        ### if you haven't found anything, replace with the GENERIC file
        if(sum(is.na(data))>0){
            generic <- read.table(file=paste(path.package("SNPolisher"),"/Perl/multiallele_generics.txt", sep=""),header=T,stringsAsFactors=FALSE)
            if(is.na(num_multi_alleles)||is.null(num_multi_alleles)){
                num_multi_alleles <- 3
            }
            generic <- generic[generic$nAlleles==num_multi_alleles,]
            d.diploid <- generic[generic$copynumber==2,]
            if(length(d.diploid$probeset_id)>0){
                d.diploid <- d.diploid[d.diploid$cluster%in%c(paste(biallele[1],biallele[1],sep=""),paste(biallele[1],biallele[2],sep=""),paste(biallele[2],biallele[1],sep=""),paste(biallele[2],biallele[2],sep="")),]
            }
            if(length(d.diploid$probeset_id)==0){
                d.diploid <- NA 
            }
            d.haploid <- generic[generic$copynumber==1,]
            if(length(d.haploid$probeset_id)>0){
                d.haploid <- d.haploid[d.haploid$cluster%in%c(paste(biallele[1],biallele[1],sep=""),paste(biallele[2],biallele[2],sep="")),]
            }
            if(length(d.haploid$probeset_id)==0){
                d.haploid <- NA 
            }
        }
        if(sum(is.na(d.diploid))==0){
            if(length(d.diploid$probeset_id)==0){
                d.diploid <- NA
            }
        }
        if(sum(is.na(d.haploid))==0){
            if(length(d.haploid$probeset_id)==0){
                d.haploid <- NA
            }
        }
        if(sum(is.na(d.diploid))==0){
            if(length(d.diploid$probeset_id)>0){
                diploid <- rbind(rep(NA,5),rep(NA,5),rep(NA,5))
                diploid.col <- rep(NA,dim(diploid)[1])
                for(r in 1:length(d.diploid$probeset_id)){
                    row <- d.diploid[r,]
                    diploid[r,1] <- as.numeric(unlist(strsplit(row$mean,","))[which(unlist(strsplit(mean.codes[mean.codes[,1]==row$nAlleles,2],","))==biallele[1])])
                    diploid[r,3] <- as.numeric(unlist(strsplit(row$mean,","))[which(unlist(strsplit(mean.codes[mean.codes[,1]==row$nAlleles,2],","))==biallele[2])])
                    diploid[r,2] <- as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("var",biallele[1],sep=""))])
                    diploid[r,4] <- as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("var",biallele[2],sep=""))])
                    if(length(as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("cov",biallele[1],biallele[2],sep=""))]))>0){
                        diploid[r,5] <-    as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("cov",biallele[1],biallele[2],sep=""))])
                    } else{
                        diploid[r,5] <-    as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("cov",biallele[2],biallele[1],sep=""))])
                    }
                    diploid.col[r] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]==row$cluster,2],2]
                }
                diploid <- data.frame(apply(diploid, 2, as.numeric))
                if(sum(is.na(diploid.col))>0){
                    diploid.col[is.na(diploid.col)] <- "white"
                }
                diploid <- cbind(diploid,diploid.col)
                ### get rid of het ovals for MT and Y
                if(sum(is.na(data))==0){
                    if(unique(data$probeset_id)!="GENERIC"){
                        if(gender_keep$chr[which(gender_keep$probeset_id==unique(data$probeset_id))]%in%c("Y","MT")){
                            diploid <- diploid[-2,]
                        }
                    }
                }
                names(diploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.posterior')
                if(prepare.priors==TRUE){
                    names(diploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.prior')
                }
                d.diploid <- 1
                xlim <- c(xlim,diploid$x[!is.na(diploid$x)]+2*sqrt(diploid$vx[!is.na(diploid$vx)]),diploid$x[!is.na(diploid$x)]-2*sqrt(diploid$vx[!is.na(diploid$vx)]))
                ylim <- c(ylim,diploid$y[!is.na(diploid$y)]+2*sqrt(diploid$vy[!is.na(diploid$vy)]),diploid$y[!is.na(diploid$y)]-2*sqrt(diploid$vy[!is.na(diploid$vy)]))
            }
        }
        if(!is.null(dim(d.haploid))){
            haploid <- rbind(rep(NA,5),rep(NA,5))
            haploid.col <- rep(NA,dim(haploid)[1])
            for(r in 1:dim(d.haploid)[1]){
                row <- d.haploid[r,]
                haploid[r,1] <- as.numeric(unlist(strsplit(row$mean,","))[which(unlist(strsplit(mean.codes[mean.codes[,1]==row$nAlleles,2],","))==biallele[1])])
                haploid[r,3] <- as.numeric(unlist(strsplit(row$mean,","))[which(unlist(strsplit(mean.codes[mean.codes[,1]==row$nAlleles,2],","))==biallele[2])])
                haploid[r,2] <- as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("var",biallele[1],sep=""))])
                haploid[r,4] <- as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("var",biallele[2],sep=""))])
                if(length(as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("cov",biallele[1],biallele[2],sep=""))]))>0){
                    haploid[r,5] <-    as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("cov",biallele[1],biallele[2],sep=""))])
                } else{
                    haploid[r,5] <-    as.numeric(unlist(strsplit(row$covariance,","))[which(unlist(strsplit(var.codes[var.codes[,1]==row$nAlleles,2],","))==paste("cov",biallele[2],biallele[1],sep=""))])
                }
                haploid.col[r] <- probeset.cols[probeset.cols[,1]==call.codes[call.codes[,1]==row$cluster,2],2]
            }
            ### finish if there's anything left in haploid, otherwise set d.haploid to NA
            if(dim(haploid)[1]>0){
                haploid <- data.frame(apply(haploid, 2, as.numeric))
                if(sum(is.na(haploid.col))>0){
                    haploid.col[is.na(haploid.col)] <- multi.col[(sum(haploid.col%in%multi.col)+1):(sum(haploid.col%in%multi.col)+sum(is.na(haploid.col)))]
                }
                haploid <- cbind(haploid,haploid.col)
                names(haploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.posterior')
                if(prepare.priors==TRUE){
                    names(haploid) <- c('x', 'vx', 'y', 'vy', 'cov','col.prior')
                }
                d.haploid <- 1
                xlim <- c(xlim,haploid$x[!is.na(haploid$x)]+2*sqrt(haploid$vx[!is.na(haploid$vx)]),haploid$x[!is.na(haploid$x)]-2*sqrt(haploid$vx[!is.na(haploid$vx)]))
                ylim <- c(ylim,haploid$y[!is.na(haploid$y)]+2*sqrt(haploid$vy[!is.na(haploid$vy)]),haploid$y[!is.na(haploid$y)]-2*sqrt(haploid$vy[!is.na(haploid$vy)]))
            } else{
                d.haploid <- NA
            }
        }
    }
    ### put the data into a list to return to the plotting function
    d <- list()
    d[["d.diploid"]] <- d.diploid
    d[["d.haploid"]] <- d.haploid
    d[["diploid"]] <- diploid
    d[["haploid"]] <- haploid
    d[["xlim"]] <- xlim
    d[["ylim"]] <- ylim
    return(d)
}


sub_plot_clusters <- function(data, type, col.by, current_id, xlim, ylim, xlab, ylab, d.prior.diploid, prior.diploid, d.prior.haploid, prior.haploid, d.post.diploid, post.diploid, d.post.haploid, post.haploid, no_plot, labels, gender.title, multi.title, batch, j, multiallele, vertical.line, cex.main, cex.lab, ...){
    ### make the main title
    main_title <- current_id
    ### check if there's a label supplied for the probeset
    if(!is.null(labels) && !is.na(labels)){
        if(current_id%in%labels$probeset_id){
            main_title <- labels$label[labels$probeset_id==current_id]
        }
    }
    if(!is.na(gender.title)){
        main_title <- paste(main_title,gender.title,sep=" ")
    } else if(batch==TRUE){
        if(is.null(labels) || is.na(labels)){
            main_title <- paste(main_title," (batch: ",j,")",sep="")
        }
        if(!is.null(labels) && !is.na(labels)){
            if(!current_id%in%labels$probeset_id){
                main_title <- paste(main_title," (batch: ",j,")",sep="")
            }
        }
    }
    if(type=="intensity"){
        sub_title <- paste('\nIntensities:', col.by, sep=" ")
    } else{
        sub_title <- paste('\nGenotypes:', col.by, sep=" ")
    }
    if(!is.na(multi.title)){
        sub_title <- paste(sub_title,multi.title,sep=" ")
    }
    if(batch==TRUE && !is.na(gender.title)){
        sub_title <- paste(sub_title," (batch: ",j,")",sep="")
    }
    ### plot any probesets that aren't non-PAR X or Y
    ### make the empty plot
    plot(0, 0, type='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=paste(main_title, sub_title, sep=' '), cex.main=cex.main, cex.lab=cex.lab)
    ### only plot the data if it's not missing, otherwise just have the empty plot
    if(dim(data)[1]>0){
        ### if the vertical line option is TRUE, plot the line before the points so that it's under the data
        if(type=="genotype" && vertical.line==TRUE){
            abline(v=0, lty=2, col="#E0E0E0")
        }
        ### plot the points for all calls or where there's data for references or confidences
        if(col.by=="calls" | (col.by=="reference" && "reference"%in%names(data)) | (col.by=="confidences" && "confidences"%in%names(data))){
            if(sum(is.na(no_plot))==0){
                points(data$x_data[-no_plot], data$y_data[-no_plot], col=data$col[-no_plot], pch=data$pch[-no_plot], bg=data$bg[-no_plot])
            } else{
                points(data$x_data, data$y_data, col=data$col, pch=data$pch, bg=data$bg)
            }
        }
        ### plot priors if requested
        if(!is.na(d.prior.haploid)){
            sapply(1:dim(prior.haploid)[1], function(j){sub_plot_cov(prior.haploid$vx[j], prior.haploid$vy[j], prior.haploid$cov[j], prior.haploid$x[j], prior.haploid$y[j], as.character(prior.haploid$col.prior[j]), lwd=2, lty=2)})
        }
        if(!is.na(d.prior.diploid)){
            sapply(1:dim(prior.diploid)[1], function(j){sub_plot_cov(prior.diploid$vx[j], prior.diploid$vy[j], prior.diploid$cov[j], prior.diploid$x[j], prior.diploid$y[j],    as.character(prior.diploid$col.prior[j]), lwd=2, lty=2)})
        }
        ### plot posteriors if requested
        if(!is.na(d.post.haploid)){
            sapply(1:dim(post.haploid)[1], function(j){sub_plot_cov(post.haploid$vx[j], post.haploid$vy[j], post.haploid$cov[j], post.haploid$x[j], post.haploid$y[j], as.character(post.haploid$col.posterior[j]), lwd=2)})
        }
        if(!is.na(d.post.diploid)){
            sapply(1:dim(post.diploid)[1], function(j){sub_plot_cov(post.diploid$vx[j], post.diploid$vy[j], post.diploid$cov[j], post.diploid$x[j], post.diploid$y[j],    as.character(post.diploid$col.posterior[j]), lwd=2)})
        }
    }
}



### make the biallelic pairs for a multiallelic probeset
sub_multiallele_pairs <- function(data, call.codes, multiallele.order, multiallele.all, biallele.order, ...){
    if(multiallele.all==TRUE){
        data <- data[data>=0]
        for(m in 1:dim(call.codes)[1]){
            if(sum(data==call.codes[m,2])>0){
                data[data==call.codes[m,2]] <- call.codes[m,1]
            }
        }
        alleles <- unique(unlist(strsplit(data,split="")))
        alleles <- alleles[order(alleles)]
        if(multiallele.order=="decreasing"){
            ### put the alleles into descreasing order before creating biallele.order
            calls <- unlist(strsplit(data,split=""))
            alleles.counts <- rep(0,length(alleles))
            for(m in 1:length(alleles)){
                alleles.counts[m] <- sum(calls==alleles[m])
            }
            alleles <- alleles[order(alleles.counts,decreasing=TRUE)]
        }
        if(length(alleles)==1){
            biallele.order[[1]][1] <- "A"
            biallele.order[[1]][2] <- ifelse(alleles=="A","B",alleles)
        } else{
            counts <- 1
            for(m in 1:(length(alleles)-1)){
                for(n in (m+1):length(alleles)){
                    ### check that any possible biallelic combination of the two alleles appears in the calls
                    if(sum(data%in%c(paste(alleles[m],alleles[m],sep=""),paste(alleles[m],alleles[n],sep=""),paste(alleles[n],alleles[m],sep=""),paste(alleles[n],alleles[n],sep="")))>0){
                        biallele.order[[counts]] <- c(alleles[m],alleles[n])[order(c(alleles[m],alleles[n]))]
                        counts <- counts+1
                    }
                }
            }
        }
    }
    ### check if there are any allele pairs with only one allele: homozogyous calls
    for(m in 1:length(biallele.order)){
        if(length(biallele.order[[m]])==1){
            ### check if there are all calls<0 i.e. all NoCalls, OTVs, etc
            if(is.na(biallele.order[[m]][1])){
                biallele.order[[m]][1] <- "A"
                biallele.order[[m]][2] <- "B"
            } else if(biallele.order[[m]][1]=="A"){
                biallele.order[[m]][2] <- "B"
            } else if(biallele.order[[m]][1]=="B"){
                biallele.order[[m]][1] <- "A"
                biallele.order[[m]][2] <- "B"
            } else{
                biallele.order[[m]][2] <- biallele.order[[m]][1]
                biallele.order[[m]][1] <- "A"
            }
        }
    }
    return(biallele.order)
}




evaluate.cluster.plots.old = function(recomb.df){
    #older, slower version of evaluate.cluster.plots
    #get recombination events which immediately change back
    recomb.df2 = recomb.df[which(recomb.df$geno.before == recomb.df$geno.after2), ]
    recomb.df2 = arrange(recomb.df2, probe.before.transition)
    
    change.mid.to.nocall = as.numeric()
    
    for(i in 1:nrow(recomb.df2)){
        #write sample file for Ps_Visualization()
        sample.name = recomb.df2$individual.trunc[[i]]
        samp.file = data.frame(sample.name, "black")
        colnames(samp.file) = c("sample", "color")
        write.table(samp.file, "samp.file.txt", sep = "\t", quote = F)
        
        
        probes = switch.affy.format(c(recomb.df2$probe.before.transition[[i]], 
                                                                    recomb.df2$probe.after.transition[[i]], recomb.df2$probe.after.transition2[[i]]))
        
        #write probe file for Ps_Visualization()
        
        genotype.positions = as.numeric()
        
        for(i2 in probes){
            pid.file = c("probeset_id", i2)
            writeLines(pid.file, "pid.file.txt")
            
            cluster.data = Ps_Visualization2(pidFile = "pid.file.txt",
                                                                             sampleFile = "samp.file.txt", output.File = "plot1.pdf",
                                                                             summaryFile = "AxiomGT1.summary.txt", callFile = "AxiomGT1.calls.txt",
                                                                             posteriorFile = "AxiomGT1.snp-posteriors.txt",
                                                                             confidenceFile = "AxiomGT1.confidences.txt", plot.prior = T,
                                                                             plot.ref = F, plot.intensity = F, plot.type = "pdf", plot.width = 15,
                                                                             plot.height = 15, num.cols = 1, num.rows = 3)
            
            sample.data1 = cluster.data$current_data[[1]][which(cluster.data$current_data[[1]]$samples == sample.name), ]
            
            call1 = sample.data1$calls
            sample.x.coord = sample.data1$x_data
            sample.y.coord = sample.data1$y_data
            
            #assign call to correct posterior oval
            if(call1 == 0) num1 = 3
            if(call1 == 1) num1 = 2
            if(call1 == 2) num1 = 1
            
            oval.data = sub_plot_cov2(cluster.data$post.diploid$vx[num1], cluster.data$post.diploid$vy[num1], cluster.data$post.diploid$cov[num1], cluster.data$post.diploid$x[num1], cluster.data$post.diploid$y[num1], "#000000", np = 100)
            
            genotype.positions = c(genotype.positions, point.in.polygon(sample.x.coord, sample.y.coord, oval.data$x, oval.data$y))
            
        }
        
        if(genotype.positions[2] == 0){
            change.mid.to.nocall = c(change.mid.to.nocall, 1)
        } else {
            change.mid.to.nocall = c(change.mid.to.nocall, 0)
        }
        
        print(paste0("done row ", i, " of ", nrow(recomb.df2)))
        
        
    }
    
    
    
    change.mid.to.nocall
    
}







get.allen.genetic.map = function(allen.geno.df, allen.rqtl.format){
    marker.column.coord = grep("marker", colnames(allen.geno.df))
    if(length(marker.column.coord) > 0) colnames(allen.geno.df)[marker.column.coord] = "marker"
    allen.geno.df[match(colnames(allen.rqtl.format), allen.geno.df$marker), 1:3][-c(1, 2), ]    
}



prepare.allen.map = function(allen.df, filter.chromosomes){
    if(missing(filter.chromosomes)) filter.chromosomes = T
    non.chromo.to.rm = which(!allen.df$chr %in% listofwheatchromosomes)
    if(length(non.chromo.to.rm) > 0) allen.df = allen.df[-non.chromo.to.rm, ]
    
    colnames(allen.df)[1] = "marker"
    allen.df.orig2 = split.df.into.list.of.dfs.by.column(allen.df, "chr")
    
    csxp.phys = lapply(allen.df.orig2, function(x){
        genetictophysical(unique(x$chr), x$marker, "-")
    })
    
    lapply(csxp.phys, function(x){
        na.omit(x)
    })
    
    
    find.max.gap = function(phys.positions){
        #examines how well distributed the SNPs are across the chromosome by physical position
        #the smaller the number the better (0 = perfect, 100 = worst)
        max(unlist(lapply(1:100, function(x){
            min(abs(phys.positions - x))
        })))
    }
    
    concordance = sapply(csxp.phys, function(y){
        find.max.gap(na.omit(y))
    })
    
    chr.to.keep = names(which(concordance <= 18))
    
    
    
    allen.df = as.data.frame(t(allen.df))
    allen.df = convert.to.character.data.frame(allen.df)
    colnames(allen.df) = allen.df[1, ]
    allen.df = allen.df[-1, ]
    allen.df = allen.df[-2, ]
    allen.df = cbind(allen.df[, 1:2], allen.df)
    allen.df[, 2] = rownames(allen.df)
    allen.df[, 1] = ""
    allen.df[2:nrow(allen.df), 1] = 2:nrow(allen.df)
    allen.df[1, 2] = ""
    allen.df[1, 1] = NA
    coln1 = colnames(allen.df)[1:2]
    colnames(allen.df)[1:2] = ""
    colnames(allen.df)[3:4] = coln1
    if(filter.chromosomes == T) allen.df = allen.df[, c(1, 2, which(allen.df[1, ] %in% chr.to.keep))]
    colnames(allen.df)[2] = "probeset_id"
    allen.df[allen.df == "AA"] = "A"
    allen.df[allen.df == "BB"] = "B"
    allen.df[allen.df == "AB"] = "H"
    
    
    allen.df
}



only.skeleton.markers = function(gen.map.df){
    # browser()
    g2 = split.df.into.list.of.dfs.by.column(gen.map.df, "chr", F)    
    g3 = lapply(g2, function(x){
        
        g4 = split.df.into.list.of.dfs.by.column(x, "cM", F)
        
        #within each cM bin, select only markers that have physical positions, and select the marker that has the highest 
        #physical position
        allg4 = lapply(g4, function(y){
            # browser()
            nonna1 = which(!is.na(y$phys))
            if(length(nonna1) > 0){
                y = y[nonna1, ]
                
                phys.sorted = sort(y$phys)
                middle.phys = phys.sorted[ceiling(length(phys.sorted) / 2)]
                
                # y = y[which(y$phys == middle.phys), ] #take the middle marker in each bin, or the smallest?
                y = y[which(y$phys == min(phys.sorted)), ]
                return(y)
            } else {
                return(y[1, ])
            }
        })
        
        bind_rows(allg4)
        
        # browser()
        # x = x[match(unique(x$cM), x$cM), ]
        # x[sort(x$cM, index.return = T)$ix, ]
        # browser()
    })
    
    bind_rows(g3)
    
}


get.recomb.freq = function(x){
    x = arrange(x, individual)
    aggregate(x$num.events, list(factor(x$individual, levels = unique(x$individual))), sum)
}



prepareqtl = function(geno.df, ind.avg.comb, recomb.freq.comb, fgen, gen.map, geno.prob.step, skip.liss, qtl.algorithm, multi.core){
    #args:
    # geno.df - genotyping dataframe in rqtl format
    # ind.avg.comb - MRD scores produced by ind.avg()
    # recomb.freq.comb - Recombination frequency produced by get.recomb.freq()
    # fgen - integer specifying filial generation 
    # gen.map - genetic map dataframe with columns "marker", "chr", "cM"
    # geno.prob.step - integer, step used for calc.genoprob()
    
    if(missing(skip.liss)) skip.liss = F
    if(missing(geno.prob.step)) geno.prob.step = 1
    if(missing(fgen)) fgen = 2
    if(missing(qtl.algorithm)) qtl.algorithm = "ehk"
    if(missing(multi.core)) multi.core = F
    
    chromos1 = unique(as.character(geno.df[1, ]))
    chromos1 = chromos1[-1]
    geno.df[1, 1] = ""
    ind.avg.comb$chromo = factor(ind.avg.comb$chromo, levels = chromos1)
    
    phenotypes1 = dcast(ind.avg.comb, individual ~ chromo, value.var = "phys.dist")
    
    phenotypes1 = phenotypes1[match(geno.df$probeset_id, phenotypes1$individual), ]
    
    if(ncol(phenotypes1) == 2){
        geno.df = add.column.at.position(geno.df, 0)
        geno.df[ , 1] = phenotypes1[, 2:ncol(phenotypes1)]
        geno.w.pheno1 = geno.df
        colnames(geno.w.pheno1)[1] = colnames(phenotypes1)[2]
    } else {
        geno.w.pheno1 = bind_cols(phenotypes1[, 2:ncol(phenotypes1)], geno.df)    
    }
    
    # geno.w.pheno1 = geno.w.pheno1[, -16]
    
    geno.w.pheno1[1, 1:(ncol(phenotypes1) - 1)] = ""
    
    r.freq.pheno = recomb.freq.comb[match(geno.w.pheno1$probeset_id, recomb.freq.comb$Group.1), ]
    
    geno.w.pheno2 = bind_cols(r.freq.pheno, geno.w.pheno1)
    geno.w.pheno2 = geno.w.pheno2[, -c(1)]
    geno.w.pheno2 = geno.w.pheno2[, !names(geno.w.pheno2) %in% c("treatment")]
    colnames(geno.w.pheno2)[1] = "Frequency"
    geno.w.pheno2[1, 1] = ""
    geno.w.pheno2[1, 2] = ""
    
    # browser()
    write.csv(geno.w.pheno2, "rotation1scripts_v4/original_data/genotype.data/rqtl.export/test.csv", row.names = F)
    
    
    if(fgen == "riself"){
        rqtl.test1 = rqtl.read("rotation1scripts_v4/original_data/genotype.data/rqtl.export/test.csv", crosstype = "riself")
    } else if(fgen == "dh"){
        rqtl.test1 = rqtl.read("rotation1scripts_v4/original_data/genotype.data/rqtl.export/test.csv", crosstype = "dh")    
    } else {
        rqtl.test1 = rqtl.read("rotation1scripts_v4/original_data/genotype.data/rqtl.export/test.csv", fgen = fgen)    
    }
    
    rqtl.test1
    
    names(unlist(lapply(rqtl.test1$geno, function(x){
        names(x[1])
    })))
    
    for(i in 1:length(rqtl.test1$geno)){
        newcm = gen.map[match(names(rqtl.test1$geno[[i]]$map), gen.map$marker), ]$cM
        names(newcm) = names(rqtl.test1$geno[[i]]$map)
        rqtl.test1$geno[[i]]$map = newcm
    }
    
    
    # lapply(rqtl.test1$geno, function(x){
    #     x$map = gen.map[match(names(x$map), gen.map$marker), ]$cM
    # })
    # 
    
    # browser()
    
    rqtl.test1$pheno[] = lapply(rqtl.test1$pheno, function(x) as.numeric(as.character(x)))
    
    rqtl.test1 = calc.genoprob(rqtl.test1, step = geno.prob.step)
    rqtl.test1 = sim.geno(rqtl.test1, step = 0)
    
    
    phenotypes.to.test = colnames(rqtl.test1$pheno)
    phenotypes.to.test[2:length(phenotypes.to.test)] = multi.str.split(phenotypes.to.test, "X", 2)[-1]
    
    if(skip.liss == T){
        phen5 = 1
    } else {
        phen5 = 1:(ncol(rqtl.test1$pheno) - 2)
    }
    # browser()
    if(multi.core == F){
        qtl.scans = lapply(phen5, function(x){
        # operm = scanone(rqtl.test1, method = "ehk", pheno.col = x, n.perm = 1000)
        operm = scanone(rqtl.test1, method = qtl.algorithm, pheno.col = x, n.perm = 1000)
        lod.threshold = unclass(summary(operm, alpha = c(0.05)))[1]
        lod.threshold2 = unclass(summary(operm, alpha = c((0.05 / (ncol(rqtl.test1$pheno) - 2)))))[1]
        # outmr = as.data.frame(scanone(rqtl.test1, method = "ehk", pheno.col = x))
        outmr = as.data.frame(scanone(rqtl.test1, method = qtl.algorithm, pheno.col = x))
        sig.markers = outmr[which(outmr$lod > lod.threshold), ]
        
        outmr.realmarkers = outmr[grep("AX", rownames(outmr)), ]
        
        sig.markers = sig.markers[grep("AX", rownames(sig.markers)), ]
        
        
        # browser()
        #examine luzie's qtl on 6A CS x P - do we get the same LOD score?
        # o1 = outmr[which(outmr$chr == "6A"), ]
        # which(o1$lod == max(o1$lod))
        # o1[70:80, ]
        # 
        # o1[which(rownames(o1) == "AX-94901219"), ]
        
        
        # if(nrow(sig.markers) > 1) browser()
        
        qtl.analysis1 = NULL
        
        if(nrow(sig.markers) > 0){
            
            qtl.analysis1 = lapply(1:nrow(sig.markers), function(z1){
                qtlanal.function1 = function(y){
                # browser()
                
                # -------- GET % PHENOTYPIC VARIANCE
                
                qtl2 = makeqtl(rqtl.test1, chr = sig.markers$chr[y], pos = sig.markers$pos[y], what = "prob")
                f1 = fitqtl(rqtl.test1, qtl = qtl2, formula = y~Q1, method = "hk", pheno.col = x)
                f2 = f1$result.full
                
                
                # -------- GET EFFECT MEANS AND SE
                
                eff = effectplot(rqtl.test1, mname1 = rownames(sig.markers)[y], pheno.col = x, draw = F)
                
                # ------ GET ADDITIVE AND DOMINANCE EFFECTS
                
                additive1 = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = y~Q1, method = "hk", pheno.col = x, get.ests = T, dropone = F))
                
                additive1 = unclass(additive1)$ests
                
                
                # ------ GET FLANKING MARKERS
                
                marker.name1 = rownames(sig.markers)[y]
                marker.coord1 = grep(marker.name1, rownames(outmr.realmarkers))
                flanking.markers = outmr.realmarkers[c((marker.coord1 - 1), (marker.coord1 + 1)), ]
                flanking.markers$flank = c("before", "after")
                flanking.markers = flanking.markers[which(as.character(flanking.markers$chr) == as.character(sig.markers[y, ]$chr)), ]
                #add physical information to flanking markers
                flanking.markers$phys = gen.map[match(rownames(flanking.markers), gen.map$marker), ]$phys
                flanking.markers$phys.bp = gen.map[match(rownames(flanking.markers), gen.map$marker), ]$phys.bp
                
                #add physical information to significant markers
                sig.markers$phys = gen.map[match(rownames(sig.markers), gen.map$marker), ]$phys
                sig.markers$phys.bp = gen.map[match(rownames(sig.markers), gen.map$marker), ]$phys.bp
                
                allqtllist = list(f2, eff, additive1, flanking.markers)
                names(allqtllist) = c("percent.var", "effect.means", "additive.effects", "flanking.markers")
                allqtllist }
                
                tryCatch(qtlanal.function1(z1), error = function(e) NULL)
                
                
            })
        }
        
        all.interaction.levels = NULL
        
        print("done 1")
        
        if(nrow(sig.markers) > 1 & nrow(sig.markers) < 5){
            qtlformula1 = function(num.markers, interaction1){
                #interaction1 - string, either "+" (no interaction) or "*" (all interactions recursively) or ":" (interactions without recursion)
                if(missing(interaction1)) interaction1 = "+"
                seq = 1:num.markers
                paste0("y~", paste("Q", seq, collapse = interaction1, sep = ""))
            }
            
            qtlformula1(nrow(sig.markers), ":")

            # qtl2 = tryCatch(makeqtl(rqtl.test1, chr = sig.markers$chr, pos = sig.markers$pos, what = "prob"), error = function(e) browser())
            qtl2 = makeqtl(rqtl.test1, chr = sig.markers$chr, pos = sig.markers$pos, what = "prob")
            #additive effects
            no.interaction1 = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = qtlformula1(nrow(sig.markers), "+"), method = "hk", pheno.col = x))
            interaction1 = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = qtlformula1(nrow(sig.markers), ":"), method = "hk", pheno.col = x))
            interaction1.recur = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = qtlformula1(nrow(sig.markers), "*"), method = "hk", pheno.col = x))
            all.interaction.levels = list(no.interaction1, interaction1, interaction1.recur)
        }
        
        outputlist1 = list(which(outmr$lod > lod.threshold),
                                             ggplot(outmr, aes(x = pos, y = lod)) + geom_line(size = 1) + 
                                                 facet_grid(cols = vars(chr), scales = "free_x", space = "free_x") + 
                                                 xlab("Chromosome") + ylab("LOD") + geom_hline(yintercept = lod.threshold, linetype = "dashed") + 
                                                 geom_hline(yintercept = lod.threshold2, linetype = "dotted") +
                                                 theme_bw() +
                                                 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                             panel.spacing = unit(0, "lines"), strip.background = element_rect(),
                                                             panel.grid = element_blank()) + ggtitle(phenotypes.to.test[x]),
                                             phenotypes.to.test[x],
                                             sig.markers, rqtl.test1, qtl.analysis1, all.interaction.levels, outmr, lod.threshold, lod.threshold2)
        
        
        names(outputlist1) = c("loci.over.threshold", "plot", "phenotype", "sig.markers", "rqtl.cross", "qtl.analysis", "multi.qtl.models", "outmr", "lod.threshold", "lod.threshold2")
        outputlist1
        # browser()
        
    })
    } else {
        qtl.scans = mclapply(phen5, function(x){
        # operm = scanone(rqtl.test1, method = "ehk", pheno.col = x, n.perm = 1000)
        operm = scanone(rqtl.test1, method = qtl.algorithm, pheno.col = x, n.perm = 1000)
        lod.threshold = unclass(summary(operm, alpha = c(0.05)))[1]
        lod.threshold2 = unclass(summary(operm, alpha = c((0.05 / (ncol(rqtl.test1$pheno) - 2)))))[1]
        # outmr = as.data.frame(scanone(rqtl.test1, method = "ehk", pheno.col = x))
        outmr = as.data.frame(scanone(rqtl.test1, method = qtl.algorithm, pheno.col = x))
        sig.markers = outmr[which(outmr$lod > lod.threshold), ]
        
        outmr.realmarkers = outmr[grep("AX", rownames(outmr)), ]
        
        sig.markers = sig.markers[grep("AX", rownames(sig.markers)), ]
        
        
        # browser()
        #examine luzie's qtl on 6A CS x P - do we get the same LOD score?
        # o1 = outmr[which(outmr$chr == "6A"), ]
        # which(o1$lod == max(o1$lod))
        # o1[70:80, ]
        # 
        # o1[which(rownames(o1) == "AX-94901219"), ]
        
        
        # if(nrow(sig.markers) > 1) browser()
        
        qtl.analysis1 = NULL
        
        if(nrow(sig.markers) > 0){
            
            qtl.analysis1 = lapply(1:nrow(sig.markers), function(z1){
                qtlanal.function1 = function(y){
                    # browser()
                    
                    # -------- GET % PHENOTYPIC VARIANCE
                    
                    qtl2 = makeqtl(rqtl.test1, chr = sig.markers$chr[y], pos = sig.markers$pos[y], what = "prob")
                    f1 = fitqtl(rqtl.test1, qtl = qtl2, formula = y~Q1, method = "hk", pheno.col = x)
                    f2 = f1$result.full
                    
                    
                    # -------- GET EFFECT MEANS AND SE
                    
                    eff = effectplot(rqtl.test1, mname1 = rownames(sig.markers)[y], pheno.col = x, draw = F)
                    
                    # ------ GET ADDITIVE AND DOMINANCE EFFECTS
                    
                    additive1 = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = y~Q1, method = "hk", pheno.col = x, get.ests = T, dropone = F))
                    
                    additive1 = unclass(additive1)$ests
                    
                    
                    # ------ GET FLANKING MARKERS
                    
                    marker.name1 = rownames(sig.markers)[y]
                    marker.coord1 = grep(marker.name1, rownames(outmr.realmarkers))
                    flanking.markers = outmr.realmarkers[c((marker.coord1 - 1), (marker.coord1 + 1)), ]
                    flanking.markers$flank = c("before", "after")
                    flanking.markers = flanking.markers[which(as.character(flanking.markers$chr) == as.character(sig.markers[y, ]$chr)), ]
                    #add physical information to flanking markers
                    flanking.markers$phys = gen.map[match(rownames(flanking.markers), gen.map$marker), ]$phys
                    flanking.markers$phys.bp = gen.map[match(rownames(flanking.markers), gen.map$marker), ]$phys.bp
                    
                    #add physical information to significant markers
                    sig.markers$phys = gen.map[match(rownames(sig.markers), gen.map$marker), ]$phys
                    sig.markers$phys.bp = gen.map[match(rownames(sig.markers), gen.map$marker), ]$phys.bp
                    
                    allqtllist = list(f2, eff, additive1, flanking.markers)
                    names(allqtllist) = c("percent.var", "effect.means", "additive.effects", "flanking.markers")
                    allqtllist }
                
                tryCatch(qtlanal.function1(z1), error = function(e) NULL)
                
                
            })
        }
        
        all.interaction.levels = NULL
        
        print("done 1")
        
        if(nrow(sig.markers) > 1 & nrow(sig.markers) < 5){
            qtlformula1 = function(num.markers, interaction1){
                #interaction1 - string, either "+" (no interaction) or "*" (all interactions recursively) or ":" (interactions without recursion)
                if(missing(interaction1)) interaction1 = "+"
                seq = 1:num.markers
                paste0("y~", paste("Q", seq, collapse = interaction1, sep = ""))
            }
            
            qtlformula1(nrow(sig.markers), ":")
            
            # qtl2 = tryCatch(makeqtl(rqtl.test1, chr = sig.markers$chr, pos = sig.markers$pos, what = "prob"), error = function(e) browser())
            qtl2 = makeqtl(rqtl.test1, chr = sig.markers$chr, pos = sig.markers$pos, what = "prob")
            #additive effects
            no.interaction1 = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = qtlformula1(nrow(sig.markers), "+"), method = "hk", pheno.col = x))
            interaction1 = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = qtlformula1(nrow(sig.markers), ":"), method = "hk", pheno.col = x))
            interaction1.recur = summary(fitqtl(rqtl.test1, qtl = qtl2, formula = qtlformula1(nrow(sig.markers), "*"), method = "hk", pheno.col = x))
            all.interaction.levels = list(no.interaction1, interaction1, interaction1.recur)
        }
        
        outputlist1 = list(which(outmr$lod > lod.threshold),
                                             ggplot(outmr, aes(x = pos, y = lod)) + geom_line(size = 1) + 
                                                 facet_grid(cols = vars(chr), scales = "free_x", space = "free_x") + 
                                                 xlab("Chromosome") + ylab("LOD") + geom_hline(yintercept = lod.threshold, linetype = "dashed") + 
                                                 geom_hline(yintercept = lod.threshold2, linetype = "dotted") +
                                                 theme_bw() +
                                                 theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                             panel.spacing = unit(0, "lines"), strip.background = element_rect(),
                                                             panel.grid = element_blank()) + ggtitle(phenotypes.to.test[x]),
                                             phenotypes.to.test[x],
                                             sig.markers, rqtl.test1, qtl.analysis1, all.interaction.levels, outmr, lod.threshold, lod.threshold2)
        
        
        names(outputlist1) = c("loci.over.threshold", "plot", "phenotype", "sig.markers", "rqtl.cross", "qtl.analysis", "multi.qtl.models", "outmr", "lod.threshold", "lod.threshold2")
        outputlist1
        # browser()
        
    }, mc.cores = 15)
    }
    
    qtl.scans
}

qtlanalysis.all = function(geno.df, genetic.map, pop.type, geno.return, marker.format, geno.prob.step, skip.liss, fgen, qtl.algorithm, 
                                                     skip.qtl.anal, multi.core){
    #get physical positions for markers without interpolating between NAs
    
    #geno.return - Boolean, if True, returns the genotyping dataframe after filtering for coverage and congruence with physical map (longest incr. subsequence)
                                #thus skipping the QTL analysis
    #skip.liss - Boolean, if TRUE, doesn't perform filtering of markers by chromosome coverage and longest increasing subsequence
    #fgen - unless fgen is 2, rQTL cross type will be set to riself, which does not recognize heterozygote genotypes
    if(missing(geno.return)) geno.return = F
    if(missing(marker.format)) marker.format = "-"
    if(missing(geno.prob.step)) geno.prob.step = 1
    if(missing(skip.liss)) skip.liss = F
    if(missing(fgen)) fgen = 6
    if(missing(qtl.algorithm)) qtl.algorithm = "ehk"
    if(missing(skip.qtl.anal)) skip.qtl.anal = F
    if(missing(multi.core)) multi.core = F
    # --------- MAIN ANALYSIS: 
    sxr.phys1 = get.phys.for.chromos(geno.df, F, T, marker.format)
    
    
    #add physical positions to genetic map dataframe
    genetic.map$phys = sxr.phys1$phys.pos[match(genetic.map$marker, sxr.phys1$marker)]
    genetic.map$phys.bp = sxr.phys1$phys.bp[match(genetic.map$marker, sxr.phys1$marker)]
    
    genetic.map = order.by.physical.dist.within.bins(genetic.map, "chr", "cM", "phys")
    
    
    # genetic.map[which(!genetic.map$marker %in% genetic.map2$marker), ]
    # genetic.map2[which(duplicated(genetic.map2$marker)), ]
    # 
    # genetic.map2[3175:3177, ]
    
    if(skip.liss == F){
        
        
        #get longest increasing subsequence of physical positions
        genetic.map3 = split.df.into.list.of.dfs.by.column(genetic.map, "chr", F)
        genetic.map4 = lapply(genetic.map3, function(x){
            # browser()
            na.coord = which(is.na(x$phys))
            if(length(na.coord) > 0) x = x[-na.coord, ]
            x = x[longest_subseq.R(na.omit(x$phys)), ]
            x
        })
        
        genetic.map5 = bind_rows(genetic.map4)
        
        #grab skeleton markers from longest increasing subsequence (take the central marker in each bin)
        #need to take only skeleton markers as QTL analysis downstream doesn't like markers in the same position
        genetic.map.ske = only.skeleton.markers(genetic.map5)
        
        #remove chromosomes with over 25% gap in chromosome physical distance
        genetic.map.ske2 = split.df.into.list.of.dfs.by.column(genetic.map.ske, "chr", do.sort = F)
        chromo.coverage1 = sapply(genetic.map.ske2, function(x){
            
            
            sequence1 = seq(0, 100, 25)
            
            distances1 = sapply(sequence1, function(z){
                min(abs(z - x$phys))
            })
            
            length(which(distances1 > 25))
        })
        
        chromo.coverage2 = which(chromo.coverage1 > 0)
        
        if(length(chromo.coverage2) > 0){
            print("removing the following chromosomes")
            print(chromo.coverage2)
            genetic.map.ske2 = genetic.map.ske2[-chromo.coverage2]
        }
        
        genetic.map.ske = bind_rows(genetic.map.ske2)
        
        

    } else {
        genetic.map = only.skeleton.markers(genetic.map)
        genetic.map.ske = genetic.map
    }
    
    #make new genotyping dataframe containing only skeleton markers in longest increasing subsequence of physical positions
    sxr.ls = geno.df[, c(1, 2, match(genetic.map.ske$marker, colnames(geno.df)))]
    
    #recalculate centimorgan distances between markers
    sxr.cm2 = extract.centimorgan.from.genotypes(sxr.ls, pop.type)
    
    genetic.map.ske$cM =    sxr.cm2[match(genetic.map.ske$marker, sxr.cm2$marker), ]$cm
    
    #make sure map and genotyping df only contain skeleton markers
    genetic.map.ske = only.skeleton.markers(genetic.map.ske)
    sxr.ls = sxr.ls[, c(1, 2, match(genetic.map.ske$marker, colnames(sxr.ls)))]
    
    
    
    if(geno.return == T){
        return(list(sxr.ls, genetic.map.ske, sxr.cm2))
    }
    
    #are any of the markers in the savannah x rialto map monomorphic?
    t.len2 = sapply(sxr.ls[2:nrow(sxr.ls), 3:ncol(sxr.ls)], function(x){
        if(length(which(x == "B")) == 0 | length(which(x == "A")) == 0){
            return(0)
        } else {
            return(1)
        }
    })
    
    which(t.len2 == 0) #no, they are not
    
    sxr.recomb = detect.recombination(sxr.ls)
    
    sxr.recomb$phys.pos = genetic.map.ske[match(sxr.recomb$probe.before.transition, genetic.map.ske$marker), ]$phys
    
    sxr.ind = ind.avg(sxr.recomb, 2)
    
    sxr.freq = get.recomb.freq(sxr.recomb)
    
    # browser()s
    if(skip.qtl.anal == T) return(list(sxr.ls, sxr.ind, sxr.freq, genetic.map.ske))
    list(prepareqtl(sxr.ls, sxr.ind, sxr.freq, fgen, genetic.map.ske, geno.prob.step, skip.liss, multi.core = multi.core), genetic.map.ske)
}



put.geno.in.order.of.chromo = function(geno.df){
    #rearranges a genotyping dataframe so that the chromosomes are in numerical and alphabetical order (e.g. 1A, 1B, 1D, 2A etc.)
    chromo.coords1 = sapply(listofwheatchromosomes, function(x){
        which(as.character(geno.df[1, ]) == x)
    })
    
    cc.tokeep = unlist(lapply(chromo.coords1, function(x){
        if(length(x) > 0){
            return(1)
        } else {
            return(0)
        }
    }))
    
    final.coords1 = c(1, 2, unlist(chromo.coords1[which(cc.tokeep == 1)]))
    
    geno.df[, final.coords1]
}




combinelist.of.genos = function(list.genos){
    #combine a list of genotyping dataframes into one
    list.genos[[1]][1, 2] = ""
    list.genos[[1]][1, 1] = ""
    
    allgcomb = do.call(rbind, list.genos)
    allgcomb = reset.rownames(allgcomb)
    
    allgcomb2 = allgcomb[-which(is.na(allgcomb[, 1])), ]
    
    allgcomb2 = convert.to.character.data.frame(allgcomb2)
    allgcomb2
}


grab.sig.phenotype.qtl = function(qtlanal.all.output){
    lapply(qtlanal.all.output[[1]], function(x){
        x[[4]]
    })    
}



parse.qtl.data = function(qtl.analysis.all.output){
    qtl.info1 = bind_rows(lapply(qtl.analysis.all.output$qtl.analysis, function(x){
        x$percent.var = as.data.frame(x$percent.var)
        x2 = data.frame(x$percent.var$LOD[1], x$percent.var$`%var`[1], x$percent.var$`Pvalue(Chi2)`[1], x$percent.var$`Pvalue(F)`[1])
        colnames(x2) = c("LOD", "% Var", "p-value Chi", "p-value F")
        x2
        
        eff1 = data.frame(x$effect.means$Means, x$effect.means$SEs)
        eff1$geno = rownames(eff1)
        colnames(eff1) = c("means", "se", "geno")
        
        if(nrow(eff1) == 3){
            eff2 = data.frame(eff1$means[1], eff1$means[2], eff1$means[3], eff1$se[1], eff1$se[2], eff1$se[3])
            colnames(eff2) = c("aa.mean", "ab.mean", "bb.mean", "aa.se", "ab.se", "bb.se")    
        } else if(nrow(eff1) == 2){
            eff2 = data.frame(eff1$means[1], eff1$means[2], eff1$se[1], eff1$se[2])
            colnames(eff2) = c("aa.mean", "bb.mean", "aa.se", "bb.se")    
        }
        
        
        
        
        add = as.data.frame(x$additive.effects)
        
        if(nrow(add) == 3){
            add2 = data.frame(add$est[2], add$est[3], add$SE[2], add$SE[3], add$t[2], add$t[3])
            colnames(add2) = c("add.est", "dom.est", "add.se", "dom.se", "add.t", "dom.t")            
        } else {
            add2 = data.frame(add$est[2], add$SE[2], add$t[2])
            colnames(add2) = c("add.est", "add.se", "add.t")            
        }
        
        
        f1 = x$flanking.markers[1, ]
        f1 = f1[, -which(colnames(f1) == "flank")]
        colnames(f1) = paste0(colnames(f1), ".before")    
        
        f2 = x$flanking.markers[2, ]
        f2 = f2[, -which(colnames(f2) == "flank")]
        colnames(f2) = paste0(colnames(f2), ".after")    
        
        f3 = as.data.frame(cbind(f1, f2))
        
        bind_cols(x2, eff2, add2, f3)
        
        
    }))
    
    
    
    returndat1 = as.data.frame(cbind(qtl.analysis.all.output$sig.markers, qtl.info1))
    #add phenotype to dataframe
    if(qtl.analysis.all.output$phenotype %in% listofwheatchromosomes) qtl.analysis.all.output$phenotype = paste0("MRD.", qtl.analysis.all.output$phenotype)
    returndat1$phenotype = qtl.analysis.all.output$phenotype
    returndat1 = convert.to.character.data.frame(returndat1)
    
    #grab only the numeric columns (from all being character)
    numeric.cols1 = which(!sapply(lapply(returndat1, as.numeric), function(x) all(is.na(x))))
    for(i in numeric.cols1){
        returndat1[, i] = round(as.numeric(returndat1[, i]), digits = 3)
    }
    
    
    returndat1
}


get.sig.qtl = function(x){
    g = lapply(x[[1]], function(y) y[[1]])
    names(g) = 1:length(g)
    # browser()
    g[which(lapply(g, length) > 0)]
}


relabel.treat.w.pop.names = function(data.df, treatment.col.name, list.of.pop.names, reverse.treatment.order){
    if(missing(reverse.treatment.order)) reverse.treatment.order = F
    if(length(unique(data.df[[treatment.col.name]])) != length(list.of.pop.names)){
        stop("number of treatments does not equal number of population names")
    }
    data.df[[treatment.col.name]] = as.character(data.df[[treatment.col.name]])
    for(i in 1:4){
        data.df[[treatment.col.name]][which(data.df[[treatment.col.name]] == i)] = list.of.pop.names[i]
    }
    data.df[[treatment.col.name]] = factor(data.df[[treatment.col.name]], levels = unique(data.df[[treatment.col.name]]))
    
    if(reverse.treatment.order == T){
        data.df = split(data.df, data.df[[treatment.col.name]])
        data.df = data.df[length(list.of.pop.names):1]
        data.df = bind_rows(data.df)
        data.df[[treatment.col.name]] = factor(data.df[[treatment.col.name]], levels = unique(data.df[[treatment.col.name]]))
    }
    
    data.df
}

