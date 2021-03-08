#### FUNCTIONS ####

# setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

load("rotation1scripts_v4/saved.objects/axc.chr1a.cm")
binned.markers1a = unique(axc.chr1a.cm)
load("rotation1scripts_v4/saved.objects/cxa.6b.cm")

full.path = "/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/"

make.founder.genotypes = function(num.markers, file.path.to.write){
    fg = c("marker\tP1_1\tP1_2\tP2_1\tP2_2", paste0("marker", 1:num.markers, "\tA\tA\tB\tB"))
    write(fg, file.path.to.write, ncolumns = 1)
}

make.pedigree.file = function(f.gen, number.individuals, file.path.to.write){
    if(f.gen < 2) return("f.gen must be at least 2")
    top.lines = list(c("Name\tParent1\tParent2", "P1\tNA\tNA", "P2\tNA\tNA", "F1_1\tP1\tP2"))
    f2.lines = list(paste0("F2_", 1:number.individuals, "\tF1_1\tF1_1"))
    
    if(f.gen > 2){
        more.lines = lapply(3:f.gen, function(x){
            paste0("F", x, "_", 1:number.individuals, "\tF", (x - 1), "_", 1:number.individuals, "\tF", (x - 1), "_", 1:number.individuals)
        })
        all.lines = c(top.lines, f2.lines, more.lines)
    } else {
        all.lines = c(top.lines, f2.lines)
    }
    
    
    all.lines2 = do.call(c, all.lines)
    write(all.lines2, file.path.to.write, ncolumns = 1)
    
}

make.map.file = function(cm.positions, file.path.to.write){
    lines1 = c("marker\tchromosome\tposition", paste0("marker", 1:length(cm.positions), "\tA\t", cm.positions))
    write(lines1, file.path.to.write, ncolumns = 1)
}

make.chrom.file = function(cm.length, file.path.to.write){
    lines1 = c("chromosome\tlength\tcentromere\tprefPairing\tquadrivalents", paste0("A\t", cm.length, "\t", (cm.length / 2), "\t1.00\t0.0"))
    write(lines1, file.path.to.write, ncolumns = 1)
}

make.pedsim.par.file = function(project.name, file.path.to.write, map.function, gen.file.path, ped.file.path){
    if(missing(map.function)) map.function = "KOSAMBI"
    if(missing(gen.file.path)) gen.file.path = p(full.path, "config.files/", project.name, ".gen")
    if(missing(ped.file.path)) ped.file.path = p(full.path, "config.files/", project.name, ".ped")
    
    full.path = "/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/"
    lines1 = c("PLOIDY = 2", p("MAPFUNCTION = ", map.function), "MISSING = NA",
                         p("CHROMFILE = ", full.path, "config.files/", project.name, ".chrom"),
                         p("MAPFILE = ", full.path, "config.files/", project.name, ".map"),
                         p("FOUNDERFILE = ", gen.file.path),
                         p("PEDFILE = ", ped.file.path),
                         p("OUTPUT = ", full.path, "simresults/", project.name, ".output"))
    write(lines1, file.path.to.write, ncolumns = 1)
}

make.pedsim.par.file.f3.selec = function(project.name, uniqueid, file.path.to.write, map.function){
    if(missing(map.function)) map.function = "KOSAMBI"
    full.path = "/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/"
    lines1 = c("PLOIDY = 2", p("MAPFUNCTION = ", map.function), "MISSING = NA",
                         p("CHROMFILE = ", full.path, "config.files/", project.name, ".chrom"),
                         p("MAPFILE = ", full.path, "config.files/", project.name, ".map"),
                         p("FOUNDERFILE = ", full.path, "config.files/", project.name, ".", uniqueid, "f3.selec.gen"),
                         p("PEDFILE = ", full.path, "config.files/", project.name, ".", uniqueid, "f3.selec.ped"),
                         p("OUTPUT = ", full.path, "simresults/", uniqueid, "_", project.name, ".output"))
    write(lines1, file.path.to.write, ncolumns = 1)
}

make.pedsim.config.files.and.run = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1){
    if(missing(run1)) run1 == F
    if(missing(map.function)) map.function = "KOSAMBI"
    make.founder.genotypes(num.markers, p(write.path.folder, "/config.files/", project.name, ".gen"))
    make.pedigree.file(f.gen, num.individuals, p(write.path.folder, "/config.files/", project.name, ".ped"))
    make.map.file(cm.positions, p(write.path.folder, "/config.files/", project.name, ".map"))
    make.chrom.file(max(cm.positions), p(write.path.folder, "/config.files/", project.name, ".chrom"))
    make.pedsim.par.file(project.name, p(write.path.folder, "/", project.name, ".par"), map.function = map.function)
    
    if(run1 == T){
        system(p("java -jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/PedigreeSim.jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/", project.name, ".par"))    
    }    
}

make.pedsim.par.files.and.run.f3.selec = function(write.path.folder, project.name, uniqueid, map.function, run1){
    if(missing(run1)) run1 == F
    if(missing(map.function)) map.function = "KOSAMBI"
    make.pedsim.par.file.f3.selec(project.name, uniqueid, p(write.path.folder, "/", uniqueid, "_", project.name, ".par"), map.function = map.function)
    
    if(run1 == T){
        system(p("java -jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/PedigreeSim.jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/", uniqueid, "_", project.name, ".par"))    
    }
}

convert.ped.sim.to.zygote = function(ped.sim){
    # ped.sim = read.delim(ped.sim.genotype.path, sep = "\t", stringsAsFactors = F)
    rownames(ped.sim) = ped.sim[, 1]
    #remove marker column
    ped.sim = ped.sim[, -1]
    ped.sim = convert.to.character.data.frame(ped.sim)
    
    #combine gametes
    ped.sim2 = lapply(seq(1, ncol(ped.sim), 2), function(q){
        combined.geno1 = unlist(Map(function(x, y){
            if(x == "A" & y == "A") z = "A"
            if(x == "A" & y == "B") z = "H"
            if(x == "B" & y == "A") z = "H"
            if(x == "B" & y == "B") z = "B"
            z
        }, ped.sim[, q], ped.sim[, (q + 1)]))    
        combined.geno1
    })
    
    ped.sim3 = do.call(cbind, ped.sim2)
    ped.sim3 = as.data.frame((ped.sim3))
    colnames(ped.sim3) = colnames(ped.sim)[seq(1, ncol(ped.sim), 2)]
    ped.sim3 = convert.to.character.data.frame(ped.sim3)
    ped.sim3
}

read.ped = function(ped.path) read.delim(ped.path, sep = "\t", stringsAsFactors = F)

f2.selection.procedure = function(project.name, selection.strength, selection.pos){
    #args:
    # selection.strength - integer, denominator of selection fraction
    # selection.pos - integer, locus at which to apply selection
    
    #if geno.10k doesn't already exists, load it and process it

    if(!exists("geno.10k")){
        geno.10k = read.ped("rotation1scripts_v4/original_data/simulation/pedigreesim/10k.sim/axc1a.pop10k.output_genotypes.dat")
        geno.10k = geno.10k[8:ncol(geno.10k)]
        geno.10k = geno.10k[, which(sapply(geno.10k, function(x) x[[selection.pos]] == "B"))]
    } 
    

    f2.genotypes = read.ped(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", project.name, ".output_genotypes.dat"))
    
    to.replace = names(which(sapply(f2.genotypes[4:ncol(f2.genotypes)], function(x){
        newdice = runif(1)
        threshold1 = 1 / selection.strength
        if(newdice < threshold1){
                dice = 1
            } else {
                dice = 2
            }
        # dice = sample(selection.strength, 1)
        if(dice == 1 & x[[selection.pos]] == "A"){
            return(T)
        } else {
            return(F)
        }
    })))

    replacements1 = geno.10k[, sample(ncol(geno.10k), length(to.replace))]
    

    if(length(to.replace) > 1){
        colnames(replacements1) = to.replace    
    }
    
    f2.genotypes[, to.replace] = replacements1
    write.table(f2.genotypes,
                            p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", project.name, ".output_genotypes.dat_wselec.f2.dat"),
                            sep = "\t",
                            row.names = F,
                            quote = F)
    return()
}

f3.selection.procedure = function(project.name, uniqueid, f.gen, selection.strength, selection.pos){
    #selection procedure for filial generations above F2
    #args:
    # output.genotypes.path - string, path to genotypes produced by PedSim
    # uniqueid - string, unique identifier for files produced by this procedure
    # f.gen - integer indicating filial generation
    # selection.strength - integer, denominator of selection pressure fraction
    # selection.pos - integer, locus at which    to apply selection
    
    # browser()
    
    f.gen2 = p("F", f.gen)
    prev.f.gen = p("F", (f.gen - 1))
    psr.path = "rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/"
    output.genotypes.path = p(psr.path, project.name, ".output_genotypes.dat")
    prev.ped.file.path = p("rotation1scripts_v4/original_data/simulation/pedigreesim/config.files/",
                                                 project.name, ".", f.gen, ".ped")
    
    
    
    ped1 = read.delim(output.genotypes.path, sep = "\t", stringsAsFactors = F)
    
    #first get individuals in which selection via gamete competition is possible (heterozygotes)
    ped.zygote = convert.ped.sim.to.zygote(ped1)
    
    ped.zygote2 = ped.zygote[, grep(f.gen2, colnames(ped.zygote))]
    
    
    selection.ind = which(sapply(ped.zygote2, function(x){
        newdice = runif(1)
        threshold1 = 1 / selection.strength
        if(newdice < threshold1){
                dice = 1
            } else {
                dice = 2
            } 
        # dice1 = sample(selection.strength, 1)
        x[[selection.pos]] == "H" & dice1 == 1
    }))
    
    if(length(selection.ind) > 0){
        
        #### SELECTION PROCEDURE ####
        
        #get columns names of all gametes in which selection is possible
        snames = strsplit(names(selection.ind), "_")
        snames2 = unlist(lapply(snames, function(x){
            c(paste0(x[[1]], "_", x[[2]], "_", 1),
                paste0(x[[1]], "_", x[[2]], "_", 2))
        }))
        
        snames.prev.filial.gen = gsub(f.gen2, prev.f.gen, snames2)
        snames.prev.filial.gen.trunc = gsub("_[1-2]$", "", snames.prev.filial.gen)
        
        #prepare founder genotypes file
        ped.selec = ped1[, snames.prev.filial.gen]
        marker = ped1$marker
        ped.selec2 = as.data.frame(cbind(marker, ped.selec))
        
        #make new founder genotypes file with lines that are to undergo selection
        write.table(ped.selec2, p("rotation1scripts_v4/original_data/simulation/pedigreesim/config.files/", project.name, ".", uniqueid, "f3.selec.gen"), sep = "\t", row.names = F, quote = F)
        
        #process pedigree file
        pedigree.file1 = read.delim(prev.ped.file.path, sep = "\t", stringsAsFactors = F)
        pedigree.file2 = pedigree.file1[which(pedigree.file1$Parent1 %in% snames.prev.filial.gen.trunc), ]
        pedigree.file3 = pedigree.file2[match(sort(rep(pedigree.file2$Parent1, 5)), pedigree.file2$Parent1), ]
        
        
        # if(nrow(pedigree.file3) == 0) browser()
        
        pedigree.file3 = reset.rownames(pedigree.file3)
        pedigree.file3$Name = paste0(pedigree.file3$Name, ".", rep(1:5, nrow(pedigree.file3) / 5))
        
        # pedigree.file3$Parent1 = gsub("_", ".", pedigree.file3$Parent1)
        # pedigree.file3$Parent2 = gsub("_", ".", pedigree.file3$Parent2)
        
        print("ped.parents")
        
        ped.parents = data.frame(unique(pedigree.file3$Parent1), rep(NA, length(unique(pedigree.file3$Parent1))), rep(NA, length(unique(pedigree.file3$Parent1))))
        # if(is.null(ncol(ped.parents))) browser()
        ped.parents = reset.colnames(ped.parents)
        # if(is.null(ncol(ped.parents))) browser()
        colnames(ped.parents) = c("Name", "Parent1", "Parent2")
        
        pedigree.file4 = rbind(ped.parents, pedigree.file3)
        
        #make new pedigree file with lines that are to undergo selection
        write.table(pedigree.file4, p("rotation1scripts_v4/original_data/simulation/pedigreesim/config.files/", project.name, ".", uniqueid, "f3.selec.ped"), sep = "\t",
                                row.names = F, quote = F)
        
        perform.again = T
        counter1 = 1
        while(perform.again == T){
            make.pedsim.par.files.and.run.f3.selec("rotation1scripts_v4/original_data/simulation/pedigreesim/", project.name, uniqueid, "KOSAMBI", T)
            
            # if(counter1 > 6) browser()
            
            #run pedigreesim with new pedigree and founder files
            # browser()
            
            selection.test = read.ped(p(psr.path, uniqueid, "_", project.name, ".output_genotypes.dat"))
            
            selection.test.f3 = selection.test[, grep(f.gen2, colnames(selection.test))]
            
            selection.test.f3.2 = selection.test.f3[, which(sapply(selection.test.f3, function(x){
                x[[selection.pos]] == "B"
            }))]
            
            # if(is.null(colnames(selection.test.f3.2))) browser()
            
            
            
            colnames(selection.test.f3.2) = gsub("\\.[1-5]", "", colnames(selection.test.f3.2))
            trunc.colnames1 = gsub("_[1-2]$", "", colnames(selection.test.f3.2))
            #only need one of each gamete 
            selection.test.f3.3 = as.data.frame(selection.test.f3.2[, match(unique(trunc.colnames1), trunc.colnames1)])
            selection.test.f3.3 = convert.to.character.data.frame(selection.test.f3.3)
            colnames(selection.test.f3.3) = unique(trunc.colnames1)
            
            # q1 = gsub("_[1-2]$", "", colnames(selection.test.f3.3))
            
            #do we have gametes with a B at locus 200 for all of the lines we performed selection on ?
            perform.again = !all(unique(gsub("_[1-2]$", "", snames2)) %in% colnames(selection.test.f3.3))
            counter1 = counter1 + 1
        }
        
        #grab columns of ped1 (original genotype file) to replace
        original.colnames.to.replace = names(which(sapply(ped1[, snames2], function(x) x[[selection.pos]] == "A")))
        
        #check columns match, if so replace original genotyping data with new selected genotyping data.
        if(all(gsub("_[1-2]$", "", colnames(ped1[, original.colnames.to.replace])) == gsub("_[1-2]$", "", colnames(selection.test.f3.3)))){
            colnames(selection.test.f3.3) = colnames(ped1[, original.colnames.to.replace])
            #perform replacement
            ped1[, original.colnames.to.replace] = selection.test.f3.3
        }
        
        write.table(ped1, p(output.genotypes.path, "_wselec.f", f.gen, ".dat"), sep = "\t",
                                row.names = F, quote = F)
        
        
        #### END OF SELECTION PROCEDURE ####
    } else {
        write.table(ped1, p(output.genotypes.path, "_wselec.f", f.gen, ".dat"), sep = "\t",
                                row.names = F, quote = F)
    }
    
    
    
    # browser()
}

make.files.next.iteration = function(project.name, f.gen){
    #this function assumes one of the selection procedure functions has been run previously
    
    prev.fgen = f.gen - 1
    
    prev.geno = read.ped(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", 
                                                 project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".dat"))
    
    
    to.keep = colnames(prev.geno)[grep(p("F", prev.fgen), colnames(prev.geno))]
    
    parents1 = unique(gsub("_[1-2]$", "", to.keep))
    pedigree1 = data.frame(parents1, NA, NA)
    
    new.gen = gsub(p("F", prev.fgen), p("F", f.gen), parents1)
    
    ped2 = data.frame(new.gen, parents1, parents1)
    # if(is.null(colnames(pedigree1))) browser()
    colnames(pedigree1) = c("Name", "Parent1", "Parent2")
    # if(is.null(colnames(ped2))) browser()
    colnames(ped2) = c("Name", "Parent1", "Parent2")
    ped2 = ped2[grep(p("F", f.gen), ped2$Name), ]
    
    pedigree1 = convert.to.character.data.frame(pedigree1)
    ped2 = convert.to.character.data.frame(ped2)
    ped3 = rbind(pedigree1, ped2)
    # print("before row reset 1")
    # tryCatch(function(g) reset.rownames(ped3), error = function(x) browser())
    ped3 = reset.rownames(ped3)
    write.table(ped3,
                            p("rotation1scripts_v4/original_data/simulation/pedigreesim/config.files/", project.name, ".", f.gen, ".ped"),
                            sep = "\t",
                            row.names = F,
                            quote = F)
    
    if(f.gen > 2){
        gen.file1 = read.ped(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/",
                                                     project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".dat"))
        
        gen.file2 = gen.file1[, c(1, grep(p("F", prev.fgen), colnames(gen.file1)))]
        
        
        
        # newcolnames1 = strsplit(colnames(gen.file2), "_")
        # colnames(gen.file2) = sapply(newcolnames1, function(x) p(x[1], ".", x[2], "_", x[3]))
        # colnames(gen.file2)[1] = "marker"
        
        write.table(gen.file2,
                                p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/",
                                    project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".trunc.gen"),
                                sep = "\t",
                                row.names = F,
                                quote = F)
        
        make.pedsim.par.file(project.name = project.name,
                                                 file.path.to.write = p("rotation1scripts_v4/original_data/simulation/pedigreesim/", project.name, ".f", f.gen, ".par"),
                                                 map.function = "KOSAMBI",
                                                 gen.file.path = p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", 
                                                                                     project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".trunc.gen"),
                                                 ped.file.path = p("rotation1scripts_v4/original_data/simulation/pedigreesim/config.files/", project.name, ".", f.gen, ".ped")
        )
        
    } else {
        make.pedsim.par.file(project.name = project.name,
                                                 file.path.to.write = p("rotation1scripts_v4/original_data/simulation/pedigreesim/", project.name, ".f", f.gen, ".par"),
                                                 map.function = "KOSAMBI",
                                                 gen.file.path = p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", 
                                                                                     project.name, ".output_genotypes.dat_wselec.f", prev.fgen, ".dat"),
                                                 ped.file.path = p("rotation1scripts_v4/original_data/simulation/pedigreesim/config.files/", project.name, ".", f.gen, ".ped")
        )
    }
    
    
    
    
    
    # browser()
    
    system(p("java -jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/PedigreeSim.jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/", project.name, ".f", f.gen, ".par"))
    
    
}

clean.files = function(){
    write.path.folder = "rotation1scripts_v4/original_data/simulation/pedigreesim/"
    files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = "_founderalleles.dat")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = "_alleledose.dat")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = ".hsa")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = ".hsb")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".gen")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".map")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".ped")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
    
    files.to.rm = list.files(p(write.path.folder, "config.files/"), pattern = ".chrom")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "config.files/", x)))
    
    files.to.rm = list.files(p(write.path.folder),    pattern = ".par")
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, x)))
}

complete.f5.selection.procedure = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder,
                                                                                     project.name, map.function, run1, selection.strength, selection.pos){

    make.pedsim.config.files.and.run(num.markers, num.individuals, 2, cm.positions, write.path.folder,
                                                                     project.name, map.function, run1)

    f2.selection.procedure(project.name, selection.strength, selection.pos)

    if(f.gen > 2){
        lapply(3:f.gen, function(x){
            make.files.next.iteration(project.name, x)
            unique.code = paste0(LETTERS[sample(26, 5)], collapse = "", sample(9, 4))
            f3.selection.procedure(project.name, unique.code, x, selection.strength, selection.pos)
            print("done 1")
            
            files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = unique.code)
            lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
            
            files.to.rm2 = list.files(p(write.path.folder, "config.files/"), pattern = unique.code)
            lapply(files.to.rm2, function(x) file.remove(p(write.path.folder, "config.files/", x)))
            
            files.to.rm3 = list.files(p(write.path.folder), pattern = unique.code)
            lapply(files.to.rm3, function(x) file.remove(p(write.path.folder, x)))
            
        })
    }
    
    
    
    
    files.to.rm = list.files(p(write.path.folder, "simresults/"), pattern = p(project.name, "\\.output_genotypes\\.dat"))
    files.to.rm = files.to.rm[-grep(p("wselec", "\\.f", f.gen, "\\.dat"), files.to.rm)]
    lapply(files.to.rm, function(x) file.remove(p(write.path.folder, "simresults/", x)))
}


complete.fgen.no.selec.procedure = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1){
    make.pedsim.config.files.and.run(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1)
}

make.pedsim.config.files.and.run = function(num.markers, num.individuals, f.gen, cm.positions, write.path.folder, project.name, map.function, run1){
    if(missing(run1)) run1 == F
    if(missing(map.function)) map.function = "KOSAMBI"
    make.founder.genotypes(num.markers, p(write.path.folder, "/config.files/", project.name, ".gen"))
    make.pedigree.file(f.gen, num.individuals, p(write.path.folder, "/config.files/", project.name, ".ped"))
    make.map.file(cm.positions, p(write.path.folder, "/config.files/", project.name, ".map"))
    make.chrom.file(max(cm.positions), p(write.path.folder, "/config.files/", project.name, ".chrom"))
    make.pedsim.par.file(project.name, p(write.path.folder, "/", project.name, ".par"), map.function = map.function)
    
    if(run1 == T){
        system(p("java -jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/PedigreeSim.jar ~/project.phd.main/rotation1scripts_v4/original_data/simulation/pedigreesim/", project.name, ".par"))    
    }    
}

run.1000.sims = function(project.name, f.gen, num.sims, w.selec, selec.str, selec.pos, pop.size, cm.pos, perform.analysis.flag){
    if(missing(num.sims)) num.sims = 1000
    if(missing(w.selec)) w.selec = F
    if(missing(selec.str)) selec.str = 10
    if(missing(selec.pos)) selec.pos = 200
    if(missing(cm.pos)) cm.pos = axc.chr1a.cm
    if(missing(perform.analysis.flag)) perform.analysis.flag = T
    
    t1 = mclapply(1:num.sims, function(x){
        if(w.selec == T){
            complete.f5.selection.procedure(length(cm.pos), pop.size, f.gen, cm.pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", paste0(project.name, x), "KOSAMBI", T, selec.str, selec.pos)    
        } else {
            complete.fgen.no.selec.procedure(length(cm.pos), pop.size, f.gen, cm.pos, "rotation1scripts_v4/original_data/simulation/pedigreesim/", paste0(project.name, x), "KOSAMBI", T)
        }
    }, mc.cores = 60)    
    
    clean.files()    
    
    # browser()
    
    if(w.selec == T){
        files.to.combine = list.files("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", 
                                                                    pattern = p(project.name, ".*output_genotypes.dat_wselec.f", f.gen, ".dat"))    
    } else {
        files.to.combine = list.files("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", 
                                                                    pattern = p(project.name, ".*output_genotypes.dat"))
    }
    
    
    all.sim.files = mclapply(files.to.combine, function(x){
        
        
        
        g = read.ped(p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", x))
        g = convert.ped.sim.to.zygote(g)
        g = reset.rownames(g)
        rownames(g) = paste0("marker", rownames(g))
        g = as.data.frame(t(g))
        g = g[grep(p("F", f.gen), rownames(g)), ]
    }, mc.cores = 60)
    
    if(class(all.sim.files[[1]]) != "data.frame"){
        all.sim.files = lapply(all.sim.files, as.data.frame)
    }
    
    # browser()
    
    if(f.gen == 2){
        f2flag = T
    } else {
        f2flag = F
    }
    
    # browser()
    
    s(all.sim.files, p("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", project.name), "NA")
    
    if(perform.analysis.flag == T){
        if(w.selec == F){
            # perform.analysis = function(geno.list, selection, selection.position, selection.strength, recombination.position, F2, fgen)
            analysis1 = perform.analysis(all.sim.files, recombination.position = "AxC1A", F2 = f2flag, fgen = f.gen)
        } else {
            analysis1 = perform.analysis(all.sim.files, selection = T, selection.position = selec.pos, selection.strength = selec.str, recombination.position = "AxC1A", F2 = f2flag, fgen = f.gen)
        }
        
        write.csv(analysis1, p("rotation1scripts_v4/original_data/simulation/pedigreesim/analyses/", project.name, ".csv"), row.names = F)
    }
    

    
    
    files.to.rm = list.files("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", pattern = "genotypes.dat")
    lapply(paste0("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/", files.to.rm), file.remove)
    
}