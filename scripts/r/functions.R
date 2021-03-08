#written by Alex Coulton - alex.coulton@bristol.ac.uk

#writing a function to parse BLAST information
#converts a list of markers into percentage of physical distance coverage

#SETUP

# setwd("C:/Users/ac14037/project.phd.main/")
# setwd("~/project.phd.main/")

options(stringsAsFactors = FALSE)
options(scipen = 999)

#FUNCTIONS

put.col.first = function(df1, col.name.to.move){
    #move a column to the first position in a dataframe
    g = which(colnames(df1) == col.name.to.move)
    if(length(g) == 0) stop("Column not found")
    g2 = 1:ncol(df1)
    
    g2 = g2[-which(g2 == g)]
    
    df1[c(g, g2)]
    
}

h2 = function(x){
	rows1 = nrow(x)
	cols1 = ncol(x)
	if(rows1 >= 10) rows1 = 10
	if(cols1 >= 10) cols1 = 10
    x[1:rows1, 1:cols1]
}

p = function(...){
    #quicker paste function. must supply vector of 
    #elements to paste together (i.e. with c())
    paste(..., collapse = "", sep = "")
}

add.comparison.bars = function(ggplot.obj, tukey.obj, y.start, treatment.names1){
    #function to add comparison bars to a barplot with geom_signif given a tukey posthoc test from an ANOVA
    tuk.comp.p = unclass(tukey.obj)[[1]][, 4]
    tuk.comps = rownames(unclass(tukey.obj)[[1]])
    tuk.comps = strsplit(tuk.comps, "-")
    
    tuk.sig.coord = which(tuk.comp.p < 0.05)
    tuk.comp.p = tuk.comp.p[tuk.sig.coord]
    tuk.comps = tuk.comps[tuk.sig.coord]
    
    tuk.annotations = unname(unlist(lapply(tuk.comp.p, function(x){
        if(x < 0.0001) return("****")
        if(x < 0.001) return("***")
        if(x < 0.01) return("**")
        if(x < 0.05) return("*")
    })))
    
    #can get treatment names from ggplot object
    # treatment.names1 = unclass(ggplot.obj)$data$treatment
    tuk.comps = lapply(tuk.comps, function(x){
        treatment.names1[as.numeric(x)]
    })
    
    # browser()
    
    for(i in 1:length(tuk.comps)){
        ggplot.obj = ggplot.obj + geom_signif(comparisons = list(tuk.comps[[i]]), y_position = y.start + (i * 5), annotation = tuk.annotations[i])
    }
    
    ggplot.obj
}


calculate.mean.diff.neighbours = function(x){
    #calculates the difference between this element and both of its neighbouring elements,
    #then returns the mean of this value
    #args:
    # x - a numeric vector
    differences = x
    for(i in 2:(length(x)-1)){
        differences[i] = (abs(x[i] - x[(i-1)]) + abs(x[i] - x[(i+1)]))/2
    }
    differences[1] = NA
    differences[length(differences)] = NA
    return(differences)
}

genetictophysical=function(chromosome, listofmarkers, markerformat, threshold, return.bp, nrgene){
    if(!("dplyr" %in% .packages())) library(dplyr)
    #args:
    # chromosome - a character string in the following format: "3B"
    # listofmarkers - a character vector of marker names
    # markerformat - a character string indicating listofmarkers contain "." or "-" to seperate AX from number 
    # threshold - an integer indicating at what percentage difference from neighbours should the function look for alternative hits
    # return.bp - Boolean flag: if T, function returns a vector of bp positions of markers instead of percentages
    
    
    if(missing(return.bp)){
        return.bp = F
    }
    if(missing(markerformat)) markerformat = "-"
    if(missing(threshold)) threshold = 30
    if(missing(nrgene)) nrgene = F
    
    if(nrgene == F){
        # allprobegenomeblast=read.csv("rotation1scripts_v4/processed_data/BLAST/820k/820k.full.unique.blast.4brev.ordered.csv")
        #create new variable in the global environment containing probe BLAST table
        if(exists("allprobegenomeblast") == F) allprobegenomeblast <<- read.csv("rotation1scripts_v4/processed_data/BLAST/820k/820k.full.unique.blast.ordered.csv", stringsAsFactors = F)
        if(exists("chromosomecounts") == F) chromosomecounts <<- read.table("rotation1scripts_v4/original_data/iwgsc.chromosome.counts.txt", stringsAsFactors = F)
    } else {
        allprobegenomeblast=read.table("bioinf/blast/nr_gene/results.blast/probesvsnrgene.blast")
        chromosomecounts=read.table("bioinf/blast/nr_gene/nrgene.chromo.lengths.txt")
        chromosomecounts = chromosomecounts[, 2:1]
        # browser()
    }
    
    allprobegenomeblast = reset.colnames(allprobegenomeblast)
    
    if(markerformat=="."){
    allprobegenomeblast$V1=gsub("-",".",allprobegenomeblast$V1)
    }
    blast_chr=filter(allprobegenomeblast, V2==paste("chr",chromosome,sep=""))
    blast_chr=filter(blast_chr, V11<1e-19)
    blast_chr_map=blast_chr[match(listofmarkers, blast_chr[,1]),]
    
    countsindex=grep(paste("chr",chromosome,sep=""), chromosomecounts[,2])
    
    blast_chr_map$percentage=(as.numeric(as.character(blast_chr_map[,9]))/as.numeric(as.character(chromosomecounts[countsindex, 1])))*100
    
    differences = calculate.mean.diff.neighbours(as.numeric(blast_chr_map$percentage))
    max.difference = max(differences, na.rm = T)
    max.diff.index = which(differences == max.difference)
    
    if(max.difference > threshold){
        #if there is an abnormally large difference between the physical distance of a marker and its neighbours,
        #have a look for alternative BLAST hits for that same marker that are a better fit.     
        big.index = max.diff.index
        grep.indexes = grep(blast_chr_map[big.index, 1], allprobegenomeblast$V1)
        possible.hits = allprobegenomeblast[grep.indexes,]
        possible.hits$percentage = (as.numeric(as.character(possible.hits[,9]))/as.numeric(as.character(chromosomecounts[countsindex, 1])))*100
        possible.hits$per.diff.comp = possible.hits$percentage - blast_chr_map$percentage[big.index]
        possible.hits$per.diff.comp = abs(as.numeric(possible.hits$per.diff.comp))
        smallest.diff = max(possible.hits$per.diff.comp)
        best.hit.index = which(possible.hits$per.diff.comp == smallest.diff)[1]
        new.hit.assignment = possible.hits[best.hit.index,1:(ncol(possible.hits)-1)]
        blast_chr_map[big.index,] = new.hit.assignment
    }
    
    return.val = blast_chr_map$percentage
    attr(return.val, "bp") = (0.01*blast_chr_map$percentage)*as.numeric(chromosomecounts[countsindex, 1])
    
    return(return.val)
    
}

grab.individual.hit = function(markername, chromosome){
    allprobegenomeblast=read.table("bioinf/blast/probe.vs.genome.blast/results.blast/820kprobes.vs.iwgsc.4b.rev.blast")
    if(missing(chromosome) == F){
        allprobegenomeblast = filter(allprobegenomeblast, V2 == paste("chr", chromosome, sep = ""))
    }
    grep.coords = grep(markername, allprobegenomeblast$V1)
    allprobegenomeblast[grep.coords,]
}

interpolate.na.values = function(x){
    lapply(x, function(p){
        if(is.na(p) == T){
            
        }
    })
}

#make a new dataframe
newdf = function(..., no.rows){
    #...: an unlimited number of column names
    #no.rows: boolean value, should the dataframe contain zero rows?
    if(missing(no.rows)) no.rows = F
    df=as.data.frame(matrix(nrow=1, ncol=length(c(...))))
    df[is.na(df)]=""
    colnames(df)=c(...)
    
    if(no.rows == T){
        df = df[-1, ]
    }
    
    return(df)
}

listofwheatchromosomes=c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D")

dataframe.replace = function(dataframe1, val1, replacement1){
    #replace all instances of val1 in a dataframe with replacement1
    #args:
    # dataframe1 - a dataframe
    # val1 - the value to replace
    # replacement1 - the replacement value
    dataframe1[] = lapply(dataframe1, function(x){
        x[which(x == val1)] = replacement1
        x
    })
    
    dataframe1
}

add.column.at.position = function(dataframe1, after_pos){
    #adds a new column to a dataframe, populated 1/2 with "A" and 1/2 "B"
    #args:
    # dataframe1 - arbitrary dataframe to add column to
    # after_pos - column coordinate to insert new column after
    if((nrow(dataframe1) %% 2) == 0){
        a1 = rep("A", (nrow(dataframe1) / 2))
        b1 = rep("B", (nrow(dataframe1) / 2))
        empty.col = as.character(c(a1, b1))
    } else {
        a1 = rep("A", ceiling(nrow(dataframe1) / 2))
        b1 = rep("B", floor(nrow(dataframe1) / 2))
        empty.col = as.character(c(a1, b1))
    }
    
    empty.col = as.character(rep("", nrow(dataframe1)))
    if(after_pos == 0){
        dataframe2 = cbind(empty.col, dataframe1[, 1:ncol(dataframe1)])
    } else {
        dataframe2 = cbind(dataframe1[, 1:after_pos], empty.col, dataframe1[, (after_pos+1):ncol(dataframe1)])    
    }
    
    dataframe2 = convert.to.character.data.frame(dataframe2)
    dataframe2
}

add.row.at.position = function(dataframe1, after_pos){
    #adds a new row to a dataframe, populated 1/2 with "A" and 1/2 "B"
    #args:
    # dataframe1 - arbitrary dataframe to add column to
    # after_pos - column coordinate to insert new column after
    if((ncol(dataframe1) %% 2) == 0){
        a1 = rep("A", (nrow(dataframe1) / 2))
        b1 = rep("B", (nrow(dataframe1) / 2))
        empty.row = as.character(c(a1, b1))
    } else {
        a1 = rep("A", ceiling(nrow(dataframe1) / 2))
        b1 = rep("B", floor(nrow(dataframe1) / 2))
        empty.row = as.character(c(a1, b1))
    }
    
    # browser()
    
    empty.row = as.character(rep("", ncol(dataframe1)))
    if(after_pos == 0){
        dataframe2 = rbind(empty.row, dataframe1[1:nrow(dataframe1), ])
    } else {
        dataframe2 = rbind(dataframe1[1:after_pos, ], empty.row, dataframe1[(after_pos+1):nrow(dataframe1), ])    
    }
    
    dataframe2 = convert.to.character.data.frame(dataframe2)
    dataframe2
}



convert.aas.to.rqtl = function(aasoutputfile, parents, call.code.numeric){
    #aasoutputfile - .txt file generated by axiom analysis suite
    #parents (optional) - character vector of parental variety names
    if(missing(call.code.numeric)) call.code.numeric = F
    if(missing(parents)) parents = c("Paragon", "Apogee")
    if(!"readr" %in% (.packages())) library(readr)
    
    
    # convert output of axiom analysis suite to a data.frame suitable for analysis with rQTL
    nd = read_delim(aasoutputfile, "\t", col_names = F, comment = "#")
    # nd = read.delim(aasoutputfile, sep = "\t", header = T, stringsAsFactors = F)
    nd = as.data.frame(nd)
    nd=t(nd)
    # colnames(nd) = nd[1, ]
    # nd = nd[-1, ]
    
    if(call.code.numeric == T){
        nd = convert.to.character.data.frame(nd)
        nd = dataframe.replace(nd, "0", "A")
        nd = dataframe.replace(nd, "1", "H")
        nd = dataframe.replace(nd, "2", "B")
        nd = dataframe.replace(nd, "-1", "-")
        
        #this code stalled for some reason with certain files, so I made a "dataframe.replace()" function instead 11/01/2019
        # nd[nd=="0"] = "A"
        # nd[nd=="1"] = "H"
        # nd[nd=="2"] = "B" 
        # nd[nd=="-1"] = "-"    
    } else {
        nd[nd=="AA"] = "A"
        nd[nd=="BB"] = "B" 
        nd[nd=="AB"] = "H"
        nd[nd=="NoCall"] = "-"    
    }
    
    
    
    colnames(nd)=nd[1,]
    nd[1,]=1
    nd=cbind(nd[,1], nd)
    nd[,1]=seq(0,nrow(nd)-1,1)
    nd[1,1]=""
    nd[1,2]=""
    
    #bring parents to top of dataframe
    if(length(grep(parents[[1]], nd[,2])) == 0){
     nd2 = nd 
    } else {
    nd2=rbind(nd[grep(parents[[1]],nd[,2]),], nd[-c(grep(parents[[1]],nd[,2])),])
    }
    
    if(length(grep(parents[[2]], nd[,2])) == 0){
        nd3 = nd2 
    } else {
        nd3=rbind2(nd2[grep(parents[[2]],nd2[,2]),], nd2[-c(grep(parents[[2]],nd2[,2])),])
    }
    
    rqtl.chromo.row.coord = which(nd3[, 3] == 1)
    
    nd3=rbind(nd3[rqtl.chromo.row.coord,],nd3[-rqtl.chromo.row.coord,])
    
    nd3[,1]=seq(0,nrow(nd)-1,1)
    nd3[1,1]=""
    nd3[1,2]=""
    
    rownames(nd3)=seq(1,nrow(nd3),1)
    nd3 = as.data.frame(nd3)
    return(nd3)
}

v = function(x){
    #Rstudio v1.1 is very slow at loading large dataframes with View(). 
    #This function should be used instead of View() for previewing dataframes
    dfname = deparse(substitute(x))
    
    if(ncol(x) > 30){
        if(nrow(x) > 50){
            View(x[1:50, 1:30], title = paste(dfname, "[1:50, 1:30]", sep = ""))
            print("Truncating to 1:50, 1:30")
        } else {
            View(x[1:nrow(x), 1:30], title = paste(dfname, "[1:", nrow(x), ", 1:30]", sep = ""))
            print("Truncating to 1:nrow(x), 1:30")
        }
    } else {
        View(x, title = dfname)
    }
}

combine.list.of.data.frames = function(list.name){
    g = newdf(colnames(list.name[[1]]), no.rows = T)
    for(i in 1:length(list.name)){
        g = rbind(g, list.name[[i]])
    }
    return(g)
}

conv.listoflists.to.df = function(x) reset.rownames(as.data.frame(t(as.data.frame(x))))

split.df.into.list.of.dfs.by.column = function(dataframe.input, colname.to.split, do.sort){
    #takes a dataframe and returns a list of dataframes, each with a unique value of a column of your choice
    #args:
    #dataframe.input - dataframe to split
    #colname.to.split - character string of the name of the column by which to split the dataframe
        if(missing(do.sort)) do.sort = T
        if(do.sort == T){
            g = sort(unique(dataframe.input[[colname.to.split]]))
        } else {
            g = unique(dataframe.input[[colname.to.split]])
        }
        
        
        count1 = make.counter()
        q = lapply(g, function(x){
            # print(count1())
            dataframe.input[which(dataframe.input[[colname.to.split]] == x), ]
        })
        
        names(q) = g
        return(q)
        
}

subapply = function(data.frame.to.edit, colname.to.split.df.by, function.to.apply){
    #Splits a dataframe D into a list of dataframes L of length n, with n being the number
    #of unique elements contained in the column specified. Each element of L therefore contains
    #all rows of D that contain a particular unique value. Then applies function.to.apply
    #to each element of L and combines all of these elements into a new dataframe, which is returned
    #
    #args: 
    #data.frame.to.edit: a dataframe
    #colname.to.split.df.by: character vector specifying the column name by which to split D
    #function.to.apply: the function to be applied to L
    g = split.df.into.list.of.dfs.by.column(data.frame.to.edit, colname.to.split.df.by)
    g2 = lapply(g, function.to.apply)
    g3 = combine.list.of.data.frames(g2)
    return(g3)
}

list.consecutive.subsequences = function(x, gap.size){
    #args:
    #gap.size: how big should the difference between consecutive elements in the vector
        #be before they are split into different elements of the output list?
    if(missing(gap.size)) gap.size = 1
    split(x, cumsum(c(1, diff(x) > gap.size)))
} 

#     ____________________________________________________________________________
#     BIOSTRINGS FUNCTIONS                                                                                            ####

list.of.stringset.to.stringset = function(x){
    #args:
    # x - a list of DNAStringSet objects
    newset = DNAStringSet()
    for(i in x){
        newset = c(newset, i)
    }
    newset
}

list.of.aastringset.to.stringset = function(x){
    #args:
    # x - a list of DNAStringSet objects
    newset = AAStringSet()
    for(i in x){
        newset = c(newset, i)
    }
    newset
}

list.of.dnastring.to.stringset = function(x){
    #args:
    # x - a list of DNAStringSet objects
    newset = DNAString()
    for(i in x){
        newset = c(newset, i)
    }
    newset
}


#     ____________________________________________________________________________
#     CHROMOSOME IDENTIFICATION                                                                                             ####
#functions for making chromosome assignments to linkage groups

grab.consensus.chromosome = function(list.of.markers, desktop){
    #args: list.of.markers - character vector of marker names in "." format 
    #examines cerealsdb nullisomic line consensus data and returns a table of the most common 
    #chromosome assignment values for this group of markers 
    if(missing(desktop)) desktop = F
    if(("rs2" %in% ls(.GlobalEnv)) == F){
        library(RMySQL)    
        if(desktop == T){
            mydb = setup.mysql.connection(T)    
        } else {
            mydb = setup.mysql.connection()    
        }
        
        
        rs = dbSendQuery(mydb, "SELECT affycode, consensus FROM axiom820") #grab Avalon genotyping data for 35k array
        rs2 <<- fetch(rs, n=-1)
    }
    
    g = rs2$consensus[match(list.of.markers, rs2$affycode)]
    g2 = na.omit(g[-c(which(g == ""), which(g == "none"))])
    g3 = gsub("L", "", g2)
    g4 = gsub("S", "", g3)
    
    g4 = g4[-grep("^\\d$", g4)]
    
    table(g4)
}

find.best.physical.match = function(list.of.markers){
    #provide with a list of probe names in "." format
    #returns the name of the chromosome for which there are the most BLAST hits
    #args: list.of.markers - character vector of marker names
    phys.lists = lapply(listofwheatchromosomes, function(x){
        genetictophysical(x, list.of.markers, ".")
    })
    
    num.na = unlist(lapply(phys.lists, function(x) length(which(is.na(x)))))
    min.na = min(num.na)
    listofwheatchromosomes[which(num.na == min.na)]
}


#     ____________________________________________________________________________
#     GENE ANNOTATION / FUNCTIONAL ANALYSIS                                                                     ####

check.gene.density = function(chromosome, start.region.bp, end.region.bp, return.df){
    if(missing(return.df)) return.df = F
    if(("genes" %in% ls(.GlobalEnv)) == F) genes = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.gff")
    genes2d = genes[genes$V1 == p("chr", chromosome), ]
    
    if(return.df == F){
        #how many genes are in this particular region from cs.x.p.seg.dist.markers2    
        return(nrow(genes2d[genes2d$V4 > start.region.bp & genes2d$V4 < end.region.bp, ]))
    } else {
        return(genes2d[genes2d$V4 > start.region.bp & genes2d$V4 < end.region.bp, ])
    }
}

grab.genes.and.functions = function(chromosome, start.bp, end.bp){
    func.temp = lapply(check.gene.density(chromosome, start.bp, end.bp, T)$gene.names, function(x) grep(x, functional.anno$`Gene-ID`))
    genes.func = functional.anno[na.omit(unlist(lapply(func.temp, function(x) x[1]))), ]
    genes.func = convert.to.character.data.frame(genes.func)
    return(genes.func)
}

generate.gff3.for.sequence.extraction = function(chromosome, bp.start, bp.end, description, file.name){
    main.path = "rotation1scripts_v4/processed_data/sequence.extraction.gff.files/"
    g = paste(p("chr", chromosome), "IWGSC_March2017", "sequence.of.interest", bp.start, bp.end, ".", "+", ".", p("ID=", description), sep = "\t")
    
    outputgff.file = file(p(main.path, file.name, ".gff3"), "wb") #open a binary connection to write UNIX line ending on Windows
    outputscript.file = file(p(main.path, file.name, ".sh"), "wb")
    
    writeLines(g, outputgff.file)
    close(outputgff.file)
    extraction.script = p("bedtools getfasta -fi ~/project.phd.main/genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta -bed /home/ac14037/project.phd.main/rotation1scripts_v4/processed_data/sequence.extraction.gff.files/", file.name, ".gff3 -fo /home/ac14037/project.phd.main/rotation1scripts_v4/processed_data/sequence.extraction.gff.files/extracted.sequences/", file.name, ".fa")
    write(extraction.script, outputscript.file)
    close(outputscript.file)
    print(p("Saved gff and extraction script to: rotation1scripts_v4/processed_data/sequence.extraction.gff.files/", file.name))
}

#     ____________________________________________________________________________
#     GENETIC MAP MANIPULATION                                                                                                ####

convert.rqtl.to.asmap.format = function(rqtl.df){
    #Convert a dataframe of genotyping data in the format of alexq1_man_cur_sk_map_genogeno.csv for use with mstmap.data.frame (ASMap function)
    rqtl.df = rqtl.df[-1, ]
    rqtl.df = rqtl.df[, -1]
    rqtl.df = as.data.frame(t(rqtl.df))
    rqtl.df = convert.to.character.data.frame(rqtl.df)
    colnames(rqtl.df) = rqtl.df[1, ]
    rqtl.df = rqtl.df[-1, ]
    rqtl.df
} 

calculate.cm.diff = function(list.of.centimorgan.dist) c(0, diff(list.of.centimorgan.dist))

check.if.inverted = function(list.of.bp.pos){
    #check if the physical positions of markers in a genetic map are inverted compared to the IWGSC refseq
    #returns T if they are inverted, F if not.
    #args:
    #list.of.bp.pos: numeric vector of physical positions of markers in base pairs
    
    start.mean = mean(na.omit(list.of.bp.pos[1:6]))
    end.mean = mean(na.omit(list.of.bp.pos[(length(list.of.bp.pos)-5):length(list.of.bp.pos)]))
    if(start.mean > end.mean){
        return(T)
    } else {
        return(F)
    }
}

check.if.inverted2 = function(gen.map, cm.colname, phys.pos.colname){
    #updated, better version of check.if.inverted() (uses differences between middle marker in bins instead of 
    #first and last 6 markers)
    #check if the physical positions of markers in a genetic map are inverted compared to the IWGSC refseq
    #returns T if they are inverted, F if not.
    #args:
    #list.of.bp.pos: numeric vector of physical positions of markers in base pairs
    gen.df4 = split.df.into.list.of.dfs.by.column(gen.map, cm.colname, F)
    q = sapply(gen.df4, function(y){
        sort(y[[phys.pos.colname]])[ceiling((nrow(y) / 2))]
        # browser()
    })
    
    # browser()
    diffs1 = diff(na.omit(q))
    len.negative = length(which(diffs1 < 0))
    len.positive = length(which(diffs1 > 0))
    
    #if the overall direction of markers on the chromosome is from high to low,
    #order markers within each bin so that they are descending.
    if(len.negative > len.positive){
        inverted = T
    } else {
        inverted = F
    }
    
    inverted
}



reverse.cm.distances = function(list.of.cm.dist){
    total.map.length = max(na.omit(list.of.cm.dist))
    list.of.cm.dist = abs(list.of.cm.dist - total.map.length)
    return(list.of.cm.dist)
}

reverse.genetic.map = function(data.frame.to.rev, cm.colname){
    #only works for one chromosomes at a time.
    #args:
    #data.frame.to.rev: data.frame in which the genetic map is kept
    #cm.col.name: character string indicating the name of the column in which cM distances are kept
    data.frame.to.rev = data.frame.to.rev[nrow(data.frame.to.rev):1, ]
    data.frame.to.rev[[cm.colname]] = reverse.cm.distances(data.frame.to.rev[[cm.colname]])
    return(data.frame.to.rev)
}

order.by.physical.dist.within.bins = function(gen.df, chr.colname, cm.colname, phys.pos.colname){
    #orders markers by physical position within bins
    #args:
    #gen.df: the dataframe containing the genetic map
    #chr.colname: name of column containing chromosome assignments
    #cm.colname: character vector containing the name of the column containing centimorgan distances
    #phys.pos.colname: character vector containing the name of the column containing physical positions in bp
    # browser()

    gen.df2 = split.df.into.list.of.dfs.by.column(gen.df, chr.colname, F)
    
    gen.df3 = lapply(gen.df2, function(x){
        # browser()
        unique.bins = unique(x[[cm.colname]])
        
        gen.df4 = split.df.into.list.of.dfs.by.column(x, cm.colname, F)
        q = sapply(gen.df4, function(y){
            sort(y[[phys.pos.colname]])[ceiling((nrow(y) / 2))]
            # browser()
        })
        
        # browser()
        diffs1 = diff(na.omit(q))
        len.negative = length(which(diffs1 < 0))
        len.positive = length(which(diffs1 > 0))
        
        #if the overall direction of markers on the chromosome is from high to low,
        #order markers within each bin so that they are descending.
        if(len.negative > len.positive){
            descend1 = T
        } else {
            descend1 = F
        }
        
        
        # sapply(gen.df4, function(z){
        #     "AX-95072696" %in% z$marker
        # })
        # which(gen.df$marker == "AX-95072696")
        # gen.df[3175:3177, ]
        # 
        # if("AX-95072696" %in% x$marker) browser()
        
        # browser()
        for(i in unique.bins){
            
            
            
            g = x[which(x[[cm.colname]] == i), ]    
            
            
            
            # if("AX-95072696" %in% g$marker) browser()
            
            sorted = sort(g[[phys.pos.colname]], na.last = F, decreasing = descend1, index.return = T)$ix
            g2 = g[sorted, ]
            
            
            
            # g2 = g[match(sorted, g[[phys.pos.colname]]), ]
            # g.na = g[which(is.na(g[[phys.pos.colname]])), ]
            # 
            # g2[which(is.na(g2[[phys.pos.colname]])), ] = g.na
            # g2 = reset.rownames(g2)
            # 
            #deal with markers that have the same physical position
            # dups1 = which(duplicated(g2$marker))
            # if(length(dups1) > 0){
            #     dups2 = g[which(g[[phys.pos.colname]] %in% g2[dups1, ][[phys.pos.colname]]), ]
            #     g2[dups1, ] = dups2[which(!dups2$marker %in% g2$marker), ]
            # }


            
            x[which(x[[cm.colname]] == i), ] = g2
        }
        
        
        
        if(x[[cm.colname]][1] > x[[cm.colname]][nrow(x)]) x = x[nrow(x):1, ]
        if(check.if.inverted2(x, cm.colname, phys.pos.colname) == T) x = reverse.genetic.map(x, cm.colname)
        x
    })
    
    
    bind_rows(gen.df3)
}

evaluate.hets = function(df1){
    #gets the % number of heterozygotes out of the entire genotyping dataset
    #	args:
    #	df1 - genotyping dataframe
    num.hets = sum(unlist(lapply(df1, function(x){
        length(which(x == "H"))
    })))
    
    (num.hets / (nrow(df1) * ncol(df1))) * 100
    
}

#     ____________________________________________________________________________
#     SEG. DIST. INVESTIGATION                                                                                                ####

grab.geno.vector.from.rqtl.format = function(markername, rqtl.df){
    #NB. Assumes the marker names are the column names of the data frame
    g = unname(unlist(rqtl.df[, which(colnames(rqtl.df) == markername)]))
    g[2:length(g)]
}

geno.count = function(genovector){
    #genovector - a character vector of genotypes "A", "B" and "H"
    A = length(which(genovector == "A"))
    B = length(which(genovector == "B"))
    H = length(which(genovector == "H"))
    all = c(A, H, B)
    names(all) = c("A", "H", "B")
    return(all)
}

geno.count.chi = function(genovector){
    #genovector - a character vector of genotypes "A", "B" and "H"
    A = length(which(genovector == "A"))
    B = length(which(genovector == "B"))
    H = length(which(genovector == "H"))
    all = c(A, B)
    c = chisq.test(all)
    res = c(c[[1]], c[[2]], c[[3]], A, B, H)
    names(res) = c("X-squared", "df", "p-value", "A", "B", "H")
    return(res)
}

geno.count.chi.f2 = function(genovector){
    #genovector - a character vector of genotypes "A", "B" and "H"
    total.genotypes = length(genovector)
    
    A = length(which(genovector == "A"))
    B = length(which(genovector == "B"))
    H = length(which(genovector == "H"))
    
    all.obs = c(A, B, H)
    
    
    exp.A = (total.genotypes / 4)
    exp.B = (total.genotypes / 4)
    exp.H = (total.genotypes / 2)
    
    all.exp = c(exp.A, exp.B, exp.H)
    
    sum = 0
    for(i in 1:3){
        sum = sum + (((all.obs[[i]] - all.exp[[i]])^2) / all.exp[[i]])
    }
    
    p.value = pchisq(sum, df = 2, lower.tail = F) #calculate p-value from chi-square value
    
    res = c(sum, 2, p.value, A, B, H)
    names(res) = c("X-squared", "df", "p-value", "A", "B", "H")
    return(res)
}

grab.seg.dist.markers = function(path.to.genotype.data, is.f2, p.threshold){
    if(missing(is.f2)) is.f2 = F
    if(missing(p.threshold)) p.threshold = 0.005
    
    if(is.f2) p.threshold = 0.0005 #my custom weighted chi-square function seems to produce much smaller p values than the inbuilt one, so a lower threshold is needed
    
    
    geno.data = read.csv(path.to.genotype.data, header = T, stringsAsFactors = F)
    
    pop.size = nrow(geno.data)
    
    if(is.f2 == F){
        g = lapply(geno.data[, 3:ncol(geno.data)], function(x){
            geno.count.chi(x[2:length(x)])
        })
    } else {
        g = lapply(geno.data[, 3:ncol(geno.data)], function(x){
            geno.count.chi.f2(x[2:length(x)])
        })
    }
    
    
    if(is.f2 == F){
        seg.dist.markers = g[which(unlist(lapply(g, function(x){
            x[[3]] < p.threshold
        })))]
    } else {
        seg.dist.markers = g[which(unlist(lapply(g, function(x){
            x[[3]] < p.threshold
        })))]
    }
    
    seg.dist.markers = t(as.data.frame(seg.dist.markers))
    
    if(nrow(seg.dist.markers) == 0){
        print("No seg. dist. markers found for this threshold")
        return()
    } else {
        if(is.f2 == F & length(which(seg.dist.markers[, 6] > 20)) != 0) seg.dist.markers = seg.dist.markers[-which(seg.dist.markers[, 6] > 20), ]
        
        #seg.dist.markers = names(seg.dist.markers)
        return(seg.dist.markers)
    }
    
}

#     ____________________________________________________________________________
#     GENERIC BLAST FILE PROCESSING FUNCTIONS                                                                 ####

split.blast.file = function(blastdf1, listn){
    #splits a BLAST dataframe (imported from tabular output format with read.blast()) into a list of dataframes. 
    #Ensures that no query hit is present in more than one than one list.
    #This is useful for parallel processing of the BLAST file
    #args:
    #blastdf1: a BLAST dataframe imported with read.blast()
    #listn: number of elements to split the dataframe into 
    
    row.num = nrow(blastdf1)
    g = round(row.num / listn)
    start.coords = seq(1, row.num, g)
    start.coords = start.coords[-length(start.coords)]
    end.coords = start.coords - 1
    end.coords = end.coords[-1]
    end.coords = c(end.coords, row.num)
    
    for(i in 1:length(start.coords)){
        if(i != length(start.coords)){
            while(blastdf1$qseqid[(end.coords[i] + 1)] == blastdf1$qseqid[end.coords[i]]){
                end.coords[i] = end.coords[i] + 1
                start.coords[(i + 1)] = start.coords[(i + 1)] + 1
            }    
        }
        
    }
    
    Map(function(s.coord, e.coord){
        blastdf1[s.coord:e.coord, ]
    }, start.coords, end.coords)
}

sort.blastdf = function(blastdf){
    #first clusters BLAST hits by chromosome, then sorts by the start location of each hit within clusters
    #blastdf - a dataframe containing BLAST output in tabular format
    sorted = newdf(colnames(blastdf), no.rows = T)
    
    #do some sorting
    for(i in unique(blastdf[, 2])){
        temp = blastdf[blastdf[, 2] == i, ]
        temp = temp[sort(as.numeric(temp[, 9]), index.return = T)$ix, ]
        sorted = rbind(sorted, temp)
    }
    return(sorted)
}

parse.scaffold.blast = function(blastdf1, dist.threshold, only.best){
    #finds all homologues (whether orthologue, paralogue or homeologue)
    #parses a BLAST dataframe of a short query sequence against a genome assembly
    #composed of scaffolds or chromosomes. If the assemblie is chromosomal, the parser will split 
    #the chromosome up into groups of hits where hits are more than dist.threshold bp apart.
    #returns a dataframe containing the best groups of hits (average bitscore higher than 200, individual hits no more than dist.threshold bp apart)
    #args:
    # blastdf1 - a BLAST dataframe imported using read.blast()
    # dist.threshold - Integer; the maximum number of bases between two hits for them to be considered part of the same group
    
    if(missing(only.best)) only.best = F

    blastdf1 = sort.blastdf(blastdf1)
    
    unique.groups = convert.to.character.data.frame(unique(blastdf1[, 1:2]))
    
    
    potential.homeologues = newdf(c("query", "scaffold", "start", "end", "length", "rev.comp", "avg.bitscore"), no.rows = T)
    
    count1 = make.counter()
    
    for(i in 1:nrow(unique.groups)){
        temp.df = filter(blastdf1, qseqid == unique.groups[i, 1], sseqid == unique.groups[i, 2])
        
        # browser()
        
        split.numeric.vector.by.difference = function(x, threshold){
            #takes a sorted numeric vector and splits into a list of vectors,
            #the number of elements in the list being the number of elements in x
            #that have a difference from adjacent elements that exceeds the threshold 
            #args:
            # x - a sorted numeric vector
            # threshold - an integer specifying the difference at which to split
            if(all(diff(x) < 0)){
                descending1 = T
                x = sort(x)
            } else {
                descending1 = F
            }
            
            # print("x:")
            # print(x)
            
            index.over.thres = which(diff(x) > threshold) + 1 #vector of coordinates of elements that need to be split
            
            #if there is more than one group of hits in the cluster, do some processing
            if(length(index.over.thres) > 0){
                #how far apart are the indices that need to be split?
                #i.e. how many elements will each split contain?
                distance.between.split.indices = c(diff(index.over.thres), 0) 
                
                new.list2 = list(x[1:(min(index.over.thres) - 1)])
                
                # print("index.over.thres")
                # print(index.over.thres)
                # print("distance.between.split.indices")
                # print(distance.between.split.indices)
                
                for(i in 1:length(index.over.thres)){
                    
                    # print("i")
                    # print(i)
                    # print("index.over.thres[i]")
                    # print(index.over.thres[i])
                    # print("distance.between.split.indices[i]")
                    # print(distance.between.split.indices[i])
                    
                    if(distance.between.split.indices[i] > 1){
                        new.list2 = c(new.list2, list(x[index.over.thres[i]:(index.over.thres[i] + distance.between.split.indices[i] - 1)]))
                    } else {
                        if(i == length(distance.between.split.indices)){ #if this is the last element in the group, add all remaining elements of input vector to the list
                            new.list2 = c(new.list2, list(x[index.over.thres[i]:length(x)]))
                        } else {
                            new.list2 = c(new.list2, list(x[index.over.thres[i]]))    
                        }
                        
                    }
                }
                if(descending1 == T){
                    new.list2 = new.list2[length(new.list2):1]
                    new.list2 = lapply(new.list2, sort, decreasing = T)
                } 
                
                new.list2
            } else {
                new.list2 = list(x)
            }
            
            
        }
        
        #determine whether there is more than one locus involved in this group of hits
        sstart.categories = split.numeric.vector.by.difference(temp.df$sstart, dist.threshold)
        
        #modify temp.df to have unique sseqids for each unique locus (loci more than 10000 bases apart)
        if(length(sstart.categories) > 1){
            temp.df$sseqid = as.character(temp.df$sseqid)
            
            for(x in 1:length(sstart.categories)){
                coords = which(temp.df$sstart %in% sstart.categories[[x]])
                temp.df$sseqid[coords] = paste(temp.df$sseqid[coords], ".!!$", x, sep = "")
            }    
        }
        
        # print("sstart.categories")
        # print(sstart.categories)
        
        
        #if two groups of hits are present in the same scaffold / chromosome,
        #these hits will no longer be ordered by bitscore due to sort.blastdf()
        #at the start of the script. Here we calculate mean bitscore for these newly identified
        #groups of hits, and sort groups of hits by this in descending order.
        temp.df.unique.scaffolds = unique(temp.df$sseqid)
        mean.bitscores = unlist(lapply(temp.df.unique.scaffolds, function(x){
            temp.df.filtered = filter(temp.df, sseqid == x)
            mean(as.numeric(temp.df.filtered$bitscore))
        }))
        
        mean.bitscores.sorted = sort(mean.bitscores, decreasing = T)
        
        transformation.coords = match(mean.bitscores.sorted, mean.bitscores)
        
        transformation.coords2 = unlist(lapply(temp.df.unique.scaffolds[transformation.coords], function(x){
            which(x == temp.df$sseqid)
        }))
        
        temp.df = temp.df[transformation.coords2, ]
        
        for(i2 in 1:length(unique(temp.df$sseqid))){
            temp.df2 = filter(temp.df, sseqid == unique(temp.df$sseqid)[i2])
            temp.df2$qseqid = as.character(temp.df2$qseqid)
            temp.df2$sseqid = as.character(temp.df2$sseqid)
            
            min.start = min(temp.df2$sstart)
            min.end = min(temp.df2$send)
            
            # print("temp.df2")
            # print(temp.df2)
            if(only.best == T) temp.df2 = temp.df2[which.max(temp.df2$bitscore), ]
            
            # browser()
            
            #populate potential.homeologues dataframe where average bitscore is higher than 200    
            if(mean(temp.df2$bitscore) > 200){
                if(min.start < min.end){
                    #if this is in normal orientation, do x...    
                    potential.homeologues = add_row(potential.homeologues)
                    potential.homeologues$query[nrow(potential.homeologues)] = temp.df2[1, 1]
                    potential.homeologues$scaffold[nrow(potential.homeologues)] = temp.df2[1, 2]
                    potential.homeologues$start[nrow(potential.homeologues)] = min(temp.df2$sstart)
                    potential.homeologues$end[nrow(potential.homeologues)] = max(temp.df2$send)
                    potential.homeologues$avg.bitscore[nrow(potential.homeologues)] = mean(temp.df2$bitscore)
                    potential.homeologues$rev.comp[nrow(potential.homeologues)] = F
                } else {
                    #else if this is a reverse complement sequence, do y...
                    potential.homeologues = add_row(potential.homeologues)
                    potential.homeologues$query[nrow(potential.homeologues)] = temp.df2[1, 1]
                    potential.homeologues$scaffold[nrow(potential.homeologues)] = temp.df2[1, 2]
                    potential.homeologues$start[nrow(potential.homeologues)] = min(temp.df2$send)
                    potential.homeologues$end[nrow(potential.homeologues)] = max(temp.df2$sstart)
                    potential.homeologues$avg.bitscore[nrow(potential.homeologues)] = mean(temp.df2$bitscore)
                    potential.homeologues$rev.comp[nrow(potential.homeologues)] = T
                }
            }
            
            # print("potential.homeologues")
            # print(potential.homeologues)
            
            # browser()
            
        }
        
    }
    
    potential.homeologues$length = as.numeric(potential.homeologues$end) - as.numeric(potential.homeologues$start)
    potential.homeologues$scaffold = multi.str.split(potential.homeologues$scaffold, "\\.\\!\\!\\$", 1)
    potential.homeologues    
}

extract.sequence = function(genome1, blast.df.parsed, row.coords, start.buffer, end.buffer){
    #extracts a related sequence from a genome assembly
    #args:
    # genome1 - DNAStringSet object containing the genome assembly of interest
    # blast.df.parsed - dataframe produced by parse.scaffold.blast()
    # row.coords - Numeric vector; the row coordinates of the sequences to used in blast.df.parsed
    # start.buffer - Integer; how much extra sequence before the start indicated in blast.df.parsed to extract
    # end.buffer - Integer; how much extra sequence after the end indicated in blast.df.parsed to extract
    
    #parse the original scaffold name from blastdf1.parsed (remove the appended position in kb)
    original.scaf.names = multi.str.split(blast.df.parsed$scaffold, ".$!", 1) 
    browser()
    genome2 = genome1[match(original.scaf.names[row.coords], names(genome1))]    
    
    for(i in 1:length(row.coords)){
        #check if start.buffer reaches before the start of the scaffold, and likewise if end.buffer extends bigger than the total length
        extract.start = (as.numeric(blast.df.parsed$start[i]) - start.buffer)
        if(extract.start < 1) extract.start = 1
        extract.end = (as.numeric(blast.df.parsed$end[i]) + end.buffer)
        if(extract.end > length(genome2[[i]])) extract.end = length(genome2[[i]])        
        
        genome2[[i]] = genome2[[i]][extract.start:extract.end]
        if(blast.df.parsed$rev.comp[i] == T) genome2[[i]] = reverseComplement(genome2[[i]])
    }
    
    #append the start coordinate (in kb) of the blast hit to the name of the sequence for identification later
    names(genome2) = paste(blast.df.parsed$scaffold[row.coords], ".$!", round((as.numeric(blast.df.parsed$start[row.coords]) / 1000)), sep = "")
    
    genome2    
}

extract.sequence_v2 = function(genome1, blast.df.parsed, row.coords, start.buffer, end.buffer){
    #extracts a related sequence from a genome assembly
    #args:
    # genome1 - DNAStringSet object containing the genome assembly of interest
    # blast.df.parsed - dataframe produced by parse.scaffold.blast()
    # row.coords - Numeric vector; the row coordinates of the sequences to used in blast.df.parsed
    # start.buffer - Integer; how much extra sequence before the start indicated in blast.df.parsed to extract
    # end.buffer - Integer; how much extra sequence after the end indicated in blast.df.parsed to extract
    
    #parse the original scaffold name from blastdf1.parsed (remove the appended position in kb)
    original.scaf.names = multi.str.split(blast.df.parsed$scaffold, ".$!", 1) 
    
    g1 = DNAStringSet()
    for(i in 1:length(row.coords)){
        #check if start.buffer reaches before the start of the scaffold, and likewise if end.buffer extends bigger than the total length
        chromo1 = blast.df.parsed$scaffold[[i]]
        
        extract.start = (as.numeric(blast.df.parsed$start[i]) - start.buffer)
        if(extract.start < 1) extract.start = 1
        extract.end = (as.numeric(blast.df.parsed$end[i]) + end.buffer)
        if(extract.end > length(genome1[[chromo1]])) extract.end = length(genome1[[chromo1]])        
        
        
        
        seq1 = genome1[[chromo1]][extract.start:extract.end]
        if(blast.df.parsed$rev.comp[i] == T) seq1 = reverseComplement(seq1)
        g1 = c(g1, DNAStringSet(seq1))
    }
    
    
    #append the start coordinate (in kb) of the blast hit to the name of the sequence for identification later
    names(g1) = paste0("tau_", blast.df.parsed$query)
    
    g1    
}

plot.blastdf = function(blastdf){
    #Plots a BLAST dataframe by chromosome. Use sort.blastdf() first.
    plot(match(blastdf$V2, unique(blastdf$V2)), blastdf$V9, xaxt = "n")
    axis(1, at=1:length(unique(blastdf$V2)), labels= unique(blastdf$V2))
}

filter.by.e.value = function(blastdf, threshold){
    blastdf$V11 = as.numeric(blastdf$V11)
    blastdf = blastdf[which(blastdf$V11 < threshold),]
    return(blastdf)
}

grab.best.hits = function(blastdf, only.unique){
    #returns the best hits, by e-value, for each query sequence. 
    #if two hits for the same query sequence share the same best e-value, 
    #they are both returned (unless only.unique is set to TRUE, in which case a random hit is returned).
    
    if(!("dplyr" %in% .packages())) library(dplyr)
    
    if(missing(only.unique) == T) only.unique = F
    
    all = newdf(colnames(blastdf))
    for(i in unique(blastdf[, 1])){
        g = filter(blastdf, blastdf[, 1] == i)
        g = g[which(g[, 11] == min(g[, 11])), ]
        if(only.unique == T){
            if(nrow(g) > 1) g = g[sample(1:nrow(g), 1), ]
        }
        all = rbind(all, g)
    }
    all = all[-1, ]
    return(all)
}

convert.to.character.data.frame = function(df){
    df = as.data.frame(df)
    g = df
    g[] = lapply(df, as.character)
    return(g)
}

grab.best.groups.of.hits = function(blastdf){
    #Edit: This function is for processing BLAST files in which you have more than one hit that contains the same query-subject pair. I.e. when the subject is very large, and BLAST has broken the query up into multiple hits. 
    #compare groups of hits to different subject sequences --- primarily used when subjects are scaffolds (e.g. EI Paragon Assembly)
    #calculate the average e-value of hits to a particular subject sequence, 
    #then compare these e-values between subjects (all for a single query / input sequence).
    #This function returns the original blastdf, but sorted and with some new columns added which will make it easier to select the best groups of hits (e.g. per.cov)
    
    #NOTE: the per.cov column produced by this function is not 100% accurate. In one very specific situation it will produce the wrong value (undershoots true coverage) - See "Problems with grab.best.group.of.hits()" in OneNote
    unique.query = unique(blastdf$qseqid)
    #the following is an index of the fasta file used as a query in the BLAST search, generated using "samtools faidx fasta.fa" on the remote server
    # genelengths = read.table("rotation1scripts_v4/processed_data/fasta/unique.seg.dist.genes.full.genonmic.sequence.v3.fa.fai")
    genelengths = read.table("bioinf/blast/genes.vs.paragon.genome/sacha.seg.dist.genes.full.genomic.sequence.2a.v2.fa.fai")
    
    
    all = newdf(c(colnames(blastdf)))
    blastdf$mean.e.value = "" #make new column for mean e-value
    #calculate the mean e-value for groups of hits within the same query-subject pair (V11 is e-value)
    for(i in unique(blastdf$qseqid)){
        for(p in unique(blastdf[blastdf$qseqid == i, ]$sseqid)){
            blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$mean.e.value = mean(blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$evalue)
        }
    }
    
    #grab lengths of genes (bp) and add them to blastdf
    blastdf$gene.length = ""
    for(i in 1:nrow(blastdf)){
        blastdf$gene.length[i] = genelengths[match(blastdf[i, 1], genelengths$V1), 2]
    }
    
    
    blastdf$coverage = ""
    blastdf$per.cov = ""
    
    for(i in unique(blastdf$qseqid)){
        for(p in unique(blastdf[blastdf$qseqid == i, ]$sseqid)){
            #make a list of the query start positions and query end positions
            g = list()
            g = c(list(blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$qstart), list(blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$qend))
            
            sort.indices = sort(g[[1]], index.return = T)[[2]]
            g[[1]] = g[[1]][sort.indices]
            g[[2]] = g[[2]][sort.indices]
            
            g[] = lapply(g, as.numeric)
            
            #calculate the amount of overlap (in bp) between hits
            overlap = 0
            for(count in as.numeric(1:(length(g[[1]])-1))){
                overlap = overlap + (g[[2]][count] - g[[1]][count + 1])
            }
            
            #calculate how much subject overlap there is between hits of the same query-subject pair
            subject.stats = blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ][, 9:10]
            subject.stats = subject.stats[sort(subject.stats$sstart, index.return = T)[[2]], ]
            
            
            
            
            overlap.subject = 0
            if(nrow(subject.stats) > 1){
                if(all(subject.stats[, 1] < subject.stats[, 2]) == T){
                    for(count in 1:(nrow(subject.stats)-1)){
                        overlap.subject = overlap.subject + (subject.stats[count+1, 1] - subject.stats[count, 2])
                    }
                } else if(all(subject.stats[, 1] > subject.stats[, 2]) == T){
                    for(count in 1:(nrow(subject.stats)-1)){
                        print(p(c(subject.stats[count+1, 2], "-", subject.stats[count, 1], "=", (subject.stats[count+1, 2] - subject.stats[count, 1]) )))
                        overlap.subject = overlap.subject + (subject.stats[count+1, 2] - subject.stats[count, 1])
                        print(p(c("overlap.subject = ", overlap.subject)))
                    }
                } else {
                    overlap.subject = "unknown"
                }
            }
            
            if(overlap.subject != "unknown") overlap.subject = -overlap.subject
            
            # if(p == "Triticum_aestivum_Paragon_EIv1.1_scaffold_035529") browser()
            
            
            #if there is overlap, remove from coverage value
            if(length(overlap) == 0){
                coverage = sum(g[[2]] - g[[1]])
            } else if(overlap > 0){
                    coverage = sum(g[[2]] - g[[1]]) - overlap
                } else {
                    coverage = sum(g[[2]] - g[[1]])
                }
            
            #if the subject hits are more than 500 bp apart (or different query-subject pairs are in different orientations; see "unknown"), set the coverage to that of the first hit only
            if(overlap.subject < -500 | overlap.subject == "unknown"){
                coverage = (g[[2]] - g[[1]])[[1]]
            }
            
            #add coverage and percentage coverage columns for groups of hits (i.e. each unique query-subject pair)
            blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$coverage = coverage
            blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$per.cov = (coverage/as.numeric(blastdf[blastdf$qseqid == i & blastdf$sseqid == p, ]$gene.length[1])*100)
            
            
        }
    }
    
    #sort blastdf by query start (V7) within each group of query-subject pair hits
    unique.q.s.pairs = unique(blastdf[, 1:2])
    for(i in 1:nrow(unique.q.s.pairs)){
        blastdf[blastdf[, 1] == unique.q.s.pairs[i, 1] & blastdf[, 2] == unique.q.s.pairs[i, 2], ] = blastdf[blastdf[, 1] == unique.q.s.pairs[i, 1] & blastdf[, 2] == unique.q.s.pairs[i, 2], ][sort(blastdf[blastdf[, 1] == unique.q.s.pairs[i, 1] & blastdf[, 2] == unique.q.s.pairs[i, 2], ]$qstart, index.return = T)[[2]], ]
    }
    
    return(blastdf)
    
}

clean.blast = function(x, threshold.diff.best.hits){
    #grab only the query-subject pairs which have a definitive best hit
    #x - blast dataframe in tabular format (-outfmt 6)
    #threshold.diff.best.hits - integer value; e-value difference between first and second best hit at which first hit will be kept, anything smaller and it will be discarded.
    if(missing(threshold.diff.best.hits)) threshold.diff.best.hits = 5.9e-26 #this value has been empirically determined for probes of around 80 bp
    diff.between.best.hits = 0
    single.hit = F
    
    
    newblast = newdf(colnames(x))
    for(i in unique(x$V1)){
        #if(i == "AX.94847111") browser()
        g = sort(as.numeric(x[x$V1 == i, ]$V11))
        if(length(g) > 1){
            diff.between.best.hits = abs(g[[1]] - g[[2]])
        } else {
            single.hit = T
        }
        
        row.coords = which(as.numeric(x[x$V1 == i, ]$V11) == min(as.numeric(x[x$V1 == i, ]$V11)))
        
        if(diff.between.best.hits > threshold.diff.best.hits | single.hit == T){
            newblast = rbind(newblast, x[x$V1 == i, ][row.coords, ])
        }
        
        diff.between.best.hits = 0
        single.hit = F
    }
    newblast = newblast[-1, ]
    return(newblast)
}

get.list.of.differences.between.best.and.second.best.e.values = function(blastdf){
    #returns a named list containing the differences between the first and second best e-values for each group of hits containing the same query
    differences = lapply(unique(blastdf$V1), function(x){
        diff = 0
        g = blastdf[blastdf$V1 == x, ]
        e.vals = sort(as.numeric(g$V11), decreasing = T)
        if(length(e.vals) > 1) diff = e.vals[[1]] - e.vals[[2]]
        names(diff) = x
        return(diff)
    })
    differences = differences[!differences == 0] #remove subject-query pairs where there is only one hit
    sort(unlist(differences))
}

remove.hits.with.same.best.e.value = function(blastdf){
    #removes query-subject pairs in which the two best hits for a single query both have the same e-value
    #blastdf - a dataframe containing a tabular blast output file
    
    hits.to.keep = newdf(colnames(blastdf), no.rows = T)
    blastdf[, 11] = as.numeric(blastdf[, 11])
    counter1 = make.counter()
    un.hits = unique(blastdf[, 1])
    blastdf = as.data.table(blastdf)
    setkey(blastdf, qseqid)
    lapply(un.hits, function(x){
        # temp.df = blastdf[blastdf[, 1] == x, ]
        # temp.df = filter(blastdf, qseqid == x)
        temp.df = blastdf[.(x), nomatch = 0L]
        best.e.value = min(temp.df[, 11])
        if(nrow(filter(temp.df, evalue == best.e.value) == 1)) hits.to.keep <<- rbind(hits.to.keep, temp.df)
        # if(length(which(temp.df[, 11] == best.e.value)) == 1) hits.to.keep <<- rbind(hits.to.keep, temp.df)
        print(counter1())
    })
    
    
    # for(i in unique(blastdf[, 1])){
    #     temp.df = blastdf[blastdf[, 1] == i, ]
    #     best.e.value = min(as.numeric(temp.df[, 11]))
    #     if(length(which(temp.df[, 11] == best.e.value)) == 1) hits.to.keep = rbind(hits.to.keep, temp.df)
    # }
    return(hits.to.keep)
}

gen.blast.script = function(input.fasta, blastdb, culling_limit, word_size){
    #generates a .sh file for performing a BLAST search
    #input.fasta: character string specifying the filename of fasta file to be used in BLAST 
    #blastdb: character string specifying the BLAST database to be used (three options)
    #culling_limit: integer
    #word_size: integer
    if(missing(culling_limit)) culling_limit = 2
    if(missing(word_size)) {
        word_size = ""
    } else {
        word_size = p("-word_size ", word_size, " ")
    }
    if(blastdb != "genes.vs.paragon.genome" & blastdb != "probe.vs.genes.blast" & blastdb != "probe.vs.genome.blast"){
        print("blastdb should be either of the following: genes.vs.paragon.genome, probe.vs.genes.blast, probe.vs.genome.blast")
    } else {
        if(blastdb == "genes.vs.paragon.genome") blastdb.path = "./blastdb/earlham.paragon.db"
        if(blastdb == "probe.vs.genes.blast") blastdb.path = "./blastdb/genedb"
        if(blastdb == "probe.vs.genome.blast") blastdb.path = "wheat_genome_blast_db"
        
        script = c("#!/bin/bash", p("blastn -db ", blastdb.path, " -query ./", input.fasta, " -outfmt 6 -out results.blast/", input.fasta, ".blast -num_threads 40 ", word_size, "-culling_limit ", culling_limit),
                                p("blastn -db ", blastdb.path, " -query ./", input.fasta, " -outfmt 5 -out results.blast/", input.fasta, ".xml.blast -num_threads 40 ", word_size, "-culling_limit ", culling_limit))
     
        outputfile = file(p("bioinf/blast/", blastdb, "/", input.fasta, ".sh"), "wb") #open a binary connection to write UNIX line ending on Windows
        writeLines(script, outputfile)
        close(outputfile)
        
    }
}

add.blast.colnames = function(blastdf){
    colnames(blastdf)[1:12] = c("qseqid", "sseqid", "percentage.identical", "alignment.length", "no.mismatch",
                                                            "no.gap.openings", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    return(blastdf)
}

reset.colnames = function(df, which.columns){
    #which.colnames: numeric vector, if only a subset of the columns are to be reset, specify which ones
    if(missing(which.columns)) which.columns = 1:ncol(df)
    g = which.columns
    g = paste("V", g, sep = "")
    colnames(df)[which.columns] = g
    return(df)
}


read.blast = function(filepath){
    g = read.table(filepath)
    g = add.blast.colnames(g)
    g[, 3:12] = lapply(g[, 3:12], as.numeric) #convert number columns to numeric
    return(g)
}
            

#     ____________________________________________________________________________
#     OTHER FUNCTIONS                                                                                                                 ####

reset.rownames = function(x){
    rownames(x) = 1:nrow(x)
    return(x)
}

find.dupes = function(x){
    #returns the elements of a vector that are present more than once
    #x - a vector
    unique.elements = unique(x)
    dupe.elements = ""
    for(i in unique.elements){
        if(length(which(x == i)) > 1) dupe.elements = c(dupe.elements, i)
    }
    dupe.elements = dupe.elements[-1]
    return(dupe.elements)
}

#reverse a string
strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, 
                                                                 collapse="") 

#write blast dataframe to tabular output
write.blast = function(blastdf, file.path){
    write.table(blastdf, file.path, sep = "\t", quote = F, row.names = F, col.names = F)
}

#return number of unique values in a vector
nunique = function(x){
    length(unique(x))
}


multi.str.split = function(x, character.split, split.id){
    #An implementation of strsplit() that allows selection of a particular element of the output value of strsplit(), 
    #for all elements of the input vector to which the split is applied.
    #--------------
    #args:
    #x = character vector
    #character.split = character specifying what character to split with 
    #split.id = integer specifying which element of the split you want
    
    unlist(lapply(x, function(q) strsplit(q, character.split)[[1]][split.id]))
}

setup.mysql.connection = function(desktop){
    if(!("RMySQL" %in% .packages())) library(RMySQL)    
    if(missing(desktop)) desktop = F
    
    if(desktop == F){
        #wilkins connection
        mydb = dbConnect(MySQL(), user = 'ac14037', password = 'mnb56ghj', db = 'alexcereals', host = '127.0.0.1')
    } else {
        # desktop computer connection
        mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')
            
    }
    return(mydb)
    
}

replace.aas.genotype.format.w.rqtl.format = function(dataframe){
    dataframe[dataframe == "AA"] = "A"
    dataframe[dataframe == "BB"] = "B"
    dataframe[dataframe == "AB"] = "H"
    dataframe[dataframe == "BA"] = "H"
    dataframe[dataframe == "NoCall"] = "-"
    return(dataframe)
}

s <- function(object, filepath, original.script.name){
    #see stackexchange post for details of how the code works: https://stackoverflow.com/questions/49557017/saving-object-attributes-using-a-function
    #modified save() function
    #stores the name of the script from which the object originates as an attribute, then saves as normal
    objectname <- deparse(substitute(object))
    attr(object, "original.script") = original.script.name
    save_envir <- new.env()
    save_envir[[objectname]] <- object
    save(list=objectname, file = filepath, envir=save_envir)
}

switch.affy.format = function(markernames){
    #changes vector of marker names (e.g. "AX.93493493") to other format (either "-" or ".")
    if(length(grep("\\.", markernames)) != 0){
        markernames = gsub("\\.", "-", markernames)
    } else {
        markernames = gsub("-", ".", markernames)
    }
    return(markernames)
}


make.counter = function(){
    counter = 0
    ret.count = function(){
        counter <<- counter + 1
        counter
    }
}



make.caption = function(prefix){
    caption.display = as.character()

    ret.count = function(caption.name, display, show.all, return.blank){
		if(missing(show.all)) show.all = F
		if(missing(display)) display = ''
		if(missing(caption.name)) caption.name = "NA"
		if(missing(return.blank)) return.blank = F

		g = unlist(lapply(caption.name, function(caption.name){
			if(caption.name %in% caption.display){
				if(display != 'num') return(paste0(prefix, '', which(caption.display == caption.name)))
				if(display == 'num') return(paste0(
                                    strsplit(prefix, '')[[1]][nchar(prefix) - 1],
                                    ".",
                                    which(caption.display == caption.name)))
			} else {
				caption.display <<- c(caption.display, caption.name)
				if(display != 'num') return(paste0(prefix, '', which(caption.display == caption.name)))
				if(display == 'num') return(paste0(
                                    strsplit(prefix, '')[[1]][nchar(prefix) - 1],
                                    ".",
                                    which(caption.display == caption.name)))
			}
        }))
	    if(return.blank == T) return() 
		if(show.all == T) return(caption.display)
	    return(g)	

    }
}

library(compiler) #function for finding longest increasing subsequence
#from https://www.r-bloggers.com/compute-longest-increasingdecreasing-subsequence-using-rcpp/
longest_subseq.R = cmpfun(function(x) {
    P = integer(length(x))
    M = integer(length(x) + 1)
    L = newL = 0
    for (i in seq_along(x) - 1) {
        lo = 1
        hi = L
        while (lo <= hi) {
            mid = (lo + hi)%/%2
            if (x[M[mid + 1] + 1] < x[i + 1]) {
                lo = mid + 1
            } else {
                hi = mid - 1
            }
        }
        newL = lo
        P[i + 1] = M[newL]
        if (newL > L) {
            M[newL + 1] = i
            L = newL
        } else if (x[i + 1] < x[M[newL + 1] + 1]) {
            M[newL + 1] = i
        }
    }
    k = M[L + 1]
    re = integer(L)
    for (i in L:1) {
        re[i] = k + 1
        k = P[k + 1]
    }
    re
})

check.if.coordinate.is.gene = function(coordinate, chromosome){
    #Checks whether a particular coordinate is located within a gene described by the IWGSC high confidence genes gff file.
    #args:
    # coordinate - an integer denoting the coordinate to check
    # chromosome - a character vector denoting the chromosome to check, e.g. "5B"
    
    if(!"IRanges" %in% (.packages())) library(IRanges)
    if(!"dplyr" %in% (.packages())) library(dplyr)
    if(!exists("genes")){
        genes = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", header = T, stringsAsFactors = F)
        genes$gene.names = multi.str.split(as.character(genes$V10), "=", 2)
        genes$gene.names = multi.str.split(genes$gene.names, ";", 1)
        genes = genes[, -1]
        genes = reset.colnames(genes)
        colnames(genes) = c(colnames(genes)[1:9], "gene.names")
        genes$V4 = as.numeric(genes$V4)
    }
    
    genes.filtered = filter(genes, V1 == p("chr", chromosome))
    
    # browser()
    
    gene.ranges = IRanges(genes.filtered$V4, genes.filtered$V5, names = genes.filtered$gene.names)
    
    coord.to.check = IRanges(coordinate, (coordinate + 1))
    
    if(mcols(distanceToNearest(coord.to.check, gene.ranges))$distance == 0){
        gene.ranges[nearest(coord.to.check, gene.ranges)]
    } else {
        "Not in a gene"
    }
    
}

make.ggplot.histogram = function(histogram.vector, num.bins, breaks, xlabel, ylabel, plot.title){
    #converts a numeric vector of values into a ggplot histogram
    #args:
    # histogram.vector - numeric vector containing numbers for histogram
    # num.bins - integer specifying the number of bins to include in the histogram
    # breaks - numeric vector specifying bin locations. overrides num.bins
    if(missing(num.bins)) num.bins = 100
    if(missing(xlabel)) xlabel = "X"
    if(missing(plot.title)) plot.title = F
    temphistdata = data.frame(histogram.vector, 1:length(histogram.vector))
    colnames(temphistdata) = c("c1", "c2")
    
    if(plot.title == F){
        ggplot(temphistdata, aes(temphistdata$c1)) + geom_histogram(bins = num.bins, breaks = breaks) + xlab(xlabel) + ylab(ylabel)
    } else {
        ggplot(temphistdata, aes(temphistdata$c1)) + geom_histogram(bins = num.bins, breaks = breaks) + xlab(xlabel) + ylab(ylabel) + ggtitle(plot.title)
    }
}

rqtl.read = function(rqtl.csv.path, fgen, crosstype){
    #easier version of read.cross
    #args:
    # rqtl.csv.path - string, path to genotype file in rQTL format
    if(missing(crosstype)) crosstype = "NONE"
    if(missing(fgen)) fgen = 2
    g = strsplit(rqtl.csv.path, "/")
    filename = g[[1]][length(g[[1]])]
    path1 = paste(g[[1]][-length(g[[1]])], collapse = "/")
    path1 = paste0(path1, "/")
    if(crosstype == "NONE"){
        cross.return = read.cross("csv", path1, filename, genotypes = c("A", "H", "B"), estimate.map = F, F.gen = fgen)    
    } else {
        if(crosstype == "dh"){
            cross.return = read.cross("csv", path1, filename, genotypes = c("A", "B"), estimate.map = F, crosstype = crosstype)    
        } else {
            cross.return = read.cross("csv", path1, filename, genotypes = c("A", "B"), estimate.map = F, crosstype = crosstype)        
        }
    }
    cross.return
}


transpose.list.of.lists = function(x){
    lapply(1:length(x[[1]]), function(y){
        lapply(x, function(z){
            z[[y]]
        })
    })
}


vector.captioner = function(cap.function, figure.names){
	#vectorized version of the captioner function for generating figure numbers 
	#in Rmarkdown
	#Args:
	#    cap.function: the captioner function e.g. fig(); tab()
	#    figure.names: vector of figure / table names
	lap.func = function(x){
		cap.function(x)
	}

	g <<- unlist(lapply(figure.names, lap.func))
}

