#BLAST PARSER FOR SCAFFOLDS
library(tibble)
library(dplyr)
library(Biostrings)
setwd("/home/ac14037/project.phd.main/rotation1scripts_v4/auto.primer.picker")
source("/home/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")
setwd("/home/ac14037/project.phd.main/rotation1scripts_v4/auto.primer.picker")

args = commandArgs(trailingOnly = T)
gene.name = args[1]

if(length(args) == 0){
    print("Arguments: gene.name")
} else {

#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####

parse.scaffold.blast = function(blastdf1, dist.threshold){
    #parses a BLAST dataframe of a short query sequence against a genome assembly
    #composed of scaffolds or chromosomes. If the assemblie is chromosomal, the parser will split 
    #the chromosome up into groups of hits where hits are more than dist.threshold bp apart.
    #returns a dataframe containing the best groups of hits (average bitscore higher than 200, individual hits no more than dist.threshold bp apart)
    #args:
    # blastdf1 - a BLAST dataframe imported using read.blast()
    # dist.threshold - Integer; the maximum number of bases between two hits for them to be considered part of the same group
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
        
        transformation.coords = match(mean.bitscores, mean.bitscores.sorted)
        
        temp.df = temp.df[match(temp.df.unique.scaffolds[transformation.coords], temp.df$sseqid), ]
        
        for(i2 in 1:length(unique(temp.df$sseqid))){
            temp.df2 = filter(temp.df, sseqid == unique(temp.df$sseqid)[i2])
            temp.df2$qseqid = as.character(temp.df2$qseqid)
            temp.df2$sseqid = as.character(temp.df2$sseqid)
            
            min.start = min(temp.df2$sstart)
            min.end = min(temp.df2$send)
            
            # print("temp.df2")
            # print(temp.df2)
            
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
    
    
    genome2 = genome1[match(blast.df.parsed$scaffold[row.coords], names(genome1))]
    
    
    
    for(i in 1:length(row.coords)){
        genome2[[i]] = genome2[[i]][(as.numeric(blast.df.parsed$start[i]) - start.buffer):(as.numeric(blast.df.parsed$end[i]) + end.buffer)]
        if(blast.df.parsed$rev.comp[i] == T) genome2[[i]] = reverseComplement(genome2[[i]])
    }
    
    #append the start coordinate (in kb) of the blast hit to the name of the sequence for identification later
    names(genome2) = paste(blast.df.parsed$scaffold[row.coords], ".", round((as.numeric(blast.df.parsed$start[row.coords]) / 1000)), sep = "")
    
    genome2
    
}

#     ____________________________________________________________________________
#     BEGIN PROCESSING                                                                                                                ####

#read the configuration file
config.file = readLines("./config.txt")
config.variables = multi.str.split(config.file, "=", 1)

#parse number of genomes in configuration file
number.genomes = max(na.omit(unique(as.numeric(multi.str.split(config.variables, "_", 2)))))

lapply(1:number.genomes, function(x){
    #begin by parsing the config file for the genome name, fasta file path and blastdb path
    config.settings = config.file[grep(x, config.variables)]
    
    config.settings.temp = config.settings[[1]]
    genome.name = strsplit(config.settings.temp[1], "=")
    genome.name = genome.name[[1]][2]
    
    config.settings.temp = config.settings[[2]]
    fa.path = strsplit(config.settings.temp[1], "=")
    fa.path = fa.path[[1]][2]
    
    config.settings.temp = config.settings[[3]]
    blastdb.path = strsplit(config.settings.temp[1], "=")
    blastdb.path = blastdb.path[[1]][2]
    
    blast.files = list.files("blast.results")
    
    blastdf1 = read.blast(p("blast.results/", blast.files[x]))
    
    #THIS NEEDS TO BE OPTIMISED WITH fa.index()
    genome.assembly = readDNAStringSet(fa.path)
    
    #for main genome (first in config.txt), take four best sequences (homologues included)
    if(x == 1){
        blastdf1.parsed = parse.scaffold.blast(blastdf1, 2000)[1:4, ]
        sequences = extract.sequence(genome.assembly, blastdf1.parsed, 1:4, 2000, 2000)
    } else {
        #for additional varietal genomes (after the first in config.txt), take only the first sequence,
        #the allele from the same locus as the query
        blastdf1.parsed = parse.scaffold.blast(blastdf1, 2000)[1, ]    
        sequences = extract.sequence(genome.assembly, blastdf1.parsed, 1, 2000, 2000)
    }
    
    writeXStringSet(sequences, p(gene.name, "seq/extended/seqs/", genome.name, ".seqs.fa"))
    
})


}