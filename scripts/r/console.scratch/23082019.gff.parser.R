getwd()
setwd("~/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")





gff = read.table("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3", sep = "\t", stringsAsFactors = F)

library(stringr)

g1 = str_extract(gff$V9, "Traes[A-Z0-9]*")
gff$names = g1

library(readxl)
library(Biostrings)
alab1 = readDNAStringSet("rotation1scripts_v4/original_data/recombination.gene.analysis/alabdullah.genes1.fa")
alab1 = multi.str.split(names(alab1), "_", 1)
alab1 = list.files("~/pipeline/jobs/alab1/all.alignments/", pattern = ".fa")
alab1 = multi.str.split(alab1, "\\.", 1)


random1 = multi.str.split(list.files("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/", pattern = ".fa"), ".fa", 1)

#add.cds("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/", random1)

genes1 = read.csv("rotation1scripts_v4/original_data/recombination.gene.analysis/meiosis.genes.csv", stringsAsFactors = F)

colnames(genes1) = "gene"

#### SET VARIABLES ####
alignments.dir = "rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/random.genes/all.alignments/"
gene.names = random1




# add.cds = function(alignments.dir, gene.names){
    #args:
    # alignments.dir: character string containing path to dir containing alignments, without trailing slash
    #gene.names: character vector of gene names
    # browser()
    coords = sapply(gene.names, function(x){
        which(gff$names == x)
    })
    
    gff2 = lapply(coords, function(x){
        gff[x, ]
    })
    
    rm(gff)
    
    #are there any exons which start after the start of the gene? (no)
    lapply(gff2, function(x){
        g = x[which(x$V3 == "gene"), ]
        g2 = x[which(x$V3 == "exon"), ]
        # browser()
        min(g2$V4) == min(g$V4)
    })
    
    gff3 = lapply(gff2, function(x){
        start1 = x$V4[1]
        x$V4 = x$V4 - start1
        x$V5 = x$V5 - start1
        x
    })
    
    #reverse exon orientation for antisense genes
    gff4 = lapply(gff3, function(x){
        x = x[which(x$V3 == "exon"), ]
        # browser()
        if(x$V7[1] == "-"){
            stop1 = x$V5[nrow(x)]
            x$V4 = abs(x$V4 - stop1)
            x$V5 = abs(x$V5 - stop1)
            end1 = x$V4
            x$V4 = x$V5
            x$V5 = end1
            x = x[sort(x$V4, index.return = T)$ix, ]
        }
        x
    })
    
    sapply(gff4, function(x){
        length(unique(x$V9))
    })
    
    library(parallel)
    gff5 = mclapply(gff4, function(x){
        #split.df.into.list.of.dfs.by.column(x, "V9", F)
        split(x, x$V9)
    }, mc.cores = 50)
    
    gff6 = lapply(gff5, function(x) x[[1]])
    print("done gff formatting")
    align1 = list.files(alignments.dir, pattern = ".fa")
    align2 = list.files(alignments.dir, pattern = ".fa", full.names = T)
    
    align3 = multi.str.split(align1, ".fa", 1)
    
    library(Biostrings)
    align4 = lapply(align2, readDNAStringSet)
    names(align4) = align3
    
    
    gff6 = gff6[match(names(align4), names(gff6))]
    gff6 = lapply(gff6, function(x){
        x$V4 = x$V4 + 1
        x$V5 = x$V5 + 1
        x
    })
    
    #function to get consecutive element ranges
    cons1 = function(x){
        diffs = c(1, diff(x))
        start_indexes = c(1, which(diffs > 1))
        end_indexes = c(start_indexes - 1, length(x))
        # coloned = paste(x[start_indexes], x[end_indexes], sep = ":")
        # coloned
        # paste0(coloned, collapse = ", ")
        g = data.frame(x[start_indexes], x[end_indexes])
        colnames(g) = c("start", "end")
        g
    }
    print("performing hyphen insertion")
    # browser()
    
    
    hyphen.function1 = function(seq1, coord1){
        
        
        
        seq2 = strsplit(as.character(seq1[[1]]), "")
        seq2 = seq2[[1]]    
        seq3 = seq2[-which(seq2 == "-")]
        seq4 = as.data.frame(seq3)
        seq4 = convert.to.character.data.frame(seq4)
        seq4$exon = ""
        for(i in 1:nrow(coord1)){
            start.coord1 = coord1$V4[i]
            end.coord1 = coord1$V5[i]
            seq4$exon[start.coord1:end.coord1] = seq4$seq3[start.coord1:end.coord1]
        }
        # print("x1")
        seq4$exon[which(seq4$exon == "")] = "-"
        
        ranges1 = cons1(which(seq2 == "-"))
        ranges1 = as.data.frame(t(ranges1))
        test2 = seq2[which(seq2 != "-")]
        
        #generate list of dataframes containined this coordinate position and length of insert
        coord.pos.df = lapply(ranges1, function(x){
            len.insert = (x[2] - x[1]) + 1
            pos.insert = x[1]
            data.frame(len.insert, pos.insert)
        })
        
        main.seq = test2
        
        insert.hyphens = function(coord.pos.df, main.seq){
            #perform insertion of hyphens are coordinate positions, special cases for start and end insertions
            lapply(coord.pos.df, function(x){
                insert1 = rep("-", x$len.insert)
                if(x$pos.insert == 1){
                    main.seq <<- c(insert1, main.seq)
                } else if(x$pos.insert > length(main.seq)){
                    main.seq <<- c(main.seq, insert1)
                } else {
                    main.seq <<- c(main.seq[1:(x$pos.insert - 1)], insert1, main.seq[x$pos.insert:length(main.seq)])
                }
                T
            })
            g2 = paste(main.seq, collapse = "")
            g2
        }
        
        
        exon.seq = insert.hyphens(coord.pos.df, seq4$exon)
        main.seq2 = insert.hyphens(coord.pos.df, seq4$seq3)
        
        all.seq = c(DNAStringSet(DNAString(exon.seq)), seq1)
        names(all.seq)[1] = "cds"    
        names(all.seq)[2] = "genomic"
        all.seq
    }
    
    hyphen.function1(align4[[1]], gff6[[1]])
    
    align.w.cds = mclapply(1:length(align4), function(x){
        tryCatch(hyphen.function1(align4[[x]], gff6[[x]]), error = function(e) NULL)
    }, mc.cores = 50)
    names(align.w.cds) = names(align4)
    
    align.w.cds = align.w.cds[which(unlist(lapply(align.w.cds, class)) == "DNAStringSet")]    
    
    
    
    
    count1 = 1
    dir.create(paste0(alignments.dir, "/wcds"))
    lapply(align.w.cds, function(x){
        writeXStringSet(x, paste0(alignments.dir, "wcds/", names(align.w.cds)[count1], ".fa"), width = 20000)
        count1 <<- count1 + 1
    })
    

    
    align.cds.no.gaps = mclapply(align.w.cds, function(x){
        g = strsplit(as.character(x[[1]]), "")
        g = g[[1]]
        
        x2 = as.matrix(x)
        x2 = x2[, which(g != "-")]
        x3 = as.data.frame(t(x2))
        x3 = convert.to.character.data.frame(x3)
        x4 = lapply(x3, function(x) DNAStringSet(DNAString(paste(x, collapse = ""))))
        x5 = c(x4[[1]], x4[[4]], x4[[5]], x4[[6]],
            x4[[7]], x4[[8]])
        names(x5) = names(x)[c(1, 4, 5, 6, 7, 8)]
        x5
    }, mc.cores = 50)
    
    count1 = 1
    dir.create(paste0(alignments.dir, "cds_nogaps"))
    lapply(align.cds.no.gaps, function(x){
        writeXStringSet(x, paste0(alignments.dir, "cds_nogaps/", names(align.cds.no.gaps)[count1], ".fa"), width = 20000)
        count1 <<- count1 + 1
    })
    
    align.cds.no.gaps[[1]]
    
    align.cds.no.gaps.random = lapply(list.files("~/pipeline/jobs/random.genes/all.alignments/cds_nogaps/", pattern = ".fa$", full.names = T), readDNAStringSet)
    names(align.cds.no.gaps.random) = list.files("~/pipeline/jobs/random.genes/all.alignments/cds_nogaps/", pattern = ".fa$")
    
    #### SUPERALIGNMENT ANALYSIS ####
    produce.superalignment = function(list.of.dnastringsets){
        
        superalignment1 = mclapply(list.of.dnastringsets, function(x){
            as.matrix(x)
        }, mc.cores = 50)
        
        superalignment1 = do.call(cbind, superalignment1)
        superalignment1 = as.data.frame(t(superalignment1))
        superalignment1 = convert.to.character.data.frame(superalignment1)    
        superalignment2 = lapply(superalignment1, function(x) DNAStringSet(DNAString(paste(x, collapse = ""))))         
        superalignment3 = superalignment2[[1]]
        for(i in 2:length(superalignment2)){
            superalignment3 = c(superalignment3, superalignment2[[i]])
        }
        
        names(superalignment3) = names(list.of.dnastringsets[[1]])
        superalignment3
    }
    
    
    
    
    writeXStringSet(superalignment1, paste0(alignments.dir, "cds_nogaps/superalignment.fa"))    
#### REMOVE BARLEY ANALYSIS ####
    

align.cds.no.gaps.no.barley = lapply(align.cds.no.gaps, function(x){
    x[-6]    
})    

    dir.create("~/pipeline/jobs/alab1/all.alignments/cds_nogaps/nobarley")
    count1 = 1
    lapply(align.cds.no.gaps.no.barley, function(x){
        writeXStringSet(x, paste0("~/pipeline/jobs/alab1/all.alignments/cds_nogaps/nobarley/", names(align.cds.no.gaps.no.barley)[count1], ".fa"))
        count1 <<- count1 + 1
    })
        
align.cds.no.gaps.random.no.barley = lapply(align.cds.no.gaps.random, function(x){
    x[-6]
})    

ran.nobar.super = produce.superalignment(align.cds.no.gaps.random.no.barley)



dir.create("~/pipeline/jobs/random.genes/all.alignments/cds_nogaps/nobarley")
count1 = 1
lapply(align.cds.no.gaps.random.no.barley, function(x){
    writeXStringSet(x, paste0("~/pipeline/jobs/random.genes/all.alignments/cds_nogaps/nobarley/", names(align.cds.no.gaps.random.no.barley)[count1]))
    count1 <<- count1 + 1
})


