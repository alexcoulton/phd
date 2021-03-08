source("/home/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")
setwd("/home/ac14037/project.phd.main/bioinf/blast/probe.vs.tauschii.blast/results.blast/")
library(dplyr)


blast1 = read.blast("chr1D.vs.tauschii.blast")

process.blast1 = function(x){
    blast1.2 = parse.scaffold.blast(x, 5000)
    
    blast1.3 = split(blast1.2, blast1.2$query)
    blast1.4 = lapply(blast1.3, function(x){
        if(nrow(x) > 1){
            x2 = x[which.max(x$avg.bitscore), ]
            # if(x2$scaffold != 7) browser()
        } else {
            x2 = x
        }
        x2
    })
    
    blast1.5 = bind_rows(blast1.4)
    g2 = multi.str.split(blast1.5$query, "_", 2)
    g3 = as.numeric(multi.str.split(g2, "-", 1))
    g4 = as.numeric(multi.str.split(g2, "-", 2))
    blast1.5$cs.pos.start = g3
    blast1.5$cs.pos.end = g4
    blast1.5$cs.length = blast1.5$cs.pos.end - blast1.5$cs.pos.start
    blast1.5
}


blast1.2 = process.blast1(blast1)
plot(blast1.2$start, blast1.2$cs.pos.start)
nrow(blast1.2)
nunique(blast1.2$query)
nunique(blast1$qseqid)


blast2 = read.blast("chr2D.vs.tauschii.blast")

blast2.2 = process.blast1(blast2)
# plot(blast2.2$start, blast2.2$cs.pos.start)

blast3 = read.blast("chr3D.vs.tauschii.blast")
blast3.2 = process.blast1(blast3)
# plot(blast3.2$start, blast3.2$cs.pos.start)

unique(blast3.2$scaffold)

blast4 = read.blast("chr4D.vs.tauschii.blast")
blast4.2 = process.blast1(blast4)
# plot(blast4.2$start, blast4.2$cs.pos.start)

blast5 = read.blast("chr5D.vs.tauschii.blast")
blast5.2 = process.blast1(blast5)
# plot(blast5.2$start, blast5.2$cs.pos.start)

blast6 = read.blast("chr6D.vs.tauschii.blast")
blast6.2 = process.blast1(blast6)
# plot(blast6.2$start, blast6.2$cs.pos.start)


blast7 = read.blast("chr7D.vs.tauschii.blast")
blast7.2 = process.blast1(blast7)
# plot(blast7.2$start, blast7.2$cs.pos.start)

library(Biostrings)

tau = readDNAStringSet("/home/ac14037/project.phd.main/genome_assemblies/aegilops.tauschii/chromosomesv2.fa")

all.blast = list(blast1.2, blast2.2, blast3.2, blast4.2, blast5.2, blast6.2, blast7.2)

source("/home/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")
allbseq = lapply(all.blast, function(x){
    extract.sequence_v2(tau, x, 1:nrow(x), 1000, 1000)    
})

# s(allbseq, "/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/allbseq", "tauschii.blast.processing.R")

cs.seqs1 = lapply(list.files("/local-scratch/alex/Watkins_exome_capture_data/all.cov.seqs", full.names = T), function(x){
    readDNAStringSet(x)
})

library(msa)
library(parallel)
count1 = 1

# all.alignments1 = Map(function(x, y){
#     g = multi.str.split(names(y), "tau_", 2)
#     x2 = x[match(g, names(x))]
#     print(paste0("doing stuff ", count1))
#     count1 <<- count1 + 1
#     mclapply(1:length(x2), function(z){
#         g2 = msa(c(x2[z], y[z]), method = "Muscle")
#         g2
#     }, mc.cores = 60)
# }, cs.seqs1, allbseq)

# s(all.alignments1, "/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.alignments1", "tauschii.blast.processing.R")
load("/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.alignments1")


#verifying that my SNP placement was correct
alignment.coords1 = as.numeric(multi.str.split(multi.str.split(lapply(all.alignments1[[1]], function(x) rownames(x)[1]), "_", 2), "-", 1))
which.min(abs(alignment.coords1 - 4868023))

load("/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.alignments1")

count1 = 1

#for each MSA, get the number of mismatching bases between sequences
msa.mismatch.sites = mclapply(all.alignments1, function(q){
    mclapply(q, function(x){    
        name1 = rownames(x)[1]
        g = as.data.frame(as.matrix(x))
        
        #get position at which to end the first cut (trimming initial "---------" sequences from alignment)
        cut.coords1 = which(sapply(g, function(z){
            z[1] != "-" & z[2] != "-"
            }))

        g = g[, -(1:(min(cut.coords1) - 1))]

        #get position at which to start the second cut
        cut.coords2 = which(sapply(g, function(z){
            z[1] != "-" & z[2] != "-"
            }))

        g = g[, -((max(cut.coords2) + 1):ncol(g))]



        # png(paste0("/home/ac14037/project.phd.main/rotation1scripts_v4/plots/msa.test/", name1, ".png"), 800, 800)
        # dev.off()
        # plot(sapply(g, nunique))
        sapply(g, nunique)
        
        # count1 <<- count1 + 1
    }, mc.cores = 50)
})

# s(msa.mismatch.sites, "~/project.phd.main/rotation1scripts_v4/saved.objects/msa.mismatch.sites", "tauschii.blast.processing.R")

load("~/project.phd.main/rotation1scripts_v4/saved.objects/msa.mismatch.sites")

#return the number of mismatches for each alignment
msa.mismatch.analysis = lapply(msa.mismatch.sites, function(y){
    sapply(y, function(x){
        length(which(x > 1))
    })
})

msa.mismatch.analysis2 = sum(unlist(lapply(msa.mismatch.sites, function(y){
    sapply(y, function(x){
        length(which(x <= 1))
    })
})))

msa.mismatch.analysis / msa.mismatch.analysis2


msa.mismatch2 = do.call(c, msa.mismatch.analysis)

length(which(msa.mismatch2 < 0.1)) / length(msa.mismatch2)

length(which(msa.mismatch2 < 0.1))

hist(msa.mismatch2)


lapply(msa.mismatch.analysis, function(x){
        (length(which(x > 0.1)) / length(x)) * 100
        # which(x > 0.3)
    })

count1 = 1
all.final.snps2 = lapply(all.alignments1, function(x){    
    processing.test1 = function(z){
        
        g3 = as.data.frame(as.matrix(z))        

        #all hyphens need to be removed from the CS sequence otherwise the coordinates in the SNP dataframe will be incorrect
        ins.pos1 = which(g3[1, ] == "-")
        if(length(ins.pos1) > 0){
                g4 = g3[, -ins.pos1]
            } else {
                g4 = g3
            }
        
        if((length(which(g4[1, ] != g4[2, ])) / ncol(g4)) > 0.4){
            snp.positions = "divergent"
        } else {        

            snp.positions = which(g4[1, ] != g4[2, ])
            snp.genotypes = g4[2, which(g4[1, ] != g4[2, ])]
            cs.genotypes = g4[1, which(g4[1, ] != g4[2, ])]
            
            chr1 = strsplit(rownames(g4)[[1]], "_")
            chr1 = chr1[[1]][1]
            start1 = strsplit(rownames(g4)[[1]], "_")
            start1 = start1[[1]][2]
            start1 = strsplit(start1, "-")
            start1 = as.integer(start1[[1]][1])                    
        }
        
        if(length(snp.positions) == 0){
            final.snps = NULL
        } else if(class(snp.positions) == "character"){
            final.snps = "divergent"
        } else {
            final.snps = data.frame(chr1, snp.positions)
            #I need to take away 1 here for the SNP positions to refer to the correct location in the CS genome
            final.snps$snp.positions = (final.snps$snp.positions + start1 - 1)
            final.snps$genotypes = as.character(convert.to.character.data.frame(snp.genotypes)[1, ])
            final.snps$cs.genotypes = as.character(convert.to.character.data.frame(cs.genotypes)[1, ])
            final.snps
        }
        
        final.snps
    }
        g = mclapply(x, function(qq){
            tryCatch(processing.test1(qq), error = function(e) "ERROR")
        }, mc.cores = 60)
        print(paste("done", count1))
        count1 <<- count1 + 1
        g
    })

# s(all.final.snps2, "/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.final.snps2", "tauschii.blast.processing.R")
load("/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.final.snps2")

#check for errors
lapply(1:7, function(z){
    unlist(all.final.snps2[[z]][which(sapply(all.final.snps2[[z]], function(x) length(x)) == 1)])    
    })


#check for divergence
g = sapply(all.final.snps2[[1]], function(x) class(x) == "character")


# load("/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/allbseq")
# load("/home/ac14037/project.phd.main/rotation1scripts_v4/saved.objects/all.final.snps2")
names(all.final.snps2) = 1:7

all.final.snps2.5 = lapply(all.final.snps2, function(x){
    null.coords1 = which(sapply(x, function(y) class(y)) == "NULL")
    divergent.coords1 = which(sapply(x, function(y) class(y)) == "character")
    all.coords1 = c(null.coords1, divergent.coords1)
    x[-all.coords1]
})

all.final.snps3 = lapply(all.final.snps2.5, bind_rows)
all.final.snps4 = bind_rows(all.final.snps3)

#have to remove Ns here rather than earlier as otherwise it would mess up the coordiantes of the bases
afs4 = all.final.snps4[-which(all.final.snps4$genotypes == "N"), ]
afs4 = afs4[-which(afs4$cs.genotypes == "N"), ]

#### BEGIN INTEGRATION OF A. TAUSCHII SNPS AND ALL WATKINS SNPS #### 

filtered.coords1 = lapply(allbseq, function(x){
    multi.str.split(names(x), "D_", 2)
})

filtered.coords2 = lapply(filtered.coords1, function(x){
    startcoord = as.numeric(multi.str.split(x, "-", 1))
    endcoord = as.numeric(multi.str.split(x, "-", 2))
    df1 = data.frame(startcoord, endcoord)
    df2 = df1[sort(df1$startcoord, index.return = T)$ix, ]
    df2
})

all.snps1 = read.csv("/local-scratch/alex/Watkins_exome_capture_data/all.snps_w_dels.csv", stringsAsFactors = F)

all.snps2 = split(all.snps1, all.snps1$chr)
all.snps3 = all.snps2[grep("D", names(all.snps2))]

#trim down all.snps1 to only the SNPs located in the BLAST regions that we identified in 
#Aegilops tauschii
all.snps4 = Map(function(x, y){
    g = lapply(x$pos, function(z){
        which(y$startcoord <= z & y$endcoord >= z)
    })
    length(which(sapply(g, length) > 0))
    g2 = x[which(sapply(g, length) > 0), ]
    g2
}, all.snps3, filtered.coords2)


all.snps5 = bind_rows(all.snps4)
# write.csv(all.snps5, "/local-scratch/alex/Watkins_exome_capture_data/all.snps.blast.filtered.csv", row.names = F)


#check which of the Watkins SNPs are also present in the Tauschii SNPs. Those that are not 
#will need to be added to the integrated SNPs dataframe. 
length(which(afs4$snp.positions %in% all.snps5$pos))


afs5 = split(afs4, afs4$chr1)
all.snps6 = split(all.snps5, all.snps5$chr)
#afs5 = tauschii snps (w/ deletions w.r.t CS)
#all.snps6 = watkins snps (w/ deletions w.r.t CS)
integrated.snps = Map(function(x, y){
    y$tauschii = ""
    x = x[sort(x$snp.positions, index.return = T)$ix, ]
    x.present = x[which(x$snp.positions %in% y$pos), ]
    y[match(x.present$snp.positions, y$pos), ]$tauschii = x.present$genotypes
    
    
    x.notpresent = x[which(!x$snp.positions %in% y$pos), ]
    
    colnames(x.notpresent)[1:2] = c("chr", "pos")
    colnames.to.add = colnames(y)
    colnames.to.add = colnames.to.add[-(1:2)]
    
    for(i in colnames.to.add){
        x.notpresent[[i]] = x.notpresent$cs.genotypes
    }
    
    x.notpresent$cov = 0
    x.notpresent$tauschii = x.notpresent$genotypes
    x.notpresent$ref = x.notpresent$cs.genotypes    
    x.notpresent = x.notpresent[, -c(3, 4)]
    
    integrated.snps = bind_rows(y, x.notpresent)
    integrated.snps$tauschii[which(integrated.snps$tauschii == "")] = integrated.snps$ref[which(integrated.snps$tauschii == "")]
    
    integrated.snps
    
}, afs5, all.snps6)

integrated.snps2 = bind_rows(integrated.snps)

all.vcf.files = c("Sample_10-23_1190092/all2.vcf", "Sample_11-24_1190110/all2.vcf", "Sample_12-25_1190126/all2.vcf", "Sample_1-3_1190034/all2.vcf", "Sample_13-26_1190128/all2.vcf", "Sample_14-27_1190139/all2.vcf", "Sample_15-28_1190145/all2.vcf", "Sample_16-29_1190149/all2.vcf", "Sample_17-30_1190160/all2.vcf", "Sample_18-31_1190166/all2.vcf", "Sample_19-32_1190181/all2.vcf", "Sample_1a-1_1190141/all2.vcf", "Sample_20-33_1190199/all2.vcf", "Sample_21-34_1190209/all2.vcf", "Sample_22-35_1190216/all2.vcf", "Sample_23-36_1190218/all2.vcf", "Sample_24-37_1190219/all2.vcf", "Sample_25-38_1190224/all2.vcf", "Sample_26-39_1190224/all2.vcf", "Sample_27-40_1190231/all2.vcf", "Sample_2-8_1190299/all2.vcf", "Sample_28-41_1190239/all2.vcf", "Sample_29-42_1190246/all2.vcf", "Sample_2a-2_1190292/all2.vcf", "Sample_30-43_1190254/all2.vcf", "Sample_3-14_1190007/all2.vcf", "Sample_31-44_1190264/all2.vcf", "Sample_32-45_1190273/all2.vcf", "Sample_33-46_1190281/all2.vcf", "Sample_34-47_1190291/all2.vcf", "Sample_35-48_1190300/all2.vcf", "Sample_36-50_1190313/all2.vcf", "Sample_37-51_1190324/all2.vcf", "Sample_38-52_1190326/all2.vcf", "Sample_39-53_1190349/all2.vcf", "Sample_3a-4_1190238/all2.vcf", "Sample_40-54_1190352/all2.vcf", "Sample_41-55_1190355/all2.vcf", "Sample_4-16_1190032/all2.vcf", "Sample_42-56_1190360/all2.vcf", "Sample_43-57_1190387/all2.vcf", "Sample_44-58_1190396/all2.vcf", "Sample_45-59_1190398/all2.vcf", "Sample_46-60_1190406/all2.vcf", "Sample_47-61_1190420/all2.vcf", "Sample_48-62_1190433/all2.vcf", "Sample_49-63_1190440/all2.vcf", "Sample_4a-5_1190308/all2.vcf", "Sample_50-65_1190451/all2.vcf", "Sample_51-66_1190460/all2.vcf", "Sample_5-17_1190040/all2.vcf", "Sample_52-67_1190468/all2.vcf", "Sample_53-68_1190471/all2.vcf", "Sample_54-69_1190474/all2.vcf", "Sample_55-70_1190475/all2.vcf", "Sample_56-71_1190481/all2.vcf", "Sample_57-72_1190483/all2.vcf", "Sample_58-73_1190496/all2.vcf", "Sample_59-74_1190507/all2.vcf", "Sample_5a-6_1190103/all2.vcf", "Sample_60-75_1190546/all2.vcf", "Sample_61-76_1190551/all2.vcf", "Sample_6-18_1190042/all2.vcf", "Sample_62-77_1190560/all2.vcf", "Sample_63-78_1190562/all2.vcf", "Sample_64-81_1190579/all2.vcf", "Sample_65-82_1190580/all2.vcf", "Sample_66-83_1190591/all2.vcf", "Sample_67-85_1190624/all2.vcf", "Sample_68-86_1190629/all2.vcf", "Sample_69-87_1190637/all2.vcf", "Sample_6a-7_1190827/all2.vcf", "Sample_70-88_1190639/all2.vcf", "Sample_71-89_1190651/all2.vcf", "Sample_7-19_1190044/all2.vcf", "Sample_72-90_1190652/all2.vcf", "Sample_73-91_1190662/all2.vcf", "Sample_74-92_1190670/all2.vcf", "Sample_75-93_1190671/all2.vcf", "Sample_76-95_1190683/all2.vcf", "Sample_77-96_1190685/all2.vcf", "Sample_78-97_1190690/all2.vcf", "Sample_79-98_1190694/all2.vcf", "Sample_7a-9_1190627/all2.vcf", "Sample_80-99_1190698/all2.vcf", "Sample_81-100_1190700/all2.vcf", "Sample_8-20_1190045/all2.vcf", "Sample_82-101_1190704/all2.vcf", "Sample_83-102_1190705/all2.vcf", "Sample_84-103_1190707/all2.vcf", "Sample_85-104_1190722/all2.vcf", "Sample_86-105_1190731/all2.vcf", "Sample_87-107_1190740/all2.vcf", "Sample_88-108_1190742/all2.vcf", "Sample_89-109_1190746/all2.vcf", "Sample_8a-10_1190811/all2.vcf", "Sample_90-110_1190747/all2.vcf", "Sample_91-111_1190749/all2.vcf", "Sample_9-21_1190079/all2.vcf", "Sample_92-112_1190753/all2.vcf", "Sample_93-113_1190771/all2.vcf", "Sample_94-114_1190777/all2.vcf", "Sample_95-115_1190784/all2.vcf", "Sample_96-116_1190788/all2.vcf")

all.vcf.files = gsub("all2.vcf", "", all.vcf.files)
colnames(integrated.snps2)[grep("allele", colnames(integrated.snps2))] = all.vcf.files





#### add geo data to integrated.snps2 df ####

library(readxl)
wilkins.data = read_excel("/home/ac14037/project.phd.main/rotation1scripts_v4/original_data/wild.relative.accessions/winfieldsup1.xlsx", sheet = 2)
wilkins.data$`Watkins Number` = gsub("Watkins_", "", wilkins.data$`Watkins Number`)
wilkins.data$`Watkins Number` = paste0("1190", wilkins.data$`Watkins Number`)
wilkins.data$`Country of Origin` = gsub(" ", "_", wilkins.data$`Country of Origin`)

for(i in 1:nrow(wilkins.data)){
    colnames(integrated.snps2) = gsub(wilkins.data$`Watkins Number`[i], 
                                                                        paste0(wilkins.data$`Country of Origin`[i],
                                                                                     "_", wilkins.data$`Watkins Number`[i]), colnames(integrated.snps2))
}

# colnames(integrated.snps2) = gsub("Sample_.*?_", "", colnames(integrated.snps2)) # have two samples of the same var so can't remove sample ID
colnames(integrated.snps2) = gsub("/$", "", colnames(integrated.snps2))

#make dataframe without deletions
deletion.coords2 = lapply(integrated.snps2, function(x){
    which(x == "-")
    })

deletion.coords3 = unique(do.call(c, deletion.coords2))
integrated.snps.no.indels = integrated.snps2[-deletion.coords3, ]
s(integrated.snps.no.indels, "rotation1scripts_v4/saved.objects/integrated.snps.no.indels", "tauschii.blast.processing.R")

sequence.names = colnames(integrated.snps.no.indels)[4:(ncol(integrated.snps.no.indels))]
sequences1 = lapply(integrated.snps.no.indels[, 4:(ncol(integrated.snps.no.indels))], function(x){
    paste0(x, collapse = "")
})

interleave <- function(v1,v2)                                                                                                                                                                        
{                                                                                                                                                                                                                            
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}            

seq2 = do.call(c, sequences1) 

sequence.names = paste0(">", sequence.names)
fasta1 = interleave(sequence.names, seq2)
writeLines(fasta1, "/local-scratch/alex/Watkins_exome_capture_data/full_sequences_w_tauschii_NO_indels.fa")






#### DELETIONS SHOULD BE INTEGRATED INTO THIS DATAFRAME, NOW NEED TO DEAL WITH INSERTIONS #### 


int.snps.w.dels = integrated.snps2
s(int.snps.w.dels, "~/project.phd.main/rotation1scripts_v4/saved.objects/int.snps.w.dels", "tauschii.blast.processing.R")
load("rotation1scripts_v4/saved.objects/int.snps.w.dels")
integrated.snps2 = int.snps.w.dels

length(which(integrated.snps2[, 109] == "-"))
length(which(integrated.snps2[, 107] == "-"))



count1 = 1
all.final.snps.insertions = lapply(all.alignments1, function(x){    
    processing.test1 = function(z){
        

        g3 = as.data.frame(as.matrix(z)) 
        g3 = convert.to.character.data.frame(g3)     
        #get position at which to end the first cut (trimming initial "---------" sequences from alignment)
        cut.coords1 = which(sapply(g3, function(z){
            z[1] != "-" & z[2] != "-"
            }))

        #trim start buffer
        g3 = g3[, -(1:(min(cut.coords1) - 1))]

        #get position at which to start the second cut
        cut.coords2 = which(sapply(g3, function(z){
            z[1] != "-" & z[2] != "-"
            }))

        #trim end buffer
        g3 = g3[, -((max(cut.coords2) + 1):ncol(g3))]

        g4 = g3        
        
        ins.coords1 = which(g4[1, ] == "-")             
        if((length(which(g4[1, ] != g4[2, ])) / ncol(g4)) > 0.4){
            ranges1 = "divergent"
        } else {
            if(length(ins.coords1) > 0){
            
            diffs <- c(1, diff(ins.coords1))
            start_indexes <- c(1, which(diffs > 1))
            end_indexes <- c(start_indexes - 1, length(ins.coords1))
            coloned <- paste(ins.coords1[start_indexes], ins.coords1[end_indexes], sep=":")
            ranges1 = data.frame(coloned, stringsAsFactors = F)
            r1 = multi.str.split(ranges1[, 1], ":", 1)
            r2 = multi.str.split(ranges1[, 1], ":", 2)
            ranges1 = data.frame(r1, r2, stringsAsFactors = F)
            ranges1$length = (as.numeric(ranges1$r2) - as.numeric(ranges1$r1) + 1)
            
            init.length = 0
            ranges1$insert.pos = ""

            #algorithm for calculating insertion coordinates
            #for each range of insertions, calculate insertion position based on the length of all the previous insertions (translates the MSA coordinates into the CS coordinates)
            for(i in 1:nrow(ranges1)){
                insert.position = as.numeric(ranges1[i, ]$r1) - as.numeric(init.length) - 1 #take away 1, the hyphens are therefore inserted after the base located at insert.position
                ranges1$insert.pos[i] = insert.position
                init.length = as.numeric(init.length) + as.numeric(ranges1[i, ]$length)
            }

            ranges1$insertion = ""
            for(i in 1:nrow(ranges1)){
                insertion1 = paste0(g4[2, ranges1[i, ]$r1:ranges1[i, ]$r2], collapse = "")
                ranges1$insertion[i] = insertion1
            }

            chr1 = strsplit(rownames(g4)[[1]], "_")
            chr1 = chr1[[1]][1]
            start1 = strsplit(rownames(g4)[[1]], "_")
            start1 = start1[[1]][2]
            start1 = strsplit(start1, "-")
            start1 = as.integer(start1[[1]][1])                    

            ranges1$cs.insert.pos = as.numeric(ranges1$insert.pos) + start1 - 1
            ranges1$chr = chr1            
            } else {
                ranges1 = "no.inserts"
            }
        } 
    
    ranges1

    }
        
        g = mclapply(x, function(qq){
            tryCatch(processing.test1(qq), error = function(e) "ERROR")
        }, mc.cores = 60)
        print(paste("done", count1))
        count1 <<- count1 + 1
        g
    })




afs.ins = bind_rows(lapply(all.final.snps.insertions, function(x){
    g = sapply(x, function(y){
            class(y)
        })

    x2 = x[-which(g == "character")]
    x3 = bind_rows(x2)
    x3 = x3[, c(7, 6, 5)]
    x3
}))

afs.ins$var = "tauschii"

load("~/project.phd.main/rotation1scripts_v4/saved.objects/all.vcf.ins4")

all.vcf.ins4$cs.insert.pos = as.numeric(all.vcf.ins4$cs.insert.pos)

all.inserts1 = bind_rows(all.vcf.ins4, afs.ins)

# s(all.inserts1, "~/project.phd.main/rotation1scripts_v4/saved.objects/all.inserts1", "tauschii.blast.processing.R")

#update variety names for inserts

for(i in 1:nrow(wilkins.data)){
    all.inserts1$var = gsub(wilkins.data$`Watkins Number`[i], 
                                                                        paste0(wilkins.data$`Country of Origin`[i],
                                                                                     "_", wilkins.data$`Watkins Number`[i]), all.inserts1$var)
}

# all.inserts1$var = gsub("Sample_.*?_", "", all.inserts1$var)
all.inserts1$var = gsub("/$", "", all.inserts1$var)


s(all.inserts1, "~/project.phd.main/rotation1scripts_v4/saved.objects/all.inserts1", "tauschii.blast.processing.R")



all.inserts2 = split(all.inserts1, all.inserts1$chr)


total.num.rows.to.add = sum(sapply(all.inserts2, function(x){
    g = split(x, x$cs.insert.pos)
    sum(sapply(g, function(y){
        max(sapply(y$insertion, nchar))
        }))
    }))


#preallocating a dataframe with rows for all of the insertions
#this needs to be done as adding rows to very large dataframes in R is very slow
new.data.frame.w.ins = as.data.frame(matrix(data = "", nrow = (nrow(integrated.snps2) + total.num.rows.to.add), ncol = ncol(integrated.snps2)))
new.data.frame.w.ins2 = convert.to.character.data.frame(new.data.frame.w.ins)

new.data.frame.w.ins2[1:nrow(integrated.snps2), 1:ncol(integrated.snps2)] = integrated.snps2
colnames(new.data.frame.w.ins2) = colnames(integrated.snps2)

new.data.frame.w.ins2$pos = as.numeric(new.data.frame.w.ins2$pos)
new.data.frame.w.ins2$pos2 = 0
#testing insertion speed - much faster than adding new rows
# new.data.frame.w.ins2[6100000:6100010, 20] = c("A", "A", "G", "T", "G", "C", "A", "T", "C", "G", "A")
# new.data.frame.w.ins2[6100000:6100010, 20] = ""

#format inserts to match integrated.snps dataframe
count1 = 1
all.inserts3 = lapply(all.inserts2, function(x){
    g = split(x, x$cs.insert.pos)
    g2 = lapply(g, function(y){                        
            max.ins.length = max(sapply(y$insertion, nchar))

            new.ins2 = sapply(y$insertion, function(z){
                    if(nchar(z) != max.ins.length){
                        z = paste0(z, paste0(rep("-", (max.ins.length - nchar(z))), collapse = ""))
                    }
                    z
                })

            y$insertion = new.ins2

            ins3 = lapply(y$insertion, function(x){
                    g = strsplit(x, "")[[1]]
                    as.data.frame(g)
                })

            ins4 = bind_cols(ins3)
            colnames(ins4) = y$var
            ins4$chr = y$chr[1]
            ins4$pos = y$cs.insert.pos[1]
            ins4$pos2 = 1:nrow(ins4)
            ins4$cov = 100

            full.ins.df.template = new.data.frame.w.ins2[1:nrow(ins4), ]
            full.ins.df.template[, ] = "-"
            col1coord = match(colnames(ins4), colnames(full.ins.df.template))
            
            #debugging
            # test.func1 = function(){
            #     full.ins.df.template[, col1coord] = ins4
            # }
            # tryCatch(test.func1(), error = function(e) browser())
            
            full.ins.df.template[, col1coord] = ins4
            
            full.ins.df.template

        })
    print(paste0("done ", count1))
    count1 <<- count1 + 1
    g2
    })

all.inserts4 = bind_rows(lapply(all.inserts3, function(x){
    bind_rows(x)
    }))


s(all.inserts4, "~/project.phd.main/rotation1scripts_v4/saved.objects/all.inserts4", "tauschii.blast.processing.R")

integrated.snps2$pos2 = 0
integrated.snps3 = bind_rows(integrated.snps2, all.inserts4)




integrated.snps4 = split(integrated.snps3, integrated.snps3$chr)

integrated.snps4 = lapply(integrated.snps4, function(x){
    x[sort(x$pos, index.return = T)$ix, ]
    })

integrated.snps5 = bind_rows(integrated.snps4)


s(integrated.snps5, "~/project.phd.main/rotation1scripts_v4/saved.objects/integrated.snps5", "tauschii.blast.processing.R")


#### remove non-snp loci ####

load("rotation1scripts_v4/saved.objects/integrated.snps.no.indels")

count1 = 1
snp.verify = apply(integrated.snps.no.indels[, c(4:107, 109)], 1, function(x){
    g = length(which(x != integrated.snps.no.indels$ref[count1]))
    count1 <<- count1 + 1
    g
})

integrated.snps.no.indels2 = integrated.snps.no.indels[which(snp.verify != 0), ]

# write.csv(integrated.snps3, "/local-scratch/alex/Watkins_exome_capture_data/all.integrated.snps.v2.csv", row.names = F)



sequence.names = colnames(integrated.snps.no.indels2)[4:(ncol(integrated.snps.no.indels2) - 1)]
sequences1 = lapply(integrated.snps.no.indels2[, 4:(ncol(integrated.snps.no.indels2) - 1)], function(x){
    paste0(x, collapse = "")
})

interleave <- function(v1,v2)                                                                                                                                                                        
{                                                                                                                                                                                                                            
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
}            

seq2 = do.call(c, sequences1) 

sequence.names = paste0(">", sequence.names)
fasta1 = interleave(sequence.names, seq2)
writeLines(fasta1, "/local-scratch/alex/Watkins_exome_capture_data/onlysnps_w_tauschii_no_indels.fa")


library(Biostrings)
bio1 = readDNAStringSet("/local-scratch/alex/Watkins_exome_capture_data/onlysnps_w_tauschii_no_indels.fa")

bio2 = as.data.frame(as.matrix(bio1))


total.number.of.snps = nrow(integrated.snps2)


count1 = 1
snp.verify = apply(integrated.snps.no.indels2[, c(4:107)], 1, function(x){
    g = length(which(x != integrated.snps.no.indels2$ref[count1]))
    count1 <<- count1 + 1
    g
})



#### SNP COMPARISONS ####
watkins.snp.comparisons = sapply(integrated.snps.no.indels[, c(4:107)], function(x){
    length(which(x != integrated.snps.no.indels$ref))
})

watkins.snp.comparisons2 = lapply(integrated.snps.no.indels[, c(4:107, 109)], function(x){
    which(x != integrated.snps.no.indels$ref)
})

nunique(unlist(watkins.snp.comparisons2))


tauschii.snp.comparison = length(which(integrated.snps.no.indels$tauschii != integrated.snps.no.indels$ref))

mean(watkins.snp.comparisons)
sd(watkins.snp.comparisons)



### generate tauschii vcf file to functional SNP annotation #### 

tauschii.snps1 = integrated.snps2[which(integrated.snps2$tauschii != integrated.snps2$ref), 
                                                                    c(1, 2, 108, 109)]
tauschii.snps1 = add.column.at.position(tauschii.snps1, 2)
colnames(tauschii.snps1) = c("CHROM", "POS", "ID", "REF", "ALT")
tauschii.snps1$QUAL = 100
tauschii.snps1$FILTER = "."
tauschii.snps1$INFO = ""

write.table(tauschii.snps1, "/local-scratch/alex/Watkins_exome_capture_data/tauschii.snps.vcf",    sep = "\t", row.names = F, quote = F)


# functional.snp1 = read.table("~/temp2/ensembl-vep/output.files/all.snps.output.txt", comment.char = "#", sep = "\t",
                                                         #stringsAsFactors = F)

functional.snp1 = read.table("~/temp2/ensembl-vep/output.files/all/all.snps.no.tauschii.txt", comment.char = "#", sep = "\t", stringsAsFactors = F)
functional.snp1.d = functional.snp1[grep("D_", functional.snp1$V1), ]

functional.snp1.d$chromo = multi.str.split(functional.snp1.d$V2, ":", 1)
functional.snp1.d$pos = multi.str.split(functional.snp1.d$V2, ":", 2)
functional.snp1.d$pos = as.numeric(functional.snp1.d$pos)


all.vcf.comb = read.table("/local-scratch/alex/Watkins_exome_capture_data/all.filtered.coverage.combined.vcf", sep = "\t", stringsAsFactors = F)
all.vcf.comb.d = all.vcf.comb[grep("D", all.vcf.comb$V1), ]

which(!functional.snp1.d$pos %in% all.vcf.comb.d$V2)


#### some basic stats ####

