gff = read.delim("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_2017Mar13.gff3", sep = "\t", header = F)
gff = gff[which(gff$V1 == "chr1A"), ]
write.csv(gff, "rotation1scripts_v4/original_data/datasonification/chr1a.gff")

library(Biostrings)
iwgsc = readDNAStringSet("genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta")
chr1a = iwgsc[1]
writeXStringSet(chr1a, "rotation1scripts_v4/original_data/datasonification/chr1a.fa")


chr1a.sub = chr1a[[1]][3000000:4000000]
chr1a.sub = DNAStringSet(chr1a.sub)
writeXStringSet(chr1a.sub, "rotation1scripts_v4/original_data/datasonification/chr1a.sub.fa")


gff.sub = gff[which(gff$V4 >= 3000000 & gff$V4 <= 4000000), ]
write.csv(gff.sub, "rotation1scripts_v4/original_data/datasonification/chr1a.sub.gff")




# generate random DNA sequence

nucl = c("A", "T", "C", "G")
nucl.seq = sapply(1:100000, function(x){
	nucl[sample(4, 1)]
	})

nucl.seq2 = paste(nucl.seq, collapse = "")


### look at gene sequences

genes = readDNAStringSet("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_HighConf_CDS_2017Mar13.fa")

library(parallel)
genes2 = mclapply(genes, function(x){
	 tryCatch(translate(DNAStringSet(x)), error = function(e) NULL)
	}, mc.cores = 40)

genes3 = genes2[1:1000]
genes3 = genes3[-which(genes3 == "NULL")]
list.of.aastringset.to.stringset = function(x){
    #args:
    # x - a list of DNAStringSet objects
    newset = AAStringSet()
    for(i in x){
        newset = c(newset, i)
    }
    newset
}

genes4 = list.of.aastringset.to.stringset(genes3)
names(genes4) = names(genes3)


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

names(genes4) = multi.str.split(names(genes4), "gene=", 2)



writeXStringSet(genes4, "rotation1scripts_v4/original_data/datasonification/genes.fa")



#### EXAMINE FUNCTIONAL ANNOTATION ####



functional.anno = read.csv("rotation1scripts_v4/original_data/IWGSC/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0.csv", stringsAsFactors = F)

gene.coords = lapply(names(genes4), function(x){
	grep(x, functional.anno$Gene.ID)
	})

gene.annotations = lapply(gene.coords, function(x){
	functional.anno[x, ]$Human.Readable.Description[1]
	})

writeLines(unlist(gene.annotations), "rotation1scripts_v4/original_data/datasonification/genes.annotations.txt")