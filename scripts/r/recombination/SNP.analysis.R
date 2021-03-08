# setwd("rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/all.alignments/")

setwd("/home/ac14037/project.phd.main/rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/genomic/all.alignments/wcds/")
fasta.files = list.files("./", pattern = "[0-9].fa")

library(Biostrings)

alignments = lapply(fasta.files, readDNAMultipleAlignment)


grab.homeologous.snps_new = function(input.row, template.row, homologue.rows, multiple.alignment){
    #gets homeologous snps when there is only one genome
    #takes a DNAMultipleAlignment object and returns a numeric vector of the column coordinates containing homeologous SNPs
    #args:
    # input.row - Integer; the row of the sequence inputted by the user (usually 1)
    #template.row - Integer; the row of the sequence to design primers from (usually 2)
    # homologue.rows - numeric vector containing the row coordinates of the homologous sequences (either paralogous or homeologous)
    # multiple.alignment - DNAMultipleAlignment class 
    
    #NB. Insertions "-" in the template sequence cannot be allowed when classifying SNPs, as these
    #will be subsequently removed by get.coordinates.after.removing.hyphens(), meaning that the
    #program maps the primers to the wrong locations in the final multiple sequence alignment output.
    
    mult.align.mat1 = convert.to.character.data.frame(as.data.frame(as.matrix(multiple.alignment)))
    g = lapply(mult.align.mat1, function(x){
        # if(!all((unique(x[homologue.rows]) == "-"))) {
        #code for not allowing "-" to be called as SNPs
        # if("-" %in% x[2:length(x)]){
        #     snp = 0
        # } else {
        	if(x[input.row] == "-"){
        	    snp = 0
        	} else {
        	    remove.ns = which(x == "N")
        	    if(length(remove.ns) > 0) x = x[-remove.ns]
        	    
        	    if(length(unique(x)) > 1){
        	        snp = 1
        	    } else {
        	        snp = 0
        	    }
        	}


	    # if(x[template.row] %in% x[homologue.rows]){
	    #     snp = 0
	    # } else {
	    #     snp = 1
	    # }    
        # }        

        snp
        })
        

	q1 = strsplit(as.character(conv.mult.align.dnastringset(multiple.alignment)[1]), "")
	min.coord = min(which(!unname(unlist(lapply(q1, function(x) x == "-")))))
	max.coord = max(which(!unname(unlist(lapply(q1, function(x) x == "-")))))

    	# browser()
    

    # print("ncol(mult.align.mat1)")
    # print(ncol(mult.align.mat1))
    # print("")
    # print("length(unlist(g))")
    # print(length(unlist(g)))
    g1 = which(unlist(g) == 1)
    g1[which(g1 >= min.coord & g1 <= max.coord)]
}

conv.mult.align.dnastringset = function(mult.align){
    #converts a multiplealignment object to dnastringset
    mult.mat = as.matrix(mult.align)
    
    mult.mat3 = apply(mult.mat, 1, function(x) paste(x, collapse = ""))
    
    #convert each object in the list to a DNAStringSet object
    mult.mat4 = lapply(mult.mat3, DNAStringSet)
    
    #make new DNAStringSet object
    mult.mat5 = DNAStringSet()
    
    #add our sequences to the new DNAStringSet object one by one
    for(i in 1:length(mult.mat4)){
        mult.mat5 = c(mult.mat5, mult.mat4[[i]])
    }
    
    names(mult.mat5) = rownames(mult.align)
    
    mult.mat5
}

# source("E:/phd.project.main/rotation1scripts_v4/scripts/r/functions.R")
source("/home/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")


show.snp = function(snp.coord, mult.alignment){	
	g = conv.mult.align.dnastringset(mult.alignment)
	as.data.frame(t(as.data.frame(lapply(g, function(x) x[(snp.coord - 1):(snp.coord + 1)]))))
}



grab.homeologous.snps_new(1, 3, c(4, 5, 6, 7), alignments[[53]])


show.snp(2093, alignments[[53]])
library(parallel)
all.snps = mclapply(alignments, function(x){
	grab.homeologous.snps_new(1, 3, c(4, 5, 6, 7), x)
	}, mc.cores = 50)


all.snps[[1]]
show.snp(2001, alignments[[1]])


sapply(all.snps, length)
all.snps[[52]]
alignments[[52]]
show.snp(21989, alignments[[52]])

fasta.files[which(sapply(all.snps, length) == max(sapply(all.snps, length)))]
fasta.files[which(sapply(all.snps, length) > 500)]

all.snps[[25]][1:500]

show.snp(7574, alignments[[25]])

alignments[[23]]
all.snps[[23]]



erroneous.alignments = fasta.files[which(sapply(all.snps, length) > 100)]
lapply(erroneous.alignments, function(x){
    file.copy(x, paste0("erroneous/", x))
})


show.snp(949, alignments[[23]])



show.snp(5937, allalignments[[34]])



t1 = conv.mult.align.dnastringset(alignments[[1]])





