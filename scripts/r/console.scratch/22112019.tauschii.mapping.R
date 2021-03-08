#the aim of this script is to split the IWGSC genome sequence into 
#smaller subsequences based on coverage of the Watkins exome capture data.
#these subsequences will then be used in a BLAST search against the Aegilops tauschii genome
#in order to find homologous sequences

cov = read.csv("all.coverage.20x.csv", stringsAsFactors = F)
cov2 = split(cov, cov$V1)



#this function takes consecutive integers. I'll be taking positions in the coverage data frame
#that are next to each other and transforming them into ranges
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


cons2 = lapply(cov2, function(x){
	g = cons1(x$V2)
	g$length= g$end - g$start
	g
})

library(Biostrings)

iwgsc = readDNAStringSet("~/project.phd.main/genome_assemblies/iwgsc/161010_Chinese_Spring_v1.0_pseudomolecules.fasta")

# count = 1

#here i will concatenate ranges that are within 7000 bp of each other.
sequence.ranges = lapply(cons2, function(x){
	end.vals = x$end[1:(nrow(x) - 1)]
	start.vals = x$start[2:nrow(x)]
	range.diffs = start.vals - end.vals
	transitions.end.coord = which(range.diffs > 7000)
	transitions.start.coord = transitions.end.coord + 1
	transitions.start.coord = c(1, transitions.start.coord)
	transitions.end.coord = c(transitions.end.coord, nrow(x))
	g1 = x$start[transitions.start.coord]
	g2 = x$end[transitions.end.coord]
	g3 = data.frame(g1, g2)
	colnames(g3) = c("start", "end")
	
	# if(count == 3) browser()
	# count <<- count + 1

	g3$length = g3$end - g3$start
	g3
	})

#extract only D genome ranges
seq.ranges.d = sequence.ranges[grep("D", names(sequence.ranges))]

list.of.dnastringset.to.stringset = function(x){
    #args:
    # x - a list of DNAStringSet objects
    newset = DNAStringSet()
    for(i in x){
        newset = c(newset, i)
    }
    newset
}

count1 = 1
all.sequences = lapply(seq.ranges.d, function(x){
	seqname = names(seq.ranges.d)[count1]
	sequences1 = apply(x, 1, function(y){		
		# g = tryCatch(DNAStringSet(iwgsc[[seqname]][y[1]:y[2]]), error = function(e) browser())
		g = DNAStringSet(iwgsc[[seqname]][y[1]:y[2]])
		names(g) = paste0(seqname, "_", y[1], "-", y[2])
		g
		})
	sequences2 = list.of.dnastringset.to.stringset(sequences1)
	count1 <<- count1 + 1
	sequences2
	})


hcgff = read.table("~/project.phd.main/rotation1scripts_v4/original_data/IWGSC/iwgsc.hc.genesonly.tab.gff", stringsAsFactors = F, header = T)

hcgff2 = split(hcgff, hcgff$V2)

hcgff3 = hcgff2[grep("D", names(hcgff2))]

seq.ranges.d

full.gene.coverage = Map(function(x, y){
	apply(y, 1, function(z){
			which(x$start < z[5] & x$end > z[6])
		})
	}, seq.ranges.d, hcgff3)

