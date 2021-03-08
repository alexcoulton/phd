read.vcf = function(x){
    read.table(x, comment.char = "#", sep = "\t", stringsAsFactors = F)
}


all.vcf.files = paste0(list.files("../samples", full.names = T), "/all.filtered.coverage.vcf")

allv2 = lapply(all.vcf.files, read.vcf)

#split by chromosome
allv3 = lapply(allv2, function(x){
	split(x, x$V1)
	})

positions1 = lapply(1:22, function(y){
	sort(unique(unlist(lapply(allv3, function(x){
		x[[y]]$V2
	}))))
	})

geno.df = as.data.frame(matrix(nrow = sum(sapply(positions1, length)), ncol = (length(all.vcf.files) + 2)))

chromos = unique(allv2[[1]]$V1)

positions2 = Map(function(x, y){
	data.frame(y, x)
	}, positions1, chromos)

library(dplyr)

positions3 = bind_rows(positions2)

geno.df[, 1:2] = positions3

#get reference alleles to fill geno.df
all.comb = bind_rows(allv2)
all.comb2 = split(all.comb, all.comb$V1)
all.comb3 = lapply(all.comb2, function(x){
	x[match(unique(sort(x$V2)), x$V2), ]
	})

all.comb4 = bind_rows(all.comb3)

geno.df[, 3:ncol(geno.df)] = all.comb4$V4



#change values of geno.df to alt alleles
count1 = 3
lapply(allv3, function(x){
	lapply(x, function(y){
		geno.df[which(geno.df$V1 == y$V1), ][match(y$V2, geno.df[which(geno.df$V1 == y$V1), ]$V2), count1] <<- y$V5		
		})
	count1 <<- count1 + 1
	})



library(reshape2)

het.codes = unique(unlist(sapply(geno.df[, 3:ncol(geno.df)], unique)))

library(reshape2)
geno.df.melt = melt(geno.df, id = c("V1", "V2"))



colnames(geno.df.melt) = c("chr", "pos", "var", "base1")
geno.df.melt$base2 = geno.df.melt$base1
het.genos = geno.df.melt[which(nchar(geno.df.melt$base1) == 3), ]$base1

source("~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")

geno.df.melt[which(nchar(geno.df.melt$base1) == 3), ]$base2 = multi.str.split(het.genos, ",", 2)
geno.df.melt[which(nchar(geno.df.melt$base2) == 3), ]$base2 = multi.str.split(het.genos, ",", 2)
geno.df.melt[which(nchar(geno.df.melt$base1) > 3), ]$base1 = 10000
geno.df.melt[which(nchar(geno.df.melt$base2) > 3), ]$base2 = 10000

array.data1 = dcast(geno.df.melt, var ~ pos, value.var = "base1")


geno.df.melt$marker = paste0(geno.df.melt$chr, geno.df.melt$pos)

geno.df.melt$base1[which(geno.df.melt$base1 == "A")] = 1
geno.df.melt$base1[which(geno.df.melt$base1 == "T")] = 2
geno.df.melt$base1[which(geno.df.melt$base1 == "C")] = 3
geno.df.melt$base1[which(geno.df.melt$base1 == "G")] = 4
geno.df.melt$base2[which(geno.df.melt$base2 == "A")] = 1
geno.df.melt$base2[which(geno.df.melt$base2 == "T")] = 2
geno.df.melt$base2[which(geno.df.melt$base2 == "C")] = 3
geno.df.melt$base2[which(geno.df.melt$base2 == "G")] = 4

array.data1 = dcast(geno.df.melt, var ~ marker, value.var = "base1")
array.data1$var = multi.str.split(all.vcf.files, "/", 3)

array.data2 = dcast(geno.df.melt, var ~ marker, value.var = "base2")
array.data2$var = multi.str.split(all.vcf.files, "/", 3)

new.array.data = as.data.frame(matrix(nrow = (nrow(array.data1) * 2), ncol = (ncol(array.data1))))
new.array.data[seq(1, nrow(new.array.data), 2), ] = array.data1
new.array.data[seq(2, nrow(new.array.data), 2), ] = array.data2

new.array.data2 = cbind(new.array.data[, 1], new.array.data)
new.array.data2 = cbind(new.array.data2[, 1], new.array.data2)

new.array.data2[, 1] = new.array.data2[, 3]
new.array.data2[, 2] = 0
new.array.data2[, 3] = 0

write.table(new.array.data2, "~/project.phd.main/rotation1scripts_v4/original_data/STRUCTURE/exome.input.diploid.str", sep = " ", row.names = F, quote = F, col.names = F)



geno.pco = geno.df
colnames(geno.pco) = c("chr", "pos", multi.str.split(all.vcf.files, "/", 3))
write.csv(geno.pco, "~/project.phd.main/rotation1scripts_v4/original_data/STRUCTURE/geno.pco.csv", row.names = F)

s(geno.pco, "~/project.phd.main/rotation1scripts_v4/saved.objects/geno.pco", "exome.popgen.R")


