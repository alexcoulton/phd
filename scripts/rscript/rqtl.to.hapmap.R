#!/usr/bin/Rscript
library(readr)
source("~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")

vcf = read_delim("~/project.phd.main/rotation1scripts_v4/original_data/vcf/820k.vcf", delim = "\t", comment = "#", col_names = F)
# ult = read_csv("~/project.phd.main/rotation1scripts_v4/original_data/genotype.data/820k.genotypes/all.geno.cons.order.barley.4b.rev.csv")
ult = read_csv("~/project.phd.main/rotation1scripts_v4/original_data/genotype.data/820k.genotypes/all.geno.bar.cons.d2.csv")
ult = ult[, 1:ncol(ult)]

vcf1 = vcf[match(ult$probe_ID, vcf$X3), ]

library(dplyr)
library(tibble) 

colnames(ult)[1:3] = c("chrom", "pos", "rs#")

ult = ult[, c(3, 1, 2, 4:ncol(ult))]

ult = add_column(ult, .before = "chrom", alleles = paste(vcf1$X4, "/", vcf1$X5, sep = ""))

ult = add_column(ult, .before = "BELLAROI", strand = ".")
ult = add_column(ult, .before = "BELLAROI", `assembly#` = NA)
ult = add_column(ult, .before = "BELLAROI", center = NA)
ult = add_column(ult, .before = "BELLAROI", protLSID = NA)
ult = add_column(ult, .before = "BELLAROI", assayLSID = NA)
ult = add_column(ult, .before = "BELLAROI", panelLSID = NA)
ult = add_column(ult, .before = "BELLAROI", QCcode = NA)

# ult$REF = vcf1$X4
# ult$ALT = vcf1$X5
#ult = as.data.frame(ult)

count1 = make.counter()

ap.res2 = apply(ult, 1, function(q){
	print(count1())
	a.allele = substr(q[2], 1, 1)
	b.allele = substr(q[2], 3, 3)
	g = lapply(q[12:length(q)], function(x){
		x[x == "AA"] = p(a.allele, a.allele)
		x[x == "AB"] = p(a.allele, b.allele)
		x[x == "BB"] = p(b.allele, b.allele)
		x
	})
	q[12:ncol(ult)] = unlist(g)
	q
	# browser()
})

ap.res2 = t(ap.res2)
ap.res2[ap.res2 == "NoCall"] = "NN"

ult = as.data.frame(ap.res2)

# for(i in 1:nrow(ult)){
# 	print(count1())
# 	a.allele = substr(ult[i, 2], 1, 1)
# 	b.allele = substr(ult[i, 2], 3, 3)
# 	g = lapply(ult[i, 12:ncol(ult)], function(x){
# 		x[x == "AA"] = p(a.allele, a.allele)
# 		x[x == "AB"] = p(a.allele, b.allele)
# 		x[x == "BB"] = p(b.allele, b.allele)
# 		x
# 	})
# 	ult[i, 12:ncol(ult)] = unlist(g)
# }

# ult[ult == "AA"] = "0/0"
# ult[ult == "BB"] = "1/1:1000"
# ult[ult == "AB"] = "0/1:1000"
# ult[ult == "NoCall"] = "./.:1000"



ult$chrom = unlist(lapply(substr(ult$chrom, 4, 5), function(x){
	match(x, listofwheatchromosomes)
}))

write_delim(ult[, 1:ncol(ult)], "~/project.phd.main/rotation1scripts_v4/original_data/genotype.data/820k.genotypes/bar.cons.d2.hapmap", delim = "\t")

