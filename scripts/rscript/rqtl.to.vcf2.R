#!/usr/bin/Rscript
library(readr)

vcf = read_delim("~/project.phd.main/rotation1scripts_v4/original_data/vcf/820k.vcf", delim = "\t", comment = "#", col_names = F)
ult = read_csv("~/project.phd.main/rotation1scripts_v4/original_data/genotype.data/820k.genotypes/all.geno.cons.order.barley.4b.rev.csv")

vcf1 = vcf[match(ult$probe_ID, vcf$X3), ]

library(dplyr)
library(tibble) 

ult = add_column(ult, .before = "BELLAROI", REF = "")
ult = add_column(ult, .before = "BELLAROI", ALT = "")
ult = add_column(ult, .before = "BELLAROI", QUAL = ".")
ult = add_column(ult, .before = "BELLAROI", FILTER = ".")
ult = add_column(ult, .before = "BELLAROI", INFO = ".")
ult = add_column(ult, .before = "BELLAROI", FORMAT = "GT:DP")

ult$REF = vcf1$X4
ult$ALT = vcf1$X5
#ult = as.data.frame(ult)



ult[ult == "AA"] = "0/0:1000"
ult[ult == "BB"] = "1/1:1000"
ult[ult == "AB"] = "0/1:1000"
ult[ult == "NoCall"] = "./.:1000"

source("~/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")

ult$chromosome = unlist(lapply(substr(ult$chromosome, 4, 5), function(x){
	match(x, listofwheatchromosomes)
}))

write_delim(ult[, 1:ncol(ult)], "~/project.phd.main/rotation1scripts_v4/original_data/genotype.data/820k.genotypes/bar.cons.vcf", delim = "\t")

