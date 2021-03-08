library(readr)
paulvcf = read_delim("FINAL_sorted_axiom820_update_for_baron.vcf", delim = "\t", comment = "#", col_names = F)
axpvcf = read_delim("vcf/axp10degreesv2.vcf", delim = "\t", comment = "#", col_names = F)


source("/home/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")

paul.markers = multi.str.split(paulvcf$X3, ";", 2)

paulvcf2 = paulvcf[na.omit(match(axpvcf$X3, paul.markers)), ]

paulvcf2$X3 = paul.markers[na.omit(match(axpvcf$X3, paul.markers))]

axpvcf2 = axpvcf[which(axpvcf$X3 %in% paulvcf2$X3), ]


paulvcf3 = paulvcf2[match(axpvcf2$X3, paulvcf2$X3), ]