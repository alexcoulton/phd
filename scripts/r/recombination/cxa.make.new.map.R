all.cxa = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/polyprobes.good.txt")
all.cxa = convert.to.character.data.frame(all.cxa)



all.cxa$probeset_id[grep("M23_AxC 61-1 G12.CEL_call_code", all.cxa$probeset_id)] = "Cadenza1"
all.cxa$probeset_id[grep("O23_AxC 61-1 H12.CEL_call_code", all.cxa$probeset_id)] = "Avalon1"
all.cxa$probeset_id[grep("M23_CxA 65-1 G12.CEL_call_code", all.cxa$probeset_id)] = "Cadenza2"
all.cxa$probeset_id[grep("O23_CxA 65-1 H12.CEL_call_code", all.cxa$probeset_id)] = "Avalon2"



all.cxa2 = all.cxa

parents1 = all.cxa2[c(grep("Avalon", all.cxa2$probeset_id), grep("Cadenza", all.cxa2$probeset_id)), ]

all.cxa2 = all.cxa2[-c(grep("Avalon", all.cxa2$probeset_id), grep("Cadenza", all.cxa2$probeset_id)), ]

all.cxa3 = bind_rows(all.cxa2[1, ], parents1, all.cxa2[3:nrow(all.cxa2), ])

all.cxa4 = all.cxa3[-c(2, 4), ]

source("rotation1scripts_v4/scripts/r/axiom.data.processing/generic.axiom.data.preparation.functions.R")

combined.geno.data3 = cleanup(all.cxa4, 2:3)

combined.geno.data3 = convert.to.character.data.frame(combined.geno.data3)

combined.geno.data4 = assign.parental.genotypes(combined.geno.data3, 2:3)

write.csv(combined.geno.data4, "rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/all.replicates.flipped.geno.data.newprobeset.csv")


c.x.a.f2 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/cadenza.x.avalon.f2/c.x.a.f2.txt", c("Cadenza", "Avalon"))




