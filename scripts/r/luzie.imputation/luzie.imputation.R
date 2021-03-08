setwd("E:/phd.project.main/")
library(readr)
vcf1 = read_delim("rotation1scripts_v4/original_data/luzie_imputation/watkins_1190141_vs_paragon.var.vcf", delim = "\t", comment = "##")
vcf2 = read_delim("rotation1scripts_v4/original_data/luzie_imputation/watkins_1190141_vs_paragon.var.vcf", delim = "\t", comment = "##")




library(readxl)
pxwat = read_excel("rotation1scripts_v4/original_data/luzie_imputation/ParagonxWatkins141_2014_MDmap.xlsx")
parents1 = read_csv("rotation1scripts_v4/original_data/luzie_imputation/Paragon_and_Watkins141.2_axiom_gts_raw.csv")


library(RMySQL)
# install.packages("RMySQL")
g = setup.mysql.connection(T)

bacodes = pxwat$ParagonxWatkins141
bacodes2 = paste0(" OR BS_code = '", bacodes[2:length(bacodes)], "'")


bacodes[1] = paste0(bacodes[1], "'", collapse = "")
bacodes3 = paste0(c(bacodes[1], bacodes2), collapse = "")

# q1 = dbSendQuery(g, "SELECT affycode FROM axiom820new WHERE BS_code = 'BS00110651'")

q1 = dbSendQuery(g, p("SELECT BS_code, affycode, Sequence FROM axiom820new WHERE BS_code = '", bacodes3))
q2 = fetch(q1)

markers.not.in.cerealsdb = pxwat$ParagonxWatkins141[which(!pxwat$ParagonxWatkins141 %in% q2$BS_code)]
writeLines(markers.not.in.cerealsdb, "rotation1scripts_v4/original_data/luzie_imputation/markers.not.in.cereals.db.txt")


#clean sequences 
q2$Sequence = gsub("\\[", "", q2$Sequence)
q2$Sequence = gsub("\\]", "", q2$Sequence)

library(Biostrings)
seq1 = DNAStringSet(q2$Sequence)
names(seq1) = q2$BS_code
# writeXStringSet(seq1, "rotation1scripts_v4/original_data/luzie_imputation/bscode.seqs.fa")


blast1 = read.blast("bioinf/blast/genes.vs.paragon.genome/results.blast/luzie.kasp.probes.vs.paragon.blast")


luzie.axiom = read.csv("rotation1scripts_v4/original_data/IWGSC/Luzie.probes.JIC_RefSeq1_axiom_order_Nov18.csv", stringsAsFactors = F)


length(which(bacodes %in% luzie.axiom$BS.code.BLAST.derived))


parents2 = parents1[na.omit(match(luzie.axiom$SNP.id, parents1$SNP.id)), ]

pxwat2 = add.column.at.position(pxwat, 1)
pxwat2$empty.col = ""
colnames(pxwat2)[2] = "affycode"

pxwat2[match(q2$BS_code, pxwat2$ParagonxWatkins141), ]$affycode = q2$affycode

#convert seperate LGs into single chromosomes
chromo = gsub("[a-z]", "", pxwat$`LG group`)
pxwat3 = as.data.frame(cbind(chromo, pxwat2))

pxwat4 = split.df.into.list.of.dfs.by.column(pxwat3, "chromo", F)

#check how order of genetic map compares to Luzie chinese spring axiom order / chromosome assignments
lapply(pxwat4, function(x){
    
    g = match(x$ParagonxWatkins141, luzie.axiom$BS.code.BLAST.derived)
    g2 = na.omit(g)
    
    
    all(g2 == cummax(g2))
    
    length(which(diff(g2) < 0)) / length(g2)
    
})

luzie.axiom2 = split.df.into.list.of.dfs.by.column(luzie.axiom, colname.to.split = "WGA.Chr")

match(unique(pxwat2$affycode), luzie.axiom$SNP.id)



