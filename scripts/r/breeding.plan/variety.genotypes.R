setwd("C:/Users/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

install.packages("readxl")
library(readxl)

axiomgeno = read.csv("rotation1scripts_v4/processed_data/genotypes/apogee.x.paragon.sacha/all.axiom.genotypes.for.selected.markers.csv", stringsAsFactors = F)


cappele = read_delim("rotation1scripts_v4/original_data/genotype.data/35k/Cappele_desprez.txt", delim = "\t")

cappele[which(cappele$`Axiom probe` %in% as.character(unique(axiomgeno$Probe_row))), ]

colnames(cappele) = c(colnames(axiomgeno), "chr")

axiomgeno = rbind(axiomgeno, cappele[, 1:3])
axiomgeno = as.data.frame(axiomgeno)
axiomgeno = convert.to.character.data.frame(axiomgeno)

sachalines = read_excel("rotation1scripts_v4/original_data/sacha.genetic.maps/sacha.kasp.varieties2013paper.xlsx")


g = sapply(sachalines$Variety, function(x){
    grep(x, axiomgeno$Var_col)
})


axiomgeno2 = axiomgeno[unlist(g), ]


lapply(unique(axiomgeno2$Probe_row), function(x){
    table(axiomgeno2[which(axiomgeno2$Probe_row == x), ]$Matrix_value)
    
})

#lines not in cerealsDB:
#Badger, Bobwhite
# does "Hereford" = "Hereward" ?

lines2 = c("Apogee", "Paragon", "Avalon", "Cadenza", "Rialto", "Savannah", "Robigus", "Skyfall", "Bacanora", "Premio", "Solstice", "Spark")

lines3 = c("Shamrock", "shango", "Keilder", "Chinese_Spring", "Xi19", "Cappele", "Oratorio", "Longbow", "Soissons")

lines4 = c(lines2, lines3)

axiomgeno[grep("Cappele", axiomgeno$Var_col), ]





g1 = lapply(lines4, function(x){
    # axiomgeno[grep(x, axiomgeno$Var_col), ]
    axiomgeno[which(axiomgeno$Var_col == x), ]
    
})

g2 = combine.list.of.data.frames(g1)

marker.table = lapply(unique(g2$Probe_row), function(x){
    table(g2[which(g2$Probe_row == x), ]$Matrix_value)
})

names(marker.table) = unique(g2$Probe_row)

v(lines2)



search1 = function(varname, markername){
    if(missing(markername)) markername = "-"
    g = axiomgeno[grep(varname, axiomgeno$Var_col), ]
    
    if(markername == "-"){
        g
    } else {
        g[which(g$Probe_row == markername), ]
    }
    
}


search1("Paragon", "AX-94501412")
search1("Apogee", "AX-94501412")

search1("Paragon", "AX-94504766")
search1("Apogee", "AX-94504766")



markername = function(markername, genotype){
    axiomgeno[which(axiomgeno$Probe_row == markername & axiomgeno$Matrix_value == genotype), ]
}


markername("AX-95232967", "BB")
markername("AX-94504766", "BB")

markername("AX-95236314", "BB")



alexf2 = read_delim("rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.alex.f2/all.snps.txt", delim = "\t")

alexf2 = as.data.frame(t(alexf2))
alexf2 = convert.to.character.data.frame(alexf2)

selected.lines = read_excel("rotation1scripts_v4/processed_data/genotypes/apogee.x.paragon.sacha/selected.lines.all.3b.xlsx")


selected.lines[, ]


markers = c("AX-94504766", "AX-94457592", "AX-95232967", "AX-94785822", "AX-94716868", "AX-94851468", "AX-94907975", "AX-94457937", "AX-95236314", "AX-94402393", "AX-94531806", "AX-94445993", "AX-94501412", "AX-94452872")

alexf2 = alexf2[, which(alexf2[1, ] %in% markers)]
colnames(alexf2) = alexf2[1, ]
alexf2 = alexf2[-1, ]

alexf2 = alexf2[markers]

write.csv(alexf2, "rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.alex.f2/selected.snps.csv")
