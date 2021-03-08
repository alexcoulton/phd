setwd("C:/Users/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")

install.packages("readxl")
library(readxl)

axiomgeno = read.csv("rotation1scripts_v4/processed_data/genotypes/apogee.x.paragon.sacha/all.axiom.genotypes.for.selected.markers.csv", stringsAsFactors = F)
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


lines3 = c("Shamrock", "Shango", "Keilder", "Chinese_Spring", "Xi19", "Capelle", "Oratorio", 
                     "Longbow", "Soissons")


lines4 = c(lines2, lines3)
g1 = lapply(lines4, function(x){
    axiomgeno[grep(x, axiomgeno$Var_col), ]
    
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


markername("AX-94501412", "BB")
markername("AX-94504766", "AB")

q1 = markername("AX-95236314", "AA")

length(marker.table)

14*24


wantmore = markername("AX-94452872", "BB") # want more


q1 = markername("AX-94452872", "AA") # to remove

sapply(lines4, function(x) grep(x, q1$Var_col))


unique(wantmore$Var_col)




