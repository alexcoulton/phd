setwd("C:/Users/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
source("rotation1scripts_v4/scripts/r/generic.functions/fasta.R")

library(readr)
rimbert.data = read_delim("rotation1scripts_v4/original_data/rimbert.array/Supplemental Table S2.txt", "\t")

rimbert.data$CONTEXT_SEQUENCE = gsub("\\[A\\/G\\]", "R", rimbert.data$CONTEXT_SEQUENCE)
rimbert.data$CONTEXT_SEQUENCE = gsub("\\[C\\/G\\]", "S", rimbert.data$CONTEXT_SEQUENCE)
rimbert.data$CONTEXT_SEQUENCE = gsub("\\[C\\/T\\]", "Y", rimbert.data$CONTEXT_SEQUENCE)
rimbert.data$CONTEXT_SEQUENCE = gsub("\\[A\\/T\\]", "W", rimbert.data$CONTEXT_SEQUENCE)
rimbert.data$CONTEXT_SEQUENCE = gsub("\\[A\\/C\\]", "M", rimbert.data$CONTEXT_SEQUENCE)
rimbert.data$CONTEXT_SEQUENCE = gsub("\\[G\\/T\\]", "K", rimbert.data$CONTEXT_SEQUENCE)

rimbert.data$SNP_ID = paste(">", rimbert.data$SNP_ID, sep = "")
g = c(rbind(rimbert.data$SNP_ID, rimbert.data$CONTEXT_SEQUENCE))
writeLines(g, "bioinf/blast/probe.vs.genome.blast/rimbert.probes.fa")

gen.blast.script("rimbert.probes.fa", "probe.vs.genome.blast", culling_limit = 2)


#     ____________________________________________________________________________
#     LOAD GENETIC MAP                                                                                                                ####

rimbert.map = read_delim("rotation1scripts_v4/original_data/rimbert.array/Supplemental Table S4.txt", "\t")

v(rimbert.map)

plot.bris.vs.rimbert = function(chromosome){
    g = filter(rimbert.map, CHROMOSOME == chromosome)
    dgen = genetictophysical(chromosome, g$SNP_ID, markerformat = "-", threshold = 30, rimbert = T)
    
    
    
    ours = filter(cs.x.p.map, chr == chromosome)
    
    dours = genetictophysical(chromosome, ours$ï..marker)
    
    pdf(p("rotation1scripts_v4/plots/bristol.array.vs.rimbert.array/", chromosome, ".pdf"))
    par(mfrow = c(1, 2))
    plot(1:length(dours), dours, xlab = "Marker number", ylab = "% Physical sequence", main = p("Bristol ", chromosome))
    plot(1:length(dgen), dgen, xlab = "Marker number", ylab = "% Physical sequence", main = p("Rimbert ", chromosome))
    dev.off()
}

lapply(listofwheatchromosomes, plot.bris.vs.rimbert)
