#correcting structure population assignment (change from 0 to 1) for CLUMPPAK

pop1 = as.numeric(as.factor(wilkins.data2$Region))

files1 = list.files("software/strauto/array.results/results_f/k1/", full.names = T)
files2 = list.files("software/strauto/array.results/results_f/k2/", full.names = T)


correct.structure.files = function(k){
    files1 = list.files(paste0("software/strauto/exome.results/results_f/k", k, "/"), full.names = T)
    g = lapply(files1, readLines)
    
    g2 = lapply(g, function(x){
        x[31] = gsub("0:", "1:", x[31])
        for(i in grep("Watkins_119", g[[1]])){
            x[i] = gsub("0 :", "1 :", x[i])
        }
        x
    })
    
    g2
    
    
}

corrections1 = lapply(1:10, correct.structure.files)

write.structure.corrections = function(k){
    if(!dir.exists("software/strauto/exome.results.corrected")) dir.create("software/strauto/exome.results.corrected")
    if(!dir.exists("software/strauto/exome.results.corrected/results_f")) dir.create("software/strauto/exome.results.corrected/results_f")
    if(!dir.exists(paste0("software/strauto/exome.results.corrected/results_f/k", k))) dir.create(paste0("software/strauto/exome.results.corrected/results_f/k", k))
    
    files2 = list.files(paste0("software/strauto/exome.results/results_f/k", k, "/"), full.names = T)
    files3 = gsub("exome\\.results", "exome.results.corrected", files2)
    
    lapply(1:5, function(x){
        writeLines(corrections1[[k]][[x]], files3[[x]])
    })
    
}

lapply(1:10, write.structure.corrections)





