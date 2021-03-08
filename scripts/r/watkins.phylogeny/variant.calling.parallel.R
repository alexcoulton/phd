#!/home/ac14037/bin/Rscript

library(parallel)

input.files = paste0("/local-scratch/alex/Watkins_exome_capture_data/", sample.names, "/all.filtered.coverage.vcf")
output.files = paste0("output.files/", sample.names, ".sift.output.txt")
