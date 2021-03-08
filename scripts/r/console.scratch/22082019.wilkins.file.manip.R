setwd("/home/ac14037/project.phd.main/rotation1scripts_v4/website/autoprimerpicker/primer/pipeline/jobs/genomic/all.alignments/all.alignments.w.cds/")


list.files()
align1 = list.files("./", pattern = ".align")
orig = list.files("./", pattern = "[0-9].fa")

source("/home/ac14037/project.phd.main/rotation1scripts_v4/scripts/r/functions.R")

align1.1 = multi.str.split(align1, "\\.", 1)
orig1.1 = multi.str.split(orig, "\\.", 1)

missing = orig[which(is.na(match(orig1.1, align1.1)))]

dir.create("missing")

file.copy(missing, "missing/")

t1 = readDNAStringSet("TraesCS2A01G208200.1.fa")

t1[[1]] = reverseComplement(t1[[1]])
writeXStringSet(t1, "rev.TraesCS2A01G208200.1.fa")
