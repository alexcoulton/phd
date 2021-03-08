alab.tree = readLines("~/pipeline/jobs/alab1/all.alignments/cds_nogaps/tree.lengths.txt")
alab.tree2 = multi.str.split(alab.tree, ": ", 2)
alab.tree2 = as.numeric(alab.tree2)

random.tree = readLines("~/pipeline/jobs/random.genes/all.alignments/cds_nogaps/tree.lengths.txt")
random.tree2 = multi.str.split(random.tree, ": ", 2)
random.tree2 = as.numeric(random.tree2)
t.test(alab.tree2, random.tree2)

