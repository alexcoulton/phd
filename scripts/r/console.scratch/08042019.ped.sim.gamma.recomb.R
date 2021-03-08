mean(rgamma(1000000, shape = 2.63, scale = (2.63 / 0.5)))


mean(rgamma(1000000, shape = 2.63, scale = (2.63 / 0.5)))

g = rexp(100000, rate = 2)
(length(which(g > 1.3)) / length(g)) * 100


haldane1 = read.delim("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axc1a.pop300.haldane.output_genotypes.dat", sep = "\t", stringsAsFactors = F)
table(sapply(haldane1, function(x) length(unique(x))))

(174 / (174 + 432)) * 100


hist(rpois(100000, 1.3), breaks = 100)

ptm <- proc.time()
gen.data(225, 300, prob.vec = prob.vec3)
proc.time() - ptm


