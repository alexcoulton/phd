#chinese spring x paragon rqtl analysis

library(qtl)

axiomdata = read.cross("csv", "./processed_data/genotypes/chinese.spring.x.paragon/", "cs.x.para.flipped.csv", estimate.map=F)

#######
#look at pariwise recombination fractions
adat = est.rf(axiomdata) #changed from markerlrt to est.rf

#do some checks for possible switched alleles
checkAlleles(adat, threshold=5)
rf=pull.rf(adat)
lod=pull.rf(adat, what="lod")

png(filename = "plots/cs.x.para.lodvsrf.png", width = 1000, height = 1000)
plot(as.numeric(rf), as.numeric(lod), xlab="recombination fraction", ylab="LOD score")
dev.off()

