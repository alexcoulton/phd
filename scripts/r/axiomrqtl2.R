#setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
setwd("~/Google Drive/R/scripts/rotation1scripts")
library(qtl)
install.packages(SNPolisher)
??snpolisher
install.packages("qtl")
####LOAD UP NEW AXIOM DATA WITH 6K SNPS (POST FILTERING W/ ANALYSIS SUITE)####
# nd = read.csv("recommendedsnpsv3.csv", stringsAsFactors = F)
# nd[nd=="AA"] = "A"
# nd[nd=="BB"] = "B"
# nd[nd=="AB"] = "H"
# nd[nd=="NoCall"] = "-"
# nd[1,]=1
# nd[1,1]=""
# nd[1,2]=""
# nd=nd[,-(6978:ncol(nd))]
# write.csv(nd, "nd.csv", row.names = F)



axiomdata = read.cross("csv", "./", "nd.csv", estimate.map=F)
axiomdata_untouched = axiomdata
?read.cross

summary(axiomdata)
plotMissing(axiomdata)

#do some plots 
?par
par(mfrow=c(1,2), las=1)
plot(ntyped(axiomdata), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(axiomdata, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

#drop individuals with low information
axiomdata = subset(axiomdata, ind=(ntyped(axiomdata)>50))

#drop markers with low information
nt.bymar = ntyped(axiomdata, "mar")
todrop = names(nt.bymar[nt.bymar < 200])
axiomdata = drop.markers(axiomdata, todrop)

#####remove duplicate individuals####
cg = comparegeno(axiomdata)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

wh = which(cg > 0.9, arr=TRUE)
wh = wh[wh[,1] < wh[,2],]
wh

g = pull.geno(axiomdata)
table(g[3,], g[11,])

for(i in 1:nrow(wh)) {
    tozero = !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],] 
    axiomdata$geno[[1]]$data[wh[i,1],tozero] = NA
}

axiomdata = subset(axiomdata, ind=-wh[,2])
#######

####remove markers with same genotypes####
dup = findDupMarkers(axiomdata, exact.only=F)

for(i in length(dup)){
 axiomdata = drop.markers(axiomdata, dup[[i]])
}

axiomdata = drop.markers(axiomdata, dup)

####inspect seg. distotion####
gt = geno.table(axiomdata)
gt[gt$P.value < 0.05/totmar(axiomdata),]

todrop = rownames(gt[gt$P.value < 1e-200,]) # drop poor markers - is p value set properly? this removes ~ 2/3 of markers
axiomdata2 = drop.markers(axiomdata, todrop)
summary(axiomdata2)

gt2 = geno.table(axiomdata2)
gt2[gt2$P.value < 0.05/totmar(axiomdata2),]

todrop = rownames(gt[gt$P.value < 1e-10,]) # drop poor markers - is p value set properly? this removes ~ 2/3 of markers
axiomdata3 = drop.markers(axiomdata, todrop)
summary(axiomdata3)

gt3 = geno.table(axiomdata3)
gt3[gt3$P.value < 0.05/totmar(axiomdata3),]

#######
#look at pariwise recombination fractions
adat = markerlrt(axiomdata)

#do some checks for possible switched alleles
checkAlleles(adat, threshold=5)
rf=pull.rf(adat)
lod=pull.rf(adat, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="recombination fraction", ylab="LOD score")
    
    