setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#setwd("~/Google Drive/R/scripts/rotation1scripts/")
library(qtl)

fildata=read.csv("flipped_data.csv", header = T)
fildata[1,1]=""
pc3=read.csv("pc3.csv")

####do some filtering to remove any parental heterozygosities and miscalls####
#remove columns in which any of the rows contain a particular value, in this case "H"
#https://stackoverflow.com/questions/38024320/remove-columns-from-data-frame-if-any-row-contains-a-specific-string
#note the one line function within Filter() is known as an anonymous function
pc3=Filter(function(x) !any(x=="H"), pc3)
#also remove columns containing "-"
pc3=Filter(function(x) !any(x=="-"), pc3)
#select the columns in axiomcsv which are also in pc3
fildata=fildata[,c(which(names(fildata) %in% names(pc3)))]
write.csv(fildata, "flipped_data2.csv", row.names=F)

#########BEGIN RQTL ANALYSES###########
axiomdata = read.cross("csv", "./", "flipped_data2.csv", estimate.map=F)
axiomdata_untouched = axiomdata
#axiomdata=convert2riself(axiomdata) #convert cross to RIL? assumes 1:1 segregation, but removes heterozygotes (i think)
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
?drop.markers
axiomdata = drop.markers(axiomdata, dup)
summary(axiomdata)

####inspect seg. distotion (here only for 1:2:1 ratio. how check deviation from 1:1??)####
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
adat = est.rf(axiomdata) #changed from markerlrt to est.rf

#do some checks for possible switched alleles
checkAlleles(adat, threshold=5)
rf=pull.rf(adat)
lod=pull.rf(adat, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="recombination fraction", ylab="LOD score")

?formLinkageGroups

#form linkage groups
lg <- formLinkageGroups(axiomdata, max.rf=0.37, min.lod=12)
table(lg[,2])

?formLinkageGroups
?est.rf

axiomdata2=axiomdata
axiomdata2=formLinkageGroups(axiomdata2, max.rf=0.37, min.lod=12, reorgMarkers=TRUE)
plotRF(axiomdata2, alternate.chrid=TRUE)

pull.map(axiomdata2, chr=8)
pull.map(axiomdata2, chr=3)
?pull.map

######DOING SOME TESTING OF MARKER ORDERING
######Without the orderMarkers function, the markers are listed in the group as ~ 10 cM apart.
######Therefore need to do orderMarkers to get proper centimorgan distances between markers,
#can then split linkage groups apart that have distances of over ~ 50 cM between markers. 

axiomdata2 = orderMarkers(axiomdata2, chr=8) #done
axiomdata2 = orderMarkers(axiomdata2, chr=1) #in progress

axiomdata2 = orderMarkers(axiomdata2, chr=7)
axiomdata2 = orderMarkers(axiomdata2, chr=6)
axiomdata2 = orderMarkers(axiomdata2, chr=5)
axiomdata2 = orderMarkers(axiomdata2, chr=4)
axiomdata2 = orderMarkers(axiomdata2, chr=3)
axiomdata2 = orderMarkers(axiomdata2, chr=2)
axiomdata2 = orderMarkers(axiomdata2, chr=9)
axiomdata2 = orderMarkers(axiomdata2, chr=10)
axiomdata2 = orderMarkers(axiomdata2, chr=11)
axiomdata2 = orderMarkers(axiomdata2, chr=12)
axiomdata2 = orderMarkers(axiomdata2, chr=13)
axiomdata2 = orderMarkers(axiomdata2, chr=14)
axiomdata2 = orderMarkers(axiomdata2, chr=15)
axiomdata2 = orderMarkers(axiomdata2, chr=16)
axiomdata2 = orderMarkers(axiomdata2, chr=17)
axiomdata2 = orderMarkers(axiomdata2, chr=18)
axiomdata2 = orderMarkers(axiomdata2, chr=19)
axiomdata2 = orderMarkers(axiomdata2, chr=20)
axiomdata2 = orderMarkers(axiomdata2, chr=21)
axiomdata2 = orderMarkers(axiomdata2, chr=22)
axiomdata2 = orderMarkers(axiomdata2, chr=23)
axiomdata2 = orderMarkers(axiomdata2, chr=24) #only one marker in this group...

#marker distribution per linkage group
#     1        2        3        4        5        6        7        8        9     10     11     12     13     14     15     16     17     18     19     20     21     22     23     24 
#2604    364    222    202    116     73     53     50     25     24     18     12     10     10        9        7        6        3        3        2        2        2        2        1 


#going to remove 70% of markers from chromosome 1 at random (the linkage group is too big to be ordered on this computer),
#then order the markers and then compare to cerealsdb data
?orderMarkers
?pull.map
chr1=pull.map(axiomdata2, chr="1", T)
nrow(chr1)
chr1$marker=rownames(chr1)
#make a vector of marker names (size = 500)    
chr1_sample=sample(chr1$marker, 500)
#remove markers in chr1_sample from chr1_marker (so these arn't removed from data_pruned2 in subsequent removal step)
chr1_pruned=chr1$marker[-which(chr1$marker %in% chr1_sample)]
length(chr1_pruned)
length(chr1$marker)
length(chr1_sample)
head(chr1_sample)
data_pruned=read.csv("flipped_data2.csv", header=T)
data_pt=colnames(data_pruned)
data_pt=data_pt[-c(1,2)]
data_pruned2=data_pruned[,-which(data_pt %in% chr1_pruned)]
ncol(data_pruned)
ncol(data_pruned2)
data_pruned2[1,1]=""
write.csv(data_pruned2, "flipped_data2_pruned.csv")
