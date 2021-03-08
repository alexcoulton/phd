setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
setwd("~/Google Drive/R/scripts/rotation1scripts/")
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

axiomcsv = read.csv("nd.csv", stringsAsFactors = F)

####stage 1 - remove conflicting parental markers####
apo=axiomcsv[c(340, 343, 345),] #select parent plants (apogee) from dataset
para=axiomcsv[c(349, 351, 353),] #select parent plants (paragon) from dataset

#define marker deletion function
#arguments: list of columns in which markers differ from all.equal, dataset from which markers are to be removed
delncm = function(acomp, parents){
    length(acomp)
    c=1
    markerstodel=list()
    i=4
    while(i<length(acomp)+1){
        t = strsplit(acomp[i], " ")[[1]][2] #parse string from acomp to get marker
        t = substring(t, 2, nchar(t)-2) #more parsing
        t
        
        markerstodel[c]=t
        c=c+1
        i=i+1
    }
    
    parents=parents[,!names(parents) %in% markerstodel]
    return(parents)
}

#setup comparisons for apogee
acomp1=all.equal(apo[1,], apo[2,]) #find markers in which parents of same cultivar differ
acomp2=all.equal(apo[2,], apo[3,]) #find markers in which parents of same cultivar differ
acomp3=all.equal(apo[1,], apo[3,]) #find markers in which parents of same cultivar differ
#setup comparisons for paragon
pcomp1=all.equal(para[1,], para[2,]) #find markers in which parents of same cultivar differ
pcomp2=all.equal(para[2,], para[3,]) #find markers in which parents of same cultivar differ
pcomp3=all.equal(para[1,], para[3,]) #find markers in which parents of same cultivar differ

#remove conflicting markers from apogee
apo=delncm(acomp1, apo)
apo=delncm(acomp2, apo)
apo=delncm(acomp3, apo)
apo=delncm(pcomp1, apo)
apo=delncm(pcomp2, apo)
apo=delncm(pcomp3, apo)

#remove conflicting markers from paragon
para=delncm(pcomp1, para)
para=delncm(pcomp2, para)
para=delncm(pcomp3, para)
para=delncm(acomp1, para)
para=delncm(acomp2, para)
para=delncm(acomp3, para)

#remove the now redundant parent individuals
apo=apo[1,]
para=para[1,]

pc = rbind(apo,para)

#remove elements of pc which are the same (parents do not differ at these loci)
pc2=pc[vapply(pc, function(x) length(unique(x)) > 1, logical(1L))] # https://stackoverflow.com/questions/30544282/how-to-remove-columns-with-same-value-in-r

#get list of markers for swapping
n = list()
for(i in 1:ncol(pc2)){
    if(pc2[2,i] == "A" & pc2[1,i] == "B")
        n=c(n,i)
}

pc3 = pc2

#swap markers according to list
for(i in 1:4854){
    if(i %in% n == TRUE){
        pc3[1,i]="A"
        pc3[2,i]="B"
    }
}

?colnames
colnames(b)
#UPDATE - NEED TO DO FLIPPING OF MAIN DATA BEFORE MORE FILTERING
#remove columns in which any of the rows contain a particular value, in this case "H"
#https://stackoverflow.com/questions/38024320/remove-columns-from-data-frame-if-any-row-contains-a-specific-string
#note the one line function within Filter() is known as an anonymous function
#pc3=Filter(function(x) !any(x=="H"), pc3)
#also remove columns containing "-"
#pc3=Filter(function(x) !any(x=="-"), pc3)
#select the columns in axiomcsv which are also in pc3
#fildata_preflip=axiomcsv[,c(which(names(axiomcsv) %in% names(pc3)))]
#fildata=axiomcsv[,c(which(names(axiomcsv) %in% names(pc3)))]
ncol(pc3)
ncol(fildata)
nrow(fildata)
#flip signs for main dataset
for(g in 1:355){
    for(i in 1:4099){
        if(i %in% n == T){
            if(fildata[g,i]=="A"){
                fildata[g,i]="B"
            }
            else if(fildata[g,i]=="B"){
                fildata[g,i]="A"
            }
        }
    }
}

write.csv(fildata, "flipped_data.csv", row.names=F)

#########BEGIN RQTL ANALYSES###########
axiomdata = read.cross("csv", "./", "flipped_data.csv", estimate.map=F)
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
lg <- formLinkageGroups(axiomdata, max.rf=0.35, min.lod=6)
table(lg[,2])
    
    