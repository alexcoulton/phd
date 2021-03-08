setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#setwd("~/Google Drive/R/scripts/rotation1scripts/")
library(qtl)
suppressPackageStartupMessages(library(gmailr)) #GMAILR - FOR SENDING EMAILS AFTER LONG COMPUTE TIMES
#install.packages("qtl")
###LOAD UP NEW AXIOM DATA WITH 6K SNPS (POST FILTERING W/ ANALYSIS SUITE)####
#nd = read.csv("recommendedsnpsv3.csv", stringsAsFactors = F)
# nd = read.table("alexaxnew(A8_incl).txt") #UPDATED 24/11/2017 - POST FIXING MISLABELLED WELLS AND INCLUSION OF B16
nd = read.table("manually curated snps.txt") #UPDATED 06/12/2017
nd=t(nd)
nd[nd=="AA"] = "A"
nd[nd=="BB"] = "B"
nd[nd=="AB"] = "H"
nd[nd=="NoCall"] = "-"
colnames(nd)=nd[1,]
nd[1,]=1
nd=cbind(nd[,1], nd)
nd[,1]=seq(0,nrow(nd)-1,1)
nd[1,1]=""
nd[1,2]=""
#bring parents to top of dataframe
nd2=rbind(nd[grep("Paragon",nd[,2]),], nd[-c(grep("Paragon",nd[,2])),])
nd3=rbind2(nd2[grep("Apogee",nd2[,2]),], nd2[-c(grep("Apogee",nd2[,2])),])
nd3=rbind(nd3[12,],nd3[-12,])

nd3[,1]=seq(0,nrow(nd)-1,1)
nd3[1,1]=""
nd3[1,2]=""

rownames(nd3)=seq(1,nrow(nd3),1)
#write out dataframe, now in format appropriate for rQTL
write.csv(nd3, "manuallycuratedsnps.csv", row.names = F)

#ROWENA DATA
#axiomcsv = read.csv("nd.csv", stringsAsFactors = F)
#MY DATA
axiomcsv = read.csv("processed_data/genotypes/allpolyhighresmarkers.csv", stringsAsFactors = F) #edited on 07/02/2018

##############################################
#GENOTYPE SIMILARITY CHECKS - CAN SKIP
##############################################
#check similarities for paragon

b=ncol(axiomcsv)
similarities=list() 
for(i in 13:nrow(axiomcsv)){
    similarities=c(similarities, (length(which(axiomcsv[12,]==axiomcsv[i,]))/b)*100) #C18_Paragon_3_samp1.CEL_call_code vs all
}

simpara=data.frame(similarities)
simpara=t(simpara)
summary(simpara)
simpara2=data.frame(simpara, list(seq(1, nrow(simpara), 1)))
subset(simpara2, simpara2[,1]>60) #id's 19, 230 and 232 in simpara2 - need to check whole probe list to see if these are mislabelled
                                                                    #adding 12: 31, 242 244

#check similarities for apogee

b=ncol(axiomcsv)
similaritiesapo=list()
for(i in 13:nrow(axiomcsv)){
    similaritiesapo=c(similaritiesapo, (length(which(axiomcsv[43,]==axiomcsv[i,]))/b)*100) #C10_Apogee_2_samp3.CEL_call_code vs all #EDIT: now checking B16 against all (was previously excluded from AAS by mistake)
}

simapo=data.frame(similaritiesapo)
simapo=t(simapo)
summary(simapo)
simapo2=data.frame(simapo, list(seq(1,nrow(simapo),1)))
subset(simapo2, simapo2[,1]>60) #ids 274, 296 and 316 v. high
                                                                #adding 12: 286, 308, 328

hist(similarities2)

#add comparison data to original dataframe (axiomcsv) - new dataframe called axiomcsv3
add=data.frame(matrix(ncol=2,nrow=12))
names(add)=names(simapo2)
simapo3=rbind(add,simapo2)
axiomcsv2=cbind(simapo3,axiomcsv)
axiomcsv2=axiomcsv2[,-2]

add=data.frame(matrix(ncol=2,nrow=12))
names(add)=names(simpara2)
simpara3=rbind(add,simpara2)
axiomcsv3=cbind(simpara3,axiomcsv2)
axiomcsv3=axiomcsv3[,-2]


#CHECK SIMILARITIES FOR ALL PROBES
b=ncol(axiomcsv)
simparaall=list() 
for(i in 13:nrow(axiomcsv)){
    simparaall=c(simparaall, (length(which(axiomcsv[12,]==axiomcsv[i,]))/b)*100) #C18_Paragon_3_samp1.CEL_call_code vs all
}

b=ncol(axiomcsv)
similaritiesapo=list()
for(i in 13:nrow(axiomcsv)){
    similaritiesapo=c(similaritiesapo, (length(which(axiomcsv[6,]==axiomcsv[i,]))/b)*100) #C10_Apogee_2_samp3.CEL_call_code vs all
}


#reorganise allprobes dataset from allprobesalexdata.R (put parents at the top - see allprobeschecker.RData for backup)
bua=grep("Paragon", allprobes[,1])
allprobes2=rbind(allprobes[c(bua),], allprobes[-c(bua),])

you=grep("Apogee", allprobes2[,1])
allprobes2=rbind(allprobes2[c(you),], allprobes2[-c(you),])

jim=allprobes2[12,]
allprobes2=allprobes2[-12,]
allprobes2=rbind(jim, allprobes2)

rownames(allprobes2)=seq(1,nrow(allprobes2),1)

gogo=ncol(allprobes2)
#allprobe comparisons for paragon - all look good
(length(which(allprobes2[8,]==allprobes2[31,]))/gogo)*100
(length(which(allprobes2[8,]==allprobes2[242,]))/gogo)*100
(length(which(allprobes2[8,]==allprobes2[244,]))/gogo)*100

#allprobe comparisons for apogee - all look good
(length(which(allprobes2[6,]==allprobes2[286,]))/gogo)*100
(length(which(allprobes2[6,]==allprobes2[308,]))/gogo)*100
(length(which(allprobes2[6,]==allprobes2[328,]))/gogo)*100

allprobes2[328,1]

grep("Apogee_2_samp1", axiomcsv3[,4])

#WAS GOING TO RENAME MISLABELLED PARENTAL WELLS IN RQTL FILE, BUT RENAMED .CEL FILES INSTEAD.
# axiomcsv[axiomcsv=='Apogee_2_samp1']="A23_Q1_10_A12"
# axiomcsv[axiomcsv=='Apogee_2_samp2']="B23_Q3_14_A12"
# axiomcsv[axiomcsv=='Apogee_1_samp1']="B24_Q4_28_A12"
# axiomcsv[axiomcsv=='Paragon_2_samp1']="C23_Q1_10_B12"
# axiomcsv[axiomcsv=='Paragon_2_samp2']="F23_Q3_14_C12"
# axiomcsv[axiomcsv=='Paragon_1_samp1']="P24_Q4_28_H12"
# 
# axiomcsv[axiomcsv=='B02_Q4_28_A1.CEL_call_code']='Paragon_misc_1'
# axiomcsv[axiomcsv=='M01_Q1_10_G1.CEL_call_code']='Paragon_misc_2'
# axiomcsv[axiomcsv=='L23_Q3_14_F12.CEL_call_code']='Paragon_misc_3'
# axiomcsv[axiomcsv=='O01_Q1_10_H1.CEL_call_code']='Apogee_misc_1'
# axiomcsv[axiomcsv=='P02_Q4_28_H1.CEL_call_code']='Apogee_misc_2'
# axiomcsv[axiomcsv=='P23_Q3_14_H12.CEL_call_code']='Apogee_misc_3'
##############################################
#axiomcsv[2:nrow(axiomcsv),1]=seq(1,nrow(axiomcsv)-1,1) #setup id column
#write.csv(axiomcsv, "alexdata4.csv", row.names=F)

####stage 1 - remove conflicting parental markers####
#ROWENA ROWS
#apo=axiomcsv[c(341,344,346),] #select parent plants (apogee) from dataset
#para=axiomcsv[c(350,352,354),] #select parent plants (paragon) from dataset

#remove excess parents
axiomcsv=axiomcsv[-c(2,5,7,10,11),]
#clean rownames and id column
rownames(axiomcsv)=seq(1,nrow(axiomcsv),1)
axiomcsv[,1]=seq(0,nrow(axiomcsv)-1,1)
axiomcsv[1,1]=""

#MY ROWS
apo=axiomcsv[c(2, 3, 4),] #select parent plants (apogee) from dataset
para=axiomcsv[c(5,6,7),] #select parent plants (paragon) from dataset

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

length(acomp1)
length(acomp2)
length(acomp3)
length(pcomp1)
length(pcomp2)
length(pcomp3)

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

####do some filtering to remove any parental heterozygosities and miscalls####
#remove columns in which any of the rows contain a particular value, in this case "H"
#https://stackoverflow.com/questions/38024320/remove-columns-from-data-frame-if-any-row-contains-a-specific-string
#note the one line function within Filter() is known as an anonymous function
pc2=Filter(function(x) !any(x=="H"), pc2)
#also remove columns containing "-"
pc2=Filter(function(x) !any(x=="-"), pc2)


#get list of markers for swapping
n = list()
for(i in 1:ncol(pc2)){
    if(pc2[2,i] == "A" & pc2[1,i] == "B")
        n=c(n,i)
}

pc3 = pc2
write.csv(pc3,"pc3.csv",row.names = F)
#swap markers according to list
for(i in 1:ncol(pc3)){
    if(i %in% n == TRUE){
        pc3[1,i]="A"
        pc3[2,i]="B"
    }
}

#UPDATE - NEED TO DO FLIPPING OF MAIN DATA BEFORE MORE FILTERING
fildata_preflip=axiomcsv[,c(which(names(axiomcsv) %in% names(pc3)))]
fildata=axiomcsv[,c(which(names(axiomcsv) %in% names(pc3)))]
#fildata=axiomcsv[,c(which(names(axiomcsv) %in% names(pc3)))]
ncol(pc3)
ncol(fildata)
nrow(fildata)

###########INITIATE FLIPPING ALGORITHM FOR MAIN DATA#########
#flip signs for main dataset #THIS TAKES A LONG TIME (~2 hours)
for(g in 1:nrow(fildata)){
    for(i in 1:ncol(fildata)){
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
test_email <- mime(
    To = "alexcoulton@gmail.com",
    From = "apodforce@gmail.com",
    Subject = "flip alg. complete again 2",
    body = "Hi Alex, the flipping algorithm is done.")
send_message(test_email)

################################################################


fildata_postflip=fildata
fildata_postflip=fildata_postflip[-c(3,4,6,7),]
fildata_preflip=fildata_preflip[-c(3,4,6,7),]


#select the columns in axiomcsv which are also in pc3
#fildata_post_fil=fildata_postflip[,c(which(names(fildata_postflip) %in% names(pc3_fil)))]
#fildata_preflip=fildata_preflip[,c(which(names(fildata_preflip) %in% names(fildata_post_fil)))]

##CHECKS TO SEE IF FLIPPING ALGORITHM WORKED##
#note: flipping algorithm will not work where H is present
#this is fine, the filtering stage (immediately above this block)
#will take care of this
length(which(fildata_postflip[2,]=="B"))
length(which(fildata_postflip[3,]=="A"))

write.csv(fildata_postflip, "processed_data/genotypes/allpolyhighresmarkers.flipped.csv", row.names=F)

ncol(fildata_postflip)
nrow(fildata_postflip)

mancur_ske=read.csv("manuallycuratedsnps_flipped_skeletononly.csv", header=T)

####SEPARATE INTO QUADRANTS####
# alexq1=rbind(fildata_postflip[c(1,2,3),], fildata_postflip[c(grep("Q1",fildata_postflip[,2])),])
# alexq2=rbind(fildata_postflip[c(1,2,3),], fildata_postflip[c(grep("Q2",fildata_postflip[,2])),])
# alexq3=rbind(fildata_postflip[c(1,2,3),], fildata_postflip[c(grep("Q3",fildata_postflip[,2])),])
# alexq4=rbind(fildata_postflip[c(1,2,3),], fildata_postflip[c(grep("Q4",fildata_postflip[,2])),])

#update 07/12/2017 for manually curated skeleton markers only
alexq1=rbind(mancur_ske[c(1,2,3),], mancur_ske[c(grep("Q1",mancur_ske[,2])),])
alexq2=rbind(mancur_ske[c(1,2,3),], mancur_ske[c(grep("Q2",mancur_ske[,2])),])
alexq3=rbind(mancur_ske[c(1,2,3),], mancur_ske[c(grep("Q3",mancur_ske[,2])),])
alexq4=rbind(mancur_ske[c(1,2,3),], mancur_ske[c(grep("Q4",mancur_ske[,2])),])

write.csv(alexq1, "alexq1_mancur_ske.csv", row.names = F)
write.csv(alexq2, "alexq2_mancur_ske.csv", row.names = F)
write.csv(alexq3, "alexq3_mancur_ske.csv", row.names = F)
write.csv(alexq4, "alexq4_mancur_ske.csv", row.names = F)

alexq1=read.csv("alexq1_mancur_ske.csv", header=T, stringsAsFactors = F)
alexq2=read.csv("alexq2_mancur_ske.csv", header=T, stringsAsFactors = F)

length(unique(colnames(alexq1)))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### INVESTIGATING MULTIPOINT BOUND FUNCTION                                                                 ####
#do these markers show exactly the same segregation pattern?

which(alexq1$AX.94771499 != alexq1$AX.94523513)

#see onebq1_comb from prelim_map_comp.R
length(which(alexq2$AX.94428799 != alexq2$AX.94970814)) #multipoint has put these markers in the same bin even though they don't have the same pattern of recombination
#is this due to NoCalls?
alexq2$AX.94428799
alexq2$AX.94970814
#yes, the only differences between these two markers are NoCalls

alexq2$AX.94647201 != alexq2$AX.94480226


alexq2$AX.95004030 != alexq2$AX.94996993

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

