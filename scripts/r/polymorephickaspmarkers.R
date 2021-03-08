#polymorphic markers script KASP apogee paragon

###########################################################################################################
#####SAVED WORKSPACE IMAGE TO C:\Users\ac14037\Google Drive\R\workspaces\rotation1workspace ###############
###########################################################################################################


setwd("C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/")


codeconversion=read.csv("AX vs BS codes.csv", header=T)

apogee1=read.csv("cerealsdb/apogee1axiomgenotypedata.csv", header=T)
apogee2=read.csv("cerealsdb/apogee2axiomgenotypedata.csv", header=T)
apogee4=read.csv("cerealsdb/apogee4axiomgenotypedata.csv", header=T)

apocomb=cbind(apogee1, apogee2)
apocomb=cbind(apocomb, apogee4)


b = apogee1$Matrix_value == apogee2$Matrix_value
c = apogee1$Matrix_value == apogee4$Matrix_value
d = b==c
#consensus between apogee genotypes
apocombcons=apocomb[which(d),]
apocombcons=apocombcons[,1:3]
nrow(apocombcons)
nrow(codeconversion)
#select rows of codeconversion for which we have genotyping data
kasp=codeconversion[which(codeconversion$ï..Affymetrix.Code %in% apocombcons$Probe_row),]
kasp2=apocombcons[which(apocombcons$Probe_row %in% kasp$ï..Affymetrix.Code),]

kasp3=cbind(kasp, kasp2)
kasp3=kasp3[,-5]
all(as.character(kasp3[,1])==as.character(kasp3[,5])) #are all affy probe codes the same?

#LOAD IN PARAGON DATA
pa=read.csv("cerealsdb/kasp Haplotypes_paragon.csv", header=T)
head(pa)

pa2=pa[which(pa$SNP.name %in% kasp3$Bristol.SNP.Code)]

pa=pa[!pa$NIAB_Paragon=="NS",]
head(pa)
pa=pa[!pa$Paragon=="NS",]
pa=pa[!pa$Paragon=="-",]
pa=pa[!pa$NIAB_Paragon=="-",]
pa$NIAB_Paragon=as.character(pa$NIAB_Paragon)
pa$Paragon=as.character(pa$Paragon)
pa2=pa[which(pa$NIAB_Paragon==pa$Paragon),]
length(which(pa2$NIAB_Paragon==pa2$Paragon))
length(pa2$Paragon)
pa2=pa2[pa2$Co.dominant.Dominant=="Co-dominant",]

kasp4=kasp3[which(kasp3$Bristol.SNP.Code %in% pa2$SNP.name),]
pa3=pa2[which(pa2$SNP.name %in% kasp4$Bristol.SNP.Code),]
kasp4[order]
kasp4=kasp4[order(kasp4$Bristol.SNP.Code),]
pa3=pa3[order(pa3$SNP.name),]
?order
nrow(kasp4)
nrow(pa3)

#are all of the values the same and in the same order? (yes)
all(as.character(kasp4$Bristol.SNP.Code)==as.character(pa3$SNP.name))


#checking to see if order function works
b=data.frame(c(1,2,3))
y=data.frame(c("c","b","a"))
ta=cbind(b,y)
ta=ta[order(ta$c..c....b....a..),]
#it does

#need to cross reference affyids of Apogee with probe SNP data
snp=read.csv("affy_820K_array_data_probe_set.csv",header=T)
snp2=snp[which(snp$Bristol.SNP.Code %in% kasp4$Bristol.SNP.Code),]
snp2=snp2[order(snp2$Bristol.SNP.Code),]
all(as.character(kasp4$Bristol.SNP.Code)==as.character(snp2$Bristol.SNP.Code))
#combine dataframes
kasp5=data.frame(kasp4,snp2$Bristol.SNP.Code, snp2$Sequence..including.SNP.ambiguity.code.)
#convert data type to character
kasp5$snp2.Sequence..including.SNP.ambiguity.code. = as.character(kasp5$snp2.Sequence..including.SNP.ambiguity.code.)

#find all amb. codes used
table(reg3)

nrow(kasp5)

regexpr("\\[.\\]", kasp5$snp2.Sequence..including.SNP.ambiguity.code.)
?regexpr

#some probes are not in "[snp amb. code]" format... 
#need to clean data
#specify data cleaning function: extract data
#w/ no square brackets
repp=function(x){
#can't get replacement function to work... #EDIT: FIXED
#see the following: functions make a local copy of the variable
#from the global environment before making changes
#https://stackoverflow.com/questions/3969852/update-data-frame-via-function-doesnt-work
#NOTE: use <<- to make assignments to global variables from inside functions
f=regexpr(x, kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi])
f
f1=which(f!=-1)
f[f1[2]]
f1

substr(kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi][f1[1]],f[f1[1]],f[f1[1]])
length(f1)

#note global assignment (<<-)
for(i in 1:length(f1)){
    kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi][f1[i]]<<-substr(kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi][f1[i]],f[f1[i]],f[f1[i]])
}

}
#isolate values using function
repp("M");repp("R");repp("S");repp("W");repp("Y");repp("K")

#check it worked
kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi]

#do some gsub: add square brackets [] to values missing them
#so subsequent replacement function works
#CHANGE VALUE IN GSUB LINE TO REPLACE LETTERS YOU NEED
kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi] = gsub("W","\\[W\\]", kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi])
#check it worked
kasp5$snp2.Sequence..including.SNP.ambiguity.code.[whi]

#info on regex matching in R https://www.regular-expressions.info/rlanguage.html
#http://www.rexegg.com/regex-quickstart.html
#https://emacs.stackexchange.com/questions/19645/escaping-brackets-in-regexp
#extracting SNP information from probe sequences
reg=regexpr("\\[.\\]",kasp5$snp2.Sequence..including.SNP.ambiguity.code.)
reg2=regmatches(kasp5$snp2.Sequence..including.SNP.ambiguity.code.,reg)
#using substr to select the 2nd character within each string
reg3=substr(reg2, 2, 2)

#reg3 not same length as no. rows in kasp5... troubleshooting
length(reg)
length(reg3)
nrow(kasp5)

kasp5=data.frame(kasp5, reg3) #ERROR: not same length 

#load in paragon axiom genotype data
paraax=read.csv("cerealsdb/paragonaxiomgenotypedata.csv", header=T)
head(paraax)

paraax2=paraax[which(paraax$Probe_row %in% kasp5$Probe_row),]
paraax2=paraax2[order(paraax2$Probe_row),]
head(paraax2)
?match
#see match function
#orders the rows in one column according to their order
#in another column
paraax2=paraax2[match(kasp5$Probe_row, paraax2$Probe_row),]
#make vector of polymorphic markers
poly=paraax2$Matrix_value == kasp5$Matrix_value
table(poly)
?which
kasp6=data.frame(kasp5,paraax2)
#filter kasp6 to only include polymorphic loci
kasp6=kasp6[-which(kasp6$Matrix_value == kasp6$Matrix_value.1),]
kasp7=kasp6
kasp7$Matrix_value=as.character(kasp7$Matrix_value)
kasp7=kasp7[!kasp7$Matrix_value=="AB",]
kasp7=kasp7[!kasp7$Matrix_value=="NoCall",]
kasp7=kasp7[!kasp7$Matrix_value.1=="AB",]
kasp7=kasp7[!kasp7$Matrix_value.1=="NoCall",]
kasp_clean=data.frame(kasp7$Bristol.SNP.Code, kasp7$Var_col, kasp7$Matrix_value, kasp7$Var_col.1, kasp7$Matrix_value.1)

write.csv(kasp_clean, "para_apo_kasp.csv")

#add in column explaining SNP ambiguity codes
for(i in 1:nrow(kasp7)){
if(kasp7$reg3[i] == "Y"){
kasp7$ambi_code[i] = "C or T"
} else if(kasp7$reg3[i] == "R"){ 
    kasp7$ambi_code[i] = "A or G"
} else if(kasp7$reg3[i] == "W"){ 
    kasp7$ambi_code[i] = "A or T"
} else if(kasp7$reg3[i] == "S"){ 
    kasp7$ambi_code[i] = "G or C"
} else if(kasp7$reg3[i] == "K"){ 
    kasp7$ambi_code[i] = "T or G"
} else if(kasp7$reg3[i] == "M"){ 
    kasp7$ambi_code[i] = "C or A"} 
}

