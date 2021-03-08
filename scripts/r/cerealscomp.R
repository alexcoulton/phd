setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#setwd("~/Google Drive/R/scripts/rotation1scripts/")
library(qtl)
load("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts/axiomrqtl5.RData") #load rqtl data
cereals=read.csv("C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/cerealsdb/axiom820_2.csv", header=T)
cereals$affycode=gsub('-','.',cereals$affycode) #fix affycode field so matches with flipped_data

#USING MULTIPOINT DATA (Y)? ELSE GO THROUGH ITERATION PROCEDURE
use_mp="Y"

if(use_mp =="Y"){
#begin section for processing multipoint output
####################################################
#aka. "Result file_Sk&Ext.txt"
mult_dat=read.table("multipointq1mapv2.txt",sep="\n", stringsAsFactors = F)

mult_dat2=cbind(mult_dat,mult_dat,mult_dat,mult_dat)

length(strsplit(mult_dat[10,], split="\\s+"))

#put data into multiple columns
for(i in 1:nrow(mult_dat)){
    b=strsplit(mult_dat[i,], split="\\s+")
        for(l in 1:length(b[[1]])){
            mult_dat2[i,l]=b[[1]][l]
    }
}

#remove blank rows
rowtodel=vector()
for(i in 1:nrow(mult_dat)){
    if(((grepl("", mult_dat2[i,1])==T)==T & (grepl("\\s+",mult_dat2[i,2])==T)==T & (grepl("\\s+",mult_dat2[i,3])==T)==T & (grepl("\\s+",mult_dat2[i,4])==T)==T)==T){
        rowtodel=c(rowtodel, i)
    }
}
rowtodel
mult_dat3=mult_dat2[-c(rowtodel),]

#populate LG column
for(i in 2:nrow(mult_dat3)){
    if((grepl("LG\\d{1,2}",mult_dat3[(i-1),1]) == T) & mult_dat3[i,1]==""){
        mult_dat3[i,1]=mult_dat3[(i-1),1]
    }
}

rowstodel=vector()
for(i in 1:nrow(mult_dat3)){
    if((mult_dat3[i,1] == mult_dat3[i,2]) & (mult_dat3[i,1] == mult_dat3[i,3])){
        rowstodel=c(rowstodel,i)
    }
}
mult_dat3=mult_dat3[-c(rowstodel),]
######################################################
#GRAB LIST OF MARKERS IN MULTIPOINT LINKAGE GROUP mlg
mlg="LG4"
listofmlg=unique(mult_dat3[,1])
#COMPARE TO CEREALSDB DATA
for(i in listofmlg){
bn=mult_dat3[mult_dat3==i,][,4]
t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group
print(table(droplevels(t$consensus)))
}
} else {


quadvec=c("2","3","4")
#quadrant="1"

for(quadrant in quadvec){

zseq=seq(0.1,0.3,0.01)
lodseq=c(12,14,16,18,20)
lodseq

axiomdata = read.cross("csv", "./", paste("alexq",quadrant,".csv", sep=""), estimate.map=F) #NOTE THE THIRD ARGUMENT — CHANGE BETWEEN flipped_data2.csv AND flipped_data2_pruned.csv ACCORDINGLY
axiomdata_untouched = axiomdata
#axiomdata=convert2riself(axiomdata) #convert cross to RIL? assumes 1:1 segregation, but removes heterozygotes (i think)
summary(axiomdata)

for(z in zseq){
for(L in lodseq){
rfthres=z
lodthres=L
#rfthres=0.18
#lodthres=16


#form linkage groups
lg <- formLinkageGroups(axiomdata, max.rf=rfthres, min.lod=lodthres)
table(lg[,2])
length(table(lg[,2]))

print("forming linkage groups")
axiomdata2=axiomdata
axiomdata2=formLinkageGroups(axiomdata2, max.rf=rfthres, min.lod=lodthres, reorgMarkers=TRUE)
print("done forming linkage groups")

#pull genetic map data for first chromosome as table


###FUNCTIONALIZE THIS CODE###
#####
#POPULATE MYMAT
csvoutput=function(listofmarkers,dformat){
    #supply the function with a character vector of marker names as well as the data source (either rQTL or multipoint)
    if(dformat=="rqtl"){
        numchromosomes=(length(table(lg[,2]))-10)
    } else if(dformat="multipoint"){
        numchromosomes=length(listofmlg)
    }
    
    #SETUP DATAFRAME FOR POPULATION WITH CHROMOSOME COUNT DATA FOR ITERATIONS ON LINKAGE GROUP, LOD AND RF
    #####
    mymat=matrix(data=NA,nrow=2000,ncol=52)
    mymat=data.frame(mymat)
    mymat[1,1]="linkage group"
    cerealstitles=data.frame(table(cereals$consensus))
    cerealstitles <- rapply(cerealstitles, as.character, classes="factor", how="replace")
    cerealstitles[1,1]="BLANK"
    cerealstitles[,1]
    length(cerealstitles[,1])
    mymat[1,2:(length(cerealstitles[,1])+1)]=cerealstitles[,1]
    mymat=mymat[,-c(ncol(mymat),(ncol(mymat)-1))]
    colnames(mymat)=mymat[1,]
    mymat[1,]=NA
    mymat[is.na(mymat)]=""
    mymat=sapply(mymat, as.numeric)
    

    for(p in 1:(numchromosomes)){
        if(dformat=="rqtl"){
            print("populating mymat")
            mymat[p,1]=p #populate linkage group column
            
            b=pull.map(axiomdata2, chr=p, as.table=T)
            bn=rownames(b)
            
            t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group
            table(droplevels(t$consensus))
        } else if(dformat=="multipoint"){
                bn=mult_dat3[mult_dat3==listofmlg[p],][,4]
                t=cereals[c(which(cereals[,1] %in% bn)),] #find markers in cerealsdb which are also in linkage group
        }
        
        #SAVE CHROMOSOME COUNT BARPLOTS FOR EACH LINKAGE GROUP
        print("creating directories")
        if(dir.exists(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, sep=""))==F){
            dir.create(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat, sep=""))
        }
        if(dir.exists(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat,"/rf_", rfthres, "_lod_", lodthres, sep=""))==F){
            dir.create(paste("chromosome_count_plots_quad",quadrant, "_form_", dformat,"/rf_", rfthres, "_lod_", lodthres, sep=""))
        }
        png(filename=paste("chromosome_count_plots_quad",quadrant, "_form_", dformat,"/rf_", rfthres, "_lod_", lodthres, "/","chromo_", p, ".png", sep=""))
        barplot(table(droplevels(t$consensus)))
        dev.off()
        
        #make dataframe containing table output
        mydat=data.frame(table(droplevels(t$consensus)), stringsAsFactors = F) #tabulate counts of consensus. NB: droplevels() removes zero counts
        mydat <- rapply(mydat, as.character, classes="factor", how="replace")
        
        if(mydat[1,1]==""){ mydat[1,1]="BLANK"}
        #populate mymat with chromosome frequency data for linkage group p
        for(x in 1:nrow(mydat)){
            coord=which(colnames(mymat)==mydat[x,1])
            mymat[p,coord]=mydat[x,2]
        }
    }
    
    mymat=cbind(mymat[,1],mymat)
    colnames(mymat)[1]="linkage group"
    colnames(mymat)[2]="most frequent"
    
    
    for(i in 1:length(table(lg[,2]))){
        gogo=data.frame(na.exclude(mymat[i,4:(ncol(mymat)-1)]))
        biggest=max(as.numeric(as.character(gogo[,1])))
        mymat[i,2]=paste(rownames(data.frame(which(mymat[i,]==biggest))), collapse=" ")
    }
    
    mymat=cbind(mymat[,1],mymat)
    colnames(mymat)[1]="linkage group"
    colnames(mymat)[2]="totalinlg"
    
    for(i in 1:length(table(lg[,2]))){
        gogo=data.frame(na.exclude(mymat[i,4:(ncol(mymat)-1)]))
        total=sum(as.numeric(as.character(gogo[,1])))
        mymat[i,2]=total
    }
    
    if(dformat="rqtl"){
        print("generating rfplot")
        png(filename=paste("chromosome_count_plots_quad",quadrant,"/rf_", rfthres, "_lod_", lodthres, "/","rf_", rfthres,"_lod_",lodthres, ".png", sep=""))
        plotRF(axiomdata2, alternate.chrid=TRUE)
        dev.off()
    }
    
    if(dir.exists(paste("linkage_grouping_it_quad",quadrant, "_form_", dformat, sep=""))==F){
        dir.create(paste("linkage_grouping_it_quad",quadrant, "_form_", dformat, sep=""))
    }
    write.csv(mymat,paste("linkage_grouping_it_quad",quadrant, "_form_", dformat,"/lgroups_rf_",rfthres,"_lod_",lodthres,".csv",sep=""), row.names=F)
    print(paste("end of iteration",z,"...", (length(z)-which(z==zseq)),"to go", sep=" "))
    }
}
}
}
}