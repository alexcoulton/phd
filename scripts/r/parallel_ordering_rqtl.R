#!/home/ac14037/R/bin/Rscript
#setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#setwd("~/Google Drive/R/scripts/rotation1scripts/")
#ARGUMENTS AS FOLLOWS --args -starting_rf_parameter_iteration -ending_rf_parameter_iteration -starting_lod_parameter_iteration -ending_lod_parameter_iteration
setwd("/home/ac14037/wheat/rotation1scripts")
library(qtl)
args=commandArgs(TRUE)
print(args)

quadrant="1"
if(dir.exists(paste("chromosomes_preordering_quad",quadrant,sep=""))==F){
    dir.create(paste("chromosomes_preordering_quad",quadrant,sep=""))
}
if(dir.exists(paste("chromosomes_postordering_quad",quadrant,sep=""))==F){
    dir.create(paste("chromosomes_postordering_quad",quadrant,sep=""))
}

#read in the chromosome data into new cross object. NB: specify the folder path in parallel
chrdata=read.cross("csv", "./", args[1], estimate.map=F, crosstype="f2")

#extract the chromosome number from the filename for saving purposes later
g=regexpr("chr_[[:digit:]]{1,2}", args[1])
gl=as.numeric(attributes(g)[1])
argchr=substring(args[1],g,(g-1+gl))
g=regexpr("[[:digit:]]{1,2}", argchr)
gl=as.numeric(attributes(g)[1])
chrnum=as.numeric(substring(argchr,g,(g-1+gl)))

print("Processing chromosome number:")
print(chrnum)

#order the markers within this particular chromosome
chrdata=orderMarkers(chrdata, chr=chrnum)

#write the order chromosome information to file
write.cross(chrdata, format="csv", filestem=paste("chromosomes_postordering_quad",quadrant,"/chr_",chrnum,"_post-ord",sep=""), chr=chrnum)
chrcsv=read.csv(paste("chromosomes_postordering_quad",quadrant,"/chr_",chrnum,"_post-ord.csv",sep=""),header=T,stringsAsFactors = F)
chrcsv[chrcsv=="AA"]="A"
chrcsv[chrcsv=="BB"]="B"
chrcsv[chrcsv=="AB"]="H"
chrcsv[1,1]=""
chrcsv[2,1]=""
write.csv(chrcsv, paste("chromosomes_postordering_quad",quadrant,"/chr_",chrnum,"_post-ord.csv", sep=""), row.names=F)

