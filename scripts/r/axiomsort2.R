#setwd("/Users/alexcoulton/Google Drive/R/scripts/rotation1scripts")
setwd("C:/Users/ac14037/Google Drive/R/scripts/rotation1scripts")
#rm(list=ls()) #clear the environment
#WINDOWS DATA
axiomdatalocation="C:/Users/ac14037/Google Drive/University/PhD/Rotation 1/Some_data_to_get_you_started/Apogee_Paragon_F5_export"
#MAC DATA
#axiomdatalocation="/Users/alexcoulton/Google Drive/University/PhD/Rotation 1/Some_data_to_get_you_started/Apogee_Paragon_F5_export"
#list all axiom files in directory
files = list.files(axiomdatalocation, full.names=T, recursive=F, pattern=".txt")
files_j_t = list.files(axiomdatalocation)
# print(files) #verify files variables are set correctly
# print(files_j_t)

#import text file to data frame, ignore first 3 lines of file
######LOAD UP FIRST FILE##################
firstfile = read.table(files[1], skip = 3)
firstfile=t(firstfile)
colnames(firstfile)=firstfile[1,] #make the column names in firstfile the same as the first row
firstfile=cbind(id=files_j_t[1],firstfile) #add a new column to firstfile at the start of the data frame
d=firstfile
d[]=NA
firstfile=rbind(d,firstfile)
firstfile=firstfile[-3,]

files=files[-1]
files_j_t=files_j_t[-1]
########################

cbindcount=1 #initiate counter for files_j_t (just title) - for some reason for loop variable wasn't working
#LOOP THROUGH FILES AND ADD DATA TO FIRSTFILE DATAFRAME
for (i in files){
    secondfile=read.table(i,skip=3)
    secondfile=t(secondfile)
    secondfile=cbind(id=files_j_t[cbindcount],secondfile)
    cbindcount=cbindcount+1
    ###NEED TO IMPLEMENT CHECK TO SEE IF MARKERS MATCH UP IN EVERY AXIOM FILE###
    firstcolnames=colnames(firstfile)
    secondcolnames=secondfile[1,]
    names(secondcolnames) = NULL
    firstcolnames[1] = NA
    if(isTRUE(all.equal(firstcolnames[c(2:length(firstcolnames))], secondcolnames[c(2:length(secondcolnames))]))){
        markersame=TRUE
    } else {
        print(paste(files[i], "not the same marker"))
    }
    if(markersame == TRUE){
        secondfile=secondfile[-1,]
        firstfile=rbind(firstfile, secondfile)
    }
}

#assign all markers to chromosome 1
head(firstfile[1,], 10)
firstfile[1,]=1

#rename first row
rownames(firstfile)[1]="chromosome"

#add unique id number to each sample (1-380)
newid=firstfile[,1]
length(newid)
newid[]=seq(1,380,1)
newid
firstfile=cbind(newid, firstfile)

#apply correction to id sequence (remove 1 & 2 from first 2 rows)
firstfile[3:380,1]=seq(1,378,1)
firstfile[2,1]=NA

#replace some values to conform with rqtl
firstfile[firstfile=="AB"]="H"
firstfile[firstfile=="NoCall"]="-"
firstfile[firstfile=="AA"]="A"
firstfile[firstfile=="BB"]="B"

#remove useless second row and useless cells
firstfile=firstfile[-2,]
firstfile[1,1]=""
firstfile[1,2]=""

write.csv(firstfile, file="alldata.csv", row.names=F)

##############NO LONGER USING LAPPLY - WOULDN'T WORK IN THIS SITUATION###########
#need to replace "function(x)" (and function(t)) with a function of your choice --- SPECIFY THE FUNCTION!!
# re: above, do I? how do functions work?
# lapply(files, function(x) {
#     t = read.table(x, skip = 3) # load file
#     # apply function
#     out = function(t)
#     # write to file
#     write.table(out, "/Users/alexcoulton/Google Drive/R/scripts/rotation1scripts/test.txt", sep="\t", quote=F, row.names=F, col.names=T, append=T)
# })
###########################
