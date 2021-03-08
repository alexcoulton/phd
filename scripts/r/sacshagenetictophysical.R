source("functions.R") #import genetic map function

allenmap=read.csv("Allen Characterization of Wheat breeders array suppl. docs/Allen genetic map AxP.csv",header=T,stringsAsFactors = F)
allen1a=filter(allenmap, chr=="1A")
allen1a$marker

listofchromosomes=c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D")
listofchr2=listofchromosomes[16:21]
listofchr2
for(i in listofchromosomes){
    chromosometolookat=i
    png(filename=paste("SaschaPhysicalChromosomeCoverage/chr",chromosometolookat,".png",sep=""))
    hist(genetictophysical(chromosometolookat, filter(allenmap, chr==chromosometolookat)$marker), main=chromosometolookat)
    dev.off()
}



hist(genetictophysical("6A", filter(allenmap, chr=="6A")$marker), main="6A")


#doing some checks on my own map. is 3D really 3D?

bigmap=read.csv("man_cur_full_asmap_w_chr.csv", header=T, stringsAsFactors = F)

length(as.numeric(na.omit(genetictophysical("3D", filter(bigmap, lg=="LG15")$marker, "."))))
genetictophysical("2D", filter(bigmap, lg=="LG20")$marker, ".")
length(as.numeric(na.omit(genetictophysical("3B", filter(bigmap, lg=="LG15")$marker, "."))))
