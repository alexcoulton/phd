#fasta functions

ReadFasta<-function(file) { #from https://stackoverflow.com/questions/26843995/r-read-fasta-files-into-data-frame-using-base-r-not-biostrings-and-the-like
    #authors note: "... be warned: the function does not currently handle, e.g., special characters, which are rather commonplace in sequence files (in context such as 5' or #5 rRNA)"
    # Read the file line by line
    fasta<-readLines(file)
    # Identify header lines
    ind<-grep(">", fasta)
    # Identify the sequence lines
    s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
    # Process sequence lines
    seqs<-rep(NA, length(ind))
    for(i in 1:length(ind)) {
        seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
    }
    # Create a data frame 
    DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
    # Return the data frame as a result object from the function
    return(DF)
}


writefasta = function(fa.df, path){
    #fa.df - dataframe containing fasta information, headers in column 1, sequences in column 2
    #path - path to save file
    fa.df[] = lapply(fa.df, as.character)
    g = character()
    
    #check if headers include header symbol ">", if not add one
    if(length(grep(">", fa.df[, 1])) == 0){
        fa.df[, 1] = paste(">", fa.df[, 1], sep = "")
    }
    
    # bob = colnames(fa.df)
    # 
    # g = c(rbind(fa.df[[bob[[1]]]], fa.df[[bob[[2]]]]))
    
    
    for(i in 1:nrow(fa.df)){
        g = c(g, fa.df[i, 1], fa.df[i, 2])
        print(i)
    }
    
    fileconn = file(path)
    writeLines(g, fileconn)
    close(fileconn)
}
