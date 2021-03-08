remove.hits.with.same.best.e.value = function(blastdf){
    #removes query-subject pairs in which the two best hits for a single query both have the same e-value
    #blastdf - a dataframe containing a tabular blast output file
    
    hits.to.keep = newdf(colnames(blastdf), no.rows = T)
    blastdf[, 11] = as.numeric(blastdf[, 11])
    counter1 = make.counter()
    un.hits = unique(blastdf[, 1])
    blastdf = as.data.table(blastdf)
    setkey(blastdf, qseqid)
    lapply(un.hits, function(x){
        # temp.df = blastdf[blastdf[, 1] == x, ]
        # temp.df = filter(blastdf, qseqid == x)
        temp.df = blastdf[.(x), nomatch = 0L]
        best.e.value = min(temp.df[, 11])
        if(nrow(filter(temp.df, evalue == best.e.value) == 1)) hits.to.keep <<- rbind(hits.to.keep, temp.df)
        # if(length(which(temp.df[, 11] == best.e.value)) == 1) hits.to.keep <<- rbind(hits.to.keep, temp.df)
        print(counter1())
    })
    
    
    # for(i in unique(blastdf[, 1])){
    #     temp.df = blastdf[blastdf[, 1] == i, ]
    #     best.e.value = min(as.numeric(temp.df[, 11]))
    #     if(length(which(temp.df[, 11] == best.e.value)) == 1) hits.to.keep = rbind(hits.to.keep, temp.df)
    # }
    return(hits.to.keep)
}