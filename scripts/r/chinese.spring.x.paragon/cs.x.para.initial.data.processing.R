#chinese spring x paragon f7 sacha data preparation

source("scripts/r/functions.R")
library(dplyr)
library(tibble)


#     ____________________________________________________________________________
#     INITIAL PROCESSING                                                                                                            ####


chinese.spring.geno = read.csv("original_data/genotype.data/chinese.spring.x.paragon/chinese.spring.cerealsdb.35k.array.genotype.csv", 
                                                             header = T, stringsAsFactors = F)

chinese.spring.geno = cbind(chinese.spring.geno[which(chinese.spring.geno$Var_col == "192_Chinese_Spring"), ], chinese.spring.geno[which(chinese.spring.geno$Var_col == "Chinese_Spring"), ])

chinese.spring.geno.concordant = chinese.spring.geno[which(chinese.spring.geno[, 3] == chinese.spring.geno[, 6]), ]

cs.x.p.geno.data = convert.aas.to.rqtl("original_data/axiom analysis suite output/CSxP.F7.Allen.PolyHighRes.txt", c("Paragon", "Chinese_Spring"))

#grab only markers that are concordant between chinese spring parental genotypes
cs.x.p.geno.data = cs.x.p.geno.data[, c(1, 2, which(colnames(cs.x.p.geno.data)[1:ncol(cs.x.p.geno.data)] %in% chinese.spring.geno.concordant[, 2]))]

v(cs.x.p.geno.data)

replace.aas.genotype.format.w.rqtl.format = function(dataframe){
    dataframe[dataframe == "AA"] = "A"
    dataframe[dataframe == "BB"] = "B"
    dataframe[dataframe == "AB"] = "H"
    dataframe[dataframe == "BA"] = "H"
    dataframe[dataframe == "NoCall"] = "-"
    return(dataframe)
}

cs.x.p.geno.data = replace.aas.genotype.format.w.rqtl.format(cs.x.p.geno.data)

grab.columns.in.which.parents.are.concordant = function(dataframe, parent){
    #supply dataframe in rqtl format and name of parent variety (only processes one parent at a time)
    parent.row.coords = grep(parent, dataframe[, 2])
    parent.df = dataframe[parent.row.coords, ]
    parent.df = as.data.frame(parent.df)
    parent.df[] = lapply(parent.df, as.character)
    
    g = lapply(parent.df, function(x){
        b = x[[1]] == x[[2]]
        c = x[[2]] == x[[3]]
        g = all(c(c, b))
        return(g)
    })
    
    which(unlist(g))
    
}

cs.x.p.geno.data2 = cs.x.p.geno.data[, c(1, 2, as.numeric(grab.columns.in.which.parents.are.concordant(cs.x.p.geno.data, "Paragon")))]

convert.cerealsdb.format.to.rqtl = function(df){
    g = df[, 1:3]
    g = t(g)
    colnames(g) = g[2, ]
    
    g = as.data.frame(g)
    g = add_column(g, probeset_id = g[1, 1], .before = colnames(g)[1])
    g = add_column(g, V1 = 1, .before = colnames(g)[1])
    g = g[3, ]
    g = convert.to.character.data.frame(g)
    g = replace.aas.genotype.format.w.rqtl.format(g)
    return(g)
}


#add chinese spring parent data
colnames(cs.x.p.geno.data2)[1] = "V1"

chinese.spring.geno.concordant.rqtl.format = convert.cerealsdb.format.to.rqtl(chinese.spring.geno.concordant)
chinese.spring.geno.concordant.rqtl.format = chinese.spring.geno.concordant.rqtl.format[, which(colnames(chinese.spring.geno.concordant.rqtl.format) %in% colnames(cs.x.p.geno.data2))]

chinese.spring.geno.concordant.rqtl.format = convert.to.character.data.frame(chinese.spring.geno.concordant.rqtl.format)
cs.x.p.geno.data2 = convert.to.character.data.frame(cs.x.p.geno.data2)


combined.geno.data = rbind(cs.x.p.geno.data2[1, ], chinese.spring.geno.concordant.rqtl.format, cs.x.p.geno.data2[2:nrow(cs.x.p.geno.data2), ])


#     ____________________________________________________________________________
#     DEFINE FUNCTIONS                                                                                                                ####


#remove columns in which parents are the same

rm.interparental.homozygosities = function(df, parent.rows){
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    
    if(length((which(df[parent.rows[[1]], 3:ncol(df)] == df[parent.rows[[2]], 3:ncol(df)]) + 2)) != 0) {
        df = df[, -(which(df[parent.rows[[1]], 3:ncol(df)] == df[parent.rows[[2]], 3:ncol(df)]) + 2)]    
    }
    return(df)
}

#remove parental heterozygosities
rm.parental.het = function(df, parent.rows){
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    if(length(which(df[parent.rows[1], ] == "H")) != 0) g = df[, -which(df[parent.rows[1], ] == "H")]
    if(length(which(g[parent.rows[2], ] == "H")) != 0) g = g[, -which(g[parent.rows[2], ] == "H")]
    return(g)
}


#remove parental nocalls
rm.parental.nocall = function(df, parent.rows){
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    if(length(which(df[parent.rows[1], ] == "-")) != 0) g = df[, -which(df[parent.rows[1], ] == "-")]
    if(length(which(g[parent.rows[2], ] == "-")) != 0) g = g[, -which(g[parent.rows[2], ] == "-")]
    return(g)
}

cleanup = function(df, parent.rows){
    df = rm.interparental.homozygosities(df, parent.rows)
    df = rm.parental.het(df, parent.rows)
    df = rm.parental.nocall(df, parent.rows)
    return(df)
}

assign.parental.genotypes = function(df, parent.rows){
    #performs "flipping algorithm" to assign parental genotypes to data. 
    #df - the dataframe in rqtl format for which to remove columns in which parents are heterozygous
    #parent.rows - a numeric vector of length 2 providing the row coordinates of each parent in the cross.
    
    #make new df w/ only parents (2 rows total)
    pc2 = df[c(parent.rows[[1]], parent.rows[[2]]), ]
    pc2 = as.data.frame(pc2)
    pc2 = convert.to.character.data.frame(pc2)
    
    df = as.data.frame(df)
    
    #make list of columns to flip
    n = which(lapply(pc2, function(x){
        if(x[1] == "A" & x[2] == "B"){
            return(T)
        } else return(F)
    }) == T)
    
    #flip columns in main dataframe
    df[, n] = lapply(df[, n], function(x){
        g = x
        
        g[g == "A"] = "X"
        g[g == "B"] = "A"
        g[g == "X"] = "B"
        
        return(g)
    })

    
    return(df)
}

remove.excess.parents = function(df, parent.name){
    #removes all but one of the parent individuals for one of the parents in the cross (only operates on one parent at a time)
    #df - df in rqtl format
    #parent.name - character string of the parent to perform operation on
    to.remove = grep(parent.name, df$probeset_id)
    df = df[-to.remove[2:length(to.remove)], ]
    return(df)
}


#     ____________________________________________________________________________
#     PROCESS DATA                                                                                                                        ####


combined.geno.data3 = cleanup(combined.geno.data, 2:3)

combined.geno.data4 = assign.parental.genotypes(combined.geno.data3, 2:3)

combined.geno.data4 = remove.excess.parents(combined.geno.data4, "Paragon")

write.csv(combined.geno.data4, "processed_data/genotypes/chinese.spring.x.paragon/cs.x.para.flipped.csv", row.names = F)



