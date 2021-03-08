#avalon x cadenza initial data preperation

setwd("C:/Users/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
library(dplyr)
library(tibble)


#     ____________________________________________________________________________
#     INITIAL PROCESSING                                                                                                            ####

s.x.s.array = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/shamrock.x.shango/shamrock.x.shango.aas.output.polyhighres.txt")


#     ____________________________________________________________________________
#     INTERFACE WITH MYSQL                                                                                                        ####

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

rs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Shamrock%'") #grab Avalon genotyping data for 35k array
rs2 = fetch(rs, n=-1)
rs2 = rs2[rs2$Var_col == "Shamrock", ]


#grab cadenza info
ca = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Shango%'") #grab Avalon genotyping data for 35k array
ca2 = fetch(ca, n=-1)

rs2 = replace.aas.genotype.format.w.rqtl.format(rs2)
ca2 = replace.aas.genotype.format.w.rqtl.format(ca2)


#     ____________________________________________________________________________
#     FURTHER PROCESSING                                                                                                            ####

merged.array.df = s.x.s.array

merged.array.df = convert.to.character.data.frame(merged.array.df)

merged.array.df = add_row(merged.array.df, .before = 1)
merged.array.df = add_row(merged.array.df, .before = 1)

g = rs2$Matrix_value[match(colnames(merged.array.df)[3:ncol(merged.array.df)], rs2$Probe_row)] #grab genotyping data for Shamrock

merged.array.df[1, 3:ncol(merged.array.df)] = g
merged.array.df[1, 2] = "Shamrock"

#prepare cadenza genotyping data
g = ca2$Matrix_value[match(colnames(merged.array.df)[3:ncol(merged.array.df)], ca2$Probe_row)] #grab genotyping data for Shango

merged.array.df[2, 3:ncol(merged.array.df)] = g
merged.array.df[2, 2] = "Shango"

merged.array.df = rbind(merged.array.df[3, ], merged.array.df[1:2, ], merged.array.df[4:nrow(merged.array.df), ])

merged.array.df[, 1] = c("", 1:(nrow(merged.array.df)-1))

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


merged.array.df2 = cleanup(merged.array.df, 2:3)

merged.array.df3 = assign.parental.genotypes(merged.array.df2, 2:3)

write.csv(merged.array.df3, "rotation1scripts_v4/processed_data/genotypes/shamrock.x.shango/shamrock.x.shango.flipped.csv", row.names = F)



