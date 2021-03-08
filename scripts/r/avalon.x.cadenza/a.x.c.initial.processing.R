#avalon x cadenza initial data preperation

setwd("C:/Users/ac14037/project.phd.main/")
source("rotation1scripts_v4/scripts/r/functions.R")
library(dplyr)
library(tibble)


#     ____________________________________________________________________________
#     INITIAL PROCESSING                                                                                                            ####

a.x.c.array.a = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza/avalon.x.cadenza.aas.output.array.a.polyhighres.txt", c("Avalon", "Cadenza"))

a.x.c.array.b = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza/avalon.x.cadenza.aas.output.array.b.polyhighres.txt", c("Avalon", "Cadenza"))

#find sample entries in array a that correspond to crosses between avalon and cadenza
matches1 = regmatches(a.x.c.array.a[, 2], regexec("AxC.*\\d*\\.", as.character(a.x.c.array.a[, 2])))

#remove samples that are not AxC cross from dataset (irrelevant samples)
entries.to.remove = unname(which(unlist(lapply(matches1, length)) == 0))
entries.to.remove = entries.to.remove[-(1:3)]

a.x.c.array.a = a.x.c.array.a[-entries.to.remove, ]

#do the same for the second array
matches2 = regmatches(a.x.c.array.b[, 2], regexec("AxC.*\\d*\\.", as.character(a.x.c.array.b[, 2])))

entries.to.remove = unname(which(unlist(lapply(matches2, length)) == 0))
entries.to.remove = entries.to.remove[-(1)]

a.x.c.array.b = a.x.c.array.b[-entries.to.remove, ]

#grab matches again for array a
matches1 = regmatches(a.x.c.array.a[, 2], regexec("AxC.*\\d*\\.", as.character(a.x.c.array.a[, 2])))

matches1.1 = regmatches(as.character(matches1), regexpr("\\d+", as.character(matches1)))

#grab matches again for array b
matches2 = regmatches(a.x.c.array.b[, 2], regexec("AxC.*\\d*\\.", as.character(a.x.c.array.b[, 2])))

matches2.1 = regmatches(as.character(matches2), regexpr("\\d+", as.character(matches2)))


a.x.c.array.a$a.x.c.index = ""
a.x.c.array.a$a.x.c.index = matches1.1

a.x.c.array.b$a.x.c.index = ""
a.x.c.array.b$a.x.c.index = matches2.1



#remove samples from array b that are not on array a
a.x.c.array.b = a.x.c.array.b[-which(!a.x.c.array.b$a.x.c.index %in% a.x.c.array.a$a.x.c.index)[3:4], ]

#reorganise array b to match array a in terms of sample order
a.x.c.array.b = a.x.c.array.b[match(a.x.c.array.a$a.x.c.index[3:nrow(a.x.c.array.a)], a.x.c.array.b$a.x.c.index), ]

all(a.x.c.array.a$a.x.c.index[3:124] == a.x.c.array.b$a.x.c.index)

#add two rows to array b dataframe
addrows = newdf(colnames(a.x.c.array.b))
addrows = rbind(addrows, addrows)

a.x.c.array.b = rbind(addrows, a.x.c.array.b)
a.x.c.array.b = rbind(a.x.c.array.b[3, ], a.x.c.array.b[1:2, ], a.x.c.array.b[4:nrow(a.x.c.array.b), ]) #do some rearragement of rows in array b

#TODO: merge array dataframes. find parent genotyping information for markers on array b

merged.array.df = cbind(a.x.c.array.a, a.x.c.array.b[3:ncol(a.x.c.array.b)])



#     ____________________________________________________________________________
#     INTERFACE WITH MYSQL                                                                                                        ####

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

rs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Avalon%'") #grab Avalon genotyping data for 35k array
rs2 = fetch(rs, n=-1)
rs2 = rs2[rs2$Var_col == "Avalon", ] 

#grab cadenza info
ca = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Cadenza%'") #grab Avalon genotyping data for 35k array
ca2 = fetch(ca, n=-1)
ca2 = ca2[ca2$Var_col == "Cadenza", ]



#     ____________________________________________________________________________
#     FURTHER PROCESSING                                                                                                            ####

merged.array.df = merged.array.df[, c(1:2, which(colnames(merged.array.df) %in% rs2$Probe_row))] #cut out non-35k array probes from merged df

g = rs2$Matrix_value[match(colnames(merged.array.df)[3:ncol(merged.array.df)], rs2$Probe_row)] #grab genotyping data for avalon
g[g=="AA"] = "A"
g[g=="BB"] = "B"
g[g=="AB"] = "H"
g[g=="NoCall"] = "-"

merged.array.df = convert.to.character.data.frame(merged.array.df)

merged.array.df[3, 3:ncol(merged.array.df)] = g

#prepare cadenza genotyping data
g = ca2$Matrix_value[match(colnames(merged.array.df)[3:ncol(merged.array.df)], ca2$Probe_row)] #grab genotyping data for avalon
g[g=="AA"] = "A"
g[g=="BB"] = "B"
g[g=="AB"] = "H"
g[g=="NoCall"] = "-"

merged.array.df[2, 3:ncol(merged.array.df)] = g



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

write.csv(merged.array.df3, "rotation1scripts_v4/processed_data/genotypes/avalon.x.cadenza/a.x.c.flipped.csv", row.names = F)



