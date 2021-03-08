#35k arrays
sacha.axp = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/apogee.x.paragon.sacha/all.markers.txt")
#remove excess parents
sacha.axp = sacha.axp[-c(3, 4, 6, 7), ]
sacha.axp[sacha.axp == "AA"] = "A"
sacha.axp[sacha.axp == "AB"] = "H"
sacha.axp[sacha.axp == "BB"] = "B"
sacha.axp[sacha.axp == "NoCall"] = "-"

sacha.csxp = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/chinese.spring.x.paragon/c.x.p.all.markers.txt")

#820k arrays

#### opata synthetic ####
#### opata synthetic processing ####
sacha.oxs.a = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/opata.x.synthetic/oxs.arraya.all.snps.txt")
sacha.oxs.b = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/opata.x.synthetic/oxs.all.snps.array.b.txt")

g = multi.str.split(as.character(sacha.oxs.a$probeset_id), ".CEL_call_code", 1)
g = multi.str.split(g, "_A$", 1)
g2 = multi.str.split(g, "(?U)^.*_", 3)
length(g2) == length(unique(g2))

sacha.oxs.a$probeset_id = g2

g = multi.str.split(as.character(sacha.oxs.b$probeset_id), ".CEL_call_code", 1)
g = multi.str.split(g, "_B$", 1)
g2 = multi.str.split(g, "(?U)^.*_", 3)
length(g2) == length(unique(g2))

sacha.oxs.b$probeset_id = g2

sacha.oxs.b = sacha.oxs.b[match(sacha.oxs.a$probeset_id, sacha.oxs.b$probeset_id), ]

sacha.oxs.a[which(!as.character(sacha.oxs.a$probeset_id) == as.character(sacha.oxs.b$probeset_id)), 2]
sacha.oxs.b[which(!as.character(sacha.oxs.a$probeset_id) == as.character(sacha.oxs.b$probeset_id)), 2]





as.character(sacha.oxs.a$probeset_id) == as.character(sacha.oxs.b$probeset_id)

sacha.oxs.a[, 1:4]
sacha.oxs.b[, 1:4]

#### oxs final ####
sacha.oxs.comb = bind_cols(sacha.oxs.a, sacha.oxs.b[, 3:ncol(sacha.oxs.b)])
sacha.oxs.35k = sacha.oxs.comb[, which(colnames(sacha.oxs.comb) %in% colnames(sacha.csxp))]

missing.snps1 = colnames(sacha.csxp)[which(!colnames(sacha.csxp) %in% colnames(sacha.oxs.comb))]
write(missing.snps1, "rotation1scripts_v4/original_data/sacha.genetic.maps/missing.snps.txt")


##### avalon x cadenza ####

#### axc processing ####
sacha.axc.a = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza/avalon.x.cad.all.snps.arraya.txt")
sacha.axc.b = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza/avalon.x.cad.all.snps.arrayb.txt")

g = multi.str.split(as.character(sacha.axc.a$probeset_id), ".CEL_call_code", 1)
g = multi.str.split(g, "_A$", 1)
g = multi.str.split(g, "^[A-Z0-9]*_", 2)
length(g) == length(unique(g))

sacha.axc.a$probeset_id = g

g = multi.str.split(as.character(sacha.axc.b$probeset_id), ".CEL_call_code", 1)
g = multi.str.split(g, "_B$", 1)
g = multi.str.split(g, "^[A-Z0-9]*_", 2)
length(g) == length(unique(g))

sacha.axc.b$probeset_id = g

sacha.axc.a[which(!as.character(sacha.axc.a$probeset_id) == as.character(sacha.axc.b$probeset_id)), 2]
sacha.axc.b[which(!as.character(sacha.axc.a$probeset_id) == as.character(sacha.axc.b$probeset_id)), 2]

sacha.axc.a$probeset_id[which(sacha.axc.a$probeset_id == "AxC62")] = "AxC_62"

sacha.axc.b = sacha.axc.b[match(sacha.axc.a$probeset_id, sacha.axc.b$probeset_id), ]

as.character(sacha.axc.a$probeset_id) == as.character(sacha.axc.b$probeset_id)

sacha.axc.a[, 1:4]
sacha.axc.b[, 1:4]

sacha.axc.comb = bind_cols(sacha.axc.a, sacha.axc.b[, 3:ncol(sacha.axc.b)])

#### final axc ####
sacha.axc.35k = sacha.axc.comb[, which(colnames(sacha.axc.comb) %in% colnames(sacha.csxp))]


#### rialto x savannah ####

#### rxs processing ####
sacha.rxs.a = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/rialto.x.savannah/rialto.x.sav.all.snps.arraya.txt")
sacha.rxs.b = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/rialto.x.savannah/rialto.x.sav.all.snps.arrayb.txt")

g = multi.str.split(as.character(sacha.rxs.a$probeset_id), ".CEL_call_code", 1)
g = multi.str.split(g, "_A$", 1)
g2 = multi.str.split(g, "^[A-Z][0-9][0-9]_", 2)
length(g2) == length(unique(g2))

sacha.rxs.a$probeset_id = g2

g = multi.str.split(as.character(sacha.rxs.b$probeset_id), ".CEL_call_code", 1)
g = multi.str.split(g, "_B$", 1)
g2 = multi.str.split(g, "^[A-Z][0-9][0-9]_", 2)
length(g2) == length(unique(g2))


sacha.rxs.b$probeset_id = g2

sacha.rxs.a[which(!as.character(sacha.rxs.a$probeset_id) == as.character(sacha.rxs.b$probeset_id)), 2]
sacha.rxs.b[which(!as.character(sacha.rxs.a$probeset_id) == as.character(sacha.rxs.b$probeset_id)), 2]

as.character(sacha.rxs.a$probeset_id) == as.character(sacha.rxs.b$probeset_id)

#### rxs final ####
sacha.rxs.comb = bind_cols(sacha.rxs.a, sacha.rxs.b[, 3:ncol(sacha.rxs.b)])
sacha.rxs.35k = sacha.rxs.comb[, which(colnames(sacha.rxs.comb) %in% colnames(sacha.csxp))]


#### get parents from cerealsdb ####

rearrange.mysql.output = function(x){
    x[x == "AA"] = "A"
    x[x == "AB"] = "H"
    x[x == "BB"] = "B"
    x
    
    # browser()
    
    x = transpose(x)
    colnames(x) = x[2, ]
    
    x = add.column.at.position(x, 0)
    x[3, 1] = x[1, 2]
    x = x[-c(1, 2), ]
    colnames(x)[1] = "probeset_id"
    x
    
}

library(RMySQL)

mydb = dbConnect(MySQL(), user = 'alex', password = 'alex', db = 'cerealsdb', host = '127.0.0.1')

ava = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Avalon%'") #grab Avalon genotyping data for 35k array
ava2 = fetch(ava, n=-1)
ava2 = ava2[ava2$Var_col == "Avalon", ] 

#grab cadenza info
ca = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Cadenza%'") #grab Avalon genotyping data for 35k array
ca2 = fetch(ca, n=-1)
ca2 = ca2[ca2$Var_col == "Cadenza", ]

op = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Opata%'") #grab Avalon genotyping data for 35k array
op2 = fetch(op, n=-1)
op2 = op2[op2$Var_col == "Opata", ]

syn = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Synthetic%'") #grab Avalon genotyping data for 35k array
syn2 = fetch(syn, n=-1)
syn2 = syn2[syn2$Var_col == "Synthetic", ]

rial = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Rialto%'") #grab Avalon genotyping data for 35k array
rial2 = fetch(rial, n=-1)
rial2 = rial2[rial2$Var_col == "Rialto", ]


sav = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Savannah%'") #grab Avalon genotyping data for 35k array
sav2 = fetch(sav, n=-1)
sav2 = sav2[sav2$Var_col == "Savannah", ]

cs = dbSendQuery(mydb, "SELECT * FROM 35k_array_oct2016_named WHERE Var_col LIKE '%Chinese%'") #grab Avalon genotyping data for 35k array
cs2 = fetch(cs, n=-1)
cs2 = cs2[cs2$Var_col == "Chinese_Spring", ]

cs3 = rearrange.mysql.output(cs2)



parents = list(ava2, ca2, op2, syn2, rial2, sav2)


parents = lapply(parents, rearrange.mysql.output)

rialto.final = parents[[5]][, na.omit(match(colnames(sacha.rxs.35k), colnames(parents[[5]])))]

savannah.final = parents[[6]][, na.omit(match(colnames(sacha.rxs.35k), colnames(parents[[6]])))]

opata.final = parents[[3]][, na.omit(match(colnames(sacha.oxs.35k), colnames(parents[[3]])))]

synthetic.final = parents[[4]][, na.omit(match(colnames(sacha.oxs.35k), colnames(parents[[4]])))]

avalon.final = parents[[1]][, na.omit(match(colnames(sacha.axc.35k), colnames(parents[[1]])))]

cadenza.final = parents[[2]][, na.omit(match(colnames(sacha.axc.35k), colnames(parents[[2]])))]


#### add parents to genotyping dataframes ####

#rxs
sacha.rxs.35k2 = add.row.at.position(sacha.rxs.35k, 1)
sacha.rxs.35k2[2, 2:ncol(sacha.rxs.35k2)] = rialto.final

sacha.rxs.35k2 = add.row.at.position(sacha.rxs.35k2, 1)
sacha.rxs.35k2[2, 2:ncol(sacha.rxs.35k2)] = savannah.final

#oxs
sacha.oxs.35k2 = add.row.at.position(sacha.oxs.35k, 1)
sacha.oxs.35k2[2, 2:ncol(sacha.oxs.35k2)] = synthetic.final

sacha.oxs.35k2 = add.row.at.position(sacha.oxs.35k2, 1)
sacha.oxs.35k2[2, 2:ncol(sacha.oxs.35k2)] = opata.final

#axc
#get parents from F2 array for axc as these were better than parents from CerealsDB (less heterozygosity)
a.x.cf2 = convert.aas.to.rqtl("rotation1scripts_v4/original_data/genotype.data/avalon.x.cadenza.f2/axc61-1.all.snps.txt")
a.x.cf2.parents = a.x.cf2[which(a.x.cf2$probeset_id == "M23_AxC 61-1 G12.CEL_call_code" | a.x.cf2$probeset_id == "O23_AxC 61-1 H12.CEL_call_code"), ]
a.x.cf2.parents$probeset_id = c("Cadenza", "Avalon")

a.x.cf2.parents = a.x.cf2.parents[, match(colnames(sacha.axc.35k), colnames(a.x.cf2.parents))]
a.x.cf2.parents = convert.to.character.data.frame(a.x.cf2.parents)

sacha.axc.35k2 = add.row.at.position(sacha.axc.35k, 1)
sacha.axc.35k2 = add.row.at.position(sacha.axc.35k, 1)
sacha.axc.35k2[2:3, ] = a.x.cf2.parents[2:1, ]

#csxp
sacha.csxp2 = add.row.at.position(sacha.csxp, 1)
sacha.csxp2[2, 2:ncol(sacha.csxp2)] = cs3[, na.omit(match(colnames(sacha.csxp2), colnames(cs3)))]
#remove excess parents
sacha.csxp2 = sacha.csxp2[-(4:5), ]

all.sacha.geno.dfs = list(sacha.csxp2, sacha.axc.35k2, sacha.oxs.35k2, sacha.rxs.35k2, sacha.axp)

#functions from generic.axiom.data.preparation.R

all.sacha.geno.dfs2 = lapply(all.sacha.geno.dfs, function(x){
    x[x == "NoCall"] = "-"
    x = cleanup(x, 2:3)
    x = convert.to.character.data.frame(x)
    assign.parental.genotypes(x, 2:3)
})

lapply(all.sacha.geno.dfs2, function(x) x[1:5, 1:5])

#read in allen genetic maps
# a.x.p.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic map AxP.csv", stringsAsFactors = F)
# s.x.r.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/sxr.allen.csv", stringsAsFactors = F)
# oxs.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/o.x.s.allen.csv", stringsAsFactors = F)
# csxp.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.cs.x.p.map.csv", stringsAsFactors = F)
# axc.allen = read.csv("rotation1scripts_v4/original_data/sacha.genetic.maps/allen.axc.csv", stringsAsFactors = F)

# a.x.p.allen[a.x.p.allen == "AA"] = "A"
# a.x.p.allen[a.x.p.allen == "AB"] = "H"
# a.x.p.allen[a.x.p.allen == "BB"] = "B"
# a.x.p.allen[a.x.p.allen == "NoCall"] = "-"

allen.maps = list(csxp.allen, axc.allen, oxs.allen, s.x.r.allen, a.x.p.allen)

# write.csv(sacha.csxp2, "rotation1scripts_v4/original_data/sacha.genetic.maps/redone/csxp.allen.csv", row.names = F)
# write.csv(sacha.axc.35k2, "rotation1scripts_v4/original_data/sacha.genetic.maps/redone/axc.allen.csv", row.names = F)
# write.csv(sacha.oxs.35k2, "rotation1scripts_v4/original_data/sacha.genetic.maps/redone/oxs.allen.csv", row.names = F)
# write.csv(sacha.rxs.35k2, "rotation1scripts_v4/original_data/sacha.genetic.maps/redone/sxr.allen.csv", row.names = F)
# write.csv(sacha.axp, "rotation1scripts_v4/original_data/sacha.genetic.maps/redone/a.x.p.allen.csv", row.names = F)

#filter genotyping dataframes to contain only markers that are common between geno.df and allen genetic maps
all.sacha.geno.dfs3 = Map(function(x, y){
    # browser()

    colnames(y)[1] = "marker"
    # to.rm1 = which(!y$marker %in% colnames(x))
    # print(paste0(length(to.rm1), " markers omitted"))
    # 
    # if(length(to.rm1) > 0){
    #     y = y[-to.rm1, ]
    # }
    
    x = x[, c(1, 2, na.omit(match(y$marker, colnames(x))))]
    x[1, 3:ncol(x)] = y$chr
    x
}, all.sacha.geno.dfs2, allen.maps)



allen.maps2 = Map(function(x, y){
    colnames(y)[1] = "marker"
    y[which(y$marker %in% colnames(x)), ]
}, all.sacha.geno.dfs3, allen.maps)

sapply(all.sacha.geno.dfs3, ncol)
sapply(allen.maps2, nrow)

#see if any of the genotypes from sacha's map match so I can identify which individual is which
lapply(1:5, function(x){
    g.test1 = as.character(all.sacha.geno.dfs3[[x]][8, ][3:ncol(all.sacha.geno.dfs3[[x]])])
    g.test2 = sapply(allen.maps2[[x]], function(y){
        # browser()
        length(which(y == g.test1))
    })
    
    # if(x == 5) browser()
    
    g.test2[which(g.test2 == max(g.test2))] / length(g.test1)
})


#try using unflipped markers
all.sacha.geno.dfs4 = Map(function(x, y){
    # browser()
    
    colnames(y)[1] = "marker"
    # to.rm1 = which(!y$marker %in% colnames(x))
    # print(paste0(length(to.rm1), " markers omitted"))
    
    # if(length(to.rm1) > 0){
    #     y = y[-to.rm1, ]
    # }
    
    # to.rm2 = which(!colnames(x)[3:ncol(x)] %in% y$marker)
    # x = x[, -(to.rm2 + 2)]
    
    x = x[, c(1, 2, na.omit(match(y$marker, colnames(x))))]
    x[1, 3:ncol(x)] = y$chr
    # browser()
    x
}, all.sacha.geno.dfs, allen.maps)



allen.maps3 = Map(function(x, y){
    colnames(y)[1] = "marker"
    y[which(y$marker %in% colnames(x)), ]
}, all.sacha.geno.dfs4, allen.maps)




lapply(all.sacha.geno.dfs4, function(x){
    head(x)[1:5, 1:5]
})

sapply(all.sacha.geno.dfs4, ncol)
sapply(allen.maps3, nrow)

#see if any of the genotypes from sacha's map match so I can identify which individual is which
g.test1 = all.sacha.geno.dfs4[[1]][16, ][2:nrow(all.sacha.geno.dfs3[[1]])]
g.test2 = sapply(allen.maps[[1]], function(x){
    # browser()
    length(which(x == g.test1))
})


lapply(1:5, function(x){
    g.test1 = as.character(all.sacha.geno.dfs4[[x]][8, ][3:ncol(all.sacha.geno.dfs4[[x]])])
    g.test2 = sapply(allen.maps[[x]], function(x){
        # browser()
        length(which(x == g.test1))
    })
    g.test2[which(g.test2 == max(g.test2))] / length(g.test1)
})



allen.maps[[1]][, 1][1:5] 

colnames(all.sacha.geno.dfs4[[1]])[1:5]

which(g.test2 > 5000)
