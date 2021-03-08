source("rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R")

axp.map = read_csv("rotation1scripts_v4/original_data/sacha.genetic.maps/Allen genetic map AxP.csv")

#convert dataframe into rQTL format
axp.map2 = as.data.frame(t(axp.map))
axp.map2 = convert.to.character.data.frame(axp.map2)
colnames(axp.map2) = axp.map2[1, ]
axp.map2 = axp.map2[-1, ]
axp.map2[axp.map2 == "AB"] = "H"
axp.map2[axp.map2 == "AA"] = "A"
axp.map2[axp.map2 == "BB"] = "B"

axp.map3 = axp.map2[-2, ]

#prepare dataframe for makegenotypelist function
axp.map3.1 = cbind(1:nrow(axp.map3), rownames(axp.map3), axp.map3)
axp.map3.1 = convert.to.character.data.frame(axp.map3.1)
axp.map3.1[1, 1:2] = ""



axp.map4 = split.geno.df.into.list.by.chromosome(axp.map3.1)
axp.map5 = lapply(axp.map4, function(x){
    x = x[, -1]
    x
})

axpchromos = unlist(lapply(axp.map4, function(x){
    g = as.character(x[1, ])
    unique(g)[2]
}))

source("rotation1scripts_v4/scripts/r/recombination_analysis_functions.R")
axpmap.genolist = makegenotypelist(axp.map3.1)


do.call(grid.arrange, make.plots(axp.map5, axpchromos))
