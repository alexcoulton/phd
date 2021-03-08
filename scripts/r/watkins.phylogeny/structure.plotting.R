#### ARRAY DATA ####

array.path1 = "~/project.phd.main/software/CLUMPAK/26_03_2015_CLUMPAK/CLUMPAK/array.output2/K=3/CLUMPP.files/ClumppIndFile.output"
exome.path1 = "~/project.phd.main/software/CLUMPAK/26_03_2015_CLUMPAK/CLUMPAK/exome.output2/K=3/CLUMPP.files/ClumppIndFile.output"

prepare.structure.df = function(path1){
    k3clump = read.table(path1, sep = ":", stringsAsFactors = F)
    k3clump2 = as.data.frame(lapply(k3clump$V2, function(x){
        strsplit(as.character(x), " ")
    }))
    
    k3clump2 = t(k3clump2)
    k3clump2 = k3clump2[, -(1:2)]
    k3clump2 = reset.rownames(k3clump2)
    k3clump2 = as.data.frame(k3clump2)
    k3clump2$pop = 1
    k3 = cbind(k3clump2[, 1], k3clump2)
    k3$`k3clump2[, 1]` = k3$pop
    k3 = reset.colnames(k3)
    k3 = k3[, -ncol(k3)]
    k3[] = lapply(k3, function(x) as.numeric(as.character(x)))
    k3
}

a.k3 = prepare.structure.df(array.path1)

a.k3$watkins = unique(new.array.data2$`new.array.data2[, 1]`)

a.k3$V1 = wilkins.data2$Region
pop.sort = sort(a.k3$V1, index.return = T)$ix
a.k3 = a.k3[pop.sort, ]


# a.k3 = a.k3[sort(wilkins.data2$Region, index.return = T)$ix, ]
#order the Ks so that they match between array and exome
a.k3 = a.k3[, c(1, sort(unlist(a.k3[1, 2:(ncol(a.k3) - 1)]), index.return = T)$ix + 1, ncol(a.k3))]

a.k3$ind2 = 1:nrow(a.k3)
a.k3.v2 = split(a.k3, a.k3$V1)

#order by K membership within regions
a.k3.v3 = lapply(a.k3.v2, function(x){
    max.component = which.max(na.omit(sapply(x[, 2:(ncol(x) - 2)], function(y){
        sum(as.numeric(y))
    })))
    print(max.component)
    x[sort(x[[(max.component + 1)]], index.return = T, decreasing = T)$ix, ]
    
})

a.k3.v4 = bind_rows(a.k3.v3)


a.k3.v4$ind = 1:nrow(a.k3.v4)
colnames(a.k3.v4)[2:(ncol(a.k3.v4) - 3)] = paste0("K", 1:length(colnames(a.k3.v4)[2:(ncol(a.k3.v4) - 3)]))

a.k3.v5 = melt(a.k3.v4, id.vars = c("ind", "ind2", "V1", "watkins"))
a.k3.v5$value = as.numeric(a.k3.v5$value)
a.k3.v5$watkins = factor(as.character(a.k3.v5$watkins), levels = as.character(unique(a.k3.v5$watkins)))


splot.labs = c("AS",  "A", "EE", "WE", "ME", "NA", "U")


a.k3.v5$abv = a.k3.v5$V1

Map(function(x, y){
a.k3.v5$abv[which(a.k3.v5$V1 == x)] <<- y
}, unique(a.k3.v5$V1), splot.labs)

a.k3.v5$abv = factor(a.k3.v5$abv, levels = c('WE', 'EE', 'U', 'ME', 'AS', 'NA', 'A'))

splot1 = ggplot(a.k3.v5, aes(x = watkins, y = value, color = variable, fill = variable)) + geom_bar(stat = 'identity', position = "fill") +
    scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T) +
    facet_grid(cols = vars(abv), scales = "free_x", space = "free_x") +
    theme_bw(base_size = 22) +
    theme(panel.grid.major = element_blank(), legend.title = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = 60, hjust = 1)) +
    xlab("Individuals") + ylab("Cluster membership") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())





splot1



e.k3 = prepare.structure.df(exome.path1)
e.k3$V1 = wilkins.data2$Region
e.k3 = e.k3[pop.sort, ]


# e.k3 = e.k3[sort(wilkins.data2$Region, index.return = T)$ix, ]
#order the Ks so that they match between array and exome
e.k3 = e.k3[, c(1, sort(unlist(e.k3[1, 2:ncol(e.k3)]), index.return = T)$ix + 1)]


e.k3$ind2 = 1:nrow(e.k3)

e.k3 = e.k3[match(a.k3.v4$ind2, e.k3$ind2), ]
e.k3$ind = 1:nrow(e.k3)

colnames(e.k3)[2:(ncol(e.k3) - 2)] = paste0("K", 1:length(colnames(e.k3)[2:(ncol(e.k3) - 2)]))

# swap.k1 = e.k3$K1
# 
# e.k3$K1 = e.k3$K2
# e.k3$K2 = swap.k1

e.k4 = melt(e.k3, id.vars = c("ind", "ind2", "V1"))
e.k4$value = as.numeric(e.k4$value)
e.k4$ind = factor(as.character(e.k4$ind), levels = unique(as.character(e.k4$ind)))

e.k4$watkins = a.k3.v5$watkins

# levels(e.k4$variable) = c("V5", "V4", "V6", "V3", "V2")
# levels(e.k4$variable) = c("V3", "V4", "V5", "V6", "V2")

e.k4$abv = a.k3.v5$abv


plotlabels1 = unique(new.array.data2$`new.array.data2[, 1]`)

splot2 = ggplot(e.k4, aes(x = watkins, y = value, color = variable, fill = variable)) + geom_bar(stat = 'identity', position = "fill") +
    scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T) + facet_grid(cols = vars(abv), scales = "free_x", space = "free_x") + theme_bw(base_size = 22) +
    theme(panel.grid.major = element_blank(), legend.title = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_text(angle = 60, hjust = 1)) +
    xlab("Individuals") + ylab("Cluster membership") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())



#main_plot
grid.arrange(splot1, splot2)


png("rotation1scripts_v4/plots/exome.array.comparison/structurek3.png", 1200, 800)
grid.arrange(splot1, splot2)
dev.off()

as.numeric(unique(new.array.data2$`new.array.data2[, 1]`))



plot.structure = function(k3){
    # coord1 = 1:nrow(k3)
    # coord2 = unlist(lapply(k3[, 2:ncol(k3)], function(x){
    #     g = which(x > 0.2)
    #     g2 = g[sort(x[g], index.return = T, decreasing = T)$ix]
    #     g2
    # }))
    #
    # coord3 = coord1[-coord2]
    #
    # coord4 = c(coord2, coord3)
    #
    # k3 = k3[coord4, ]

    # k3 = k3[sort(k3$V2, index.return = T)$ix, ]
    k3$ind = 1:nrow(k3)
    
    k4 = melt(k3, id.vars = c("ind", "V1"))
    k4$value = as.numeric(k4$value)
    ggplot(k4, aes(x = ind, y = value, color = variable, fill = variable)) + geom_bar(stat = 'identity', position = "fill") +
        scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T)
}




structureplot1 = plot.structure(a.k3.v4)
structureplot2 = plot.structure(e.k3)

structureplot1
grid.arrange(structureplot1, structureplot2)

#### array k3 processing ####

coord1 = 1:nrow(a.k3)
coord2 = unlist(lapply(a.k3[, 2:ncol(a.k3)], function(x){
    g = which(x > 0.5)
    g2 = g[sort(x[g], index.return = T, decreasing = T)$ix]
    g2
}))

coord3 = coord1[-coord2]

a.coord4 = c(coord2, coord3)

# a.k3 = a.k3[a.coord4, ]    

# a.k3 = a.k3[sort(a.k3$V2, index.return = T)$ix, ]
a.k3$ind = 1:nrow(a.k3)

a.k4 = melt(a.k3, id.vars = c("ind", "V1"))
a.k4$value = as.numeric(a.k4$value)

ggplot(a.k4, aes(x = ind, y = value, color = variable, fill = variable)) + geom_bar(stat = 'identity', position = "fill") +
             scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T)



#### exome k3 processing #### 

coord1 = 1:nrow(e.k3)
coord2 = unlist(lapply(e.k3[, 2:ncol(e.k3)], function(x){
    g = which(x > 0.5)
    g2 = g[sort(x[g], index.return = T, decreasing = T)$ix]
    g2
}))

coord3 = coord1[-coord2]

e.coord4 = c(coord2, coord3)

# e.k3 = e.k3[e.coord4, ]    

# e.k3 = e.k3[sort(e.k3$V2, index.return = T)$ix, ]
e.k3$ind = 1:nrow(e.k3)

e.k4 = melt(e.k3, id.vars = c("ind", "V1"))
e.k4$value = as.numeric(e.k4$value)

grid.arrange(
ggplot(a.k4, aes(x = ind, y = value, color = variable, fill = variable)) + geom_bar(stat = 'identity', position = "fill") +
    scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T),
ggplot(e.k4, aes(x = ind, y = value, color = variable, fill = variable)) + geom_bar(stat = 'identity', position = "fill") +
             scale_color_viridis(discrete = T) + scale_fill_viridis(discrete = T)
)
