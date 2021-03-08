c.x.a.f2 = read.csv("rotation1scripts_v4/processed_data/genotypes/c.x.a.f2/c.x.a.f2.csv")

low.h = which(unlist(lapply(c.x.a.f2, function(x){
    length(which(x == "H"))
})) < 30)[-c(1, 2)]

c.x.a.f2 = c.x.a.f2[, -low.h]


c.x.a.651 = c.x.a.f2[grep("65-1", c.x.a.f2$probeset_id), ]

c.x.a.653 = c.x.a.f2[grep("65-3", c.x.a.f2$probeset_id), ]


c.x.a.651.seg = unlist(lapply(c.x.a.651, function(x){
        a = length(which(x == "A"))
        b = length(which(x == "B"))
        h = length(which(x == "H"))
        a / b
}))

c.x.a.653.seg = unlist(lapply(c.x.a.653, function(x){
    a = length(which(x == "A"))
    b = length(which(x == "B"))
    h = length(which(x == "H"))
    a / b
}))

g = which(c.x.a.651.seg > 1.7)

g1 = which(c.x.a.653.seg > 1.7)

a.x.c.sacha = combine.list.of.data.frames(all.m.8.iwgsc.4b.rev[[1]])

na.omit(match(a.x.c.sacha$marker, switch.affy.format(colnames(c.x.a.653))))

g[which(names(g) %in% names(g1))]


