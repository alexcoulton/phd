#do pairwise comparison of all genotypes

library(gtools)

read.genotype.data = function(q) read.csv(paste("original_data/genotype.data/alexq", q, "_man_cur_sk_map_genogeno.csv", sep = ""))
allgenotypes = lapply(1:4, read.genotype.data)

allgenotypes = lapply(allgenotypes, function(x){ 
    x = x[-c(1:3),]
    return(x)
}
)

all.genotypes.comb = newdf(colnames(allgenotypes[[1]]))
for(i in 1:4){
    all.genotypes.comb = rbind(all.genotypes.comb, allgenotypes[[i]])
}

all.genotypes.comb = all.genotypes.comb[-1,]

options(expressions = 40000) #change R-wide limit on the number of nested expressions that can be performed
pairwise.comparisons.df = as.data.frame(combinations(315, 2, repeats.allowed = F))


percentages = unlist(Map(function(x, y) (length(which(all.genotypes.comb[x,] == all.genotypes.comb[y,])) / 1298)*100, 
pairwise.comparisons.df[,1], pairwise.comparisons.df[,2]))

pairwise.comparisons.df$percentages = percentages

write.csv(pairwise.comparisons.df, "processed_data/pairwise.genotype.comparisons.csv", row.names = F)





stats = function(x){
    c(mean(x), sd(x),    max(x),    min(x))
}

stats(pairwise.comparisons.df$percentages)
hist(pairwise.comparisons.df$percentages)

print.hist.pdf = function(vector, filename, xlabel){
    #args:
    #vector - a numeric vector
    #filename - a character string for the filename
    #xlabel - a character string for the x axis label
    pdf(file = paste("plots/", filename, sep = ""))
    hist(vector, xlab = xlabel)
    dev.off()
}

print.hist.pdf(pairwise.comparisons.df$percentages, "pairwise.genotype.comparisons.pdf", "Genotype similarity (%)")

