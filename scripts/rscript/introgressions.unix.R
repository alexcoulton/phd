#!/usr/bin/Rscript
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")

introgressions = newdf(c("elite", "wild.relative", "intro.start.bp", "intro.end.bp", "intro.length"))

# g4.small2 = split.df.into.list.of.dfs.by.column(g4.small, "chromosome")

args = commandArgs(trailingOnly = T)
file1 = args[2]
var1 = args[1]

geno.df1 = read.csv(p("rotation1scripts_v4/original_data/genotype.data/820k.genotypes/ordered.split/", file1), header = T, stringsAsFactors = F)

geno.df1 = geno.df1[, 2:ncol(geno.df1)]

geno.df1 = list(geno.df1)


comparison = function(df1, col1, col2){
  comp1 = which(df1[, col1] == df1[, col2])
  list.consecutive.subsequences(comp1, 1)
}



global.comparison = function(variety.to.compare, list.of.genotype.dfs1){
  #args:
  #variety.to.compare: character string indicating the exact name of the elite variety to compare to all
    #other varieties, as listed in the column name of list.of.genotype.dfs1
  #list.of.genotype.dfs1: a list of genotype dataframes split by chromosome using split.df.into.lsit.of.dfs.by.column
  
  count1 = make.counter()
  #run a loop for each chromosome
  main.analysis = lapply(list.of.genotype.dfs1, function(chromo){
    print(p("chromosome ", count1()))
    current.chromo = chromo$chromosome[[1]]
    var.coord = which(colnames(chromo) == variety.to.compare)
    
    comparison.list = 4:120 #4:120 is the range of columns containing wild relatives 
    
    
    #run a loop for each variety to compare the intitial variety to
    count2 = make.counter()
    lapply(comparison.list, function(q){
      print(p("comparison ", count2()))
      consec.list1 = comparison(chromo, var.coord, q)
      # browser()
      to.rm1 = which(unlist(lapply(consec.list1, function(q1) length(q1) == 1)))
      if(length(to.rm1) > 0) consec.list1 = consec.list1[-to.rm1]
      
      #generate dataframe containing physical info. etc
      parse.intros = lapply(consec.list1, function(x){
        intro.start = chromo$phys.pos[min(x)]
        intro.end = chromo$phys.pos[max(x)]
        intro.length = intro.end - intro.start
        elite = variety.to.compare
        wild.relative = colnames(chromo)[q]
        chromo1 = current.chromo
        
        all.in = c(elite, wild.relative, chromo1, intro.start, intro.end, intro.length)
        names(all.in) = c("elite", "wild.relative", "chromosome", "intro.start", "intro.end", "intro.length")
        return(all.in)
      })
      
      as.data.frame(t(as.data.frame(parse.intros)))
      
      
      
    })
    
    
    
  })
 
  main2 = lapply(main.analysis, combine.list.of.data.frames)
  main3 = combine.list.of.data.frames(main2)
  return(main3)
}


comp1 = global.comparison(var1, geno.df1)

write.csv(comp1, p("rotation1scripts_v4/processed_data/introgressions/int.", var1, ".", file1), row.names = F)




