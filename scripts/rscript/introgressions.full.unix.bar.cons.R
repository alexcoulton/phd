#!/usr/bin/Rscript
setwd("/home/ac14037/project.phd.main")
source("rotation1scripts_v4/scripts/r/functions.R")

library(parallel)

# install.packages("readr")
# library(readr)

introgressions = newdf(c("elite", "wild.relative", "intro.start.bp", "intro.end.bp", "intro.length"))

# g4.small2 = split.df.into.list.of.dfs.by.column(g4.small, "chromosome")

args = commandArgs(trailingOnly = T)

genotypes.root.path = args[1]
file1 = args[2]
# i2 = as.numeric(args[3])

# print(i2)

#i2 is the column containing the hexaploid variety to find introgressions for, this can be done in parallel

#TEMPORARY ASSIGNMENT
file1 = "data_analysis2.csv"

geno.df1 = read.csv(p(genotypes.root.path, file1), header = T, stringsAsFactors = F)

#extract sample type information from dataframe
sample.type.info = geno.df1[1, ]
#remove sample type information from dataframes
geno.df1 = geno.df1[-1, ]    

wild.rel.col.coords = which(sample.type.info == "relative")
hexaploid.col.coords = which(sample.type.info == "wheat")

all.comp = mclapply(hexaploid.col.coords, function(i2){





  geno.df1 = read.csv(p(genotypes.root.path, file1), header = T, stringsAsFactors = F)    



  # geno.df1 = geno.df1[, 2:ncol(geno.df1)]

  var1 = colnames(geno.df1)[i2]

  # print(var1)

  print(colnames(geno.df1)[1:10])
  # print(geno.df1$chromosome)

  geno.df1 = split.df.into.list.of.dfs.by.column(geno.df1, "chromosome")
  if(class(geno.df1) != "list") geno.df1 = list(geno.df1)


  # print(head(geno.df1))

  comparison = function(df1, col1, col2){
    comp1 = which(df1[, col1] == df1[, col2])
    list.consecutive.subsequences(comp1, 1)
  }

  # print(geno.df1)

  global.comparison = function(variety.to.compare, list.of.genotype.dfs1){
    #args:
    #variety.to.compare: character string indicating the exact name of the elite variety to compare to all
      #other varieties, as listed in the column name of list.of.genotype.dfs1
    #list.of.genotype.dfs1: a list of genotype dataframes split by chromosome using split.df.into.list.of.dfs.by.column
    
    count1 = make.counter()
    #run a loop for each chromosome
    main.analysis = lapply(list.of.genotype.dfs1, function(chromo){
      print(p("chromosome ", count1()))
      current.chromo = chromo$chromosome[[1]]
      chromo$phys.pos = as.numeric(chromo$phys.pos)
      var.coord = which(colnames(chromo) == variety.to.compare)
      
      comparison.list = wild.rel.col.coords #4:120 is the range of columns containing wild relatives 
      
      
      #run a loop for each variety to compare the intitial variety to
      count2 = make.counter()
      lapply(comparison.list, function(q){
        print(p("comparison ", count2()))
        consec.list1 = comparison(chromo, var.coord, q)

        # browser()
        to.rm1 = which(unlist(lapply(consec.list1, function(q1) length(q1) == 1 | length(q1) == 0)))
        if(length(to.rm1) > 0) consec.list1 = consec.list1[-to.rm1]
        
        #debug code
        # g.7 = unlist((lapply(consec.list1, length)))
        # if(0 %in% g.7) print(g.7)

        #generate dataframe containing physical info. etc
      

        parse.intros = lapply(consec.list1, function(x){
          intro.start = chromo$phys.pos[min(x)]
          intro.end = chromo$phys.pos[max(x)]
          intro.length = intro.end - intro.start
          number.snps = length(x)
          snp.start = min(x)
          snp.end = max(x)
          snp.density = number.snps / intro.length
          elite = variety.to.compare
          wild.relative = colnames(chromo)[q]
          chromo1 = current.chromo
          
          all.in = c(elite, wild.relative, chromo1, intro.start, intro.end, intro.length, number.snps, snp.start, snp.end, snp.density)
          names(all.in) = c("elite", "wild.relative", "chromosome", "intro.start", "intro.end", "intro.length", "number.snps", "snp.start", "snp.end", "snp.density")
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

  # write.csv(comp1, p("rotation1scripts_v4/processed_data/introgressions/sacha/int.", var1, ".", file1), row.names = F)

  return(comp1)

}, mc.cores = 50)

all.comp2 = combine.list.of.data.frames(all.comp)

#remove introgressions with only 3 or below SNPs
all.comp2$number.snps = as.numeric(all.comp2$number.snps)
all.comp2 = all.comp2[which(all.comp2$number.snps > 3), ]

write.csv(all.comp2, "rotation1scripts_v4/processed_data/introgressions/sacha/all.introgressions_anal2_11012019.csv", row.names = F)


####################### TEMPORARY --- CAN DELETE ############################
file1 = "data_analysis3.csv"

geno.df1 = read.csv(p(genotypes.root.path, file1), header = T, stringsAsFactors = F)

#extract sample type information from dataframe
sample.type.info = geno.df1[1, ]
#remove sample type information from dataframes
geno.df1 = geno.df1[-1, ]    

wild.rel.col.coords = which(sample.type.info == "relative")
hexaploid.col.coords = which(sample.type.info == "wheat")

all.comp = mclapply(hexaploid.col.coords, function(i2){





  geno.df1 = read.csv(p(genotypes.root.path, file1), header = T, stringsAsFactors = F)    



  # geno.df1 = geno.df1[, 2:ncol(geno.df1)]

  var1 = colnames(geno.df1)[i2]

  # print(var1)

  print(colnames(geno.df1)[1:10])
  # print(geno.df1$chromosome)

  geno.df1 = split.df.into.list.of.dfs.by.column(geno.df1, "chromosome")
  if(class(geno.df1) != "list") geno.df1 = list(geno.df1)


  # print(head(geno.df1))

  comparison = function(df1, col1, col2){
    comp1 = which(df1[, col1] == df1[, col2])
    list.consecutive.subsequences(comp1, 1)
  }

  # print(geno.df1)

  global.comparison = function(variety.to.compare, list.of.genotype.dfs1){
    #args:
    #variety.to.compare: character string indicating the exact name of the elite variety to compare to all
      #other varieties, as listed in the column name of list.of.genotype.dfs1
    #list.of.genotype.dfs1: a list of genotype dataframes split by chromosome using split.df.into.list.of.dfs.by.column
    
    count1 = make.counter()
    #run a loop for each chromosome
    main.analysis = lapply(list.of.genotype.dfs1, function(chromo){
      print(p("chromosome ", count1()))
      current.chromo = chromo$chromosome[[1]]
      chromo$phys.pos = as.numeric(chromo$phys.pos)
      var.coord = which(colnames(chromo) == variety.to.compare)
      
      comparison.list = wild.rel.col.coords #4:120 is the range of columns containing wild relatives 
      
      
      #run a loop for each variety to compare the intitial variety to
      count2 = make.counter()
      lapply(comparison.list, function(q){
        print(p("comparison ", count2()))
        consec.list1 = comparison(chromo, var.coord, q)

        # browser()
        to.rm1 = which(unlist(lapply(consec.list1, function(q1) length(q1) == 1 | length(q1) == 0)))
        if(length(to.rm1) > 0) consec.list1 = consec.list1[-to.rm1]
        
        #debug code
        # g.7 = unlist((lapply(consec.list1, length)))
        # if(0 %in% g.7) print(g.7)

        #generate dataframe containing physical info. etc
      

        parse.intros = lapply(consec.list1, function(x){
          intro.start = chromo$phys.pos[min(x)]
          intro.end = chromo$phys.pos[max(x)]
          intro.length = intro.end - intro.start
          number.snps = length(x)
          snp.start = min(x)
          snp.end = max(x)
          snp.density = number.snps / intro.length
          elite = variety.to.compare
          wild.relative = colnames(chromo)[q]
          chromo1 = current.chromo
          
          all.in = c(elite, wild.relative, chromo1, intro.start, intro.end, intro.length, number.snps, snp.start, snp.end, snp.density)
          names(all.in) = c("elite", "wild.relative", "chromosome", "intro.start", "intro.end", "intro.length", "number.snps", "snp.start", "snp.end", "snp.density")
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

  # write.csv(comp1, p("rotation1scripts_v4/processed_data/introgressions/sacha/int.", var1, ".", file1), row.names = F)

  return(comp1)

}, mc.cores = 50)

all.comp2 = combine.list.of.data.frames(all.comp)

#remove introgressions with only 3 or below SNPs
all.comp2$number.snps = as.numeric(all.comp2$number.snps)
all.comp2 = all.comp2[which(all.comp2$number.snps > 3), ]

write.csv(all.comp2, "rotation1scripts_v4/processed_data/introgressions/sacha/all.introgressions_anal3_11012019.csv", row.names = F)

file1 = "data_analysis4.csv"

geno.df1 = read.csv(p(genotypes.root.path, file1), header = T, stringsAsFactors = F)

#extract sample type information from dataframe
sample.type.info = geno.df1[1, ]
#remove sample type information from dataframes
geno.df1 = geno.df1[-1, ]    

wild.rel.col.coords = which(sample.type.info == "relative")
hexaploid.col.coords = which(sample.type.info == "wheat")

all.comp = mclapply(hexaploid.col.coords, function(i2){





  geno.df1 = read.csv(p(genotypes.root.path, file1), header = T, stringsAsFactors = F)    



  # geno.df1 = geno.df1[, 2:ncol(geno.df1)]

  var1 = colnames(geno.df1)[i2]

  # print(var1)

  print(colnames(geno.df1)[1:10])
  # print(geno.df1$chromosome)

  geno.df1 = split.df.into.list.of.dfs.by.column(geno.df1, "chromosome")
  if(class(geno.df1) != "list") geno.df1 = list(geno.df1)


  # print(head(geno.df1))

  comparison = function(df1, col1, col2){
    comp1 = which(df1[, col1] == df1[, col2])
    list.consecutive.subsequences(comp1, 1)
  }

  # print(geno.df1)

  global.comparison = function(variety.to.compare, list.of.genotype.dfs1){
    #args:
    #variety.to.compare: character string indicating the exact name of the elite variety to compare to all
      #other varieties, as listed in the column name of list.of.genotype.dfs1
    #list.of.genotype.dfs1: a list of genotype dataframes split by chromosome using split.df.into.list.of.dfs.by.column
    
    count1 = make.counter()
    #run a loop for each chromosome
    main.analysis = lapply(list.of.genotype.dfs1, function(chromo){
      print(p("chromosome ", count1()))
      current.chromo = chromo$chromosome[[1]]
      chromo$phys.pos = as.numeric(chromo$phys.pos)
      var.coord = which(colnames(chromo) == variety.to.compare)
      
      comparison.list = wild.rel.col.coords #4:120 is the range of columns containing wild relatives 
      
      
      #run a loop for each variety to compare the intitial variety to
      count2 = make.counter()
      lapply(comparison.list, function(q){
        print(p("comparison ", count2()))
        consec.list1 = comparison(chromo, var.coord, q)

        # browser()
        to.rm1 = which(unlist(lapply(consec.list1, function(q1) length(q1) == 1 | length(q1) == 0)))
        if(length(to.rm1) > 0) consec.list1 = consec.list1[-to.rm1]
        
        #debug code
        # g.7 = unlist((lapply(consec.list1, length)))
        # if(0 %in% g.7) print(g.7)

        #generate dataframe containing physical info. etc
      

        parse.intros = lapply(consec.list1, function(x){
          intro.start = chromo$phys.pos[min(x)]
          intro.end = chromo$phys.pos[max(x)]
          intro.length = intro.end - intro.start
          number.snps = length(x)
          snp.start = min(x)
          snp.end = max(x)
          snp.density = number.snps / intro.length
          elite = variety.to.compare
          wild.relative = colnames(chromo)[q]
          chromo1 = current.chromo
          
          all.in = c(elite, wild.relative, chromo1, intro.start, intro.end, intro.length, number.snps, snp.start, snp.end, snp.density)
          names(all.in) = c("elite", "wild.relative", "chromosome", "intro.start", "intro.end", "intro.length", "number.snps", "snp.start", "snp.end", "snp.density")
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

  # write.csv(comp1, p("rotation1scripts_v4/processed_data/introgressions/sacha/int.", var1, ".", file1), row.names = F)

  return(comp1)

}, mc.cores = 50)

all.comp2 = combine.list.of.data.frames(all.comp)

#remove introgressions with only 3 or below SNPs
all.comp2$number.snps = as.numeric(all.comp2$number.snps)
all.comp2 = all.comp2[which(all.comp2$number.snps > 3), ]

write.csv(all.comp2, "rotation1scripts_v4/processed_data/introgressions/sacha/all.introgressions_anal4_11012019.csv", row.names = F)