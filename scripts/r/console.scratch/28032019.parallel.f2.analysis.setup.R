loadstatements = readLines("rotation1scripts_v4/temp/load.statements.txt")
performanalysisstatements = readLines("rotation1scripts_v4/temp/perform.analysis.statements.txt")

i = 1
g = Map(function(x, y){
    q = c("#!/home/ac14037/bin/Rscript",
                "setwd('/home/ac14037/project.phd.main/')",
                "source('rotation1scripts_v4/scripts/r/functions.R')",
                "source('rotation1scripts_v4/scripts/r/modelling/recombination.modelling.R')",
                "load('rotation1scripts_v4/saved.objects/recombination.profiles1a')", 
                x, p("g = ", y),
                "write.csv(g, p('rotation1scripts_v4/temp/parallel/analysis/analysis.', g$R.object.name, '.csv'))")
    
    outputfile = file(p("rotation1scripts_v4/temp/parallel/scripts/f2.anal.", i, ".txt"), "wb")
    writeLines(q, con = outputfile)
    close(outputfile)
    i <<- i + 1
}, loadstatements, performanalysisstatements)



writeLines(g, "rotation1scripts_v4/temp/load.w.anal.txt")

load("rotation1scripts_v4/saved.objects/list.of.no.selec2")


