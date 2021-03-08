f6.no.selec.96pop = perform.1000.simulations.fgen(96, fgen = 6, corestouse = 60, recombination.profile = recombination.profiles1a)
f6.no.selec.96pop.analysis = perform.analysis(f6.no.selec.96pop, recombination.position = "recombination.profiles1a", fgen = 6)
write.csv(f6.no.selec.96pop.analysis, "rotation1scripts_v4/temp/f6.96anal.csv", row.names = F)
s(f6.no.selec.96pop, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/f6.no.selec.96pop", "NA")
rm(f6.no.selec.96pop)
gc()



f6.no.selec.1000pop = perform.1000.simulations.fgen(1000, fgen = 6, corestouse = 60, recombination.profile = recombination.profiles1a)
f6.no.selec.1000pop.analysis = perform.analysis(f6.no.selec.1000pop, recombination.position = "recombination.profiles1a", fgen = 6)
write.csv(f6.no.selec.96pop.analysis, "rotation1scripts_v4/temp/f6.1000anal.csv", row.names = F)
s(f6.no.selec.1000pop, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/f6.no.selec.1000pop", "NA")
rm(f6.no.selec.1000pop)
gc()

f6.no.selec.10000pop = perform.1000.simulations.fgen(10000, fgen = 6, corestouse = 60, recombination.profile = recombination.profiles1a)
f6.no.selec.10000pop.analysis = perform.analysis(f6.no.selec.10000pop, recombination.position = "recombination.profiles1a", fgen = 6)
write.csv(f6.no.selec.10000pop.analysis, "rotation1scripts_v4/temp/f6.10000anal.csv", row.names = F)
s(f6.no.selec.10000pop, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/f6.no.selec.10000pop", "NA")
rm(f6.no.selec.10000pop)
gc()


f6.no.selec.1000pop = perform.1000.simulations.fgen(1000, fgen = 6, corestouse = 60, recombination.profile = recombination.profiles1a)
f6.no.selec.1000pop.analysis = perform.analysis(f6.no.selec.1000pop, recombination.position = "recombination.profiles1a", fgen = 6)





f6.no.selec.96pop = perform.1000.simulations.fgen(96, fgen = 6, corestouse = 60, recombination.profile = recombination.profiles1a)
f6.no.selec.96pop.analysis = perform.analysis(f6.no.selec.96pop, recombination.position = "recombination.profiles1a", fgen = 6)
write.csv(f6.no.selec.96pop.analysis, "rotation1scripts_v4/temp/f6.96anal.csv", row.names = F)
s(f6.no.selec.96pop, "rotation1scripts_v4/saved.objects/recomb.sim/fgen/f6.no.selec.96pop", "NA")
rm(f6.no.selec.96pop)
gc()
