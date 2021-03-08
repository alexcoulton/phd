# load original f2s

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec12pop96")

axcf2selec12pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec14pop96")

axcf2selec14pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec16pop96")

axcf2selec16pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec18pop96")

axcf2selec18pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf2.sims/axcf2selec20pop96")

axcf2selec20pop96 = all.sim.files

#load f2 batches


load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf2selec12pop96.batch2")

axcf2selec12pop96.batch2 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf2selec14pop96.batch2")

axcf2selec14pop96.batch2 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf2selec16pop96.batch2")

axcf2selec16pop96.batch2 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf2selec18pop96.batch2")

axcf2selec18pop96.batch2 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf2selec20pop96.batch2")

axcf2selec20pop96.batch2 = all.sim.files



axcf2selec12pop96.comb = c(axcf2selec12pop96, axcf2selec12pop96.batch2)
axcf2selec14pop96.comb = c(axcf2selec14pop96, axcf2selec14pop96.batch2)
axcf2selec16pop96.comb = c(axcf2selec16pop96, axcf2selec16pop96.batch2)
axcf2selec18pop96.comb = c(axcf2selec18pop96, axcf2selec18pop96.batch2)
axcf2selec20pop96.comb = c(axcf2selec20pop96, axcf2selec20pop96.batch2)


axcf2selec12pop96.comb = axcf2selec12pop96.comb[1:1000]
axcf2selec14pop96.comb = axcf2selec14pop96.comb[1:1000]
axcf2selec16pop96.comb = axcf2selec16pop96.comb[1:1000]
axcf2selec18pop96.comb = axcf2selec18pop96.comb[1:1000]
axcf2selec20pop96.comb = axcf2selec20pop96.comb[1:1000]

allf2.redone = list(axcf2selec12pop96.comb, axcf2selec14pop96.comb, axcf2selec16pop96.comb, axcf2selec18pop96.comb, axcf2selec20pop96.comb)

allf2.anal = Map(function(x, y){
	perform.analysis(x, T, 200, y, "AxC1A", T, 2)
}, allf2.redone, seq(12, 20, 2))


# load original f6 files



load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6selec20.pop96")

axcf6selec20.pop96 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6selec20.pop300")

axcf6selec20.pop300 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6selec20.pop1000")

axcf6selec20.pop1000 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/axcf6.sims/axcf6selec20.pop10000")

axcf6selec20.pop10000 = all.sim.files

#load f6 batch 2 files


load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf6selec20.pop96.batch2")

axcf6selec20.pop96.batch2 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/batch2.sims/axcf6selec20.pop1000.batch2")

axcf6selec20.pop1000.batch2 = all.sim.files



axcf6selec20.pop96.comb = c(axcf6selec20.pop96, axcf6selec20.pop96.batch2)
axcf6selec20.pop1000.comb = c(axcf6selec20.pop1000, axcf6selec20.pop1000.batch2)

axcf6selec20.pop96.comb = axcf6selec20.pop96.comb[1:1000]

s(axcf6selec20.pop96.comb, "rotation1scripts_v4/saved.objects/axcf6selec20.pop96.comb", "15042019.combine.pedsim.batches.R")
s(axcf6selec20.pop1000.comb, "rotation1scripts_v4/saved.objects/axcf6selec20.pop1000.comb", "15042019.combine.pedsim.batches.R")

axcf6selec20.pop1000.comb = axcf6selec20.pop1000.comb[1:1000]
axcf6selec20.pop300 = axcf6selec20.pop300[1:1000]
axcf6selec20.pop10000 = axcf6selec20.pop10000[1:1000]





all.data1 = list(axcf6selec20.pop96.comb, axcf6selec20.pop1000.comb, axcf6selec20.pop300, axcf6selec20.pop10000)

allf6.anal = lapply(all.data1, function(x){
	perform.analysis(x, T, 200, 20, "AxC1A", F, 6)
})


load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.dense.noselec.pop1000")

axcf2.dense.noselec.pop1000 = all.sim.files

load("rotation1scripts_v4/original_data/simulation/pedigreesim/simresults/axcf2.dense.noselec.pop1000.batch2")

axcf2.dense.noselec.pop1000.batch2 = all.sim.files


axcf2.dense.noselec.pop1000 = c(axcf2.dense.noselec.pop1000, axcf2.dense.noselec.pop1000.batch2)
axcf2.dense.noselec.pop1000 = axcf2.dense.noselec.pop1000[1:1000]

dense2 = perform.analysis(axcf2.dense.noselec.pop1000, recombination.position = "AxC1A", F2 = T, fgen = 2)




axcf2.sparse.2markers_2.noselec.pop1000
axcf2.sparse.2markers_2.noselec.pop300
axcf2.sparse.2markers_2.noselec.pop96










