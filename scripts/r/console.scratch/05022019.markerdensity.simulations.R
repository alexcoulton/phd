marker.sim1 = perform.1000.simulations(300, number.of.markers = 10, recombination.profile = recombination.profiles1a.reduced, corestouse = 50)
marker.sim2 = perform.1000.simulations(1000, number.of.markers = 10, recombination.profile = recombination.profiles1a.reduced, corestouse = 50)
marker.sim3 = perform.1000.simulations(10000, number.of.markers = 10, recombination.profile = recombination.profiles1a.reduced, corestouse = 50)

number.dist.markers(marker.sim1)
number.dist.markers(marker.sim2)
number.dist.markers(marker.sim3)





marker.sim.further1 = perform.1000.simulations(300, number.of.markers = 5, recombination.profile = recombination.profiles1a.reduced.further, corestouse = 50)
marker.sim.further2 = perform.1000.simulations(1000, number.of.markers = 5, recombination.profile = recombination.profiles1a.reduced.further, corestouse = 50)
marker.sim.further3 = perform.1000.simulations(10000, number.of.markers = 5, recombination.profile = recombination.profiles1a.reduced.further, corestouse = 50)

number.dist.markers(marker.sim.further1)
number.dist.markers(marker.sim.further2)
number.dist.markers(marker.sim.further3)




marker.sim.any1 = perform.1000.simulations(300, number.of.markers = 5, corestouse = 50)
marker.sim.any2 = perform.1000.simulations(1000, number.of.markers = 5, corestouse = 50)
marker.sim.any3 = perform.1000.simulations(10000, number.of.markers = 5, corestouse = 50)

number.dist.markers(marker.sim.any1)
number.dist.markers(marker.sim.any2)
number.dist.markers(marker.sim.any3)

pop.sizes = c(300, 1000, 10000)

sims1kmarker = lapply(pop.sizes, function(x) perform.1000.simulations(x, number.of.markers = 1000, corestouse = 50))






recomb.profile1a.no.selec.sims = perform.1000.simulations(300, number.of.markers = 225, recombination.profile = recombination.profiles1a)


list.of.first.gen.w.prob.w.selec.10000pop7 = mclapply(1:1000, function(x){
    g = gen.data(225, 10000, recombination.profiles = recombination.profiles1a, selection = 3, selection.pos = 200)
    g
}, mc.cores = corestouse)


