
lapply(all.g.2, function(x){
	g = lapply(x, function(q){
	max(q$phys.pos) / (1000000 * nrow(q))
	})
	mean(unlist(g))
})