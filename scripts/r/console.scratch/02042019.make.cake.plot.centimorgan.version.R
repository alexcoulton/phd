
#centimorgan version
make.cake.plot = function(sim.data){
    sim.data.f2.bar = as.data.frame(t(sapply(sim.data, geno.count)))
    sim.data.f2.bar = as.data.frame(cbind(rownames(sim.data.f2.bar), sim.data.f2.bar))
    
    cm1 = generate.cm.for.sim.data(sim.data)$pos
    browser()
    sim.data.f2.bar$cm = cm1
    sim.data.f2.bar = sim.data.f2.bar[match(unique(sim.data.f2.bar$cm), sim.data.f2.bar$cm), ]
    num.marker1 = length(unique(cm1))
    
    diff(unique(cm1))
    length(unique(cm1))
    
    
    
    s2 = melt(sim.data.f2.bar, c("cm", "rownames(sim.data.f2.bar)"))
    colnames(s2) = c("cm", "marker", "Genotype", "count")
    s2$percent = (s2$count / nrow(sim.data)) * 100
    s2$marker = rep(1:length(unique(cm1)), 3)
    s2$seg = ""
    s2$sig.text = ""
    
    s2$num.marker1 = ""
    
    s2$num.marker1 = rep(diff(c(unique(cm1), 0)), 3)
    
    seg.ratios1 = sapply(convert.recomb.to.seg.f2(sim.data), function(x) x[[3]])
    seg.marker1 = which(seg.ratios1 < 0.05)
    s2$seg[which(s2$marker %in% seg.marker1)] = T
    s2$sig.text[which(s2$marker[1:length(unique(cm1))] %in% seg.marker1)] = "*"
    
    ggplot(s2, aes(x = marker, y = percent, fill = Genotype, alpha = seg)) + geom_bar(stat = "identity", width = 1) +
        theme_classic() + scale_fill_manual(values = c("#EE6C4D", "#231942", "#98C1D9")) + 
        scale_alpha_discrete(range = c(0.6, 1)) + guides(alpha = F) + 
        geom_text(aes(label = sig.text), nudge_y = 70) + coord_cartesian(ylim = c(0, 110)) + 
        scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
        geom_hline(yintercept = 25) + geom_hline(yintercept = 75)
    #todo: need to add centimorgan positions to x-axis / keep only skeleton markers
}
