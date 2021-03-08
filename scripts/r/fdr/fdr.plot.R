#illustrating FDR geometrically 
#https://en.wikipedia.org/wiki/False_discovery_rate
#https://community.jmp.com/t5/JMP-Blog/Not-just-filtering-coincidences-False-discovery-rate/ba-p/30325#

g = abs(rnorm(100, 0.5, 0.2))

g = sort(g)
plot(g)
abline(0, (0.05 / 100))
abline(0.05, 0)
abline((0.05 / 100), 0)


g2 = sort(rexp(100, 1000000)^2.1)
g2 = g2 / max(g2)

plot(g2)

