source('C:/Liang/Code/Pro2/R_p2/prunetree.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
library(DDD)

pars=c(0.8,0.4,200)
result = DDD::dd_sim(pars = pars,age = 15,ddmodel = 3)
# original L table
L1 = result$L
L1 = L1[,1:4]
# compare the original L table with the L table resulted from converting the original L to a phylo class and to L by 
# the phylo2L function again.

# Test when not droping extinct species
dropextinct = F

phylo_L1=L2phylo(L1,dropextinct=dropextinct)
plot(phylo_L1,show.tip.label = F)

L2 = phylo2L(phylo_L1)
phylo_L2 = DDD::L2phylo(L2,dropextinct = dropextinct)
plot(phylo_L2,show.tip.label = FALSE)

T1 = all.equal(L1,L2)
T2 = all.equal(phylo_L1,phylo_L2)

# Test when droping extinct species
dropextinct = T

phylo_L1=L2phylo(L1,dropextinct=dropextinct)
plot(phylo_L1,show.tip.label = F)

L2 = phylo2L(phylo_L1)
phylo_L2 = DDD::L2phylo(L2,dropextinct = dropextinct)
plot(phylo_L2,show.tip.label = FALSE)

L_prune = prunetree(L1)
phylo_prune = DDD::L2phylo(L_prune)
plot(phylo_prune,show.tip.label = FALSE)

T3 = all.equal(L_prune,L2)

if(T1){
  print('L = phylo2L(L2phylo(L)) for full trees')
}
if(T3){
  print('L_prune = phylo2L(L2phylo(L,dropextinct = TRUE)) for pruning trees')
}
  