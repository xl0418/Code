source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/Nindex.R', echo=TRUE)
source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/sddsim.R', echo=TRUE)
source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/event_matrix.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/Plottree_single_Pro1.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/prunetree.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
library(DDD)
library(MASS)
library(rgl)
library(stringr)
library("reshape2")
library('Matrix')
library(plyr) 

pars=c(0.8,0.3,200)
seed_fun=29
prune=1
# result = sddsim(n=2,parsN=c(2,0),age=15,pars=pars , seed_fun = seed_fun, lambda_allo0 = 0.2, M0=0,K_fix = 1)
result = DDD::dd_sim(pars = pars,age = 15,ddmodel = 3)
L = result$L
L = L[,1:4]
dropextinct = F
phylo_a=L2phylo(L,dropextinct=dropextinct)
plot(phylo_a,show.tip.label = F)

L2 = phylo2L(phylo_a)

phylo_test = DDD::L2phylo(L2,dropextinct = dropextinct)
plot(phylo_test,show.tip.label = FALSE)

all.equal(L,L2,check.attribute = FALSE)
all.equal(phylo_a,phylo_test)
