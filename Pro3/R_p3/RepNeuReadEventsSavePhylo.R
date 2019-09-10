library(DDD)
library(ape)
library(parallel)
n.cores <- 4 #detectCores()

event2L.func <- 'C:/Liang/Code/Pro3/R_p3/event2L.R'
source(event2L.func)

dir <- 'C:/Liang/Googlebox/Research/Project3/replicate_neutral/sce1/'

turnover <- 1e9+1e7

ev2l <- function(eventtable){
  L <- event2L(eventtable, turnover=1e9)
  phy = DDD::L2phylo(L,dropextinct = T)
  return(phy)
}

colnames(neu.eventtable) = c('T','ns','x','y','sp','ancestor')

sce = "NEevent"
for(i in c(1)){
  for(j in c(1)){
    sim.event.list=list()
    for(rep in c(1:100)){
      comb=paste0(i,j)
      rname = paste0(dir,'1e+09/spatialpara1e+09NE11/',sce,comb,'rep',rep,'.csv')
      events = read.csv(rname,header = FALSE)
      colnames(events) = c('T','ns','x','y','sp','ancestor')
      # sim.events = rbind(neu.eventtable,events)
      sim.event.list <- c(sim.event.list,list(events))
      
    }
    multitreefile <- paste0(dir,'1e+09/spatialpara1e+09NE11/NEmultitree.tre')
    clust <- makeCluster(n.cores)
    clusterExport(clust, list("sim.event.list","event2L","L2phylo"))
    multitree <- parLapply(clust, sim.event.list, ev2l)
    stopCluster(clust)
    class(multitree) <- "multiPhylo"
    write.tree(multitree,file=multitreefile)
  }
}

phy.list = list()
for(rep in c(1:16)){
  phy.ele = ev2l(sim.event.list[[rep]])
  phy.list = c(phy.list,list(phy.ele))
}
class(phy.list) <- "multiPhylo"

