library(DDD)
library(ape)
library(parallel)
n.cores <- detectCores()

event2L.func <- 'C:/Liang/Code/Pro3/R_p3/event2L.R'
source(event2L.func)

dir <- 'C:/Liang/Googlebox/Research/Project3/batchsim_results/'

neu.eventtablefile <- paste0(dir,'neutraldata/neuevent.csv')
neu.eventtable <- read.csv(neu.eventtablefile,header = FALSE)
turnover <- 1e9+1e7

ev2l <- function(eventtable){
  L <- event2L(eventtable, turnover=1e9+1e7)
  phy = DDD::L2phylo(L,dropextinct = T)
  return(phy)
}

colnames(neu.eventtable) = c('T','ns','x','y','sp','ancestor')

scenario = c('Levent','Mevent','Hevent')
sce.short = c('L','M','H')


for(i_n in c(1:3)){
  sce = scenario[i_n]
  f.name = sce.short[i_n]
  for(i in c(1:3)){
    for(j in c(1:6)){
      sim.event.list=list()
      for(rep in c(61:100)){
        comb=paste0(i,j)
        rname = paste0(dir,'1e+07/spatialpara1e+07',f.name,comb,'/',sce,comb,'rep',rep,'.csv')
        events = read.csv(rname,header = FALSE)
        colnames(events) = c('T','ns','x','y','sp','ancestor')
        sim.events = rbind(neu.eventtable,events)
        sim.event.list <- c(sim.event.list,list(sim.events))
        
      }
      multitreefile <- paste0(dir,'1e+07/spatialpara1e+07',f.name,comb,'/multitree',f.name,comb,'.tre')
      clust <- makeCluster(n.cores)
      clusterExport(clust, list("sim.event.list","event2L","L2phylo"))
      multitree <- parLapply(clust, sim.event.list, ev2l)
      stopCluster(clust)
      class(multitree) <- "multiPhylo"
      write.tree(multitree,file=multitreefile)
    }
  }
}

