library(DDD)
library(ape)
library(parallel)
n.cores <- detectCores()

event2L.func <- 'C:/Liang/Code/Pro3/R_p3/event2L.R'
source(event2L.func)

dir <- 'C:/Liang/Googlebox/Research/Project3/replicate_sim_final1/'
dir = 'C:/Liang/Googlebox/Research/Project3/replicate_sim_9sces/'
sce.short = c('H','M','L')
scenario = NULL
sce.short.comb.vec = NULL
for(i.letter in sce.short){
  for(j.letter in sce.short){
    sce.folder = paste0('sce',i.letter,j.letter)
    scenario = c(scenario,sce.folder)
    sce.short.comb = paste0(i.letter,j.letter)
    sce.short.comb.vec = c(sce.short.comb.vec,sce.short.comb)
  }
}
scenario = c('Levent','Mevent','Hevent')
ev2l <- function(ltable){
  ltable[,2] = ltable[,2]+1
  ltable[,3] = ltable[,3]+1
  ltable = rbind(c(10^7,0,1,-1),ltable)
  phy = DDD::L2phylo(ltable)
  return(phy)
}

for(i_n in c(2)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  sce = scenario[i_n]
  
  for(i in c(1)){
    for(j in c(1)){
      sim.event.list=list()
      comb=paste0(i,j)
      
      for(rep in c(1:100)){
        rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'psi',i,'s_phi',j,'rep',rep,'Ltable.csv')
        L.table = read.csv(rname,header = FALSE)
        sim.events = rbind(neu.eventtable,events)
        sim.event.list <- c(sim.event.list,list(sim.events))
        
      }
      multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',f.name,comb,'.tre')
      clust <- makeCluster(n.cores)
      clusterExport(clust, list("sim.event.list","L2phylo"))
      multitree <- parLapply(clust, sim.event.list, ev2l)
      stopCluster(clust)
      class(multitree) <- "multiPhylo"
      write.tree(multitree,file=multitreefile)
    }
  }
}

