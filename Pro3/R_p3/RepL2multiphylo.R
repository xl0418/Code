library(DDD)
library(ape)

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
ev2l <- function(ltable){
  if(ltable[1,2] == 0 && ltable[1,3] == 1){
    ltable[,2] = ltable[,2]+1
    ltable[,3] = ltable[,3]+1
    ltable = rbind(c(10^7,0,1,-1),ltable)
    phy = DDD::L2phylo(ltable)
    return(phy)
  }else if(ltable[1,2] == 1 && ltable[1,3] == 1){
    ltable[1,1] = 10^7
    ltable[1,2] = 0
    ltable[,2] = ltable[,2]+1
    ltable[,3] = ltable[,3]+1
    ltable = rbind(c(10^7,0,1,-1),ltable)
    phy = DDD::L2phylo(ltable)
    return(phy)
  }else{
    print(ltable[1,])
    
  }
  
}

for(i_n in c(2:9)){
  scefolder = scenario[i_n]
  letter.comb = sce.short.comb.vec[i_n]
  for(i in c(1:5)){
    for(j in c(1:6)){
      print(paste(scefolder,i,j))
      multitree=list()
      comb=paste0(i,j)
      
      for(rep in c(1:100)){

        rname = paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/',letter.comb,'psi',i,'s_phi',j,'rep',rep,'Ltable.csv')
        L.table = read.csv(rname,header = FALSE)
        # print(nrow(L.table))
        single.phylo = ev2l(L.table)
        multitree = c(multitree,list(single.phylo))
      }
      multitreefile <- paste0(dir,scefolder,'/results/1e+07/spatialpara1e+07',letter.comb,comb,'/','multitree',letter.comb,comb,'.tre')

      class(multitree) <- "multiPhylo"
      write.tree(multitree,file=multitreefile)
    }
  }
}

